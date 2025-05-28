// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#![allow(dead_code)] // TODO: Remove once all functions are used
#![allow(clippy::too_many_arguments)]
#![allow(clippy::missing_safety_doc)] // TODO: Add safety docs
#![allow(clippy::upper_case_acronyms)] // For FFI function names
// #![allow(clippy::collapsible_else_if)] 
// #![allow(clippy::needless_return)] 

use libc::{c_double, c_int, c_void};
use std::ptr;

use crate::status;
use crate::helpers::{self, WIGGLE}; 
use crate::curve;   
use crate::triangle; 
// use crate::curve_intersection; // Might be needed later

// Constants for IntersectionClassification
pub const INTERSECTION_CLASSIFICATION_UNSET: c_int = -99;
pub const INTERSECTION_CLASSIFICATION_FIRST: c_int = 0;
pub const INTERSECTION_CLASSIFICATION_SECOND: c_int = 1;
pub const INTERSECTION_CLASSIFICATION_OPPOSED: c_int = 2;
pub const INTERSECTION_CLASSIFICATION_TANGENT_FIRST: c_int = 3;
pub const INTERSECTION_CLASSIFICATION_TANGENT_SECOND: c_int = 4;
pub const INTERSECTION_CLASSIFICATION_IGNORED_CORNER: c_int = 5;
pub const INTERSECTION_CLASSIFICATION_TANGENT_BOTH: c_int = 6;
pub const INTERSECTION_CLASSIFICATION_COINCIDENT: c_int = 7;
pub const INTERSECTION_CLASSIFICATION_COINCIDENT_UNUSED: c_int = 8;

// Constants for TriangleContained
pub const TRIANGLE_CONTAINED_NEITHER: c_int = 0;
pub const TRIANGLE_CONTAINED_FIRST: c_int = 1;
pub const TRIANGLE_CONTAINED_SECOND: c_int = 2;

// Other constants
pub const MAX_LOCATE_SUBDIVISIONS_TRIANGLE: c_int = 20; 
pub const LOCATE_EPS_TRIANGLE: f64 = 7.105427357601002e-15; // 0.5^47
pub const MAX_EDGES: c_int = 10; // Used in interior_combine for loop bound
pub const ALMOST_TANGENT: f64 = 8.881784197001252e-16; // 0.5^50


#[repr(C)]
#[derive(Debug, Clone, Copy)] 
pub struct CurvedPolygonSegment {
    pub start: c_double,
    pub end_: c_double, 
    pub edge_index: c_int,
}

// Internal (non-FFI) struct, equivalent to Fortran's `Intersection` type
#[derive(Debug, Clone, Copy)]
struct Intersection {
    s: f64,
    t: f64,
    index_first: i32,  
    index_second: i32,
    interior_curve: i32, 
}

impl Default for Intersection {
    fn default() -> Self {
        Intersection {
            s: -1.0,
            t: -1.0,
            index_first: -1,
            index_second: -1,
            interior_curve: INTERSECTION_CLASSIFICATION_UNSET,
        }
    }
}

// Internal (non-FFI) struct, equivalent to Fortran's `LocateCandidate` type
#[derive(Clone)]
struct InternalTriangleCandidate {
    centroid_x: f64, 
    centroid_y: f64, 
    width: f64,      
    nodes: Vec<f64>, 
}

impl InternalTriangleCandidate {
    fn new_empty(num_nodes: usize, dimension: usize) -> Self {
        InternalTriangleCandidate {
            centroid_x: 0.0, 
            centroid_y: 0.0, 
            width: 0.0,      
            nodes: vec![0.0; num_nodes * dimension],
        }
    }
     fn new_from_nodes(
        cx: f64, cy: f64, w: f64,
        nodes_slice: &[f64] 
    ) -> Self {
        InternalTriangleCandidate {
            centroid_x: cx,
            centroid_y: cy,
            width: w,
            nodes: nodes_slice.to_vec(),
        }
    }
}

// Helper for newton_refine_triangle
unsafe fn newton_refine_solve_internal(
    jac_both_eval_ptr: *const c_double, // Pointer to [dBds_x, dBds_y, dBdt_x, dBdt_y]
    x_val: f64,
    triangle_x: f64,
    y_val: f64,
    triangle_y: f64,
    delta_s: &mut f64,
    delta_t: &mut f64,
    singular: &mut bool,
) {
    let e_val = x_val - triangle_x;
    let f_val = y_val - triangle_y;

    let j11 = *jac_both_eval_ptr.add(0); // dBds_x
    let j21 = *jac_both_eval_ptr.add(1); // dBds_y
    let j12 = *jac_both_eval_ptr.add(2); // dBdt_x
    let j22 = *jac_both_eval_ptr.add(3); // dBdt_y
    
    let denominator = j11 * j22 - j21 * j12;

    if denominator.abs() < WIGGLE {
        *singular = true;
        // delta_s and delta_t remain as they were (e.g., 0.0 if initialized)
    } else {
        *singular = false;
        *delta_s = (j22 * e_val - j12 * f_val) / denominator;
        *delta_t = (j11 * f_val - j21 * e_val) / denominator;
    }
}


// FFI Functions
#[no_mangle]
pub unsafe extern "C" fn BEZ_newton_refine_triangle(
    num_nodes_c: c_int,
    nodes: *const c_double, // (2, num_nodes)
    degree_c: c_int,
    x_val: c_double,
    y_val: c_double,
    s: c_double,
    t: c_double,
    updated_s: *mut c_double,
    updated_t: *mut c_double,
) {
    let dim: usize = 2; 
    let mut point_eval_vec: Vec<f64> = vec![0.0; dim]; 
    let mut jac_both_eval_vec: Vec<f64> = vec![0.0; 2 * dim]; 
    
    triangle::BEZ_evaluate_barycentric(
        num_nodes_c, dim as c_int, nodes, degree_c,
        1.0 - s - t, s, t,
        point_eval_vec.as_mut_ptr()
    );

    if (point_eval_vec[0] - x_val).abs() < WIGGLE && 
       (point_eval_vec[1] - y_val).abs() < WIGGLE {
        *updated_s = s;
        *updated_t = t;
        return;
    }

    let jac_degree = degree_c - 1;
    let num_deriv_nodes = degree_c * (degree_c + 1) / 2; // Num nodes for triangle of degree `degree_c - 1`

    let mut jac_nodes_vec: Vec<f64> = vec![0.0; (2 * dim) * num_deriv_nodes as usize]; 

    triangle::BEZ_jacobian_both(
        num_nodes_c, dim as c_int, nodes, degree_c,
        jac_nodes_vec.as_mut_ptr()
    );

    triangle::BEZ_evaluate_barycentric(
        num_deriv_nodes, 
        (2 * dim) as c_int, 
        jac_nodes_vec.as_ptr(),
        jac_degree, 
        1.0 - s - t, s, t,
        jac_both_eval_vec.as_mut_ptr() 
    );
    
    let mut delta_s: f64 = 0.0;
    let mut delta_t: f64 = 0.0;
    let mut singular: bool = false;

    newton_refine_solve_internal(
        jac_both_eval_vec.as_ptr(),
        x_val, point_eval_vec[0],
        y_val, point_eval_vec[1],
        &mut delta_s, &mut delta_t,
        &mut singular
    );

    if singular {
        *updated_s = s;
        *updated_t = t;
    } else {
        *updated_s = s + delta_s;
        *updated_t = t + delta_t;
    }
}


#[no_mangle]
pub unsafe extern "C" fn BEZ_free_triangle_intersections_workspace() {
    // No-op in Rust if all workspace memory is managed by Vecs within function scopes.
}

// TODO: Implement BEZ_locate_point_triangle and its helpers
// TODO: Implement BEZ_triangle_intersections and its helpers
