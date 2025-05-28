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

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::missing_safety_doc)] // TODO: Add safety docs
#![allow(clippy::upper_case_acronyms)] // For FFI function names
#![allow(clippy::collapsible_else_if)] // For bbox_intersect
#![allow(clippy::needless_return)] // To match Fortran structure more easily in some cases

use libc::{c_double, c_int, c_void}; 
use std::ptr;

use crate::status;
use crate::helpers::{self, VECTOR_CLOSE_EPS, WIGGLE}; 
use crate::curve;   

// Constants
pub const BOX_INTERSECTION_TYPE_INTERSECTION: c_int = 0;
pub const BOX_INTERSECTION_TYPE_TANGENT: c_int = 1;
pub const BOX_INTERSECTION_TYPE_DISJOINT: c_int = 2;

pub const SUBDIVIDE_FIRST: c_int = 0;
pub const SUBDIVIDE_SECOND: c_int = 1;
pub const SUBDIVIDE_BOTH: c_int = 2;
pub const SUBDIVIDE_NEITHER: c_int = -1;

pub const LINEARIZATION_THRESHOLD: f64 = 1.4901161193847656e-08; // 0.5_f64.powi(26)
pub const MAX_INTERSECT_SUBDIVISIONS: c_int = 20;
pub const MIN_INTERVAL_WIDTH: f64 = 9.094947017729282e-13; // 0.5_f64.powi(40)
pub const MAX_CANDIDATES: c_int = 64; 
pub const ZERO_THRESHOLD: f64 = 0.0009765625; // 0.5_f64.powi(10)
pub const NEWTON_ERROR_RATIO: f64 = 1.4551915228366852e-11; // 0.5_f64.powi(36)

// Helper for column-major access: get an element nodes(dim_idx, node_idx)
#[inline]
unsafe fn get_node_val(nodes: *const c_double, dimension: usize, node_idx: usize, dim_idx: usize) -> f64 {
    *nodes.add(node_idx * dimension + dim_idx)
}

// Helper for column-major access: set an element nodes(dim_idx, node_idx) = val
#[inline]
unsafe fn set_node_val(nodes: *mut c_double, dimension: usize, node_idx: usize, dim_idx: usize, val: f64) {
    *nodes.add(node_idx * dimension + dim_idx) = val;
}

// Helper for L2 norm of a 2D vector (components are contiguous)
#[inline]
fn norm2_2d(vec_ptr: *const c_double) -> f64 {
    unsafe {
        let x = *vec_ptr.add(0);
        let y = *vec_ptr.add(1);
        (x * x + y * y).sqrt()
    }
}


unsafe fn linearization_error_internal(
    num_nodes: usize,
    dimension: usize, 
    nodes: *const c_double, 
    error: &mut f64,
) {
    if num_nodes <= 2 { 
        *error = 0.0;
        return;
    }

    let num_second_deriv_nodes = num_nodes - 2;
    let mut second_deriv_vec: Vec<f64> = vec![0.0; dimension * num_second_deriv_nodes];
    let sd_ptr = second_deriv_vec.as_mut_ptr();

    for i in 0..num_second_deriv_nodes { 
        for d in 0..dimension {
            let val = get_node_val(nodes, dimension, i, d) -
                      2.0 * get_node_val(nodes, dimension, i + 1, d) +
                      get_node_val(nodes, dimension, i + 2, d);
            set_node_val(sd_ptr, dimension, i, d, val);
        }
    }
    
    let mut worst_case_vec: Vec<f64> = vec![0.0; dimension];
    for d in 0..dimension {
        let mut max_abs_val_dim = 0.0;
        for i in 0..num_second_deriv_nodes {
            let current_abs_val = get_node_val(sd_ptr, dimension, i, d).abs();
            if current_abs_val > max_abs_val_dim {
                max_abs_val_dim = current_abs_val;
            }
        }
        worst_case_vec[d] = max_abs_val_dim;
    }
    
    let norm_worst_case = norm2_2d(worst_case_vec.as_ptr()); 

    *error = 0.125 * (num_nodes - 1) as f64 * (num_nodes - 2) as f64 * norm_worst_case;
}

unsafe fn segment_intersection_internal(
    start0_ptr: *const c_double, 
    end0_ptr: *const c_double,   
    start1_ptr: *const c_double, 
    end1_ptr: *const c_double,   
    s: &mut f64,
    t: &mut f64,
) -> bool { 
    let delta0_x = *end0_ptr.add(0) - *start0_ptr.add(0);
    let delta0_y = *end0_ptr.add(1) - *start0_ptr.add(1);
    let delta0 = [delta0_x, delta0_y];

    let delta1_x = *end1_ptr.add(0) - *start1_ptr.add(0);
    let delta1_y = *end1_ptr.add(1) - *start1_ptr.add(1);
    let delta1 = [delta1_x, delta1_y];
    
    let mut cross_d0_d1 = 0.0;
    helpers::BEZ_cross_product(delta0.as_ptr(), delta1.as_ptr(), &mut cross_d0_d1);

    if cross_d0_d1.abs() < WIGGLE { 
        return false; 
    }
    
    let start_delta_x = *start1_ptr.add(0) - *start0_ptr.add(0);
    let start_delta_y = *start1_ptr.add(1) - *start0_ptr.add(1);
    let start_delta = [start_delta_x, start_delta_y];
    
    let mut other_cross_s = 0.0;
    helpers::BEZ_cross_product(start_delta.as_ptr(), delta1.as_ptr(), &mut other_cross_s);
    *s = other_cross_s / cross_d0_d1;
    
    let mut other_cross_t = 0.0;
    helpers::BEZ_cross_product(start_delta.as_ptr(), delta0.as_ptr(), &mut other_cross_t);
    *t = other_cross_t / cross_d0_d1;
    return true;
}

unsafe fn parallel_lines_parameters_internal(
    start0_ptr: *const c_double,
    end0_ptr: *const c_double,
    start1_ptr: *const c_double,
    end1_ptr: *const c_double,
    s_vals: *mut c_double, 
    t_vals: *mut c_double, 
    num_written: &mut c_int,
) {
    *num_written = 0;
    let dim: usize = 2;

    let mut start_diff_vec = [0.0; 2];
    let mut seg0_vec = [0.0; 2];
    for i in 0..dim {
        start_diff_vec[i] = *start1_ptr.add(i) - *start0_ptr.add(i);
        seg0_vec[i] = *end0_ptr.add(i) - *start0_ptr.add(i);
    }

    let mut cross_collinear = 0.0;
    helpers::BEZ_cross_product(start_diff_vec.as_ptr(), seg0_vec.as_ptr(), &mut cross_collinear);

    if cross_collinear.abs() >= WIGGLE { 
        return;
    }

    let dot_seg0_seg0 = seg0_vec[0] * seg0_vec[0] + seg0_vec[1] * seg0_vec[1];
    if dot_seg0_seg0 < WIGGLE { 
        if start_diff_vec[0].abs() < WIGGLE && start_diff_vec[1].abs() < WIGGLE {
            let seg1_dx = *end1_ptr.add(0) - *start1_ptr.add(0);
            let seg1_dy = *end1_ptr.add(1) - *start1_ptr.add(1);
            if seg1_dx.abs() < WIGGLE && seg1_dy.abs() < WIGGLE {
                *s_vals.add(0) = 0.0; 
                *t_vals.add(0) = 0.0; 
                *num_written = 1;
            }
        }
        return;
    }

    let s_val_start1 = (start_diff_vec[0] * seg0_vec[0] + start_diff_vec[1] * seg0_vec[1]) / dot_seg0_seg0;
    
    let mut end1_start0_diff_vec = [0.0; 2];
    for i in 0..dim {
        end1_start0_diff_vec[i] = *end1_ptr.add(i) - *start0_ptr.add(i);
    }
    let s_val_end1 = (end1_start0_diff_vec[0] * seg0_vec[0] + end1_start0_diff_vec[1] * seg0_vec[1]) / dot_seg0_seg0;

    let s_min = s_val_start1.min(s_val_end1);
    let s_max = s_val_start1.max(s_val_end1);

    let overlap_min_s = (0.0f64).max(s_min);
    let overlap_max_s = (1.0f64).min(s_max);

    if overlap_min_s > overlap_max_s + WIGGLE { 
        return;
    }

    let seg1_dx = *end1_ptr.add(0) - *start1_ptr.add(0);
    let seg1_dy = *end1_ptr.add(1) - *start1_ptr.add(1);
    let seg1_len_sq = seg1_dx*seg1_dx + seg1_dy*seg1_dy;

    if seg1_len_sq < WIGGLE { 
        if overlap_min_s <= s_val_start1 && s_val_start1 <= overlap_max_s { 
            *s_vals.add(0) = s_val_start1;
            *t_vals.add(0) = 0.0; 
            *num_written = 1;
        }
        return;
    }

    *s_vals.add(0) = overlap_min_s;
    if (s_val_end1 - s_val_start1).abs() < WIGGLE { 
        *t_vals.add(0) = 0.0;
    } else {
        let t0 = (overlap_min_s - s_val_start1) / (s_val_end1 - s_val_start1);
        *t_vals.add(0) = t0.clamp(0.0, 1.0);
    }
    *num_written = 1;

    if (overlap_max_s - overlap_min_s).abs() > WIGGLE {
        *s_vals.add(1) = overlap_max_s;
         if (s_val_end1 - s_val_start1).abs() < WIGGLE {
            *t_vals.add(1) = 0.0;
        } else {
            let t1 = (overlap_max_s - s_val_start1) / (s_val_end1 - s_val_start1);
            *t_vals.add(1) = t1.clamp(0.0, 1.0);
        }
        *num_written = 2;
    }
}

unsafe fn line_line_collide_internal(
    start0_ptr: *const c_double,
    end0_ptr: *const c_double,
    start1_ptr: *const c_double,
    end1_ptr: *const c_double,
    s_vals: *mut c_double, 
    t_vals: *mut c_double, 
    num_written: &mut c_int,
) {
    let mut s: f64 = 0.0;
    let mut t: f64 = 0.0;

    if segment_intersection_internal(start0_ptr, end0_ptr, start1_ptr, end1_ptr, &mut s, &mut t) {
        if s >= -WIGGLE && s <= 1.0 + WIGGLE && t >= -WIGGLE && t <= 1.0 + WIGGLE {
            *s_vals.add(0) = s.clamp(0.0, 1.0);
            *t_vals.add(0) = t.clamp(0.0, 1.0);
            *num_written = 1;
        } else {
            *num_written = 0;
        }
    } else { 
        parallel_lines_parameters_internal(start0_ptr, end0_ptr, start1_ptr, end1_ptr, s_vals, t_vals, num_written);
    }
}


#[no_mangle]
pub unsafe extern "C" fn BEZ_newton_refine_curve_intersect(
    s_in: c_double,
    num_nodes1: c_int,
    nodes1: *const c_double, 
    t_in: c_double,
    num_nodes2: c_int,
    nodes2: *const c_double, 
    new_s: *mut c_double,
    new_t: *mut c_double,
    status_code: *mut c_int, 
) {
    let dim: usize = 2; 

    let mut func_val_vec: Vec<f64> = vec![0.0; dim]; 
    let mut b1_s_vec: Vec<f64> = vec![0.0; dim];

    let s_arr = [s_in];
    let t_arr = [t_in];

    curve::BEZ_evaluate_multi(
        num_nodes2, dim as c_int, nodes2, 1, t_arr.as_ptr(), func_val_vec.as_mut_ptr()
    );
    curve::BEZ_evaluate_multi(
        num_nodes1, dim as c_int, nodes1, 1, s_arr.as_ptr(), b1_s_vec.as_mut_ptr()
    );

    let mut all_zero = true;
    for i in 0..dim {
        func_val_vec[i] -= b1_s_vec[i]; 
        if func_val_vec[i].abs() > helpers::WIGGLE { 
            all_zero = false;
        }
    }

    if all_zero {
        *new_s = s_in;
        *new_t = t_in;
        *status_code = status::STATUS_SUCCESS;
        return;
    }

    let mut jac_mat_vec: Vec<f64> = vec![0.0; dim * dim]; 
    let jac_mat_ptr = jac_mat_vec.as_mut_ptr();
    let mut b2_prime_t_vec: Vec<f64> = vec![0.0; dim];


    curve::BEZ_evaluate_hodograph(s_in, num_nodes1, dim as c_int, nodes1, jac_mat_ptr); 
    curve::BEZ_evaluate_hodograph(t_in, num_nodes2, dim as c_int, nodes2, b2_prime_t_vec.as_mut_ptr()); 

    for i in 0..dim {
        *jac_mat_ptr.add(dim + i) = -b2_prime_t_vec[i];
    }
    
    let mut delta_s: f64 = 0.0;
    let mut delta_t: f64 = 0.0;
    let mut singular_bool: bool = false;

    helpers::solve2x2_rs(
        jac_mat_ptr, 
        func_val_vec.as_ptr(), 
        &mut singular_bool, 
        &mut delta_s, 
        &mut delta_t
    );
    
    if singular_bool {
        *status_code = status::STATUS_SINGULAR;
        *new_s = s_in; 
        *new_t = t_in;
    } else {
        *status_code = status::STATUS_SUCCESS;
        *new_s = s_in + delta_s;
        *new_t = t_in + delta_t;
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_bbox_intersect(
    num_nodes1: c_int,
    nodes1: *const c_double, 
    num_nodes2: c_int,
    nodes2: *const c_double, 
    enum_: *mut c_int,       
) {
    let mut left1: f64 = 0.0; let mut right1: f64 = 0.0; 
    let mut bottom1: f64 = 0.0; let mut top1: f64 = 0.0;
    let mut left2: f64 = 0.0; let mut right2: f64 = 0.0;
    let mut bottom2: f64 = 0.0; let mut top2: f64 = 0.0;

    helpers::BEZ_bbox(num_nodes1, nodes1, &mut left1, &mut right1, &mut bottom1, &mut top1);
    helpers::BEZ_bbox(num_nodes2, nodes2, &mut left2, &mut right2, &mut bottom2, &mut top2);

    if right2 < left1 - WIGGLE || right1 < left2 - WIGGLE || 
       top2 < bottom1 - WIGGLE || top1 < bottom2 - WIGGLE {
        *enum_ = BOX_INTERSECTION_TYPE_DISJOINT;
    } else {
        let x_overlap = left1 < right2 - WIGGLE && left2 < right1 - WIGGLE;
        let y_overlap = bottom1 < top2 - WIGGLE && bottom2 < top1 - WIGGLE;
        if x_overlap && y_overlap {
            *enum_ = BOX_INTERSECTION_TYPE_INTERSECTION;
        } else {
            *enum_ = BOX_INTERSECTION_TYPE_TANGENT;
        }
    }
}


#[no_mangle]
pub unsafe extern "C" fn BEZ_free_curve_intersections_workspace() {
    // No-op for now.
}

// TODO: Implement other FFI functions and internal helpers.
// convex_hull_collide, 
// newton_simple_root, newton_double_root, newton_iterate, full_newton_nonzero, full_newton,
// from_linearized, bbox_line_intersect, check_lines, add_intersection,
// add_from_linearized, endpoint_check, tangent_bbox_intersection,
// add_candidates, intersect_one_round, make_same_degree,
// add_coincident_parameters, all_intersections, all_intersections_abi.
