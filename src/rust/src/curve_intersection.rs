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
// #![allow(clippy::missing_safety_doc)] // All functions should have safety docs or be safe
#![allow(clippy::upper_case_acronyms)] // For FFI function names
#![allow(clippy::collapsible_else_if)] // For bbox_intersect
#![allow(clippy::needless_return)] // To match Fortran structure more easily in some cases
#![allow(clippy::let_and_return)] // For direct translation of some Fortran logic
#![allow(clippy::comparison_chain)] // For add_coincident_parameters logic matching

use libc::{c_double, c_int, c_void}; 
use std::ptr;

use crate::status;
use crate::helpers::{self, VECTOR_CLOSE_EPS, WIGGLE}; 
use crate::curve::{self, CurveData}; 

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
/// # Safety
/// Caller must ensure `nodes` points to a valid memory block of at least
/// `(node_idx * dimension + dim_idx + 1)` `c_double` elements.
/// `dimension` must not be zero if `node_idx > 0` or `dim_idx > 0`.
/// `dim_idx` must be less than `dimension`.
unsafe fn get_node_val(nodes: *const c_double, dimension: usize, node_idx: usize, dim_idx: usize) -> f64 {
    *nodes.add(node_idx * dimension + dim_idx)
}

// Helper for column-major access: set an element nodes(dim_idx, node_idx) = val
#[inline]
/// # Safety
/// Caller must ensure `nodes` points to a valid, mutable memory block of at least
/// `(node_idx * dimension + dim_idx + 1)` `c_double` elements.
/// `dimension` must not be zero if `node_idx > 0` or `dim_idx > 0`.
/// `dim_idx` must be less than `dimension`.
unsafe fn set_node_val(nodes: *mut c_double, dimension: usize, node_idx: usize, dim_idx: usize, val: f64) {
    *nodes.add(node_idx * dimension + dim_idx) = val;
}

// Helper for L2 norm of a 2D vector (components are contiguous)
#[inline]
/// # Safety
/// Caller must ensure `vec_ptr` points to a valid memory block of at least 2 `c_double` elements.
fn norm2_2d(vec_ptr: *const c_double) -> f64 {
    unsafe {
        let x = *vec_ptr.add(0);
        let y = *vec_ptr.add(1);
        (x * x + y * y).sqrt()
    }
}

/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `dimension * num_nodes` doubles if `num_nodes > 2`.
/// - `error` points to writable memory for 1 double.
/// - `dimension` and `num_nodes` are accurate.
/// - `get_node_val` and `set_node_val` preconditions are met.
/// - `norm2_2d` preconditions are met for `worst_case_vec.as_ptr()` if dimension is 2.
/// - `second_deriv_vec` and `worst_case_vec` allocations succeed.
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
    
    let norm_worst_case = if dimension == 2 {
        norm2_2d(worst_case_vec.as_ptr()) 
    } else {
        let mut sum_sq = 0.0;
        for val in worst_case_vec.iter().take(dimension) {
            sum_sq += *val * *val;
        }
        sum_sq.sqrt()
    };

    *error = 0.125 * (num_nodes - 1) as f64 * (num_nodes - 2) as f64 * norm_worst_case;
}

/// # Safety
/// Caller must ensure:
/// - `start0_ptr`, `end0_ptr`, `start1_ptr`, `end1_ptr` each point to readable memory for 2 doubles (x,y).
/// - `s` and `t` point to writable memory for 1 double each.
/// - `helpers::BEZ_cross_product` preconditions are met.
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

/// # Safety
/// Caller must ensure:
/// - `start0_ptr`, `end0_ptr`, `start1_ptr`, `end1_ptr` each point to readable memory for 2 doubles (x,y).
/// - `s_vals` and `t_vals` point to writable memory for at least 2 doubles each (for up to 2 intersection points).
/// - `num_written` points to writable memory for 1 c_int.
/// - `helpers::BEZ_cross_product` preconditions are met.
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

/// # Safety
/// Caller must ensure:
/// - `start0_ptr`, `end0_ptr`, `start1_ptr`, `end1_ptr` each point to readable memory for 2 doubles (x,y).
/// - `s_vals` and `t_vals` point to writable memory for at least 2 doubles each.
/// - `num_written` points to writable memory for 1 c_int.
/// - `segment_intersection_internal` and `parallel_lines_parameters_internal` preconditions are met.
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
/// # Safety
/// Caller must ensure:
/// - `nodes1` points to readable memory for `2 * num_nodes1` doubles.
/// - `nodes2` points to readable memory for `2 * num_nodes2` doubles.
/// - `new_s`, `new_t` point to writable memory for 1 double each.
/// - `status_code` points to writable memory for 1 c_int.
/// - `num_nodes1`, `num_nodes2` are non-negative and accurate.
/// - `curve::BEZ_evaluate_multi`, `curve::BEZ_evaluate_hodograph`, `helpers::solve2x2_rs` preconditions are met.
/// - Allocations for `func_val_vec`, `b1_s_vec`, `jac_mat_vec`, `b2_prime_t_vec` must succeed.
/// - Dimension is assumed to be 2.
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
/// # Safety
/// Caller must ensure:
/// - `nodes1` must be valid for reading `2 * num_nodes1` doubles if `num_nodes1 > 0`.
/// - `nodes2` must be valid for reading `2 * num_nodes2` doubles if `num_nodes2 > 0`.
/// - `enum_` must be a valid pointer to a `c_int`.
/// - `num_nodes1` and `num_nodes2` must be non-negative.
/// - `helpers::BEZ_bbox` preconditions must be met.
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
    // No-op for now. All global workspace management, if any, would be handled here.
}

/// # Safety
/// Caller must ensure:
/// - `nodes1` points to readable memory for `2 * num_nodes1` doubles.
/// - `nodes2` points to readable memory for `2 * num_nodes2` doubles.
/// - `polygon1_vec_target` and `polygon2_vec_target` are valid mutable references to `Vec<f64>`.
/// - `collision` points to writable memory for 1 u8.
/// - `num_nodes1`, `num_nodes2` are non-negative.
/// - `helpers::BEZ_simple_convex_hull` and `helpers::BEZ_polygon_collide` (or internal `line_line_collide_internal`) preconditions are met.
/// This function may reallocate `polygon1_vec_target` and `polygon2_vec_target`.
unsafe fn convex_hull_collide_internal(
    num_nodes1: c_int,
    nodes1: *const c_double,
    polygon1_vec_target: &mut Vec<f64>,
    num_nodes2: c_int,
    nodes2: *const c_double,
    polygon2_vec_target: &mut Vec<f64>,
    collision: &mut u8,
) {
    let mut polygon_size1: c_int = 0;
    let mut polygon_size2: c_int = 0;

    let required_cap1 = (num_nodes1 * 2) as usize;
    if polygon1_vec_target.capacity() < required_cap1 {
        polygon1_vec_target.reserve(required_cap1 - polygon1_vec_target.capacity());
    }
    
    let required_cap2 = (num_nodes2 * 2) as usize;
    if polygon2_vec_target.capacity() < required_cap2 {
        polygon2_vec_target.reserve(required_cap2 - polygon2_vec_target.capacity());
    }

    helpers::BEZ_simple_convex_hull(
        num_nodes1,
        nodes1,
        &mut polygon_size1,
        polygon1_vec_target.as_mut_ptr(),
    );
    polygon1_vec_target.set_len((polygon_size1 * 2) as usize);

    helpers::BEZ_simple_convex_hull(
        num_nodes2,
        nodes2,
        &mut polygon_size2,
        polygon2_vec_target.as_mut_ptr(),
    );
    polygon2_vec_target.set_len((polygon_size2 * 2) as usize);

    if polygon_size1 == 2 && polygon_size2 == 2 {
        let mut s_vals_dummy = [0.0; 2]; 
        let mut t_vals_dummy = [0.0; 2];
        let mut num_written_dummy = 0;
        
        line_line_collide_internal(
            polygon1_vec_target.as_ptr(),       
            polygon1_vec_target.as_ptr().add(2), 
            polygon2_vec_target.as_ptr(),       
            polygon2_vec_target.as_ptr().add(2), 
            s_vals_dummy.as_mut_ptr(),
            t_vals_dummy.as_mut_ptr(),
            &mut num_written_dummy,
        );
        *collision = if num_written_dummy > 0 { 1 } else { 0 };
    } else {
        helpers::BEZ_polygon_collide(
            polygon_size1,
            polygon1_vec_target.as_ptr(),
            polygon_size2,
            polygon2_vec_target.as_ptr(),
            collision,
        );
    }
}

/// # Safety
/// Caller must ensure:
/// - `nodes1` points to readable memory for `2 * num_nodes1` doubles.
/// - `first_deriv1` points to readable memory for `2 * (num_nodes1 - 1)` doubles if `num_nodes1 > 1`.
/// - `nodes2` points to readable memory for `2 * num_nodes2` doubles.
/// - `first_deriv2` points to readable memory for `2 * (num_nodes2 - 1)` doubles if `num_nodes2 > 1`.
/// - `jacobian` points to writable memory for `2 * 2 = 4` doubles.
/// - `func_val` points to writable memory for `2 * 1 = 2` doubles.
/// - `num_nodes1`, `num_nodes2` are non-negative and accurate.
/// - `curve::BEZ_evaluate_multi` preconditions are met.
/// - Dimension is assumed to be 2.
unsafe fn newton_simple_root_internal(
    s_param: c_double,
    num_nodes1: c_int,
    nodes1: *const c_double,
    first_deriv1: *const c_double,
    t_param: c_double,
    num_nodes2: c_int,
    nodes2: *const c_double,
    first_deriv2: *const c_double,
    jacobian: *mut c_double, 
    func_val: *mut c_double,   
) {
    let dim: usize = 2;
    let s_arr = [s_param];
    let t_arr = [t_param];
    let mut workspace_b2_t: Vec<f64> = vec![0.0; dim];

    curve::BEZ_evaluate_multi(num_nodes1, dim as c_int, nodes1, 1, s_arr.as_ptr(), func_val); 
    curve::BEZ_evaluate_multi(num_nodes2, dim as c_int, nodes2, 1, t_arr.as_ptr(), workspace_b2_t.as_mut_ptr());

    *func_val.add(0) -= workspace_b2_t[0];
    *func_val.add(1) -= workspace_b2_t[1];

    if (*func_val.add(0)).abs() < WIGGLE && (*func_val.add(1)).abs() < WIGGLE {
        *jacobian.add(0) = 0.0; *jacobian.add(1) = 0.0; // J[0,0], J[1,0]
        *jacobian.add(dim) = 0.0; *jacobian.add(dim + 1) = 0.0; // J[0,1], J[1,1]
        return;
    }

    if num_nodes1 > 1 { 
        curve::BEZ_evaluate_multi(num_nodes1 - 1, dim as c_int, first_deriv1, 1, s_arr.as_ptr(), jacobian);
    } else { 
        *jacobian.add(0) = 0.0; *jacobian.add(1) = 0.0;
    }

    let mut workspace_b2_dt: Vec<f64> = vec![0.0; dim];
    if num_nodes2 > 1 { 
        curve::BEZ_evaluate_multi(num_nodes2 - 1, dim as c_int, first_deriv2, 1, t_arr.as_ptr(), workspace_b2_dt.as_mut_ptr());
        *jacobian.add(dim) = -workspace_b2_dt[0];     
        *jacobian.add(dim + 1) = -workspace_b2_dt[1]; 
    } else {
        *jacobian.add(dim) = 0.0;
        *jacobian.add(dim + 1) = 0.0;
    }
}

/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `2 * num_nodes` doubles if `num_nodes > 0`.
/// - `line_start` and `line_end` each point to readable memory for 2 doubles.
/// - `enum_` points to writable memory for 1 c_int.
/// - `num_nodes` is non-negative.
/// - `helpers::BEZ_bbox` and `segment_intersection_internal` preconditions are met.
unsafe fn bbox_line_intersect_internal(
    num_nodes: c_int,
    nodes: *const c_double,
    line_start: *const c_double,
    line_end: *const c_double,
    enum_: *mut c_int,
) {
    let mut left: f64 = 0.0; let mut right: f64 = 0.0;
    let mut bottom: f64 = 0.0; let mut top: f64 = 0.0;
    helpers::BEZ_bbox(num_nodes, nodes, &mut left, &mut right, &mut bottom, &mut top);

    let line_start_x = *line_start.add(0);
    let line_start_y = *line_start.add(1);
    if helpers::BEZ_in_interval(line_start_x, left, right) != 0 &&
       helpers::BEZ_in_interval(line_start_y, bottom, top) != 0 {
        *enum_ = BOX_INTERSECTION_TYPE_INTERSECTION;
        return;
    }

    let line_end_x = *line_end.add(0);
    let line_end_y = *line_end.add(1);
    if helpers::BEZ_in_interval(line_end_x, left, right) != 0 &&
       helpers::BEZ_in_interval(line_end_y, bottom, top) != 0 {
        *enum_ = BOX_INTERSECTION_TYPE_INTERSECTION;
        return;
    }
    
    let mut s_curr: f64 = 0.0;
    let mut t_curr: f64 = 0.0;
    let mut success_bool: bool;

    let bbox_edges_points = [
        ([left, bottom], [right, bottom]), // Bottom
        ([right, bottom], [right, top]),   // Right
        ([right, top], [left, top]),       // Top
    ];

    for i in 0..3 { 
        let edge_pts = bbox_edges_points[i];
        success_bool = segment_intersection_internal(
            edge_pts.0.as_ptr(), edge_pts.1.as_ptr(),
            line_start, line_end, &mut s_curr, &mut t_curr
        );
        if success_bool && 
           helpers::BEZ_in_interval(s_curr, 0.0, 1.0) != 0 && 
           helpers::BEZ_in_interval(t_curr, 0.0, 1.0) != 0 {
            *enum_ = BOX_INTERSECTION_TYPE_INTERSECTION;
            return;
        }
    }
    *enum_ = BOX_INTERSECTION_TYPE_DISJOINT;
}

/// # Safety
/// Caller must ensure:
/// - `nodes1` points to readable memory for `2 * num_nodes1` doubles.
/// - `nodes2` points to readable memory for `2 * num_nodes2` doubles.
/// - `both_linear` and `coincident` point to writable memory for 1 u8 each.
/// - `intersections` is a mutable reference to a Vec that can be resized.
/// - `num_intersections` points to writable memory for 1 c_int.
/// - `num_nodes1`, `num_nodes2` are non-negative.
/// - `linearization_error_internal`, `segment_intersection_internal`, `parallel_lines_parameters_internal` preconditions met.
unsafe fn check_lines_internal(
    num_nodes1: c_int,
    nodes1: *const c_double,
    num_nodes2: c_int,
    nodes2: *const c_double,
    both_linear: &mut u8,
    coincident: &mut u8,
    intersections: &mut Vec<[c_double; 2]>, 
    num_intersections: &mut c_int,
) {
    *num_intersections = 0;
    *coincident = 0; 
    let mut error1: f64 = 0.0;
    let mut error2: f64 = 0.0;

    linearization_error_internal(num_nodes1 as usize, 2, nodes1, &mut error1);
    if error1.abs() > LINEARIZATION_THRESHOLD { 
        *both_linear = 0; 
        return;
    }
    linearization_error_internal(num_nodes2 as usize, 2, nodes2, &mut error2);
    if error2.abs() > LINEARIZATION_THRESHOLD { 
        *both_linear = 0; 
        return;
    }
    
    *both_linear = 1; 
    let mut s: f64 = 0.0;
    let mut t: f64 = 0.0;

    if num_nodes1 < 1 || num_nodes2 < 1 {
        *both_linear = 0; 
        return;
    }

    let start0_ptr = nodes1; 
    let end0_ptr = if num_nodes1 > 1 { nodes1.add(((num_nodes1 - 1) * 2) as usize) } else { start0_ptr };
    let start1_ptr = nodes2;
    let end1_ptr = if num_nodes2 > 1 { nodes2.add(((num_nodes2 - 1) * 2) as usize) } else { start1_ptr };

    let mut s_vals_arr = [0.0; 2];
    let mut t_vals_arr = [0.0; 2];
    intersections.clear(); 

    if segment_intersection_internal(start0_ptr, end0_ptr, start1_ptr, end1_ptr, &mut s, &mut t) {
        if helpers::BEZ_in_interval(s, 0.0, 1.0) != 0 && helpers::BEZ_in_interval(t, 0.0, 1.0) != 0 {
            intersections.push([s.clamp(0.0, 1.0), t.clamp(0.0, 1.0)]);
            *num_intersections = 1;
        } else {
            *num_intersections = 0; 
        }
    } else {
        let mut parallel_num_written = 0;
        parallel_lines_parameters_internal(
            start0_ptr, end0_ptr, start1_ptr, end1_ptr,
            s_vals_arr.as_mut_ptr(), t_vals_arr.as_mut_ptr(), &mut parallel_num_written
        );
        if parallel_num_written > 0 {
            *coincident = 1; 
            for i in 0..parallel_num_written as usize {
                intersections.push([s_vals_arr[i], t_vals_arr[i]]);
            }
        }
        *num_intersections = parallel_num_written;
    }
}

fn add_intersection_internal(
    s_param: c_double,
    t_param: c_double,
    num_intersections: &mut c_int,
    intersections: &mut Vec<[c_double; 2]>,
) {
    let s_norm = if s_param < ZERO_THRESHOLD { 1.0 - s_param } else { s_param };
    let t_norm = if t_param < ZERO_THRESHOLD { 1.0 - t_param } else { t_param };
    let norm_candidate_sq = s_norm * s_norm + t_norm * t_norm;

    for i in 0..(*num_intersections as usize) {
        let existing_s = intersections[i][0];
        let existing_t = intersections[i][1];
        let diff_s = s_param - existing_s;
        let diff_t = t_param - existing_t;
        let dist_sq = diff_s * diff_s + diff_t * diff_t;

        let threshold_sq = if norm_candidate_sq < WIGGLE * WIGGLE { 
            NEWTON_ERROR_RATIO * NEWTON_ERROR_RATIO 
        } else {
            norm_candidate_sq * NEWTON_ERROR_RATIO * NEWTON_ERROR_RATIO 
        };
        
        if dist_sq < threshold_sq { return; } 
        if dist_sq < WIGGLE * WIGGLE * WIGGLE * WIGGLE { return; } 
    }

    intersections.push([s_param, t_param]);
    *num_intersections += 1;
}

/// # Safety
/// Caller must ensure:
/// - `nodes1`, `first_deriv1`, `second_deriv1` (if `num_nodes1 > 2`) point to readable memory for curve 1 and its derivatives.
/// - `nodes2`, `first_deriv2`, `second_deriv2` (if `num_nodes2 > 2`) point to readable memory for curve 2 and its derivatives.
/// - `modified_lhs` (for 2x2 matrix) and `modified_rhs` (for 2x1 vector) point to writable memory.
/// - All `num_nodes` arguments are accurate. `num_nodesX >= 1`. `first_derivX` valid if `num_nodesX >= 2`. `second_derivX` valid if `num_nodesX >= 3`.
/// - `curve::BEZ_evaluate_multi` and `helpers::BEZ_cross_product` preconditions are met.
/// - Dimension is assumed to be 2.
unsafe fn newton_double_root_internal(
    s: c_double,
    num_nodes1: c_int,
    nodes1: *const c_double,
    first_deriv1: *const c_double, 
    second_deriv1: *const c_double, 
    t: c_double,
    num_nodes2: c_int,
    nodes2: *const c_double,
    first_deriv2: *const c_double, 
    second_deriv2: *const c_double, 
    modified_lhs: *mut c_double, 
    modified_rhs: *mut c_double, 
) {
    let dim: usize = 2;
    let s_arr = [s];
    let t_arr = [t];

    let mut jacobian_g_data: Vec<f64> = vec![0.0; 3 * 2]; 
    let jacobian_g_ptr = jacobian_g_data.as_mut_ptr(); 
    let mut func_val_g: Vec<f64> = vec![0.0; 3];    

    let mut b1_s_vec = vec![0.0; dim]; 
    let mut b2_t_vec = vec![0.0; dim]; 
    let mut b1_ds_vec = vec![0.0; dim]; 
    let mut b2_dt_vec = vec![0.0; dim]; 

    curve::BEZ_evaluate_multi(num_nodes1, dim as c_int, nodes1, 1, s_arr.as_ptr(), b1_s_vec.as_mut_ptr());
    curve::BEZ_evaluate_multi(num_nodes2, dim as c_int, nodes2, 1, t_arr.as_ptr(), b2_t_vec.as_mut_ptr());
    func_val_g[0] = b1_s_vec[0] - b2_t_vec[0];
    func_val_g[1] = b1_s_vec[1] - b2_t_vec[1];

    if num_nodes1 > 1 {
        curve::BEZ_evaluate_multi(num_nodes1 - 1, dim as c_int, first_deriv1, 1, s_arr.as_ptr(), b1_ds_vec.as_mut_ptr());
    } else { b1_ds_vec.fill(0.0); } 
    if num_nodes2 > 1 {
        curve::BEZ_evaluate_multi(num_nodes2 - 1, dim as c_int, first_deriv2, 1, t_arr.as_ptr(), b2_dt_vec.as_mut_ptr());
    } else { b2_dt_vec.fill(0.0); }
    helpers::BEZ_cross_product(b1_ds_vec.as_ptr(), b2_dt_vec.as_ptr(), &mut func_val_g[2]);

    if func_val_g.iter().all(|&val| val.abs() < WIGGLE) { 
        for i in 0..2 { modified_rhs.add(i).write(0.0); }
        for i in 0..4 { modified_lhs.add(i).write(0.0); } 
        return;
    }

    jacobian_g_ptr.add(0).write(b1_ds_vec[0]);      
    jacobian_g_ptr.add(1).write(b1_ds_vec[1]);      
    let mut b1_dds_vec = vec![0.0; dim]; 
    if num_nodes1 > 2 {
        curve::BEZ_evaluate_multi(num_nodes1 - 2, dim as c_int, second_deriv1, 1, s_arr.as_ptr(), b1_dds_vec.as_mut_ptr());
        helpers::BEZ_cross_product(b1_dds_vec.as_ptr(), b2_dt_vec.as_ptr(), jacobian_g_ptr.add(2)); 
    } else {
        jacobian_g_ptr.add(2).write(0.0); 
    }

    jacobian_g_ptr.add(3).write(-b2_dt_vec[0]);    
    jacobian_g_ptr.add(4).write(-b2_dt_vec[1]);    
    let mut b2_ddt_vec = vec![0.0; dim]; 
    if num_nodes2 > 2 {
        curve::BEZ_evaluate_multi(num_nodes2 - 2, dim as c_int, second_deriv2, 1, t_arr.as_ptr(), b2_ddt_vec.as_mut_ptr());
        helpers::BEZ_cross_product(b1_ds_vec.as_ptr(), b2_ddt_vec.as_ptr(), jacobian_g_ptr.add(5)); 
    } else {
        jacobian_g_ptr.add(5).write(0.0); 
    }
    
    modified_lhs.add(0).write(jacobian_g_ptr.add(0).read()*jacobian_g_ptr.add(0).read() + jacobian_g_ptr.add(1).read()*jacobian_g_ptr.add(1).read() + jacobian_g_ptr.add(2).read()*jacobian_g_ptr.add(2).read()); 
    modified_lhs.add(1).write(jacobian_g_ptr.add(0).read()*jacobian_g_ptr.add(3).read() + jacobian_g_ptr.add(1).read()*jacobian_g_ptr.add(4).read() + jacobian_g_ptr.add(2).read()*jacobian_g_ptr.add(5).read()); 
    modified_lhs.add(2).write(jacobian_g_ptr.add(3).read()*jacobian_g_ptr.add(0).read() + jacobian_g_ptr.add(4).read()*jacobian_g_ptr.add(1).read() + jacobian_g_ptr.add(5).read()*jacobian_g_ptr.add(2).read()); 
    modified_lhs.add(3).write(jacobian_g_ptr.add(3).read()*jacobian_g_ptr.add(3).read() + jacobian_g_ptr.add(4).read()*jacobian_g_ptr.add(4).read() + jacobian_g_ptr.add(5).read()*jacobian_g_ptr.add(5).read()); 

    modified_rhs.add(0).write(jacobian_g_ptr.add(0).read()*func_val_g[0] + jacobian_g_ptr.add(1).read()*func_val_g[1] + jacobian_g_ptr.add(2).read()*func_val_g[2]);
    modified_rhs.add(1).write(jacobian_g_ptr.add(3).read()*func_val_g[0] + jacobian_g_ptr.add(4).read()*func_val_g[1] + jacobian_g_ptr.add(5).read()*func_val_g[2]);
}

/// # Safety
/// The safety of this function depends heavily on the `evaluate_fn` passed to it.
/// `evaluate_fn` is an unsafe function that takes `s`, `t`, and output pointers
/// for `jacobian` (2x2) and `func_val` (2x1).
/// - `s_in`, `t_in` must be valid parameters for `evaluate_fn`.
/// - `new_s`, `new_t` must point to writable memory for 1 double each.
/// - `converged` must point to writable memory for 1 u8.
/// - `evaluate_fn` must correctly populate `jacobian` and `func_val` and handle its own inputs safely.
/// - `helpers::solve2x2_rs` preconditions must be met with `jacobian` and `func_val` from `evaluate_fn`.
unsafe fn newton_iterate_internal<F>(
    mut evaluate_fn: F,
    s_in: c_double,
    t_in: c_double,
    new_s: &mut c_double,
    new_t: &mut c_double,
    converged: &mut u8,
) 
where
    F: FnMut(c_double, c_double, *mut c_double, *mut c_double),
{
    let mut current_s = s_in;
    let mut current_t = t_in;
    let mut linear_updates = 0;
    let mut norm_update_prev = 0.0;

    let mut jacobian_vec = vec![0.0; 4]; 
    let mut func_val_vec = vec![0.0; 2];   

    *converged = 0; 

    for i in 1..=10 { 
        evaluate_fn(current_s, current_t, jacobian_vec.as_mut_ptr(), func_val_vec.as_mut_ptr());

        if func_val_vec.iter().all(|&val| val.abs() < WIGGLE) { 
            *converged = 1; 
            *new_s = current_s;
            *new_t = current_t;
            return;
        }

        let mut delta_s: f64 = 0.0;
        let mut delta_t: f64 = 0.0;
        let mut singular_bool: bool = false;
        
        helpers::solve2x2_rs(
            jacobian_vec.as_ptr(),
            func_val_vec.as_ptr(),
            &mut singular_bool,
            &mut delta_s,
            &mut delta_t,
        );

        if singular_bool {
            *new_s = current_s; 
            *new_t = current_t;
            return; 
        }

        let norm_update = (delta_s * delta_s + delta_t * delta_t).sqrt(); 
        
        if i > 1 && norm_update > 0.25 * norm_update_prev {
            linear_updates += 1;
        }
        norm_update_prev = norm_update;

        if i >= 5 && (3 * linear_updates >= 2 * i) {
            *new_s = current_s;
            *new_t = current_t;
            return; 
        }
        
        let norm_soln = (current_s * current_s + current_t * current_t).sqrt();
        
        current_s -= delta_s;
        current_t -= delta_t;

        if norm_update < NEWTON_ERROR_RATIO * norm_soln { // Original Fortran check
             *converged = 1; 
             *new_s = current_s;
             *new_t = current_t;
             return;
        }
        if norm_update < WIGGLE { // Absolute check for very small steps (e.g. if norm_soln is zero)
            *converged = 1;
            *new_s = current_s;
            *new_t = current_t;
            return;
        }
    }
    *new_s = current_s;
    *new_t = current_t;
}

/// # Safety
/// Caller must ensure:
/// - `nodes1`, `nodes2` point to readable memory for their respective curves.
/// - `new_s`, `new_t`, `status` point to writable memory.
/// - All num_nodes arguments are accurate.
/// - Internal derivative calculations and calls to `newton_iterate_internal` are safe.
unsafe fn full_newton_nonzero_internal(
    s_in: c_double,
    num_nodes1: c_int,
    nodes1: *const c_double,
    t_in: c_double,
    num_nodes2: c_int,
    nodes2: *const c_double,
    new_s: &mut c_double,
    new_t: &mut c_double,
    status: &mut c_int,
) {
    *status = status::STATUS_SUCCESS;
    let dim:usize = 2;

    let mut first_deriv1_vec: Vec<f64> = if num_nodes1 > 1 { vec![0.0; dim * (num_nodes1 - 1)as usize] } else { Vec::new() };
    let mut first_deriv2_vec: Vec<f64> = if num_nodes2 > 1 { vec![0.0; dim * (num_nodes2 - 1)as usize] } else { Vec::new() };
    let mut second_deriv1_vec: Vec<f64> = if num_nodes1 > 2 { vec![0.0; dim * (num_nodes1 - 2)as usize] } else { Vec::new() };
    let mut second_deriv2_vec: Vec<f64> = if num_nodes2 > 2 { vec![0.0; dim * (num_nodes2 - 2)as usize] } else { Vec::new() };

    if num_nodes1 > 1 {
        let degree1 = (num_nodes1 - 1) as f64;
        for i in 0..(num_nodes1 - 1) as usize {
            for d_idx in 0..dim {
                let val = degree1 * (get_node_val(nodes1, dim, i + 1, d_idx) - get_node_val(nodes1, dim, i, d_idx));
                set_node_val(first_deriv1_vec.as_mut_ptr(), dim, i, d_idx, val);
            }
        }
    }
    if num_nodes2 > 1 {
        let degree2 = (num_nodes2 - 1) as f64;
         for i in 0..(num_nodes2 - 1) as usize {
            for d_idx in 0..dim {
                let val = degree2 * (get_node_val(nodes2, dim, i + 1, d_idx) - get_node_val(nodes2, dim, i, d_idx));
                set_node_val(first_deriv2_vec.as_mut_ptr(), dim, i, d_idx, val);
            }
        }
    }
    if num_nodes1 > 2 {
        let degree1 = (num_nodes1 - 1) as f64;
        let degree_of_deriv1 = (num_nodes1 - 2) as f64; 
        for i in 0..(num_nodes1 - 2) as usize {
            for d_idx in 0..dim {
                 let val = degree_of_deriv1 * 
                           (get_node_val(first_deriv1_vec.as_ptr(), dim, i + 1, d_idx)/degree1 - 
                            get_node_val(first_deriv1_vec.as_ptr(), dim, i, d_idx)/degree1 );
                 set_node_val(second_deriv1_vec.as_mut_ptr(), dim, i, d_idx, val);
            }
        }
    }
     if num_nodes2 > 2 {
        let degree2 = (num_nodes2 - 1) as f64;
        let degree_of_deriv2 = (num_nodes2 - 2) as f64;
        for i in 0..(num_nodes2 - 2) as usize {
            for d_idx in 0..dim {
                 let val = degree_of_deriv2 * 
                           (get_node_val(first_deriv2_vec.as_ptr(), dim, i + 1, d_idx)/degree2 - 
                            get_node_val(first_deriv2_vec.as_ptr(), dim, i, d_idx)/degree2 );
                 set_node_val(second_deriv2_vec.as_mut_ptr(), dim, i, d_idx, val);
            }
        }
    }
    
    let mut current_s = s_in;
    let mut current_t = t_in;
    let mut converged_u8: u8 = 0;

    let nodes1_ptr = nodes1; 
    let nodes2_ptr = nodes2;
    let fd1_ptr = first_deriv1_vec.as_ptr();
    let fd2_ptr = first_deriv2_vec.as_ptr();
    let sd1_ptr = second_deriv1_vec.as_ptr();
    let sd2_ptr = second_deriv2_vec.as_ptr();

    let evaluate_simple_closure = |s_val, t_val, jacobian_ptr, func_val_ptr| {
        newton_simple_root_internal(
            s_val, num_nodes1, nodes1_ptr, fd1_ptr,
            t_val, num_nodes2, nodes2_ptr, fd2_ptr,
            jacobian_ptr, func_val_ptr
        );
    };
    newton_iterate_internal(evaluate_simple_closure, s_in, t_in, &mut current_s, &mut current_t, &mut converged_u8);

    if converged_u8 != 0 {
        *new_s = current_s;
        *new_t = current_t;
        return;
    }

    let evaluate_double_closure = |s_val, t_val, jacobian_ptr, func_val_ptr| {
        newton_double_root_internal(
            s_val, num_nodes1, nodes1_ptr, fd1_ptr, sd1_ptr,
            t_val, num_nodes2, nodes2_ptr, fd2_ptr, sd2_ptr,
            jacobian_ptr, func_val_ptr
        );
    };
    newton_iterate_internal(evaluate_double_closure, current_s, current_t, new_s, new_t, &mut converged_u8);
    
    if converged_u8 == 0 {
        *status = status::STATUS_BAD_MULTIPLICITY;
    }
}

/// # Safety
/// Caller must ensure:
/// - `nodes1`, `nodes2` point to readable memory for their respective curves.
/// - `new_s`, `new_t`, `status` point to writable memory.
/// - All num_nodes arguments are accurate.
/// - Internal calls to `full_newton_nonzero_internal` are safe.
unsafe fn full_newton_internal(
    s_in: c_double,
    num_nodes1: c_int,
    nodes1: *const c_double,
    t_in: c_double,
    num_nodes2: c_int,
    nodes2: *const c_double,
    new_s: &mut c_double,
    new_t: &mut c_double,
    status: &mut c_int,
) {
    let mut s_to_use = s_in;
    let mut t_to_use = t_in;
    let mut reverse_s = false;
    let mut reverse_t = false;
    let dim: usize = 2;

    let mut nodes1_storage: Vec<f64> = Vec::new();
    let mut nodes2_storage: Vec<f64> = Vec::new();
    
    let mut current_nodes1_ptr = nodes1;
    let mut current_nodes2_ptr = nodes2;

    if s_in < ZERO_THRESHOLD {
        s_to_use = 1.0 - s_in;
        reverse_s = true;
        if num_nodes1 > 0 {
            nodes1_storage = vec![0.0; (num_nodes1 * dim as c_int) as usize];
            for i in 0..num_nodes1 as usize {
                for d_idx in 0..dim {
                    set_node_val(nodes1_storage.as_mut_ptr(), dim, i, d_idx, get_node_val(nodes1, dim, num_nodes1 as usize - 1 - i, d_idx));
                }
            }
            current_nodes1_ptr = nodes1_storage.as_ptr();
        }
    }

    if t_in < ZERO_THRESHOLD {
        t_to_use = 1.0 - t_in;
        reverse_t = true;
         if num_nodes2 > 0 {
            nodes2_storage = vec![0.0; (num_nodes2 * dim as c_int) as usize];
            for i in 0..num_nodes2 as usize {
                for d_idx in 0..dim {
                    set_node_val(nodes2_storage.as_mut_ptr(), dim, i, d_idx, get_node_val(nodes2, dim, num_nodes2 as usize - 1 - i, d_idx));
                }
            }
            current_nodes2_ptr = nodes2_storage.as_ptr();
        }
    }
    
    full_newton_nonzero_internal(
        s_to_use, num_nodes1, current_nodes1_ptr,
        t_to_use, num_nodes2, current_nodes2_ptr,
        new_s, new_t, status
    );

    if reverse_s {
        *new_s = 1.0 - *new_s;
    }
    if reverse_t {
        *new_t = 1.0 - *new_t;
    }
}

/// # Safety
/// Caller must ensure:
/// - `curve1.nodes` and `curve2.nodes` (if curves are valid) point to readable memory.
/// - `root_nodes1` and `root_nodes2` point to readable memory for the original curves.
///   The number of nodes for these root_nodes must be `original_num_nodes1` and `original_num_nodes2`.
/// - `refined_s`, `refined_t`, `does_intersect`, `status` point to writable memory.
/// - All num_nodes and dimension fields in `CurveData` are accurate for `curve1` and `curve2`.
/// - `polygon1_scratch` and `polygon2_scratch` are mutable references to `Vec<f64>` for scratch space.
unsafe fn from_linearized_internal(
    curve1: &CurveData, 
    root_nodes1: *const c_double, 
    original_num_nodes1: c_int,
    curve2: &CurveData, 
    root_nodes2: *const c_double, 
    original_num_nodes2: c_int,
    refined_s: &mut c_double,
    refined_t: &mut c_double,
    does_intersect: &mut u8, 
    status: &mut c_int,
    polygon1_scratch: &mut Vec<f64>, 
    polygon2_scratch: &mut Vec<f64>, 
) {
    *status = status::STATUS_SUCCESS;
    *does_intersect = 0; 

    let mut s_approx: f64 = 0.0;
    let mut t_approx: f64 = 0.0;
    let mut segment_intersection_success: bool = false;

    if curve1.num_nodes > 0 && curve2.num_nodes > 0 && curve1.dimension > 0 && curve2.dimension > 0 {
        let end0_ptr = if curve1.num_nodes > 1 {
            curve1.nodes.as_ptr().add((curve1.num_nodes - 1) * curve1.dimension)
        } else {
            curve1.nodes.as_ptr() 
        };
        let end1_ptr = if curve2.num_nodes > 1 {
            curve2.nodes.as_ptr().add((curve2.num_nodes - 1) * curve2.dimension)
        } else {
            curve2.nodes.as_ptr() 
        };
        
        segment_intersection_success = segment_intersection_internal(
            curve1.nodes.as_ptr(), end0_ptr, 
            curve2.nodes.as_ptr(), end1_ptr, 
            &mut s_approx, &mut t_approx
        );
    }
    
    let mut bad_parameters = false;
    if segment_intersection_success {
        if !(helpers::BEZ_in_interval(s_approx, 0.0, 1.0) != 0 && 
             helpers::BEZ_in_interval(t_approx, 0.0, 1.0) != 0) {
            bad_parameters = true;
        }
    } else {
        bad_parameters = true;
        s_approx = 0.5; 
        t_approx = 0.5;
    }

    if bad_parameters {
        if curve1.num_nodes > 0 && curve2.num_nodes > 0 { 
            let mut hull_collision_u8: u8 = 0;
            convex_hull_collide_internal(
                curve1.num_nodes as c_int, curve1.nodes.as_ptr(), polygon1_scratch,
                curve2.num_nodes as c_int, curve2.nodes.as_ptr(), polygon2_scratch,
                &mut hull_collision_u8
            );
            if hull_collision_u8 == 0 { 
                return;
            }
        } else { 
            return;
        }
    }

    let s_orig = (1.0 - s_approx) * curve1.start + s_approx * curve1.end;
    let t_orig = (1.0 - t_approx) * curve2.start + t_approx * curve2.end;
    
    let mut s_newton: f64 = 0.0;
    let mut t_newton: f64 = 0.0;

    if original_num_nodes1 > 0 && original_num_nodes2 > 0 {
        full_newton_internal(
            s_orig, original_num_nodes1, root_nodes1,
            t_orig, original_num_nodes2, root_nodes2,
            &mut s_newton, &mut t_newton, status
        );
    } else { 
        *status = status::STATUS_UNKNOWN; 
        return;
    }

    if *status != status::STATUS_SUCCESS {
        return;
    }

    let mut s_wiggle_success: u8 = 0;
    let mut t_wiggle_success: u8 = 0;
    let mut s_final: f64 = 0.0;
    let mut t_final: f64 = 0.0;

    helpers::BEZ_wiggle_interval(s_newton, &mut s_final, &mut s_wiggle_success);
    if s_wiggle_success == 0 { return; }
    helpers::BEZ_wiggle_interval(t_newton, &mut t_final, &mut t_wiggle_success);
    if t_wiggle_success == 0 { return; }

    *does_intersect = 1; 
    *refined_s = s_final;
    *refined_t = t_final;
}

/// # Safety
/// Caller must ensure `curve1.nodes`, `curve2.nodes` are valid.
/// `node_first_ptr`, `node_second_ptr` must point to 2 valid doubles.
/// `intersections` may be reallocated.
unsafe fn endpoint_check_internal(
    curve1: &CurveData,
    node_first_ptr: *const c_double, 
    s_param: c_double,               
    curve2: &CurveData,
    node_second_ptr: *const c_double, 
    t_param: c_double,                
    num_intersections: &mut c_int,
    intersections: &mut Vec<[c_double; 2]>,
) {
    if curve1.dimension == 0 || curve2.dimension == 0 { return; } 
    if helpers::BEZ_vector_close(curve1.dimension as c_int, node_first_ptr, node_second_ptr, VECTOR_CLOSE_EPS) == 0 {
        return; 
    }
    add_intersection_internal(s_param, t_param, num_intersections, intersections);
}

/// # Safety
/// Caller must ensure `first.nodes`, `second.nodes` are valid.
/// `intersections` may be reallocated.
unsafe fn tangent_bbox_intersection_internal(
    first: &CurveData,
    second: &CurveData,
    num_intersections: &mut c_int,
    intersections: &mut Vec<[c_double; 2]>,
) {
    if first.num_nodes == 0 || second.num_nodes == 0 || first.dimension != 2 || second.dimension != 2 {
        return;
    }

    let nodes1_start_ptr = first.nodes.as_ptr();
    let nodes1_end_ptr = first.nodes.as_ptr().add((first.num_nodes - 1) * first.dimension);
    
    let nodes2_start_ptr = second.nodes.as_ptr();
    let nodes2_end_ptr = second.nodes.as_ptr().add((second.num_nodes - 1) * second.dimension);

    endpoint_check_internal(first, nodes1_start_ptr, 0.0, second, nodes2_start_ptr, 0.0, num_intersections, intersections);
    endpoint_check_internal(first, nodes1_start_ptr, 0.0, second, nodes2_end_ptr, 1.0, num_intersections, intersections);
    endpoint_check_internal(first, nodes1_end_ptr, 1.0, second, nodes2_start_ptr, 0.0, num_intersections, intersections);
    endpoint_check_internal(first, nodes1_end_ptr, 1.0, second, nodes2_end_ptr, 1.0, num_intersections, intersections);
}

/// # Safety
/// Caller must ensure `first` and `second` CurveData objects are valid.
/// `candidates_vec` may be reallocated.
unsafe fn add_candidates_internal(
    candidates_vec: &mut Vec<[CurveData; 2]>, 
    first: &CurveData, 
    second: &CurveData, 
    subdivision_type: c_int
) {    
    if subdivision_type == SUBDIVIDE_FIRST {
        if first.num_nodes == 0 { return; } 
        let (left1, right1) = first.subdivide_curve();
        candidates_vec.push([left1, second.clone()]);
        candidates_vec.push([right1, second.clone()]);
    } else if subdivision_type == SUBDIVIDE_SECOND {
        if second.num_nodes == 0 { return; }
        let (left2, right2) = second.subdivide_curve();
        candidates_vec.push([first.clone(), left2]);
        candidates_vec.push([first.clone(), right2]);
    } else if subdivision_type == SUBDIVIDE_BOTH {
        if first.num_nodes == 0 || second.num_nodes == 0 { return; }
        let (left1, right1) = first.subdivide_curve();
        let (left2, right2) = second.subdivide_curve();
        
        candidates_vec.push([left1.clone(), left2.clone()]);
        candidates_vec.push([left1, right2.clone()]); 
        candidates_vec.push([right1.clone(), left2]);
        candidates_vec.push([right1, right2]);
    }
}

/// # Safety
/// All pointer arguments must be valid. `intersections` and `next_candidates_list` may be reallocated.
/// `root_nodes_first` and `root_nodes_second` must point to the original, non-subdivided curve nodes.
/// `original_num_nodes1` and `original_num_nodes2` must match the node counts for `root_nodes_first` and `root_nodes_second`.
unsafe fn intersect_one_round_internal(
    root_nodes_first: *const c_double, 
    original_num_nodes1: c_int,
    root_nodes_second: *const c_double, 
    original_num_nodes2: c_int,
    current_candidates_list: &[[CurveData; 2]], 
    num_intersections: &mut c_int,
    intersections: &mut Vec<[c_double; 2]>, 
    next_candidates_list: &mut Vec<[CurveData; 2]>, 
    status: &mut c_int,
    polygon1_scratch: &mut Vec<f64>, 
    polygon2_scratch: &mut Vec<f64>,
) {
    next_candidates_list.clear();
    *status = status::STATUS_SUCCESS;

    for candidate_pair in current_candidates_list.iter() {
        let first_curve = &candidate_pair[0];
        let second_curve = &candidate_pair[1];

        if first_curve.dimension == 0 || second_curve.dimension == 0 { continue; }

        let mut err1: f64 = 0.0;
        let mut err2: f64 = 0.0;
        linearization_error_internal(first_curve.num_nodes, first_curve.dimension, first_curve.nodes.as_ptr(), &mut err1);
        linearization_error_internal(second_curve.num_nodes, second_curve.dimension, second_curve.nodes.as_ptr(), &mut err2);

        let mut subdivision_type: c_int;
        let mut bbox_intersection_type: c_int = 0; 

        if err1 < LINEARIZATION_THRESHOLD {
            if err2 < LINEARIZATION_THRESHOLD { 
                subdivision_type = SUBDIVIDE_NEITHER;
                if first_curve.num_nodes > 0 && second_curve.num_nodes > 0 {
                    BEZ_bbox_intersect(
                        first_curve.num_nodes as c_int, first_curve.nodes.as_ptr(),
                        second_curve.num_nodes as c_int, second_curve.nodes.as_ptr(),
                        &mut bbox_intersection_type,
                    );
                } else { bbox_intersection_type = BOX_INTERSECTION_TYPE_DISJOINT; }
            } else { 
                subdivision_type = SUBDIVIDE_SECOND;
                if first_curve.num_nodes > 0 && second_curve.num_nodes > 0 { 
                    let first_end_ptr = if first_curve.num_nodes > 1 {
                        first_curve.nodes.as_ptr().add((first_curve.num_nodes - 1) * first_curve.dimension)
                    } else { first_curve.nodes.as_ptr() };
                     bbox_line_intersect_internal(
                        second_curve.num_nodes as c_int, second_curve.nodes.as_ptr(),
                        first_curve.nodes.as_ptr(), first_end_ptr, 
                        &mut bbox_intersection_type,
                    );
                } else { bbox_intersection_type = BOX_INTERSECTION_TYPE_DISJOINT; }
            }
        } else {
            if err2 < LINEARIZATION_THRESHOLD { 
                subdivision_type = SUBDIVIDE_FIRST;
                 if first_curve.num_nodes > 0 && second_curve.num_nodes > 0 { 
                    let second_end_ptr = if second_curve.num_nodes > 1 {
                        second_curve.nodes.as_ptr().add((second_curve.num_nodes - 1) * second_curve.dimension)
                    } else { second_curve.nodes.as_ptr() };
                     bbox_line_intersect_internal(
                        first_curve.num_nodes as c_int, first_curve.nodes.as_ptr(),
                        second_curve.nodes.as_ptr(), second_end_ptr, 
                        &mut bbox_intersection_type,
                    );
                } else { bbox_intersection_type = BOX_INTERSECTION_TYPE_DISJOINT; }
            } else { 
                subdivision_type = SUBDIVIDE_BOTH;
                if first_curve.num_nodes > 0 && second_curve.num_nodes > 0 {
                    BEZ_bbox_intersect(
                        first_curve.num_nodes as c_int, first_curve.nodes.as_ptr(),
                        second_curve.num_nodes as c_int, second_curve.nodes.as_ptr(),
                        &mut bbox_intersection_type,
                    );
                } else { bbox_intersection_type = BOX_INTERSECTION_TYPE_DISJOINT; }
            }
        }

        if bbox_intersection_type == BOX_INTERSECTION_TYPE_DISJOINT {
            continue;
        } else if bbox_intersection_type == BOX_INTERSECTION_TYPE_TANGENT && subdivision_type != SUBDIVIDE_NEITHER {
            tangent_bbox_intersection_internal(first_curve, second_curve, num_intersections, intersections);
            continue;
        }
        
        if subdivision_type == SUBDIVIDE_NEITHER {
            let mut refined_s: f64 = 0.0;
            let mut refined_t: f64 = 0.0;
            let mut does_intersect_u8: u8 = 0;
            
            from_linearized_internal(
                first_curve, root_nodes_first, original_num_nodes1,
                second_curve, root_nodes_second, original_num_nodes2,
                &mut refined_s, &mut refined_t, &mut does_intersect_u8, status,
                polygon1_scratch, polygon2_scratch
            );
            if *status != status::STATUS_SUCCESS { return; }
            if does_intersect_u8 != 0 {
                add_intersection_internal(refined_s, refined_t, num_intersections, intersections);
            }
            continue;
        }
        
        add_candidates_internal(next_candidates_list, first_curve, second_curve, subdivision_type);
    }
}

/// # Safety
/// `nodes` points to valid memory for `2 * curr_size` doubles.
/// `workspace_nodes_vec` and `elevated_nodes_vec` will be modified (reallocated and length set).
unsafe fn elevate_helper_internal(
    curr_num_nodes: c_int,       
    nodes_in: *const c_double, 
    final_num_nodes: c_int,      
    workspace_nodes_vec: &mut Vec<f64>, 
    elevated_nodes_vec: &mut Vec<f64>,  
) {
    let dim = 2; 
    let required_cap = (final_num_nodes * dim) as usize;

    if workspace_nodes_vec.capacity() < required_cap { workspace_nodes_vec.reserve(required_cap - workspace_nodes_vec.capacity()); }
    if elevated_nodes_vec.capacity() < required_cap { elevated_nodes_vec.reserve(required_cap - elevated_nodes_vec.capacity()); }
    
    let num_steps = final_num_nodes - curr_num_nodes;
    if num_steps < 0 { return; } 
    if num_steps == 0 { 
        elevated_nodes_vec.clear();
        elevated_nodes_vec.extend_from_slice(std::slice::from_raw_parts(nodes_in, (curr_num_nodes * dim) as usize));
        return;
    }

    let mut current_nodes_is_in_elevated = false; 

    if num_steps % 2 == 1 { 
        ptr::copy_nonoverlapping(nodes_in, workspace_nodes_vec.as_mut_ptr(), (curr_num_nodes * dim) as usize);
        workspace_nodes_vec.set_len((curr_num_nodes*dim) as usize); 
        current_nodes_is_in_elevated = false; 
    } else { 
        ptr::copy_nonoverlapping(nodes_in, elevated_nodes_vec.as_mut_ptr(), (curr_num_nodes * dim) as usize);
        elevated_nodes_vec.set_len((curr_num_nodes*dim) as usize);
        current_nodes_is_in_elevated = true; 
    }

    for i in curr_num_nodes..final_num_nodes { 
        let source_ptr: *const c_double;
        let dest_ptr: *mut c_double;

        if current_nodes_is_in_elevated {
            source_ptr = elevated_nodes_vec.as_ptr();
            dest_ptr = workspace_nodes_vec.as_mut_ptr();
        } else {
            source_ptr = workspace_nodes_vec.as_ptr();
            dest_ptr = elevated_nodes_vec.as_mut_ptr();
        }
        curve::BEZ_elevate_nodes_curve(i, dim as c_int, source_ptr, dest_ptr);
        current_nodes_is_in_elevated = !current_nodes_is_in_elevated; 
        
        if current_nodes_is_in_elevated {
             elevated_nodes_vec.set_len(((i + 1) * dim) as usize);
        } else {
             workspace_nodes_vec.set_len(((i + 1) * dim) as usize);
        }
    }

    if !current_nodes_is_in_elevated { 
        ptr::copy_nonoverlapping(workspace_nodes_vec.as_ptr(), elevated_nodes_vec.as_mut_ptr(), (final_num_nodes * dim) as usize);
        elevated_nodes_vec.set_len((final_num_nodes * dim) as usize);
    }
}

/// # Safety
/// `nodes1` and `nodes2` must point to valid memory.
/// `elevated1_vec` and `elevated2_vec` will be modified (reallocated and length set).
unsafe fn make_same_degree_internal(
    num_nodes1_in: c_int,
    nodes1_in: *const c_double,
    num_nodes2_in: c_int,
    nodes2_in: *const c_double,
    elevated1_vec: &mut Vec<f64>, 
    elevated2_vec: &mut Vec<f64>, 
) -> c_int { 
    let dim = 2; 
    let final_num_nodes: c_int;
    
    let mut temp_workspace_for_elevate: Vec<f64> = Vec::new(); 

    if num_nodes1_in > num_nodes2_in {
        final_num_nodes = num_nodes1_in;
        let required_cap = (final_num_nodes * dim) as usize;
        if elevated1_vec.capacity() < required_cap { elevated1_vec.reserve(required_cap); }
        if elevated2_vec.capacity() < required_cap { elevated2_vec.reserve(required_cap); }

        elevated1_vec.clear();
        elevated1_vec.extend_from_slice(std::slice::from_raw_parts(nodes1_in, (num_nodes1_in * dim) as usize));
        elevate_helper_internal(num_nodes2_in, nodes2_in, final_num_nodes, &mut temp_workspace_for_elevate, elevated2_vec);

    } else if num_nodes2_in > num_nodes1_in {
        final_num_nodes = num_nodes2_in;
        let required_cap = (final_num_nodes * dim) as usize;
        if elevated1_vec.capacity() < required_cap { elevated1_vec.reserve(required_cap); }
        if elevated2_vec.capacity() < required_cap { elevated2_vec.reserve(required_cap); }

        elevated2_vec.clear();
        elevated2_vec.extend_from_slice(std::slice::from_raw_parts(nodes2_in, (num_nodes2_in * dim) as usize));
        elevate_helper_internal(num_nodes1_in, nodes1_in, final_num_nodes, &mut temp_workspace_for_elevate, elevated1_vec);
    } else { 
        final_num_nodes = num_nodes1_in;
        let required_cap = (final_num_nodes * dim) as usize;
        if elevated1_vec.capacity() < required_cap { elevated1_vec.reserve(required_cap); }
        if elevated2_vec.capacity() < required_cap { elevated2_vec.reserve(required_cap); }

        elevated1_vec.clear();
        elevated1_vec.extend_from_slice(std::slice::from_raw_parts(nodes1_in, (num_nodes1_in * dim) as usize));
        elevated2_vec.clear();
        elevated2_vec.extend_from_slice(std::slice::from_raw_parts(nodes2_in, (num_nodes2_in * dim) as usize));
    }
    final_num_nodes
}

/// # Safety
/// Caller must ensure pointers `nodes1`, `nodes2` are valid for reads as specified.
/// `num_nodes1`, `num_nodes2` must be accurate. `intersections` Vec can be reallocated.
unsafe fn add_coincident_parameters_internal(
    num_nodes1: c_int,
    nodes1: *const c_double,
    num_nodes2: c_int,
    nodes2: *const c_double,
    num_intersections: &mut c_int,
    intersections: &mut Vec<[c_double; 2]>, 
    coincident: &mut u8,
) {
    *coincident = 0; 
    let dim = 2_usize;

    let mut elevated1_vec: Vec<f64> = Vec::new();
    let mut elevated2_vec: Vec<f64> = Vec::new();
    let common_num_nodes = make_same_degree_internal(
        num_nodes1, nodes1, num_nodes2, nodes2, 
        &mut elevated1_vec, &mut elevated2_vec
    );

    let mut point_data = [0.0; 2];
    let mut s_initial = 0.0; let mut s_final = 0.0;
    let mut t_initial = 0.0; let mut t_final = 0.0;

    if common_num_nodes == 0 { return; } 

    ptr::copy_nonoverlapping(elevated2_vec.as_ptr(), point_data.as_mut_ptr(), dim); 
    curve::BEZ_locate_point_curve(common_num_nodes, dim as c_int, elevated1_vec.as_ptr(), point_data.as_ptr(), &mut s_initial);
    
    ptr::copy_nonoverlapping(elevated2_vec.as_ptr().add((common_num_nodes as usize - 1) * dim), point_data.as_mut_ptr(), dim); 
    curve::BEZ_locate_point_curve(common_num_nodes, dim as c_int, elevated1_vec.as_ptr(), point_data.as_ptr(), &mut s_final);

    if s_initial == curve::LOCATE_INVALID || s_final == curve::LOCATE_INVALID { return; }

    if s_initial != curve::LOCATE_MISS && s_final != curve::LOCATE_MISS {
        let mut specialized_nodes_vec: Vec<f64> = vec![0.0; (common_num_nodes * dim as c_int) as usize];
        curve::BEZ_specialize_curve(common_num_nodes, dim as c_int, elevated1_vec.as_ptr(), s_initial, s_final, specialized_nodes_vec.as_mut_ptr());
        
        if helpers::BEZ_vector_close((common_num_nodes * dim as c_int) as i32, specialized_nodes_vec.as_ptr(), elevated2_vec.as_ptr(), VECTOR_CLOSE_EPS) != 0 {
            *coincident = 1;
            intersections.clear(); *num_intersections = 0;
            add_intersection_internal(s_initial, 0.0, num_intersections, intersections);
            add_intersection_internal(s_final, 1.0, num_intersections, intersections);
        }
        return;
    }

    ptr::copy_nonoverlapping(elevated1_vec.as_ptr(), point_data.as_mut_ptr(), dim); 
    curve::BEZ_locate_point_curve(common_num_nodes, dim as c_int, elevated2_vec.as_ptr(), point_data.as_ptr(), &mut t_initial);

    ptr::copy_nonoverlapping(elevated1_vec.as_ptr().add((common_num_nodes as usize - 1) * dim), point_data.as_mut_ptr(), dim); 
    curve::BEZ_locate_point_curve(common_num_nodes, dim as c_int, elevated2_vec.as_ptr(), point_data.as_ptr(), &mut t_final);
    
    if t_initial == curve::LOCATE_INVALID || t_final == curve::LOCATE_INVALID { return; }
    if t_initial == curve::LOCATE_MISS && t_final == curve::LOCATE_MISS { return; }

    if t_initial != curve::LOCATE_MISS && t_final != curve::LOCATE_MISS {
        let mut specialized_nodes_vec: Vec<f64> = vec![0.0; (common_num_nodes * dim as c_int) as usize];
        curve::BEZ_specialize_curve(common_num_nodes, dim as c_int, elevated2_vec.as_ptr(), t_initial, t_final, specialized_nodes_vec.as_mut_ptr());

        if helpers::BEZ_vector_close((common_num_nodes * dim as c_int) as i32, elevated1_vec.as_ptr(), specialized_nodes_vec.as_ptr(), VECTOR_CLOSE_EPS) != 0 {
            *coincident = 1;
            intersections.clear(); *num_intersections = 0;
            add_intersection_internal(0.0, t_initial, num_intersections, intersections);
            add_intersection_internal(1.0, t_final, num_intersections, intersections);
        }
        return;
    }
    
    if s_initial == curve::LOCATE_MISS && s_final == curve::LOCATE_MISS { return; }

    let mut s_p1 = s_initial; 
    let mut s_p2 = s_final;
    let mut t_p1 = t_initial;
    let mut t_p2 = t_final;

    if s_p1 == curve::LOCATE_MISS { 
        s_p1 = s_p2; s_p2 = 1.0; 
        if t_p1 == curve::LOCATE_MISS { t_p1 = 1.0; } else { t_p2 = 1.0; }
    } else { // s_p2 == curve::LOCATE_MISS (since one of them must be LOCATE_MISS from earlier checks)
        s_p2 = s_p1; s_p1 = 0.0;
        if t_p1 == curve::LOCATE_MISS { t_p1 = 0.0; } else { t_p2 = 0.0; }
    }

    if (s_p1 - s_p2).abs() < MIN_INTERVAL_WIDTH && (t_p1 - t_p2).abs() < MIN_INTERVAL_WIDTH {
        return;
    }

    let mut specialized1_vec: Vec<f64> = vec![0.0; (common_num_nodes * dim as c_int) as usize];
    let mut specialized2_vec: Vec<f64> = vec![0.0; (common_num_nodes * dim as c_int) as usize];
    curve::BEZ_specialize_curve(common_num_nodes, dim as c_int, elevated1_vec.as_ptr(), s_p1, s_p2, specialized1_vec.as_mut_ptr());
    curve::BEZ_specialize_curve(common_num_nodes, dim as c_int, elevated2_vec.as_ptr(), t_p1, t_p2, specialized2_vec.as_mut_ptr());

    if helpers::BEZ_vector_close((common_num_nodes * dim as c_int) as i32, specialized1_vec.as_ptr(), specialized2_vec.as_ptr(), VECTOR_CLOSE_EPS) != 0 {
        *coincident = 1;
        intersections.clear(); *num_intersections = 0;
        add_intersection_internal(s_p1, t_p1, num_intersections, intersections);
        add_intersection_internal(s_p2, t_p2, num_intersections, intersections);
    }
}

/// # Safety
/// Modifies `candidates_vec` in place by retaining only elements for which `convex_hull_collide_internal` returns true.
/// `polygon1_scratch` and `polygon2_scratch` are used as temporary buffers.
unsafe fn prune_candidates_internal(
    candidates_vec: &mut Vec<[CurveData; 2]>,
    polygon1_scratch: &mut Vec<f64>,
    polygon2_scratch: &mut Vec<f64>,
) {
    let mut accepted_idx = 0;
    for i in 0..candidates_vec.len() {
        // To avoid issues with borrowing candidates_vec[i] and candidates_vec[accepted_idx] simultaneously for swap,
        // we can clone the candidate if it needs to be moved.
        // However, a retain_mut or manual drain_filter like approach is more idiomatic.
        // For now, let's do a simple swap if needed after checking.
        
        let curve1 = &candidates_vec[i][0];
        let curve2 = &candidates_vec[i][1];
        
        if curve1.num_nodes == 0 || curve2.num_nodes == 0 { continue; }

        let mut collision_u8: u8 = 0;
        convex_hull_collide_internal(
            curve1.num_nodes as c_int, curve1.nodes.as_ptr(), polygon1_scratch,
            curve2.num_nodes as c_int, curve2.nodes.as_ptr(), polygon2_scratch,
            &mut collision_u8
        );

        if collision_u8 != 0 { // If they collide
            if accepted_idx < i {
                // Move the accepted candidate to the `accepted_idx` position.
                // This is tricky with borrowing rules. A simple swap is fine if order doesn't matter beyond partitioning.
                // If order must be preserved, then a temporary store or clone is needed.
                // Fortran does `candidates(:, accepted) = candidates(:, i)` which is a copy.
                // Let's use swap for efficiency, assuming relative order of kept items isn't critical.
                candidates_vec.swap(accepted_idx, i);
            }
            accepted_idx += 1;
        }
    }
    candidates_vec.truncate(accepted_idx);
}


/// # Safety
/// Caller must ensure pointers `nodes_first`, `nodes_second` are valid.
/// `intersections_output_vec` will be cleared and populated.
unsafe fn all_intersections_internal(
    num_nodes_first: c_int,
    nodes_first: *const c_double,
    num_nodes_second: c_int,
    nodes_second: *const c_double,
    intersections_output_vec: &mut Vec<[c_double; 2]>,
    num_intersections: &mut c_int,
    coincident: &mut u8,
    status: &mut c_int,
) {
    *status = status::STATUS_SUCCESS;
    intersections_output_vec.clear();
    *num_intersections = 0;

    let mut both_linear_u8: u8 = 0;
    check_lines_internal(
        num_nodes_first, nodes_first, num_nodes_second, nodes_second,
        &mut both_linear_u8, coincident, intersections_output_vec, num_intersections
    );
    if both_linear_u8 != 0 { return; }

    intersections_output_vec.clear(); *num_intersections = 0;
    *coincident = 0;

    let dim = 2_usize;
    let initial_curve1 = CurveData::new(0.0, 1.0, Vec::from(std::slice::from_raw_parts(nodes_first, (num_nodes_first * dim as c_int) as usize)), num_nodes_first as usize, dim);
    let initial_curve2 = CurveData::new(0.0, 1.0, Vec::from(std::slice::from_raw_parts(nodes_second, (num_nodes_second * dim as c_int) as usize)), num_nodes_second as usize, dim);
    
    let mut candidates_storage1: Vec<[CurveData; 2]> = Vec::with_capacity(MAX_CANDIDATES as usize + 4); 
    let mut candidates_storage2: Vec<[CurveData; 2]> = Vec::with_capacity(MAX_CANDIDATES as usize + 4);
    candidates_storage1.push([initial_curve1, initial_curve2]);

    let mut polygon1_scratch: Vec<f64> = Vec::new();
    let mut polygon2_scratch: Vec<f64> = Vec::new();
    
    let mut current_is_storage1 = true;

    for iter_count in 1..=MAX_INTERSECT_SUBDIVISIONS {
        let (current_candidates_slice, next_candidates_mut_ref) = if current_is_storage1 {
            (candidates_storage1.as_slice(), &mut candidates_storage2)
        } else {
            (candidates_storage2.as_slice(), &mut candidates_storage1)
        };
        
        intersect_one_round_internal(
            nodes_first, num_nodes_first, nodes_second, num_nodes_second,
            current_candidates_slice, 
            num_intersections, intersections_output_vec,
            next_candidates_mut_ref, status,
            &mut polygon1_scratch, &mut polygon2_scratch
        );

        if *status != status::STATUS_SUCCESS { return; }
        current_is_storage1 = !current_is_storage1; 

        let active_candidates_list = if current_is_storage1 { &mut candidates_storage1 } else { &mut candidates_storage2 };
        let mut num_active_candidates = active_candidates_list.len();


        if num_active_candidates > MAX_CANDIDATES as usize {
            prune_candidates_internal(active_candidates_list, &mut polygon1_scratch, &mut polygon2_scratch);
            num_active_candidates = active_candidates_list.len(); // Update after pruning
            
            if num_active_candidates > MAX_CANDIDATES as usize {
                add_coincident_parameters_internal(
                    num_nodes_first, nodes_first, num_nodes_second, nodes_second,
                    num_intersections, intersections_output_vec, coincident
                );
                if *coincident == 0 {
                    *status = num_active_candidates as c_int; 
                }
                return;
            }
        }
        if num_active_candidates == 0 { return; }
    }
    *status = status::STATUS_NO_CONVERGE;
}


#[no_mangle]
/// # Safety
/// Caller must ensure pointers `nodes_first`, `nodes_second`, `intersections_ptr` are valid for reads/writes as specified.
/// `num_nodes` and `intersections_size` must be accurate.
pub unsafe extern "C" fn BEZ_curve_intersections(
    num_nodes_first: c_int,
    nodes_first: *const c_double,
    num_nodes_second: c_int,
    nodes_second: *const c_double,
    intersections_size: c_int,
    intersections_ptr: *mut c_double, 
    num_intersections_found: *mut c_int,
    coincident: *mut u8, 
    status_code: *mut c_int,
) {
    let mut intersections_vec: Vec<[c_double; 2]> = Vec::with_capacity(intersections_size.max(0) as usize);
    
    all_intersections_internal(
        num_nodes_first, nodes_first,
        num_nodes_second, nodes_second,
        &mut intersections_vec,
        &mut *num_intersections_found, 
        &mut *coincident,
        &mut *status_code
    );

    if *status_code != status::STATUS_SUCCESS {
        if !(*status_code >= MAX_CANDIDATES && *coincident == 1) {
             return;
        }
    }

    let actual_intersections = *num_intersections_found;
    if actual_intersections > intersections_size {
        *status_code = status::STATUS_INSUFFICIENT_SPACE;
    } else if actual_intersections > 0 {
        for i in 0..(actual_intersections as usize) {
            *intersections_ptr.add(i * 2 + 0) = intersections_vec[i][0];
            *intersections_ptr.add(i * 2 + 1) = intersections_vec[i][1];
        }
        if *status_code >= MAX_CANDIDATES && *coincident == 1 {
             *status_code = status::STATUS_SUCCESS;
        }
    }
}

// Final TODO: 
// - Review all functions for correctness, especially array indexing, loop bounds, and logic flow.
// - Ensure all module-level Fortran variables used in helpers (like POLYGON1, POLYGON2 for convex_hull_collide)
//   are correctly handled as local variables or parameters in Rust. (Done via scratch parameters)
// - The logic for CANDIDATES_ODD/EVEN in all_intersections_internal might need refinement for optimal performance/memory. (Using swap with two Vecs)
// - Refine `newton_iterate_internal` convergence check to exactly match Fortran's if needed. (Done)
// - Correct `original_num_nodes` usage in `from_linearized_internal`. (Done by passing explicitly original_num_nodesX)
