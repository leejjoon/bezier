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
#![allow(clippy::upper_case_acronyms)] // For FFI function names like BEZ_...
#![allow(clippy::manual_range_contains)] // For s_approx checks

use libc::{c_double, c_int, c_void};
use std::ptr;
use std::alloc::{alloc, dealloc, Layout}; // Required for manual allocation if used for CurveData


use crate::status;
use crate::helpers::{self, WIGGLE}; 
use crate::quadpack_replacement::{self, rust_dqagse};


// Constants from Fortran module `curve`
pub const MAX_LOCATE_SUBDIVISIONS: c_int = 20;
pub const LOCATE_STD_CAP: f64 = 9.5367431640625e-07; // 0.5_f64.powi(20)
pub const SQRT_PREC: f64 = 1.4901161193847656e-08; // 0.5_f64.powi(26)
pub const REDUCE_THRESHOLD: f64 = SQRT_PREC;
pub const LOCATE_MISS: f64 = -1.0;
pub const LOCATE_INVALID: f64 = -2.0;

// Internal struct for locate_point logic
#[derive(Clone)] 
struct InternalCurveCandidate {
    start_param: f64,
    end_param: f64,
    nodes: Vec<f64>, 
}

impl InternalCurveCandidate {
    fn new_empty(num_nodes: usize, dimension: usize) -> Self {
        InternalCurveCandidate {
            start_param: 0.0,
            end_param: 0.0,
            nodes: vec![0.0; num_nodes * dimension],
        }
    }
    fn new_from_nodes(
        start_p: f64, 
        end_p: f64, 
        nodes_slice: &[f64] // Assumes nodes_slice is dimension * num_nodes
    ) -> Self {
        InternalCurveCandidate {
            start_param: start_p,
            end_param: end_p,
            nodes: nodes_slice.to_vec(),
        }
    }
}


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

// Helper to copy a column (node) from src_nodes at src_node_idx to dest_col_ptr
#[inline]
unsafe fn get_node_col(
    src_nodes: *const c_double, 
    dimension: usize, 
    src_node_idx: usize, 
    dest_col_ptr: *mut c_double
) {
    ptr::copy_nonoverlapping(
        src_nodes.add(src_node_idx * dimension), 
        dest_col_ptr, 
        dimension
    );
}

// Helper to copy src_col_ptr into a column (node) at dest_node_idx in dest_nodes
#[inline]
unsafe fn set_node_col(
    dest_nodes: *mut c_double,
    dimension: usize,
    dest_node_idx: usize,
    src_col_ptr: *const c_double
) {
    ptr::copy_nonoverlapping(
        src_col_ptr,
        dest_nodes.add(dest_node_idx * dimension),
        dimension
    );
}

// Internal Rust equivalent of Fortran's `evaluate_curve_vs`
unsafe fn evaluate_curve_vs_internal(
    num_nodes: usize,
    dimension: usize,
    nodes: *const c_double,
    num_vals: usize,
    lambda1: *const c_double,
    lambda2: *const c_double,
    evaluated: *mut c_double,
) {
    if num_nodes == 0 || num_vals == 0 {
        return;
    }

    let mut lambda2_pow_vec: Vec<f64> = vec![1.0; num_vals];
    let mut binom_val: f64 = 1.0;

    for val_idx in 0..num_vals {
        let l1_val = *lambda1.add(val_idx);
        for dim_idx in 0..dimension {
            let node_0_dim_val = get_node_val(nodes, dimension, 0, dim_idx);
            set_node_val(evaluated, dimension, val_idx, dim_idx, l1_val * node_0_dim_val);
        }
    }

    for node_idx_1based in 2..num_nodes { 
        for val_idx_pow in 0..num_vals {
            lambda2_pow_vec[val_idx_pow] *= *lambda2.add(val_idx_pow);
        }
        binom_val = (binom_val * (num_nodes - node_idx_1based + 1) as f64) / (node_idx_1based - 1) as f64;

        let current_node_0idx = node_idx_1based - 1;
        for val_idx_eval in 0..num_vals {
            let l1_val = *lambda1.add(val_idx_eval);
            let l2_pow_val = lambda2_pow_vec[val_idx_eval];
            for dim_idx in 0..dimension {
                let current_eval_val = get_node_val(evaluated, dimension, val_idx_eval, dim_idx);
                let node_i_dim_val = get_node_val(nodes, dimension, current_node_0idx, dim_idx);
                let updated_val = (current_eval_val + binom_val * l2_pow_val * node_i_dim_val) * l1_val;
                set_node_val(evaluated, dimension, val_idx_eval, dim_idx, updated_val);
            }
        }
    }
    
    for val_idx in 0..num_vals {
        let l2_val = *lambda2.add(val_idx);
        let final_l2_pow = lambda2_pow_vec[val_idx] * l2_val;
        for dim_idx in 0..dimension {
            let current_eval_val = get_node_val(evaluated, dimension, val_idx, dim_idx);
            let last_node_dim_val = get_node_val(nodes, dimension, num_nodes - 1, dim_idx);
            let updated_val = current_eval_val + final_l2_pow * last_node_dim_val;
            set_node_val(evaluated, dimension, val_idx, dim_idx, updated_val);
        }
    }
}

// Internal Rust equivalent of Fortran's `evaluate_curve_de_casteljau`
unsafe fn evaluate_curve_de_casteljau_internal(
    num_nodes: usize,
    dimension: usize,
    nodes: *const c_double,
    num_vals: usize,
    lambda1: *const c_double,
    lambda2: *const c_double,
    evaluated: *mut c_double, 
) {
    if num_nodes == 0 || num_vals == 0 { return; }
    if num_nodes == 1 {
        for val_idx in 0..num_vals {
            for dim_idx in 0..dimension {
                set_node_val(evaluated, dimension, val_idx, dim_idx, get_node_val(nodes, dimension, 0, dim_idx));
            }
        }
        return;
    }

    let workspace_size = dimension * num_vals * (num_nodes - 1);
    let mut workspace_vec: Vec<f64> = vec![0.0; workspace_size];
    let workspace_ptr = workspace_vec.as_mut_ptr();

    let get_ws = |d: usize, v: usize, n_minus_1_idx: usize| -> f64 {
        *workspace_ptr.add(n_minus_1_idx * (dimension * num_vals) + v * dimension + d)
    };
    let set_ws = |d: usize, v: usize, n_minus_1_idx: usize, val: f64| {
        *workspace_ptr.add(n_minus_1_idx * (dimension * num_vals) + v * dimension + d) = val;
    };
    
    for k_node_stage in 0..(num_nodes - 1) { 
        for val_idx in 0..num_vals {
            let l1v = *lambda1.add(val_idx);
            let l2v = *lambda2.add(val_idx);
            for dim_idx in 0..dimension {
                let val = l1v * get_node_val(nodes, dimension, k_node_stage, dim_idx) +
                          l2v * get_node_val(nodes, dimension, k_node_stage + 1, dim_idx);
                set_ws(dim_idx, val_idx, k_node_stage, val);
            }
        }
    }

    for i_len_of_segment_array in (1..(num_nodes - 1)).rev() { 
        for k_node_stage in 0..i_len_of_segment_array { 
            for val_idx in 0..num_vals {
                let l1v = *lambda1.add(val_idx);
                let l2v = *lambda2.add(val_idx);
                for dim_idx in 0..dimension {
                    let val = l1v * get_ws(dim_idx, val_idx, k_node_stage) +
                              l2v * get_ws(dim_idx, val_idx, k_node_stage + 1);
                    set_ws(dim_idx, val_idx, k_node_stage, val);
                }
            }
        }
    }

    for val_idx in 0..num_vals {
        for dim_idx in 0..dimension {
            set_node_val(evaluated, dimension, val_idx, dim_idx, get_ws(dim_idx, val_idx, 0));
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_evaluate_curve_barycentric(
    num_nodes: c_int,
    dimension: c_int,
    nodes: *const c_double,
    num_vals: c_int,
    lambda1: *const c_double,
    lambda2: *const c_double,
    evaluated: *mut c_double,
) {
    let nn = num_nodes as usize;
    let dim = dimension as usize;
    let nv = num_vals as usize;

    if num_nodes > 55 {
        evaluate_curve_de_casteljau_internal(nn, dim, nodes, nv, lambda1, lambda2, evaluated);
    } else {
        evaluate_curve_vs_internal(nn, dim, nodes, nv, lambda1, lambda2, evaluated);
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_evaluate_multi(
    num_nodes: c_int,
    dimension: c_int,
    nodes: *const c_double,
    num_vals: c_int,
    s_vals: *const c_double, 
    evaluated: *mut c_double, 
) {
    let nv = num_vals as usize;
    if nv == 0 { return; }

    let mut one_less_s_vec: Vec<f64> = Vec::with_capacity(nv);
    for i in 0..nv {
        one_less_s_vec.push(1.0 - *s_vals.add(i));
    }
    
    BEZ_evaluate_curve_barycentric(
        num_nodes, 
        dimension, 
        nodes, 
        num_vals, 
        one_less_s_vec.as_ptr(), 
        s_vals, 
        evaluated
    );
}

unsafe fn specialize_curve_generic_internal(
    num_nodes: usize,
    dimension: usize,
    nodes: *const c_double,
    start_s: f64, 
    end_s: f64,   
    new_nodes: *mut c_double, 
) {
    if num_nodes == 0 { return; }

    let block_size = dimension * (num_nodes - 1); 
    let mut workspace_vec: Vec<f64> = vec![0.0; block_size * num_nodes];
    let ws_ptr = workspace_vec.as_mut_ptr();

    let get_ws = |d: usize, r_idx: usize, c_block_idx: usize| -> f64 {
        *ws_ptr.add(c_block_idx * block_size + r_idx * dimension + d)
    };
    let set_ws = |d: usize, r_idx: usize, c_block_idx: usize, val: f64| {
        *ws_ptr.add(c_block_idx * block_size + r_idx * dimension + d) = val;
    };

    let minus_start = 1.0 - start_s;
    let minus_end = 1.0 - end_s;

    for r_idx in 0..(num_nodes - 1) {
        for d_idx in 0..dimension {
            let val1 = minus_start * get_node_val(nodes, dimension, r_idx, d_idx) +
                       start_s * get_node_val(nodes, dimension, r_idx + 1, d_idx);
            set_ws(d_idx, r_idx, 0, val1);

            let val2 = minus_end * get_node_val(nodes, dimension, r_idx, d_idx) +
                       end_s * get_node_val(nodes, dimension, r_idx + 1, d_idx);
            set_ws(d_idx, r_idx, 1, val2);
        }
    }
    
    let mut curr_size = num_nodes - 1; 
    for index_1based in 3..=num_nodes {
        curr_size -= 1; 
        let index_0based = index_1based - 1; 

        for r_idx in 0..curr_size {
            for d_idx in 0..dimension {
                let val = minus_end * get_ws(d_idx, r_idx, index_0based - 1) +
                          end_s * get_ws(d_idx, r_idx + 1, index_0based - 1);
                set_ws(d_idx, r_idx, index_0based, val);
            }
        }

        for j_0based_block_idx in 0..index_0based { 
             for r_idx in 0..curr_size {
                for d_idx in 0..dimension {
                    let val = minus_start * get_ws(d_idx, r_idx, j_0based_block_idx) +
                              start_s * get_ws(d_idx, r_idx + 1, j_0based_block_idx);
                    set_ws(d_idx, r_idx, j_0based_block_idx, val);
                }
            }
        }
    }

    for i_node_idx in 0..num_nodes {
        for d_idx in 0..dimension {
            set_node_val(new_nodes, dimension, i_node_idx, d_idx, get_ws(d_idx, 0, i_node_idx));
        }
    }
}

unsafe fn specialize_curve_quadratic_internal(
    dimension: usize,
    nodes: *const c_double, 
    start_s: f64,
    end_s: f64,
    new_nodes: *mut c_double, 
) {
    let minus_start = 1.0 - start_s;
    let minus_end = 1.0 - end_s;
    let prod_both = start_s * end_s;

    for d_idx in 0..dimension {
        let n0 = get_node_val(nodes, dimension, 0, d_idx); 
        let n1 = get_node_val(nodes, dimension, 1, d_idx); 
        let n2 = get_node_val(nodes, dimension, 2, d_idx); 

        set_node_val(new_nodes, dimension, 0, d_idx, 
            minus_start * minus_start * n0 +
            2.0 * start_s * minus_start * n1 +
            start_s * start_s * n2
        );
        set_node_val(new_nodes, dimension, 1, d_idx,
            minus_start * minus_end * n0 +
            (end_s + start_s - 2.0 * prod_both) * n1 +
            prod_both * n2
        );
        set_node_val(new_nodes, dimension, 2, d_idx,
            minus_end * minus_end * n0 +
            2.0 * end_s * minus_end * n1 +
            end_s * end_s * n2
        );
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_specialize_curve(
    num_nodes: c_int,
    dimension: c_int,
    nodes: *const c_double,
    start_s: c_double, 
    end_s: c_double,   
    new_nodes: *mut c_double,
) {
    let nn = num_nodes as usize;
    let dim = dimension as usize;

    if nn == 0 { return; } 
    if nn == 1 { 
        for d_idx in 0..dim {
            set_node_val(new_nodes, dim, 0, d_idx, get_node_val(nodes, dim, 0, d_idx));
        }
    } else if nn == 2 { 
        for d_idx in 0..dim {
            let n0 = get_node_val(nodes, dim, 0, d_idx);
            let n1 = get_node_val(nodes, dim, 1, d_idx);
            set_node_val(new_nodes, dim, 0, d_idx, (1.0 - start_s) * n0 + start_s * n1);
            set_node_val(new_nodes, dim, 1, d_idx, (1.0 - end_s) * n0 + end_s * n1);
        }
    } else if nn == 3 { 
        specialize_curve_quadratic_internal(dim, nodes, start_s, end_s, new_nodes);
    } else {
        specialize_curve_generic_internal(nn, dim, nodes, start_s, end_s, new_nodes);
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_evaluate_hodograph(
    s: c_double,
    num_nodes: c_int,
    dimension: c_int,
    nodes: *const c_double,
    hodograph: *mut c_double, 
) {
    let nn = num_nodes as usize;
    let dim = dimension as usize;

    if nn < 2 { 
        for d_idx in 0..dim {
            set_node_val(hodograph, dim, 0, d_idx, 0.0);
        }
        return;
    }

    let deriv_num_nodes = nn - 1;
    let mut first_deriv_vec: Vec<f64> = vec![0.0; dim * deriv_num_nodes];
    let first_deriv_ptr = first_deriv_vec.as_mut_ptr();

    for deriv_node_idx in 0..deriv_num_nodes { 
        for d_idx in 0..dim {
            let val = get_node_val(nodes, dim, deriv_node_idx + 1, d_idx) - 
                      get_node_val(nodes, dim, deriv_node_idx, d_idx);
            set_node_val(first_deriv_ptr, dim, deriv_node_idx, d_idx, val);
        }
    }
    
    let s_val_arr: [f64;1] = [s];
    BEZ_evaluate_multi(
        deriv_num_nodes as c_int, 
        dimension,
        first_deriv_ptr,
        1, 
        s_val_arr.as_ptr(),
        hodograph, 
    );

    let degree_factor = (num_nodes - 1) as f64;
    for d_idx in 0..dim {
        let current_val = get_node_val(hodograph, dim, 0, d_idx);
        set_node_val(hodograph, dim, 0, d_idx, degree_factor * current_val);
    }
}

unsafe fn subdivide_nodes_generic_internal(
    num_nodes: usize,
    dimension: usize,
    nodes: *const c_double,
    left_nodes: *mut c_double,
    right_nodes: *mut c_double,
) {
    if num_nodes == 0 { return; }

    let mut pascals_triangle: Vec<f64> = vec![0.0; num_nodes];
    pascals_triangle[0] = 1.0;

    for elt_idx_1based in 1..=num_nodes { 
        let elt_idx_0based = elt_idx_1based - 1;

        if elt_idx_1based > 1 {
            let mut temp_pt_slice_reversed: Vec<f64> = pascals_triangle[0..elt_idx_1based].iter().cloned().rev().collect();
            for j_idx_0based in 0..elt_idx_1based {
                pascals_triangle[j_idx_0based] = 0.5 * (pascals_triangle[j_idx_0based] + temp_pt_slice_reversed[j_idx_0based]);
            }
        }

        for d_idx in 0..dimension {
            set_node_val(left_nodes, dimension, elt_idx_0based, d_idx, 0.0);
            set_node_val(right_nodes, dimension, num_nodes - 1 - elt_idx_0based, d_idx, 0.0);
        }
        
        for pascal_idx_1based in 1..=elt_idx_1based { 
            let pascal_idx_0based = pascal_idx_1based - 1;
            let weight = pascals_triangle[pascal_idx_0based];
            
            for d_idx in 0..dimension {
                let current_left_val = get_node_val(left_nodes, dimension, elt_idx_0based, d_idx);
                let node_to_add_left = get_node_val(nodes, dimension, pascal_idx_0based, d_idx);
                set_node_val(left_nodes, dimension, elt_idx_0based, d_idx, current_left_val + weight * node_to_add_left);

                let current_right_val = get_node_val(right_nodes, dimension, num_nodes - 1 - elt_idx_0based, d_idx);
                let node_to_add_right = get_node_val(nodes, dimension, num_nodes - 1 - pascal_idx_0based, d_idx);
                set_node_val(right_nodes, dimension, num_nodes - 1 - elt_idx_0based, d_idx, current_right_val + weight * node_to_add_right);
            }
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_subdivide_nodes_curve(
    num_nodes: c_int,
    dimension: c_int,
    nodes: *const c_double,
    left_nodes: *mut c_double,
    right_nodes: *mut c_double,
) {
    let nn = num_nodes as usize;
    let dim = dimension as usize;

    if nn == 0 { return; }

    if nn == 2 { 
        for d_idx in 0..dim {
            let n0 = get_node_val(nodes, dim, 0, d_idx);
            let n1 = get_node_val(nodes, dim, 1, d_idx);
            set_node_val(left_nodes, dim, 0, d_idx, n0);
            let mid = 0.5 * (n0 + n1);
            set_node_val(left_nodes, dim, 1, d_idx, mid);
            set_node_val(right_nodes, dim, 0, d_idx, mid);
            set_node_val(right_nodes, dim, 1, d_idx, n1);
        }
    } else if nn == 3 { 
        for d_idx in 0..dim {
            let n0 = get_node_val(nodes, dim, 0, d_idx);
            let n1 = get_node_val(nodes, dim, 1, d_idx);
            let n2 = get_node_val(nodes, dim, 2, d_idx);
            set_node_val(left_nodes, dim, 0, d_idx, n0);
            set_node_val(left_nodes, dim, 1, d_idx, 0.5 * (n0 + n1));
            let mid = 0.25 * (n0 + 2.0 * n1 + n2);
            set_node_val(left_nodes, dim, 2, d_idx, mid);
            set_node_val(right_nodes, dim, 0, d_idx, mid);
            set_node_val(right_nodes, dim, 1, d_idx, 0.5 * (n1 + n2));
            set_node_val(right_nodes, dim, 2, d_idx, n2);
        }
    } else if nn == 4 { 
        for d_idx in 0..dim {
            let n0 = get_node_val(nodes, dim, 0, d_idx);
            let n1 = get_node_val(nodes, dim, 1, d_idx);
            let n2 = get_node_val(nodes, dim, 2, d_idx);
            let n3 = get_node_val(nodes, dim, 3, d_idx);
            set_node_val(left_nodes, dim, 0, d_idx, n0);
            set_node_val(left_nodes, dim, 1, d_idx, 0.5 * (n0 + n1));
            set_node_val(left_nodes, dim, 2, d_idx, 0.25 * (n0 + 2.0 * n1 + n2));
            let mid = 0.125 * (n0 + 3.0 * n1 + 3.0 * n2 + n3);
            set_node_val(left_nodes, dim, 3, d_idx, mid);
            set_node_val(right_nodes, dim, 0, d_idx, mid);
            set_node_val(right_nodes, dim, 1, d_idx, 0.25 * (n1 + 2.0 * n2 + n3));
            set_node_val(right_nodes, dim, 2, d_idx, 0.5 * (n2 + n3));
            set_node_val(right_nodes, dim, 3, d_idx, n3);
        }
    } else {
        subdivide_nodes_generic_internal(nn, dim, nodes, left_nodes, right_nodes);
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_newton_refine_curve(
    num_nodes: c_int,
    dimension: c_int,
    nodes: *const c_double,
    point: *const c_double, 
    s: c_double,
    updated_s: *mut c_double,
) {
    let dim = dimension as usize;

    let mut pt_delta_vec: Vec<f64> = vec![0.0; dim]; 
    let pt_delta_ptr = pt_delta_vec.as_mut_ptr();
    let mut derivative_vec: Vec<f64> = vec![0.0; dim];
    let derivative_ptr = derivative_vec.as_mut_ptr();
    
    let s_val_arr: [f64;1] = [s];
    BEZ_evaluate_multi(
        num_nodes,
        dimension,
        nodes,
        1, 
        s_val_arr.as_ptr(),
        pt_delta_ptr, 
    );

    for d_idx in 0..dim {
        let p_val = get_node_val(point, dim, 0, d_idx); 
        let bs_val = get_node_val(pt_delta_ptr, dim, 0, d_idx);
        set_node_val(pt_delta_ptr, dim, 0, d_idx, p_val - bs_val);
    }

    BEZ_evaluate_hodograph(
        s,
        num_nodes,
        dimension,
        nodes,
        derivative_ptr, 
    );

    let mut dot_pt_deriv = 0.0;
    let mut dot_deriv_deriv = 0.0;
    for d_idx in 0..dim {
        let pt_d_val = get_node_val(pt_delta_ptr, dim, 0, d_idx);
        let deriv_d_val = get_node_val(derivative_ptr, dim, 0, d_idx);
        dot_pt_deriv += pt_d_val * deriv_d_val;
        dot_deriv_deriv += deriv_d_val * deriv_d_val;
    }

    if dot_deriv_deriv == 0.0 {
        *updated_s = s; 
    } else {
        *updated_s = s + dot_pt_deriv / dot_deriv_deriv;
    }
}

unsafe fn split_candidate_internal(
    num_nodes: usize,
    dimension: usize,
    candidate: &InternalCurveCandidate, 
    next_candidate1: &mut InternalCurveCandidate,
    next_candidate2: &mut InternalCurveCandidate,
) {
    BEZ_subdivide_nodes_curve(
        num_nodes as c_int,
        dimension as c_int,
        candidate.nodes.as_ptr(),
        next_candidate1.nodes.as_mut_ptr(),
        next_candidate2.nodes.as_mut_ptr(),
    );

    next_candidate1.start_param = candidate.start_param;
    next_candidate1.end_param = 0.5 * (candidate.start_param + candidate.end_param);
    next_candidate2.start_param = next_candidate1.end_param;
    next_candidate2.end_param = candidate.end_param;
}

unsafe fn update_candidates_internal(
    num_nodes: usize,
    dimension: usize,
    point_ptr: *const c_double, 
    current_candidates: &[InternalCurveCandidate], 
    next_candidates_vec: &mut Vec<InternalCurveCandidate>, 
) {
    next_candidates_vec.clear(); 

    for candidate in current_candidates {
        let mut predicate_val: u8 = 0;
        helpers::BEZ_contains_nd(
            num_nodes as c_int,
            dimension as c_int,
            candidate.nodes.as_ptr(),
            point_ptr,
            &mut predicate_val,
        );

        if predicate_val != 0 { 
            if next_candidates_vec.capacity() < next_candidates_vec.len() + 2 {
                next_candidates_vec.reserve(next_candidates_vec.len().max(2)); 
            }

            let mut next_cand1 = InternalCurveCandidate::new_empty(num_nodes, dimension);
            let mut next_cand2 = InternalCurveCandidate::new_empty(num_nodes, dimension);
            
            split_candidate_internal(
                num_nodes,
                dimension,
                candidate,
                &mut next_cand1,
                &mut next_cand2,
            );
            next_candidates_vec.push(next_cand1);
            next_candidates_vec.push(next_cand2);
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_locate_point_curve(
    num_nodes_c: c_int,
    dimension_c: c_int,
    nodes: *const c_double,
    point: *const c_double, 
    s_approx: *mut c_double,
) {
    let num_nodes = num_nodes_c as usize;
    let dimension = dimension_c as usize;

    if num_nodes == 0 { 
        *s_approx = LOCATE_MISS;
        return;
    }

    let mut candidates1: Vec<InternalCurveCandidate> = Vec::with_capacity(2 * MAX_LOCATE_SUBDIVISIONS as usize);
    let mut candidates2: Vec<InternalCurveCandidate> = Vec::with_capacity(2 * MAX_LOCATE_SUBDIVISIONS as usize);

    let initial_nodes_slice = std::slice::from_raw_parts(nodes, num_nodes * dimension);
    candidates1.push(InternalCurveCandidate::new_from_nodes(0.0, 1.0, initial_nodes_slice));
    
    let mut current_candidates_list = &mut candidates1;
    let mut next_candidates_list = &mut candidates2;

    *s_approx = LOCATE_MISS; 

    for _sub_index in 0..(MAX_LOCATE_SUBDIVISIONS + 1) {
        update_candidates_internal(
            num_nodes, 
            dimension, 
            point, 
            current_candidates_list, 
            next_candidates_list
        );

        if next_candidates_list.is_empty() {
            return; 
        }
        std::mem::swap(&mut current_candidates_list, &mut next_candidates_list);
    }
    
    let final_candidates = current_candidates_list; 
    if final_candidates.is_empty() { 
        return;
    }
    
    let num_final_candidates = final_candidates.len();
    let mut s_params_sum = 0.0;
    let mut s_params_vec: Vec<f64> = Vec::with_capacity(2 * num_final_candidates);

    for candidate in final_candidates.iter() {
        s_params_vec.push(candidate.start_param);
        s_params_vec.push(candidate.end_param);
        s_params_sum += candidate.start_param + candidate.end_param;
    }
    
    let mean_s = s_params_sum / (2.0 * num_final_candidates as f64);

    let mut sum_sq_diff = 0.0;
    for s_p in s_params_vec.iter() {
        sum_sq_diff += (*s_p - mean_s) * (*s_p - mean_s);
    }
    let std_dev = (sum_sq_diff / (2.0 * num_final_candidates as f64)).sqrt();

    if std_dev > LOCATE_STD_CAP {
        *s_approx = LOCATE_INVALID;
        return;
    }

    let mut refined_s = 0.0;
    BEZ_newton_refine_curve(
        num_nodes_c, 
        dimension_c, 
        nodes, 
        point, 
        mean_s, 
        &mut refined_s
    );

    if refined_s < 0.0 {
        *s_approx = 0.0;
    } else if refined_s > 1.0 {
        *s_approx = 1.0;
    } else {
        *s_approx = refined_s;
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_elevate_nodes_curve(
    num_nodes_c: c_int,
    dimension_c: c_int,
    nodes: *const c_double,
    elevated: *mut c_double, 
) {
    let num_nodes = num_nodes_c as usize;
    let dimension = dimension_c as usize;

    if num_nodes == 0 { return; }

    for d_idx in 0..dimension {
        set_node_val(elevated, dimension, 0, d_idx, get_node_val(nodes, dimension, 0, d_idx));
    }

    for i_1based in 1..num_nodes { 
        let elevated_node_idx = i_1based;         
        let weight1 = i_1based as f64 / (num_nodes_c as f64); 
        let weight0 = 1.0 - weight1; 

        for d_idx in 0..dimension {
            let p_i = get_node_val(nodes, dimension, i_1based, d_idx);
            let p_i_minus_1 = get_node_val(nodes, dimension, i_1based - 1, d_idx);
            let val_corrected = weight0 * p_i + weight1 * p_i_minus_1;
            set_node_val(elevated, dimension, elevated_node_idx, d_idx, val_corrected);
        }
    }
    
    for d_idx in 0..dimension {
        set_node_val(elevated, dimension, num_nodes, d_idx, get_node_val(nodes, dimension, num_nodes - 1, d_idx));
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_get_curvature(
    num_nodes_c: c_int,
    nodes: *const c_double, 
    tangent_vec_ptr: *const c_double, 
    s: c_double,
    curvature: *mut c_double,
) {
    let num_nodes = num_nodes_c as usize;
    let dimension: usize = 2; 

    if num_nodes < 3 { 
        *curvature = 0.0;
        return;
    }

    let first_deriv_num_nodes = num_nodes - 1;
    let mut work_vec: Vec<f64> = vec![0.0; dimension * first_deriv_num_nodes];
    let work_ptr = work_vec.as_mut_ptr();

    for node_idx in 0..first_deriv_num_nodes {
        for d_idx in 0..dimension {
            let val = get_node_val(nodes, dimension, node_idx + 1, d_idx) -
                      get_node_val(nodes, dimension, node_idx, d_idx);
            set_node_val(work_ptr, dimension, node_idx, d_idx, val);
        }
    }

    let second_deriv_num_nodes = first_deriv_num_nodes - 1; 
    
    for node_idx in 0..second_deriv_num_nodes {
        for d_idx in 0..dimension {
            let val = get_node_val(work_ptr, dimension, node_idx + 1, d_idx) -
                      get_node_val(work_ptr, dimension, node_idx, d_idx);
            set_node_val(work_ptr, dimension, node_idx, d_idx, val);
        }
    }
    
    let mut concavity_vec: Vec<f64> = vec![0.0; dimension]; 
    let concavity_ptr = concavity_vec.as_mut_ptr();
    let s_val_arr: [f64;1] = [s];

    BEZ_evaluate_multi(
        second_deriv_num_nodes as c_int, 
        dimension as c_int, 
        work_ptr,    
        1,           
        s_val_arr.as_ptr(),
        concavity_ptr,
    );

    let scale_factor = (num_nodes - 1) as f64 * (num_nodes - 2) as f64;
    for d_idx in 0..dimension {
        concavity_vec[d_idx] *= scale_factor; 
    }

    let mut cp_result = 0.0;
    helpers::BEZ_cross_product(tangent_vec_ptr, concavity_ptr, &mut cp_result);

    let mut norm_tangent_sq = 0.0;
    for d_idx in 0..dimension {
        let val = get_node_val(tangent_vec_ptr, dimension, 0, d_idx);
        norm_tangent_sq += val * val;
    }
    
    if norm_tangent_sq.abs() < 1e-16 { 
        *curvature = 0.0; 
    } else {
        *curvature = cp_result / (norm_tangent_sq * norm_tangent_sq.sqrt());
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_reduce_pseudo_inverse(
    num_nodes_c: c_int,
    dimension_c: c_int,
    nodes: *const c_double,
    reduced: *mut c_double, 
    not_implemented: *mut u8, 
) {
    let num_nodes = num_nodes_c as usize;
    let dimension = dimension_c as usize;
    *not_implemented = 0u8; 

    if num_nodes == 2 { 
        for d_idx in 0..dimension {
            let val = 0.5 * (get_node_val(nodes, dimension, 0, d_idx) + get_node_val(nodes, dimension, 1, d_idx));
            set_node_val(reduced, dimension, 0, d_idx, val);
        }
    } else if num_nodes == 3 { 
        for d_idx in 0..dimension {
            let n0 = get_node_val(nodes, dimension, 0, d_idx);
            let n1 = get_node_val(nodes, dimension, 1, d_idx);
            let n2 = get_node_val(nodes, dimension, 2, d_idx);
            set_node_val(reduced, dimension, 0, d_idx, (5.0 * n0 + 2.0 * n1 - n2) / 6.0);
            set_node_val(reduced, dimension, 1, d_idx, (-n0 + 2.0 * n1 + 5.0 * n2) / 6.0);
        }
    } else if num_nodes == 4 { 
         for d_idx in 0..dimension {
            let n0 = get_node_val(nodes, dimension, 0, d_idx);
            let n1 = get_node_val(nodes, dimension, 1, d_idx);
            let n2 = get_node_val(nodes, dimension, 2, d_idx);
            let n3 = get_node_val(nodes, dimension, 3, d_idx);
            set_node_val(reduced, dimension, 0, d_idx, (19.0 * n0 + 3.0 * n1 - 3.0 * n2 + n3) / 20.0);
            set_node_val(reduced, dimension, 1, d_idx, (-n0 + 3.0 * n1 + 3.0 * n2 - n3) / 4.0);
            set_node_val(reduced, dimension, 2, d_idx, (n0 - 3.0 * n1 + 3.0 * n2 + 19.0 * n3) / 20.0);
        }
    } else if num_nodes == 5 { 
        for d_idx in 0..dimension {
            let n0 = get_node_val(nodes, dimension, 0, d_idx);
            let n1 = get_node_val(nodes, dimension, 1, d_idx);
            let n2 = get_node_val(nodes, dimension, 2, d_idx);
            let n3 = get_node_val(nodes, dimension, 3, d_idx);
            let n4 = get_node_val(nodes, dimension, 4, d_idx);
            set_node_val(reduced, dimension, 0, d_idx, (69.0*n0 +  4.0*n1 -  6.0*n2 +  4.0*n3 - n4) / 70.0);
            set_node_val(reduced, dimension, 1, d_idx, (-53.0*n0 + 212.0*n1 + 102.0*n2 - 68.0*n3 + 17.0*n4) / 210.0);
            set_node_val(reduced, dimension, 2, d_idx, (17.0*n0 -  68.0*n1 + 102.0*n2 + 212.0*n3 - 53.0*n4) / 210.0);
            set_node_val(reduced, dimension, 3, d_idx, (-n0 +  4.0*n1 -  6.0*n2 +  4.0*n3 + 69.0*n4) / 70.0);
        }
    } else {
        *not_implemented = 1u8;
    }
}

unsafe fn projection_error_internal(
    num_nodes: usize,
    dimension: usize,
    nodes_ptr: *const c_double,
    projected_ptr: *const c_double, 
    error_val: &mut f64,
) {
    let mut diff_sq_sum = 0.0;
    let mut nodes_sq_sum = 0.0;

    for node_idx in 0..num_nodes {
        for dim_idx in 0..dimension {
            let n_val = get_node_val(nodes_ptr, dimension, node_idx, dim_idx);
            let p_val = get_node_val(projected_ptr, dimension, node_idx, dim_idx);
            let diff = n_val - p_val;
            diff_sq_sum += diff * diff;
            nodes_sq_sum += n_val * n_val;
        }
    }

    if nodes_sq_sum == 0.0 {
        *error_val = if diff_sq_sum == 0.0 { 0.0 } else { 1.0 }; 
        return;
    }
    
    *error_val = (diff_sq_sum / nodes_sq_sum).sqrt();
}

unsafe fn can_reduce_internal(
    num_nodes: usize, 
    dimension: usize,
    nodes_ptr: *const c_double, 
    projected_work_ptr: *mut c_double, 
) -> i32 { 
    if num_nodes < 2 || num_nodes > 5 { 
        return -1; 
    }

    let reduced_degree_num_nodes = num_nodes - 1; 
    
    let mut pb_nodes_vec: Vec<f64> = vec![0.0; dimension * reduced_degree_num_nodes.max(1)];
    let pb_nodes_ptr = pb_nodes_vec.as_mut_ptr();
    let mut not_implemented_flag: u8 = 0;

    BEZ_reduce_pseudo_inverse(
        num_nodes as c_int, 
        dimension as c_int, 
        nodes_ptr, 
        pb_nodes_ptr, 
        &mut not_implemented_flag
    );
    if not_implemented_flag != 0 { return -1; } 

    if reduced_degree_num_nodes == 0 { 
        if num_nodes == 1 { 
            for d_idx in 0..dimension { 
                 set_node_val(projected_work_ptr, dimension, 0, d_idx, get_node_val(pb_nodes_ptr, dimension, 0, d_idx));
            }
        } else { 
            for d_idx in 0..dimension {
                let val = get_node_val(pb_nodes_ptr, dimension, 0, d_idx);
                set_node_val(projected_work_ptr, dimension, 0, d_idx, val);
                if num_nodes > 1 { 
                    set_node_val(projected_work_ptr, dimension, 1, d_idx, val);
                }
            }
        }
    } else {
         BEZ_elevate_nodes_curve(
            reduced_degree_num_nodes as c_int, 
            dimension as c_int,
            pb_nodes_ptr,
            projected_work_ptr 
        );
    }

    let mut relative_err = 0.0;
    projection_error_internal(num_nodes, dimension, nodes_ptr, projected_work_ptr, &mut relative_err);

    if relative_err < REDUCE_THRESHOLD {
        1 
    } else {
        0 
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_full_reduce(
    num_nodes_c: c_int,
    dimension_c: c_int,
    nodes: *const c_double,
    num_reduced_nodes_c: *mut c_int,
    reduced_nodes_output: *mut c_double, 
    not_implemented: *mut u8,
) {
    let original_num_nodes = num_nodes_c as usize;
    let dimension = dimension_c as usize;

    if original_num_nodes == 0 {
        *num_reduced_nodes_c = 0;
        *not_implemented = 0u8;
        return;
    }

    ptr::copy_nonoverlapping(nodes, reduced_nodes_output, original_num_nodes * dimension);
    *num_reduced_nodes_c = num_nodes_c;
    *not_implemented = 0u8;

    for _ in 0..(original_num_nodes - 1) { 
        let current_num_nodes_val = *num_reduced_nodes_c as usize;
        if current_num_nodes_val <= 1 { break; } 

        let mut projected_work_vec: Vec<f64> = vec![0.0; dimension * current_num_nodes_val];

        let cr_success = can_reduce_internal(
            current_num_nodes_val,
            dimension,
            reduced_nodes_output, 
            projected_work_vec.as_mut_ptr(),
        );

        if cr_success == 1 { 
            let mut work_for_reduce_vec: Vec<f64> = vec![0.0; dimension * (current_num_nodes_val - 1)];
            let mut reduce_not_impl_flag: u8 = 0;

            BEZ_reduce_pseudo_inverse(
                current_num_nodes_val as c_int,
                dimension_c,
                reduced_nodes_output, 
                work_for_reduce_vec.as_mut_ptr(), 
                &mut reduce_not_impl_flag,
            );

            if reduce_not_impl_flag != 0 {
                *not_implemented = 1u8;
                return;
            }

            *num_reduced_nodes_c -= 1;
            ptr::copy_nonoverlapping(
                work_for_reduce_vec.as_ptr(),
                reduced_nodes_output,
                (*num_reduced_nodes_c as usize) * dimension,
            );

        } else if cr_success == 0 { 
            return;
        } else { 
            *not_implemented = 1u8;
            return;
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_compute_length(
    num_nodes_c: c_int,
    dimension_c: c_int,
    nodes: *const c_double,
    length: *mut c_double,
    error_val: *mut c_int, 
) {
    let num_nodes = num_nodes_c as usize;
    let dimension = dimension_c as usize;

    if num_nodes == 0 {
        *length = 0.0;
        *error_val = status::STATUS_UNKNOWN; 
        return;
    }
    if num_nodes == 1 { 
        *length = 0.0;
        *error_val = status::STATUS_SUCCESS;
        return;
    }

    let hodograph_num_nodes = num_nodes - 1;
    let mut hodograph_nodes_vec: Vec<f64> = vec![0.0; dimension * hodograph_num_nodes];
    
    let degree_factor = (num_nodes - 1) as f64;
    for node_idx in 0..hodograph_num_nodes {
        for d_idx in 0..dimension {
            let diff = get_node_val(nodes, dimension, node_idx + 1, d_idx) -
                       get_node_val(nodes, dimension, node_idx, d_idx);
            set_node_val(hodograph_nodes_vec.as_mut_ptr(), dimension, node_idx, d_idx, degree_factor * diff);
        }
    }

    if num_nodes == 2 { 
        let mut norm_sq = 0.0;
        for d_idx in 0..dimension {
            // For num_nodes=2, hodograph_num_nodes=1. The only node is at index 0.
            let val = get_node_val(hodograph_nodes_vec.as_ptr(), dimension, 0, d_idx);
            norm_sq += val * val;
        }
        *length = norm_sq.sqrt();
        *error_val = status::STATUS_SUCCESS;
        return;
    }

    // Clone data for the closure.
    // These are the control points for B'(s), which is what we need to integrate.
    let integrand_nodes_data = hodograph_nodes_vec; 
    let integrand_num_nodes = hodograph_num_nodes;
    let integrand_dim = dimension;

    let vec_size_closure = move |s_val: f64| -> f64 {
        let mut evaluated_hodograph_vec: Vec<f64> = vec![0.0; integrand_dim];
        let s_arr = [s_val];
        
        unsafe { // BEZ_evaluate_multi is unsafe
            BEZ_evaluate_multi(
                integrand_num_nodes as c_int,
                integrand_dim as c_int,
                integrand_nodes_data.as_ptr(),
                1, 
                s_arr.as_ptr(),
                evaluated_hodograph_vec.as_mut_ptr(),
            );
        }

        let mut norm_sq = 0.0;
        for d_idx in 0..integrand_dim {
            norm_sq += evaluated_hodograph_vec[d_idx] * evaluated_hodograph_vec[d_idx];
        }
        norm_sq.sqrt()
    };

    let mut abserr = 0.0;
    let limit_dqagse = 50; 

    rust_dqagse(
        vec_size_closure,
        0.0, 1.0,
        SQRT_PREC, SQRT_PREC, 
        limit_dqagse,
        length, 
        &mut abserr, 
        error_val, 
    );
}

// Note: `curves_equal` and `subdivide_curve` from Fortran are not C-bound
// and operate on the Fortran `CurveData` type. They are not directly
// translated here as they are not part of the FFI interface.
// If their logic is needed by other FFI functions, equivalent internal
// Rust helpers would be created. For now, they appear to be for Fortran-side
// usage or tests.
