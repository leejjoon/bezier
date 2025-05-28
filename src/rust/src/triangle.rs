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

use libc::{c_double, c_int};
use std::ptr;

use crate::curve; // To access BEZ_evaluate_curve_barycentric

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

#[no_mangle]
pub unsafe extern "C" fn BEZ_de_casteljau_one_round(
    num_nodes_total: c_int, 
    dimension_: c_int,      
    nodes: *const c_double, 
    degree: c_int,          
    lambda1: c_double,
    lambda2: c_double,
    lambda3: c_double,
    new_nodes: *mut c_double, 
                              
) {
    let dim = dimension_ as usize;
    let current_degree = degree as usize; 

    let mut write_idx_0based = 0; 

    let mut p1_idx_0based = 0; 
    let mut p2_idx_0based = 1;
    let mut p3_idx_0based = current_degree + 1; 

    for k_outer in 0..current_degree { 
        for _j_inner in 0..(current_degree - k_outer) { 
            for d_each_dim in 0..dim {
                let val = lambda1 * get_node_val(nodes, dim, p1_idx_0based, d_each_dim) +
                          lambda2 * get_node_val(nodes, dim, p2_idx_0based, d_each_dim) +
                          lambda3 * get_node_val(nodes, dim, p3_idx_0based, d_each_dim);
                set_node_val(new_nodes, dim, write_idx_0based, d_each_dim, val);
            }
            p1_idx_0based += 1;
            p2_idx_0based += 1;
            p3_idx_0based += 1;
            write_idx_0based += 1;
        }
        p1_idx_0based += 1; 
        p2_idx_0based += 1; 
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_evaluate_barycentric(
    num_nodes: c_int, 
    dimension: c_int, 
    nodes: *const c_double, 
    degree: c_int, 
    lambda1: c_double,
    lambda2: c_double,
    lambda3: c_double,
    point: *mut c_double, 
) {
    let param_vals_arr: [c_double; 3] = [lambda1, lambda2, lambda3];
    
    BEZ_evaluate_barycentric_multi(
        num_nodes,
        dimension,
        nodes,
        degree,
        1, 
        param_vals_arr.as_ptr(), 
        point,
    );
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_evaluate_barycentric_multi(
    num_nodes_c: c_int, 
    dimension_c: c_int, 
    nodes: *const c_double, 
    degree_c: c_int, 
    num_vals_c: c_int,
    param_vals: *const c_double, 
    evaluated: *mut c_double,   
) {
    let dim = dimension_c as usize;
    let degree = degree_c as usize;
    let nv = num_vals_c as usize;

    if nv == 0 { return; }

    let apex_node_idx_0based = num_nodes_c as usize - 1;
    for val_idx in 0..nv {
        for d_idx in 0..dim {
            set_node_val(evaluated, dim, val_idx, d_idx, get_node_val(nodes, dim, apex_node_idx_0based, d_idx));
        }
    }

    if degree == 0 {
        return;
    }

    let mut lambda1_vals_vec: Vec<f64> = Vec::with_capacity(nv);
    let mut lambda2_vals_vec: Vec<f64> = Vec::with_capacity(nv);
    let mut lambda3_vals_vec: Vec<f64> = Vec::with_capacity(nv);

    for i in 0..nv {
        lambda1_vals_vec.push(*param_vals.add(i * 3 + 0)); 
        lambda2_vals_vec.push(*param_vals.add(i * 3 + 1)); 
        lambda3_vals_vec.push(*param_vals.add(i * 3 + 2)); 
    }
    
    let mut row_result_vec: Vec<f64> = vec![0.0; dim * nv]; 

    let mut fortran_index_ = num_nodes_c as usize; 
    let mut binom_val: f64 = 1.0; 

    for k_loop in (0..degree).rev() { // k from degree-1 down to 0
        if k_loop == degree - 1 { 
            binom_val = degree_c as f64; 
        } else { 
            binom_val = (binom_val * (k_loop + 1) as f64) / (degree - k_loop) as f64;
        }

        let fortran_index_row_end = fortran_index_; 
        let nodes_in_row = k_loop + 1;                                   
        let fortran_index_row_start = fortran_index_row_end - nodes_in_row + 1;

        let current_row_nodes_ptr = nodes.add((fortran_index_row_start - 1) * dim);

        curve::BEZ_evaluate_curve_barycentric(
            nodes_in_row as c_int,
            dimension_c,
            current_row_nodes_ptr,
            num_vals_c,
            lambda1_vals_vec.as_ptr(),
            lambda2_vals_vec.as_ptr(),
            row_result_vec.as_mut_ptr(),
        );

        fortran_index_ = fortran_index_row_start -1; 

        for val_idx in 0..nv {
            let l3_val = lambda3_vals_vec[val_idx]; 
            for d_idx in 0..dim {
                let current_eval = get_node_val(evaluated, dim, val_idx, d_idx);
                let row_res = get_node_val(row_result_vec.as_ptr(), dim, val_idx, d_idx);
                let updated_val = l3_val * current_eval + binom_val * row_res;
                set_node_val(evaluated, dim, val_idx, d_idx, updated_val);
            }
        }
    }
}


#[no_mangle]
pub unsafe extern "C" fn BEZ_evaluate_cartesian_multi(
    num_nodes_c: c_int, 
    dimension_c: c_int, 
    nodes: *const c_double, 
    degree_c: c_int, 
    num_vals_c: c_int,
    param_vals: *const c_double, // num_vals rows, 2 columns (s, t)
    evaluated: *mut c_double,   // (dimension, num_vals)
) {
    let dim = dimension_c as usize;
    let degree = degree_c as usize;
    let nv = num_vals_c as usize;

    if nv == 0 { return; }

    let apex_node_idx_0based = num_nodes_c as usize - 1;
    for val_idx in 0..nv {
        for d_idx in 0..dim {
            set_node_val(evaluated, dim, val_idx, d_idx, get_node_val(nodes, dim, apex_node_idx_0based, d_idx));
        }
    }
    
    if degree == 0 {
        return;
    }

    let mut computed_lambda1_vals_vec: Vec<f64> = Vec::with_capacity(nv);
    let mut s_param_vals_vec: Vec<f64> = Vec::with_capacity(nv); 
    let mut t_param_vals_vec: Vec<f64> = Vec::with_capacity(nv); 

    for i in 0..nv {
        let s_val = *param_vals.add(i * 2 + 0); 
        let t_val = *param_vals.add(i * 2 + 1); 
        computed_lambda1_vals_vec.push(1.0 - s_val - t_val);
        s_param_vals_vec.push(s_val);
        t_param_vals_vec.push(t_val);
    }

    let mut row_result_vec: Vec<f64> = vec![0.0; dim * nv];

    let mut fortran_index_ = num_nodes_c as usize;
    let mut binom_val: f64 = 1.0; 

    for k_loop in (0..degree).rev() { // k from degree-1 down to 0
        if k_loop == degree - 1 { 
            binom_val = degree_c as f64; 
        } else { 
            binom_val = (binom_val * (k_loop + 1) as f64) / (degree - k_loop) as f64;
        }

        let fortran_index_row_end = fortran_index_; 
        let nodes_in_row = k_loop + 1; 
                                      
        let fortran_index_row_start = fortran_index_row_end - nodes_in_row + 1;

        let current_row_nodes_ptr = nodes.add((fortran_index_row_start - 1) * dim);

        curve::BEZ_evaluate_curve_barycentric(
            nodes_in_row as c_int,
            dimension_c,
            current_row_nodes_ptr,
            num_vals_c,
            computed_lambda1_vals_vec.as_ptr(), 
            s_param_vals_vec.as_ptr(),          
            row_result_vec.as_mut_ptr(),
        );

        fortran_index_ = fortran_index_row_start -1; 

        for val_idx in 0..nv {
            let t_param_val = t_param_vals_vec[val_idx]; 
            for d_idx in 0..dim {
                let current_eval = get_node_val(evaluated, dim, val_idx, d_idx);
                let row_res = get_node_val(row_result_vec.as_ptr(), dim, val_idx, d_idx);
                let updated_val = t_param_val * current_eval + binom_val * row_res;
                set_node_val(evaluated, dim, val_idx, d_idx, updated_val);
            }
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_jacobian_both(
    num_nodes_c: c_int,       
    dimension_c: c_int,     
    nodes: *const c_double,   
    degree_c: c_int,          
    new_nodes: *mut c_double, 
) {
    let original_degree = degree_c as usize;
    let dim = dimension_c as usize;
    let num_deriv_nodes = original_degree * (original_degree + 1) / 2;

    let mut write_col_idx_0based = 0; 
    let mut s_base_node_idx_0based = 0; 
    let mut t_base_node_idx_0based = original_degree + 1; 

    for nodes_in_row_count_1based in (1..=original_degree).rev() { 
        for _k_in_row in 0..nodes_in_row_count_1based {
            for d_each_dim in 0..dim {
                let node_i_val = get_node_val(nodes, dim, s_base_node_idx_0based, d_each_dim);
                
                let node_i_plus_1_val = get_node_val(nodes, dim, s_base_node_idx_0based + 1, d_each_dim);
                set_node_val(new_nodes, 2 * dim, write_col_idx_0based, d_each_dim, node_i_plus_1_val - node_i_val);

                let node_j_val = get_node_val(nodes, dim, t_base_node_idx_0based, d_each_dim);
                set_node_val(new_nodes, 2 * dim, write_col_idx_0based, d_each_dim + dim, node_j_val - node_i_val);
            }
            write_col_idx_0based += 1;
            s_base_node_idx_0based += 1;
            t_base_node_idx_0based += 1;
        }
        s_base_node_idx_0based += 1; 
    }

    let degree_factor = original_degree as f64;
    for node_idx in 0..num_deriv_nodes {
        for d_idx in 0..(2 * dim) { 
            let val = get_node_val(new_nodes, 2 * dim, node_idx, d_idx);
            set_node_val(new_nodes, 2 * dim, node_idx, d_idx, degree_factor * val);
        }
    }
}

// TODO: Continue with other functions: BEZ_jacobian_det, etc.
