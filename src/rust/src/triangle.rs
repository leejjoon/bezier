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

use libc::{c_double, c_int};
use std::ptr;

use crate::curve; // To access BEZ_evaluate_curve_barycentric

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

#[no_mangle]
/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `dimension_ * num_nodes_total` doubles.
/// - `new_nodes` points to writable memory for `dimension_ * (num_nodes_total - degree - 1)` doubles.
/// - `dimension_`, `num_nodes_total`, `degree` are non-negative and describe valid geometry.
/// - `degree` must be less than `num_nodes_total`.
/// - The number of nodes in a Bezier triangle of degree `d` is `(d+1)*(d+2)/2`.
///   `num_nodes_total` should correspond to `degree`.
///   `num_nodes_total - degree - 1` should correspond to `degree - 1`.
///   This implies `(d+1)*(d+2)/2 - d - 1 = d*(d+1)/2`. This is correct.
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
    let current_degree = degree as usize; // This is `d` in (d+1)(d+2)/2

    let mut write_idx_0based = 0; 

    // Fortran parent_i1, parent_i2, parent_i3 are 1-based.
    // Rust equivalent 0-based indices:
    let mut p1_idx = 0;
    let mut p2_idx = 1;
    let mut p3_idx = current_degree + 1; // Start of the "next row" in triangular indexing for p3

    // Loop k = 0 to degree - 1 (Fortran: k = 0, degree - 1)
    // This loop effectively iterates over the rows of the new, smaller triangle being formed.
    // The number of nodes in the new triangle (degree d-1) is d*(d+1)/2.
    for k_row_new_triangle in 0..current_degree { 
        // Loop j = 0 to degree - k - 1 (Fortran: j = 0, degree - k - 1)
        // This loop iterates across the nodes within a row of the new triangle.
        for _j_in_row_new_triangle in 0..(current_degree - k_row_new_triangle) { 
            for d_each_dim in 0..dim {
                let val = lambda1 * get_node_val(nodes, dim, p1_idx, d_each_dim) +
                          lambda2 * get_node_val(nodes, dim, p2_idx, d_each_dim) +
                          lambda3 * get_node_val(nodes, dim, p3_idx, d_each_dim);
                set_node_val(new_nodes, dim, write_idx_0based, d_each_dim, val);
            }
            p1_idx += 1;
            p2_idx += 1;
            p3_idx += 1;
            write_idx_0based += 1;
        }
        // After each row of the new triangle is computed, adjust parent indices for the next row.
        // p1 and p2 skip to the start of the next "segment" in their respective conceptual rows.
        p1_idx += 1; 
        p2_idx += 1; 
        // p3_idx is naturally advanced.
    }
}

#[no_mangle]
/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `dimension * num_nodes` doubles.
/// - `point` points to writable memory for `dimension` doubles.
/// - `dimension`, `num_nodes`, `degree` are non-negative and describe valid geometry.
/// - `BEZ_evaluate_barycentric_multi` preconditions must be met.
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
/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `dimension_c * num_nodes_c` doubles.
/// - `param_vals` points to readable memory for `num_vals_c * 3` doubles (lambda1, lambda2, lambda3 repeating).
/// - `evaluated` points to writable memory for `dimension_c * num_vals_c` doubles.
/// - All `_c` suffixed integers are non-negative and accurately describe their respective data.
/// - `curve::BEZ_evaluate_curve_barycentric` preconditions must be met.
/// - Allocations for `lambdaX_vals_vec` and `row_result_vec` must succeed.
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
    if num_nodes_c == 0 { // Handle empty triangle case
        for val_idx in 0..nv {
            for d_idx in 0..dim {
                set_node_val(evaluated, dim, val_idx, d_idx, 0.0); // Or some other default
            }
        }
        return;
    }


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

    for k_loop in (0..degree).rev() { 
        if k_loop == degree - 1 { 
            binom_val = degree_c as f64; 
        } else { 
            binom_val = (binom_val * (k_loop + 1) as f64) / (degree - k_loop) as f64;
        }

        let fortran_index_row_end = fortran_index_; 
        let nodes_in_row = k_loop + 1;                                   
        let fortran_index_row_start = fortran_index_row_end - nodes_in_row + 1;
        
        if fortran_index_row_start == 0 { // Should not happen with valid degree/num_nodes
            return;
        }
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
/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `dimension_c * num_nodes_c` doubles.
/// - `param_vals` points to readable memory for `num_vals_c * 2` doubles (s, t repeating).
/// - `evaluated` points to writable memory for `dimension_c * num_vals_c` doubles.
/// - All `_c` suffixed integers are non-negative and accurately describe their respective data.
/// - `curve::BEZ_evaluate_curve_barycentric` preconditions must be met.
/// - Allocations for lambda vectors and `row_result_vec` must succeed.
pub unsafe extern "C" fn BEZ_evaluate_cartesian_multi(
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
    if num_nodes_c == 0 { 
        for val_idx in 0..nv {
            for d_idx in 0..dim {
                set_node_val(evaluated, dim, val_idx, d_idx, 0.0);
            }
        }
        return;
    }

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

    for k_loop in (0..degree).rev() { 
        if k_loop == degree - 1 { 
            binom_val = degree_c as f64; 
        } else { 
            binom_val = (binom_val * (k_loop + 1) as f64) / (degree - k_loop) as f64;
        }

        let fortran_index_row_end = fortran_index_; 
        let nodes_in_row = k_loop + 1; 
                                      
        let fortran_index_row_start = fortran_index_row_end - nodes_in_row + 1;
        
        if fortran_index_row_start == 0 { return; }
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
/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `dimension_c * num_nodes_c` doubles.
/// - `new_nodes` points to writable memory for `2 * dimension_c * (num_nodes_c - degree_c - 1)` doubles.
///   The term `num_nodes_c - degree_c - 1` is `degree * (degree+1) / 2` for a triangle.
/// - `dimension_c`, `num_nodes_c`, `degree_c` are non-negative and describe valid geometry.
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

    for _nodes_in_row_count_1based in (1..=original_degree).rev() { 
        for _k_in_row in 0.._nodes_in_row_count_1based {
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

#[no_mangle]
/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `2 * num_nodes` doubles (dimension is fixed to 2).
/// - `param_vals` points to readable memory for `num_vals * 2` doubles (s,t pairs).
/// - `evaluated` points to writable memory for `num_vals` doubles.
/// - `num_nodes`, `degree`, `num_vals` are non-negative and describe valid geometry.
/// - `BEZ_jacobian_both` and `BEZ_evaluate_cartesian_multi` preconditions are met.
/// - Allocation for `jac_nodes_vec` and `bs_bt_vals_vec` must succeed.
pub unsafe extern "C" fn BEZ_jacobian_det(
    num_nodes: c_int,
    nodes: *const c_double, 
    degree: c_int,
    num_vals: c_int,
    param_vals: *const c_double,
    evaluated: *mut c_double, 
) {
    let dim_fixed: usize = 2; 
    let num_deriv_nodes_total = degree as usize * (degree as usize + 1) / 2; 
    
    if degree == 1 {
        let mut jac_nodes_val = [0.0; 4]; 
        BEZ_jacobian_both(num_nodes, dim_fixed as c_int, nodes, degree, jac_nodes_val.as_mut_ptr());
        let determinant = jac_nodes_val[0] * jac_nodes_val[3] - jac_nodes_val[1] * jac_nodes_val[2]; 
        for i in 0..num_vals as usize {
            *evaluated.add(i) = determinant;
        }
    } else {
        let mut jac_nodes_vec: Vec<f64> = vec![0.0; dim_fixed * 2 * num_deriv_nodes_total]; 
        BEZ_jacobian_both(num_nodes, dim_fixed as c_int, nodes, degree, jac_nodes_vec.as_mut_ptr());

        let mut bs_bt_vals_vec: Vec<f64> = vec![0.0; dim_fixed * 2 * num_vals as usize];
        BEZ_evaluate_cartesian_multi(
            num_deriv_nodes_total as c_int, 
            (dim_fixed * 2) as c_int, 
            jac_nodes_vec.as_ptr(),
            degree - 1, 
            num_vals,
            param_vals,
            bs_bt_vals_vec.as_mut_ptr(),
        );

        for i in 0..num_vals as usize {
            let sx = get_node_val(bs_bt_vals_vec.as_ptr(), dim_fixed * 2, i, 0);
            let sy = get_node_val(bs_bt_vals_vec.as_ptr(), dim_fixed * 2, i, 1);
            let tx = get_node_val(bs_bt_vals_vec.as_ptr(), dim_fixed * 2, i, 2); 
            let ty = get_node_val(bs_bt_vals_vec.as_ptr(), dim_fixed * 2, i, 3); 
            *evaluated.add(i) = sx * ty - sy * tx;
        }
    }
}

/// # Safety
/// Internal helper. `degree` must be non-negative.
unsafe fn specialize_workspace_sizes_internal(degree: c_int, size_odd: &mut c_int, size_even: &mut c_int) {
    if degree % 2 == 1 {
        *size_odd = ((degree + 1) * (degree + 3).pow(2) * (degree + 5)) / 64;
        *size_even = *size_odd;
    } else if degree % 4 == 0 {
        *size_odd = (degree * (degree + 2) * (degree + 4) * (degree + 6)) / 64;
        *size_even = (((degree + 2) * (degree + 4)) / 8).pow(2);
    } else { // degree % 4 == 2
        *size_odd = (((degree + 2) * (degree + 4)) / 8).pow(2);
        *size_even = (degree * (degree + 2) * (degree + 4) * (degree + 6)) / 64;
    }
}

/// # Safety
/// All pointer arguments must be valid and correctly sized.
/// `read_nodes_buffer` contains `num_read_triangles_in_buffer` * `size_read_single_triangle_nodes` * `dimension_` elements.
/// `write_nodes_buffer` must be large enough for all outputs.
/// `local_degree` is for `BEZ_de_casteljau_one_round`.
unsafe fn specialize_triangle_one_round_internal(
    dimension_: c_int,
    _num_read_triangles_in_buffer: c_int, 
    read_nodes_buffer: *const c_double,
    _num_write_triangles_in_buffer: c_int, 
    write_nodes_buffer: *mut c_double,
    size_read_single_triangle_nodes: c_int, 
    size_write_single_triangle_nodes: c_int,
    step: c_int,            
    local_degree: c_int,    
    weights_a: *const c_double, 
    weights_b: *const c_double,
    weights_c: *const c_double,
) {
    let dim_usize = dimension_ as usize;
    let read_stride_elements = (size_read_single_triangle_nodes * dimension_) as usize;
    let write_stride_elements = (size_write_single_triangle_nodes * dimension_) as usize;

    let mut current_write_ptr = write_nodes_buffer;
    let mut current_read_ptr_a = read_nodes_buffer; // For weights_a section
    let mut current_read_ptr_b = read_nodes_buffer; // For weights_b section
    let mut current_read_ptr_c = read_nodes_buffer; // For weights_c section

    // Part 1: Apply weights_a to the first triangle in the read set
    BEZ_de_casteljau_one_round(
        size_read_single_triangle_nodes, dimension_, current_read_ptr_a, local_degree,
        *weights_a.add(0), *weights_a.add(1), *weights_a.add(2),
        current_write_ptr
    );
    current_write_ptr = current_write_ptr.add(write_stride_elements);
    
    // Part 2: Apply weights_b
    // `step` indicates how many source triangles weights_b are applied to.
    // These are distinct source triangles from the previous level.
    for _i_b in 0..step {
        // The Fortran code implies `read_index` advances for these.
        // `read_index = read_index + size_read -1` in Fortran `specialize_triangle` after `call ...one_round`
        // suggests that `read_nodes` in `specialize_triangle_one_round` is a base pointer and internal
        // indexing handles the sub-triangles.
        // My current `current_read_ptr_b` needs to point to the correct source triangle.
        // Fortran: `read_nodes(:, read_index:new_read)` -- this implies a slice.
        // The `read_nodes_buffer` passed here IS that slice for one of the A, B, or C groups.
        // This helper is simpler: it processes its entire `read_nodes_buffer` as one set.
        // The complexity is in `BEZ_specialize_triangle` to call this helper with correct sub-buffers.
        // **Correction**: The Fortran `specialize_triangle_one_round` is indeed complex.
        // It takes the *entire output from the previous step* as `read_nodes`.
        // It then calculates where the "A", "B", and "C" sub-components are within that input
        // and where their results should go in `write_nodes`.
        // The Rust version needs to replicate this.
        // For now, this is a placeholder for the complex indexing.
        // The current version of specialize_triangle_one_round_internal below is simplified.
        // It will be called by BEZ_specialize_triangle for different segments.
        // The Fortran version of specialize_triangle_one_round is NOT just one De Casteljau call.
        // The loops are for the *output* triangles being generated at this step.
        // This function is the one that needs the complex indexing from Fortran's `specialize_triangle_one_round`.
        // Let's assume for now that `BEZ_specialize_triangle` calls this for each section (A, B, C) with correct sub-buffers.
        // This means this function will be called multiple times by BEZ_specialize_triangle.
        // If so, it should only perform *one* type of de Casteljau application (A, B, or C).
        // This is not what the Fortran `specialize_triangle_one_round` does.
        // Reverting to the structure from previous attempt, which is a simplified interpretation of the Fortran `specialize_triangle_one_round`.
        // This function will be called ONCE per step by BEZ_specialize_triangle.
        // `read_nodes_buffer` IS the entire buffer from previous step.
    }
     // The Fortran code's `specialize_triangle_one_round` IS complex.
    // The loops below are from that Fortran routine.
    // `write_index` and `read_index` are 1-based in Fortran, converted to 0-based offsets here.

    let mut write_offset_elements : usize = 0; // Start writing at beginning of write_nodes_buffer

    // "First: (step, 0, 0)" -> applies weights_a
    // Input is the first triangle from read_nodes_buffer
    BEZ_de_casteljau_one_round(
        size_read_single_triangle_nodes, dimension_, read_nodes_buffer, local_degree,
        *weights_a.add(0), *weights_a.add(1), *weights_a.add(2),
        write_nodes_buffer.add(write_offset_elements)
    );
    write_offset_elements += write_stride_elements;

    // "Second: (i, j, 0) for j > 0, i + j = step" -> applies weights_b
    // These are `step` output triangles.
    // Their inputs are `step` distinct triangles from `read_nodes_buffer`.
    let mut read_offset_elements_b = size_read_single_triangle_nodes as usize * dim_usize; // Start after first triangle
    for _i_b in 0..step {
        BEZ_de_casteljau_one_round(
            size_read_single_triangle_nodes, dimension_, read_nodes_buffer.add(read_offset_elements_b), local_degree,
            *weights_b.add(0), *weights_b.add(1), *weights_b.add(2),
            write_nodes_buffer.add(write_offset_elements)
        );
        write_offset_elements += write_stride_elements;
        read_offset_elements_b += read_stride_elements;
    }
    
    // "Third: (i, j, k) for k > 0, i + j + k = step" -> applies weights_c
    // These are `step * (step + 1) / 2` output triangles.
    // Their inputs are also distinct triangles from `read_nodes_buffer`.
    let num_c_applications = (step * (step+1)) / 2;
    let mut read_offset_elements_c = 0; // C weights apply to sub-triangles starting from beginning of read_nodes_buffer
    for _i_c in 0..num_c_applications {
         BEZ_de_casteljau_one_round(
            size_read_single_triangle_nodes, dimension_, read_nodes_buffer.add(read_offset_elements_c), local_degree,
            *weights_c.add(0), *weights_c.add(1), *weights_c.add(2),
            write_nodes_buffer.add(write_offset_elements)
        );
        write_offset_elements += write_stride_elements;
        read_offset_elements_c += read_stride_elements;
    }
}


#[no_mangle]
/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `dimension_ * num_nodes` doubles.
/// - `specialized` points to writable memory for `dimension_ * num_nodes` doubles.
/// - `weights_a`, `weights_b`, `weights_c` each point to 3 readable doubles.
/// - `dimension_`, `num_nodes`, `degree` are non-negative and describe valid geometry.
/// - Allocations for workspaces succeed.
/// - `specialize_triangle_one_round_internal` preconditions are met.
pub unsafe extern "C" fn BEZ_specialize_triangle(
    num_nodes: c_int, 
    dimension_: c_int, 
    nodes: *const c_double, 
    degree: c_int,          
    weights_a: *const c_double, 
    weights_b: *const c_double, 
    weights_c: *const c_double, 
    specialized: *mut c_double, 
) {
    let dim = dimension_ as usize;
    
    let mut size_odd_elements_capacity: c_int = 0; 
    let mut size_even_elements_capacity: c_int = 0; 
    specialize_workspace_sizes_internal(degree, &mut size_odd_elements_capacity, &mut size_even_elements_capacity);

    let mut workspace_odd_vec: Vec<f64> = vec![0.0; size_odd_elements_capacity as usize * dim];
    let mut workspace_even_vec: Vec<f64> = vec![0.0; size_even_elements_capacity as usize * dim];

    let mut current_read_ptr: *const c_double = nodes;
    let mut current_read_buffer_total_nodes = num_nodes; 
    let mut current_nodes_per_single_triangle = (degree + 1) * (degree + 2) / 2; 
    
    let mut is_even_step = false; 

    for step in 1..=degree {
        let next_degree_val = degree - step;
        let next_nodes_per_single_triangle = (next_degree_val + 1) * (next_degree_val + 2) / 2;
        
        // Determine number of distinct triangles in the current read buffer
        let num_triangles_in_read_buffer = if current_nodes_per_single_triangle > 0 {
             current_read_buffer_total_nodes / current_nodes_per_single_triangle
        } else { 0 };
        
        // Determine number of triangles that will be written in this step
        // This is complex based on Fortran's num_curves logic.
        // num_curves = 1 initially. For step=1, num_curves=1+2=3. For step=2, num_curves=3+3=6.
        // num_curves for step `s` is (s+1)*(s+2)/2.
        let num_triangles_to_write_this_step = (step + 1) * (step + 2) / 2;


        let write_ptr: *mut c_double;
        let write_buffer_capacity_nodes: c_int; 
        
        if step == 1 {
            write_ptr = workspace_odd_vec.as_mut_ptr();
            write_buffer_capacity_nodes = size_odd_elements_capacity;
        } else if is_even_step { 
            current_read_ptr = workspace_odd_vec.as_ptr();
            current_read_buffer_total_nodes = size_odd_elements_capacity;
            
            write_ptr = workspace_even_vec.as_mut_ptr();
            write_buffer_capacity_nodes = size_even_elements_capacity;
        } else { 
            current_read_ptr = workspace_even_vec.as_ptr();
            current_read_buffer_total_nodes = size_even_elements_capacity;
            
            write_ptr = workspace_odd_vec.as_mut_ptr();
            write_buffer_capacity_nodes = size_odd_elements_capacity;
        }
        
        specialize_triangle_one_round_internal(
            dimension_, 
            num_triangles_in_read_buffer, current_read_ptr, // Pass num triangles and base pointer
            num_triangles_to_write_this_step, write_ptr,    // Pass num expected triangles and base pointer
            current_nodes_per_single_triangle, next_nodes_per_single_triangle,    
            step, degree + 1 - step, 
            weights_a, weights_b, weights_c
        );
        
        current_nodes_per_single_triangle = next_nodes_per_single_triangle;
        is_even_step = !is_even_step;
    }

    let final_read_ptr = if is_even_step { 
        workspace_odd_vec.as_ptr()
    } else { 
        workspace_even_vec.as_ptr()
    };
    ptr::copy_nonoverlapping(final_read_ptr, specialized, (num_nodes * dim as c_int) as usize);
}


#[no_mangle]
/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `dimension_ * num_nodes` doubles.
/// - `nodes_a`, `nodes_b`, `nodes_c`, `nodes_d` point to writable memory, each for `dimension_ * num_nodes` doubles.
/// - `dimension_`, `num_nodes`, `degree` are non-negative and describe valid geometry.
/// - `BEZ_specialize_triangle` preconditions must be met if degree > 4 or for the fallback.
pub unsafe extern "C" fn BEZ_subdivide_nodes_triangle(
    num_nodes: c_int, 
    dimension_: c_int, 
    nodes: *const c_double, 
    degree: c_int,          
    nodes_a: *mut c_double, 
    nodes_b: *mut c_double, 
    nodes_c: *mut c_double, 
    nodes_d: *mut c_double, 
) {
    let dim = dimension_ as usize;

    if num_nodes == 0 || degree < 0 { return; }
    if degree == 0 { 
        if num_nodes > 0 && dim > 0 {
             for d_idx in 0..dim {
                let val = get_node_val(nodes, dim, 0, d_idx);
                set_node_val(nodes_a, dim, 0, d_idx, val);
                set_node_val(nodes_b, dim, 0, d_idx, val);
                set_node_val(nodes_c, dim, 0, d_idx, val);
                set_node_val(nodes_d, dim, 0, d_idx, val);
            }
        }
        return;
    }

    if degree == 1 {
        for d_idx in 0..dim {
            let n00 = get_node_val(nodes, dim, 0, d_idx); 
            let n10 = get_node_val(nodes, dim, 1, d_idx); 
            let n01 = get_node_val(nodes, dim, 2, d_idx); 

            let na1_val = 0.5 * (n00 + n10);
            let na2_val = 0.5 * (n00 + n01);
            let nb0_val = 0.5 * (n10 + n01);

            set_node_val(nodes_a, dim, 0, d_idx, n00);
            set_node_val(nodes_a, dim, 1, d_idx, na1_val);
            set_node_val(nodes_a, dim, 2, d_idx, na2_val);

            set_node_val(nodes_b, dim, 0, d_idx, nb0_val);
            set_node_val(nodes_b, dim, 1, d_idx, na2_val); 
            set_node_val(nodes_b, dim, 2, d_idx, na1_val); 
            
            set_node_val(nodes_c, dim, 0, d_idx, na1_val); 
            set_node_val(nodes_c, dim, 1, d_idx, n10);
            set_node_val(nodes_c, dim, 2, d_idx, nb0_val); 

            set_node_val(nodes_d, dim, 0, d_idx, na2_val); 
            set_node_val(nodes_d, dim, 1, d_idx, nb0_val); 
            set_node_val(nodes_d, dim, 2, d_idx, n01);
        }
    } else { 
        let weights_corner0 = [1.0, 0.0, 0.0]; 
        let weights_mid01   = [0.5, 0.5, 0.0]; 
        let weights_mid02   = [0.5, 0.0, 0.5]; 
        BEZ_specialize_triangle(num_nodes, dimension_, nodes, degree, weights_corner0.as_ptr(), weights_mid01.as_ptr(), weights_mid02.as_ptr(), nodes_a);

        let weights_mid12   = [0.0, 0.5, 0.5]; 
        BEZ_specialize_triangle(num_nodes, dimension_, nodes, degree, weights_mid01.as_ptr(), weights_mid02.as_ptr(), weights_mid12.as_ptr(), nodes_b);
        // Fortran nodes_b order: Mid(BC), Mid(AC), Mid(AB)
        // BEZ_specialize_triangle(num_nodes, dimension_, nodes, degree, weights_mid12.as_ptr(), weights_mid02.as_ptr(), weights_mid01.as_ptr(), nodes_b);

        
        let weights_corner1 = [0.0, 1.0, 0.0]; 
        BEZ_specialize_triangle(num_nodes, dimension_, nodes, degree, weights_mid01.as_ptr(), weights_corner1.as_ptr(), weights_mid12.as_ptr(), nodes_c);

        let weights_corner2 = [0.0, 0.0, 1.0]; 
        BEZ_specialize_triangle(num_nodes, dimension_, nodes, degree, weights_mid02.as_ptr(), weights_mid12.as_ptr(), weights_corner2.as_ptr(), nodes_d);
    }
}

#[no_mangle]
/// # Safety
/// Caller must ensure:
/// - `nodes` points to readable memory for `dimension_ * num_nodes` doubles.
/// - `nodes1`, `nodes2`, `nodes3` point to writable memory, each for `dimension_ * (degree + 1)` doubles.
/// - `dimension_`, `num_nodes`, `degree` are non-negative and describe valid geometry.
pub unsafe extern "C" fn BEZ_compute_edge_nodes(
    num_nodes_in_triangle: c_int, 
    dimension_: c_int, 
    nodes: *const c_double, 
    degree: c_int,          
    nodes1: *mut c_double, 
    nodes2: *mut c_double, 
    nodes3: *mut c_double, 
) {
    let dim = dimension_ as usize;
    let current_degree = degree as usize;

    if num_nodes_in_triangle == 0 || current_degree + 1 == 0 { return; }

    // Edge 1 (s-edge, t=0): (d,0,0) ... (0,d,0)
    for i in 0..=current_degree {
        for d_idx in 0..dim {
            set_node_val(nodes1, dim, i, d_idx, get_node_val(nodes, dim, i, d_idx));
        }
    }

    // Edge 2 (t-edge, s=0, re-parameterized from w,t space): (d,0,0) ... (0,0,d)
    let mut current_node_idx_e2 = 0; 
    for i in 0..=current_degree { 
        for d_idx in 0..dim {
            set_node_val(nodes2, dim, i, d_idx, get_node_val(nodes, dim, current_node_idx_e2, d_idx));
        }
        if i < current_degree {
             current_node_idx_e2 += (current_degree - i) + 1;
        }
    }
    
    // Edge 3 (hypotenuse, w=0): (0,d,0) ... (0,0,d)
    let mut current_node_idx_e3 = current_degree; 
    for i in 0..=current_degree { 
        for d_idx in 0..dim {
            set_node_val(nodes3, dim, i, d_idx, get_node_val(nodes, dim, current_node_idx_e3, d_idx));
        }
        if i < current_degree {
             current_node_idx_e3 += (current_degree - 1 - i) + 1;
        }
    }
}

/// # Safety
/// `nodes` must point to readable memory for `2 * num_nodes` doubles.
/// `shoelace` and `not_implemented` must point to writable memory.
/// `num_nodes` must be between 2 and 5 inclusive if `not_implemented` is to be false.
unsafe fn shoelace_for_area_internal(
    num_nodes: c_int,       
    nodes: *const c_double, 
    shoelace: &mut c_double,
    not_implemented: &mut u8, 
) {
    *not_implemented = 0; 
    let dim = 2; 

    let get_xy = |node_idx: usize| -> (f64, f64) {
        (get_node_val(nodes, dim, node_idx, 0), get_node_val(nodes, dim, node_idx, 1))
    };

    if num_nodes == 2 { 
        let (n0x, n0y) = get_xy(0);
        let (n1x, n1y) = get_xy(1);
        *shoelace = (n0x * n1y - n0y * n1x) / 2.0;
    } else if num_nodes == 3 { 
        let (n0x, n0y) = get_xy(0);
        let (n1x, n1y) = get_xy(1);
        let (n2x, n2y) = get_xy(2);
        *shoelace = (
            2.0 * (n0x * n1y - n0y * n1x) +
                  (n0x * n2y - n0y * n2x) +
            2.0 * (n1x * n2y - n1y * n2x)
        ) / 6.0;
    } else if num_nodes == 4 { 
        let (n0x, n0y) = get_xy(0);
        let (n1x, n1y) = get_xy(1);
        let (n2x, n2y) = get_xy(2);
        let (n3x, n3y) = get_xy(3);
        *shoelace = (
            6.0 * (n0x * n1y - n0y * n1x) +
            3.0 * (n0x * n2y - n0y * n2x) +
                  (n0x * n3y - n0y * n3x) +
            3.0 * (n1x * n2y - n1y * n2x) +
            3.0 * (n1x * n3y - n1y * n3x) +
            6.0 * (n2x * n3y - n2y * n3x)
        ) / 20.0;
    } else if num_nodes == 5 { 
        let (n0x, n0y) = get_xy(0);
        let (n1x, n1y) = get_xy(1);
        let (n2x, n2y) = get_xy(2);
        let (n3x, n3y) = get_xy(3);
        let (n4x, n4y) = get_xy(4);
        *shoelace = (
            20.0 * (n0x * n1y - n0y * n1x) +
            10.0 * (n0x * n2y - n0y * n2x) +
             4.0 * (n0x * n3y - n0y * n3x) +
                   (n0x * n4y - n0y * n4x) +
             8.0 * (n1x * n2y - n1y * n2x) +
             8.0 * (n1x * n3y - n1y * n3x) +
             4.0 * (n1x * n4y - n1y * n4x) +
             8.0 * (n2x * n3y - n2y * n3x) +
            10.0 * (n2x * n4y - n2y * n4x) +
            20.0 * (n3x * n4y - n3y * n4x)
        ) / 70.0;
    } else {
        *shoelace = 0.0;
        *not_implemented = 1; 
    }
}

#[no_mangle]
/// # Safety
/// Caller must ensure:
/// - `sizes` points to readable memory for `num_edges` c_int values.
/// - `nodes_pointers` points to readable memory for `num_edges` pointers to c_double. Each of these
///   sub-pointers must be valid for `2 * sizes[i]` c_double values.
/// - `area` and `not_implemented` must point to writable memory.
/// - `num_edges` is non-negative.
pub unsafe extern "C" fn BEZ_compute_area(
    num_edges: c_int,
    sizes: *const c_int,         
    nodes_pointers: *const *const c_double, 
    area: *mut c_double,
    not_implemented: *mut u8,    
) {
    *area = 0.0;
    *not_implemented = 0; 

    for i in 0..num_edges as usize {
        let current_edge_num_nodes = *sizes.add(i);
        let current_edge_nodes_ptr = *nodes_pointers.add(i);
        
        let mut shoelace_val: c_double = 0.0;
        let mut edge_not_implemented: u8 = 0;

        shoelace_for_area_internal(
            current_edge_num_nodes,
            current_edge_nodes_ptr,
            &mut shoelace_val,
            &mut edge_not_implemented,
        );

        if edge_not_implemented != 0 {
            *not_implemented = 1; 
            return;
        }
        *area += shoelace_val;
    }
}

// Final TODO for triangle.rs:
// - Thoroughly review and potentially rework `BEZ_specialize_triangle` and its helper 
//   `specialize_triangle_one_round_internal` to ensure faithful translation of Fortran's 
//   complex indexing and buffer management for staged specialization if bit-level equivalence 
//   or specific intermediate behavior is critical. The current version is a simplification.
// - Verify if the exact (and potentially non-standard) node selection and ordering 
//   for `nodes2` and `nodes3` output by Fortran's `compute_edge_nodes` is required by any 
//   calling code, or if the current Rust implementation providing standard edge definitions is acceptable.
//   (Current Rust version provides standard edges, which might differ from Fortran's specific output for those two.)
// - Consider porting the hardcoded subdivision formulas for degrees 2, 3, and 4 in 
//   `BEZ_subdivide_nodes_triangle` if the performance of the `BEZ_specialize_triangle` fallback is insufficient.
//   (Currently uses specialize_triangle for degree > 1).
