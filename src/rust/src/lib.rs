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

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use libc;

pub mod status;
pub mod helpers;
pub mod quadpack_replacement;
pub mod curve;
pub mod triangle;
pub mod triangle_intersection;
pub mod curve_intersection; // Added this line

#[pyfunction]
fn hello_from_rust_py() -> PyResult<String> {
    Ok("Hello from Rust!".to_string())
}

use std::os::raw::c_double;

#[pyfunction]
fn evaluate_curve_multi_py(
    num_nodes: i32,
    dimension: i32,
    nodes: Vec<f64>,
    s_vals: Vec<f64>,
) -> PyResult<Vec<f64>> {
    if num_nodes < 1 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "num_nodes must be at least 1 for evaluate_curve_multi.",
        ));
    }
    if dimension <= 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "dimension must be positive.",
        ));
    }
    if nodes.len() != (num_nodes * dimension) as usize {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Nodes vector length must be num_nodes * dimension.",
        ));
    }

    let num_s_vals = s_vals.len() as i32;
    // If s_vals is empty, num_s_vals will be 0.
    // The output `evaluated` vector will have length 0 if num_s_vals is 0 or dimension is 0.
    // The C function BEZ_evaluate_multi handles num_s_vals = 0 by doing nothing.
    let mut evaluated: Vec<f64> =
        vec![0.0; (num_s_vals * dimension) as usize];

    // The underlying FFI function curve::BEZ_evaluate_multi does not return an error code
    // nor does it take an error_val out-parameter according to its C signature.
    // Thus, we call it and assume success if inputs are valid according to above checks.
    unsafe {
        curve::BEZ_evaluate_multi(
            num_nodes,
            dimension,
            nodes.as_ptr() as *const c_double,
            num_s_vals,
            s_vals.as_ptr() as *const c_double,
            evaluated.as_mut_ptr() as *mut c_double,
        );
    }

    Ok(evaluated)
}

#[pyfunction]
fn compute_curve_length_py(
    num_nodes: i32,
    dimension: i32,
    nodes: Vec<f64>,
) -> PyResult<f64> {
    if num_nodes < 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "num_nodes cannot be negative.",
        ));
    }
    if num_nodes == 0 || num_nodes == 1 {
        return Ok(0.0);
    }

    // num_nodes > 1 from here
    if dimension <= 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "dimension must be positive for curve length computation.",
        ));
    }
    if nodes.len() != (num_nodes * dimension) as usize {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Nodes vector length must be num_nodes * dimension.",
        ));
    }

    let mut length: f64 = 0.0;
    let mut error_val_raw: i32 = status::STATUS_SUCCESS;

    unsafe {
        curve::BEZ_compute_length(
            num_nodes,
            dimension,
            nodes.as_ptr() as *const c_double,
            &mut length as *mut c_double,
            &mut error_val_raw as *mut i32,
        );
    }

    if error_val_raw != status::STATUS_SUCCESS as i32 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Failed to compute curve length. Error code: {}", error_val_raw),
        ));
    }

    Ok(length)
}

#[pymodule]
fn bezier_rust_ffi(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hello_from_rust_py, m)?)?;
    m.add_function(wrap_pyfunction!(evaluate_curve_multi_py, m)?)?;
    m.add_function(wrap_pyfunction!(compute_curve_length_py, m)?)?;
    m.add_function(wrap_pyfunction!(evaluate_triangle_barycentric_multi_py, m)?)?;
    m.add_function(wrap_pyfunction!(subdivide_triangle_nodes_py, m)?)?;
    m.add_function(wrap_pyfunction!(compute_triangle_area_from_edges_py, m)?)?;
    m.add_function(wrap_pyfunction!(curve_intersections_py, m)?)?;
    m.add_function(wrap_pyfunction!(newton_refine_triangle_py, m)?)?;
    // m.add_function(wrap_pyfunction!(all_triangle_intersections_py, m)?)?; // Still commented out
    Ok(())
}

#[pyfunction]
fn evaluate_triangle_barycentric_multi_py(
    num_nodes: i32,
    dimension: i32,
    nodes: Vec<f64>,
    degree: i32,
    param_vals: Vec<f64>, // lambda1, lambda2, lambda3 triplets
) -> PyResult<Vec<f64>> {
    if degree < 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Degree cannot be negative.",
        ));
    }
    if num_nodes <= 0 && degree > 0 { // Or based on how num_nodes relates to degree
         return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "num_nodes must be positive if degree is positive.",
        ));
    }
    // TODO: Add more robust check for consistency between num_nodes and degree,
    // e.g., num_nodes == (degree + 1) * (degree + 2) / 2 for triangles.
    // For now, basic checks.

    let num_params = param_vals.len() / 3;
    if param_vals.len() % 3 != 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "param_vals length must be a multiple of 3",
        ));
    }
    let mut evaluated: Vec<f64> =
        vec![0.0; num_params * dimension as usize];
    
    // The FFI function triangle::BEZ_evaluate_barycentric_multi returns void and has no error parameter.
    // Thus, no error code to check.
    unsafe {
        triangle::BEZ_evaluate_barycentric_multi(
            num_nodes,
            dimension,
            nodes.as_ptr() as *const c_double,
            degree,
            num_params as i32, 
            param_vals.as_ptr() as *const c_double,
            evaluated.as_mut_ptr() as *mut c_double,
        )
    };

    Ok(evaluated)
}

#[pyfunction]
fn subdivide_triangle_nodes_py(
    num_nodes: i32,
    dimension: i32,
    nodes: Vec<f64>,
    degree: i32,
) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    if degree < 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Degree cannot be negative for subdivide_triangle_nodes.",
        ));
    }
    // TODO: Add more robust checks for num_nodes, dimension, and consistency with degree.

    let mut nodes_a_vec: Vec<f64> = vec![0.0; (num_nodes * dimension) as usize];
    let mut nodes_b_vec: Vec<f64> = vec![0.0; (num_nodes * dimension) as usize];
    let mut nodes_c_vec: Vec<f64> = vec![0.0; (num_nodes * dimension) as usize];
    let mut nodes_d_vec: Vec<f64> = vec![0.0; (num_nodes * dimension) as usize];

    // The FFI function triangle::BEZ_subdivide_nodes_triangle returns void and has no error parameter.
    // Thus, no error code to check.
    unsafe {
        triangle::BEZ_subdivide_nodes_triangle(
            num_nodes,
            dimension,
            nodes.as_ptr() as *const c_double,
            degree,
            nodes_a_vec.as_mut_ptr() as *mut c_double,
            nodes_b_vec.as_mut_ptr() as *mut c_double,
            nodes_c_vec.as_mut_ptr() as *mut c_double,
            nodes_d_vec.as_mut_ptr() as *mut c_double,
        )
    };

    Ok((nodes_a_vec, nodes_b_vec, nodes_c_vec, nodes_d_vec))
}

use libc::c_int; // For triangle::BEZ_compute_area

#[pyfunction]
fn compute_triangle_area_from_edges_py(
    edges: Vec<Vec<f64>>, // Each inner Vec<f64> is [x1,y1,x2,y2,...]
) -> PyResult<f64> {
    let dimension: i32 = 2; // Assuming 2D for this wrapper
    let num_edges = edges.len() as i32;
    if num_edges == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Edges list cannot be empty.",
        ));
    }

    let mut sizes: Vec<c_int> = Vec::with_capacity(num_edges as usize);
    let mut nodes_pointers_vec: Vec<*const c_double> = Vec::with_capacity(num_edges as usize);

    // Keep the edge data alive
    let mut edge_data_flattened: Vec<Vec<f64>> = Vec::with_capacity(num_edges as usize);

    for edge in edges {
        if edge.is_empty() || edge.len() % dimension as usize != 0 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Each edge must have an even number of coordinates and not be empty. Found edge with length {}", edge.len()),
            ));
        }
        sizes.push((edge.len() / dimension as usize) as c_int);
        // Store the edge data in a way that it outlives the pointers
        edge_data_flattened.push(edge); 
    }
    
    // Now that edge_data_flattened owns the data, create pointers
    for edge_ref in &edge_data_flattened {
        nodes_pointers_vec.push(edge_ref.as_ptr() as *const c_double);
    }


    let mut area: f64 = 0.0;
    let mut not_implemented: u8 = 0;

    unsafe {
        // Corrected order and removed dimension argument for BEZ_compute_area
        triangle::BEZ_compute_area(
            num_edges,
            sizes.as_ptr(), // sizes comes before nodes_pointers
            nodes_pointers_vec.as_ptr(),
            &mut area as *mut c_double,
            &mut not_implemented as *mut u8,
        );
    }

    if not_implemented != 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyNotImplementedError, _>(
            "Triangle area computation for the given edge configuration is not implemented.",
        ));
    }

    Ok(area)
}

#[pyfunction]
fn curve_intersections_py(
    num_nodes1: i32,
    // dimension1: i32, // BEZ_curve_intersections assumes 2D
    nodes1: Vec<f64>,
    num_nodes2: i32,
    // dimension2: i32, // BEZ_curve_intersections assumes 2D
    nodes2: Vec<f64>,
    intersections_capacity: i32, // Max number of s,t pairs to find (s_max from test)
) -> PyResult<(Vec<f64>, i32, bool, i32)> { // Returns: intersections_buffer, num_intersections_found, coincident_flag, status_code
    if num_nodes1 < 1 || num_nodes2 < 1 {
         return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Number of nodes for each curve must be at least 1.",
        ));
    }
    let dimension = 2; // FFI function is hardcoded for 2D
    if nodes1.len() != (num_nodes1 * dimension) as usize {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Nodes1 vector length must be num_nodes1 * {} (for 2D). Found len {}.", dimension, nodes1.len()),
        ));
    }
    if nodes2.len() != (num_nodes2 * dimension) as usize {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Nodes2 vector length must be num_nodes2 * {} (for 2D). Found len {}.", dimension, nodes2.len()),
        ));
    }
    if intersections_capacity < 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "intersections_capacity cannot be negative.",
        ));
    }

    // Buffer for s,t pairs. Each intersection is 2 f64 values.
    let mut intersections_buffer: Vec<f64> = vec![0.0; (intersections_capacity * 2) as usize];
    let mut num_intersections_found_raw: c_int = 0;
    let mut coincident_raw: u8 = 0;
    let mut status_code_raw: c_int = status::STATUS_SUCCESS;

    unsafe {
        curve_intersection::BEZ_curve_intersections(
            num_nodes1,
            nodes1.as_ptr() as *const c_double,
            num_nodes2,
            nodes2.as_ptr() as *const c_double,
            intersections_capacity, 
            intersections_buffer.as_mut_ptr() as *mut c_double,
            &mut num_intersections_found_raw as *mut c_int,
            &mut coincident_raw as *mut u8,
            &mut status_code_raw as *mut c_int,
        );
    }

    if status_code_raw != status::STATUS_SUCCESS {
         // For now, just returning the status code. More specific error mapping could be done.
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Failed to compute curve intersections: FFI status code {}", status_code_raw),
        ));
    }
    
    intersections_buffer.truncate((num_intersections_found_raw * 2) as usize);
    let coincident_flag = coincident_raw != 0;

    Ok((intersections_buffer, num_intersections_found_raw, coincident_flag, status_code_raw))
}


#[pyfunction]
fn newton_refine_triangle_py(
    num_nodes: i32,
    nodes: Vec<f64>,
    degree: i32,
    x_val: f64,
    y_val: f64,
    s_initial: f64,
    t_initial: f64,
) -> PyResult<(f64, f64)> {
    let mut updated_s_raw: c_double = 0.0;
    let mut updated_t_raw: c_double = 0.0;

    // Ensure nodes vector is not empty if num_nodes > 0
    if num_nodes > 0 && nodes.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Nodes vector cannot be empty if num_nodes > 0.",
        ));
    }
     // Dimension is implicitly 2 for this function as per its C signature and usage (x_val, y_val)
    if num_nodes > 0 && nodes.len() != (num_nodes * 2) as usize {
         return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Nodes vector length must be num_nodes * 2 (for 2D).",
        ));
    }


    unsafe {
        triangle_intersection::BEZ_newton_refine_triangle(
            num_nodes,
            nodes.as_ptr() as *const c_double,
            degree,
            x_val,
            y_val,
            s_initial,
            t_initial,
            &mut updated_s_raw as *mut c_double,
            &mut updated_t_raw as *mut c_double,
        );
    }

    Ok((updated_s_raw, updated_t_raw))
}

// Commenting out because BEZ_all_triangle_intersections and helpers::get_num_nodes are missing/problematic
// #[pyfunction]
// fn all_triangle_intersections_py(
//     num_triangles_in_batch: i32, 
//     nodes_first_batch: Vec<Vec<f64>>, 
//     nodes_second_batch: Vec<Vec<f64>>,
//     degrees_first: Vec<i32>,
//     degrees_second: Vec<i32>,
// ) -> PyResult<Vec<i32>> {
//     let dimension: i32 = 2; 
// 
//     if nodes_first_batch.len() != num_triangles_in_batch as usize ||
//        nodes_second_batch.len() != num_triangles_in_batch as usize ||
//        degrees_first.len() != num_triangles_in_batch as usize ||
//        degrees_second.len() != num_triangles_in_batch as usize {
//         return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
//             "Input Vec lengths must match num_triangles_in_batch",
//         ));
//     }
// 
//     let mut nodes_first_ptrs: Vec<*const c_double> = Vec::with_capacity(num_triangles_in_batch as usize);
//     let mut nodes_second_ptrs: Vec<*const c_double> = Vec::with_capacity(num_triangles_in_batch as usize);
//     
//     let mut _nodes_first_storage: Vec<Vec<f64>> = Vec::with_capacity(num_triangles_in_batch as usize);
//     let mut _nodes_second_storage: Vec<Vec<f64>> = Vec::with_capacity(num_triangles_in_batch as usize);
// 
//     for i in 0..num_triangles_in_batch as usize {
//         // let num_nodes_first = helpers::get_num_nodes(degrees_first[i], dimension); // This function is missing
//         // For now, assume degree implies num_nodes, but this needs actual get_num_nodes
//         let num_nodes_first = (degrees_first[i] + 1) * (degrees_first[i] + 2) / 2; // Placeholder
//         if nodes_first_batch[i].len() != (num_nodes_first * dimension) as usize {
//             return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
//                 format!("Incorrect number of nodes for first triangle at index {}", i),
//             ));
//         }
//         _nodes_first_storage.push(nodes_first_batch[i].clone()); 
//         nodes_first_ptrs.push(_nodes_first_storage.last().unwrap().as_ptr() as *const c_double);
// 
//         // let num_nodes_second = helpers::get_num_nodes(degrees_second[i], dimension); // This function is missing
//         let num_nodes_second = (degrees_second[i] + 1) * (degrees_second[i] + 2) / 2; // Placeholder
//         if nodes_second_batch[i].len() != (num_nodes_second * dimension) as usize {
//             return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
//                 format!("Incorrect number of nodes for second triangle at index {}", i),
//             ));
//         }
//         _nodes_second_storage.push(nodes_second_batch[i].clone()); 
//         nodes_second_ptrs.push(_nodes_second_storage.last().unwrap().as_ptr() as *const c_double);
//     }
// 
//     let mut statuses_vec: Vec<i32> = vec![0; num_triangles_in_batch as usize];
//     let mut actual_error_val: i32 = status::STATUS_SUCCESS;
// 
//     unsafe {
//         // triangle_intersection::BEZ_all_triangle_intersections(
//         //     num_triangles_in_batch,
//         //     nodes_first_ptrs.as_ptr(),
//         //     nodes_second_ptrs.as_ptr(),
//         //     degrees_first.as_ptr(),
//         //     degrees_second.as_ptr(),
//         //     dimension,
//         //     statuses_vec.as_mut_ptr() as *mut i32,
//         //     &mut actual_error_val as *mut i32,
//         // );
//     }
// 
//     if actual_error_val != status::STATUS_SUCCESS {
//         // return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
//         //     format!("Failed to compute all triangle intersections: error code {}", actual_error_val),
//         // ));
//     }
// 
//     Ok(statuses_vec)
// }
