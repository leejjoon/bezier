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
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::comparison_chain)] // To match Fortran structure more easily

use libc::{c_double, c_int}; // Removed c_bool
use std::ptr;

pub const WIGGLE: f64 = 5.6843418860808015e-14; // 0.5f64.powi(44);
pub const VECTOR_CLOSE_EPS: f64 = 9.094947017729282e-13; // 0.5f64.powi(40);

#[no_mangle]
pub unsafe extern "C" fn BEZ_cross_product(
    vec0: *const c_double,
    vec1: *const c_double,
    result: *mut c_double,
) {
    // Safety:
    // - `vec0` must be a valid pointer to at least 2 `c_double` elements.
    // - `vec1` must be a valid pointer to at least 2 `c_double` elements.
    // - `result` must be a valid pointer to 1 `c_double` element.
    // Dereferencing `vec0.offset(0)`, `vec0.offset(1)`, `vec1.offset(0)`, `vec1.offset(1)`
    // is safe if `vec0` and `vec1` point to arrays of at least 2 doubles.
    // Writing to `*result` is safe if `result` is a valid pointer.
    *result = (*vec0.offset(0) * *vec1.offset(1)) - (*vec0.offset(1) * *vec1.offset(0));
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_bbox(
    num_nodes: c_int,
    nodes: *const c_double, 
    left: *mut c_double,
    right: *mut c_double,
    bottom: *mut c_double,
    top: *mut c_double,
) {
    // Safety:
    // - If `num_nodes > 0`, `nodes` must be a valid pointer to at least `2 * num_nodes` `c_double` elements,
    //   representing (x,y) coordinates.
    // - `left`, `right`, `bottom`, `top` must be valid pointers to 1 `c_double` element each.
    // Accessing `nodes.offset((2 * i) as isize)` and `nodes.offset((2 * i + 1) as isize)`
    // is safe if `nodes` points to `2 * num_nodes` elements.
    // Writing to `*left`, `*right`, `*bottom`, `*top` is safe if they are valid pointers.
    if num_nodes == 0 {
        *left = 0.0;
        *right = 0.0;
        *bottom = 0.0;
        *top = 0.0;
        return;
    }

    *left = *nodes.offset(0); 
    *bottom = *nodes.offset(1); 
    *right = *nodes.offset(0); 
    *top = *nodes.offset(1); 

    for i in 1..num_nodes {
        let x_i = *nodes.offset((2 * i) as isize);
        let y_i = *nodes.offset((2 * i + 1) as isize);

        if x_i < *left {
            *left = x_i;
        }
        if x_i > *right {
            *right = x_i;
        }
        if y_i < *bottom {
            *bottom = y_i;
        }
        if y_i > *top {
            *top = y_i;
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_wiggle_interval(
    value: c_double,
    result: *mut c_double,
    success: *mut u8, // c_bool -> u8
) {
    // Safety:
    // - `result` must be a valid pointer to 1 `c_double` element.
    // - `success` must be a valid pointer to 1 `u8` element.
    // Writing to `*result` and `*success` is safe if they are valid pointers.
    *success = 1u8; // true as c_bool -> 1u8
    if -WIGGLE < value && value < WIGGLE {
        *result = 0.0;
    } else if WIGGLE <= value && value <= 1.0 - WIGGLE {
        *result = value;
    } else if 1.0 - WIGGLE < value && value < 1.0 + WIGGLE {
        *result = 1.0;
    } else {
        *success = 0u8; // false as c_bool -> 0u8
    }
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_contains_nd(
    num_nodes: c_int,
    dimension: c_int,
    nodes: *const c_double, 
    point: *const c_double, 
    predicate: *mut u8, // c_bool -> u8
) {
    // Safety:
    // - If `num_nodes > 0`, `nodes` must be a valid pointer to at least `dimension * num_nodes` `c_double` elements.
    // - `point` must be a valid pointer to at least `dimension` `c_double` elements.
    // - `predicate` must be a valid pointer to 1 `u8` element.
    // - `dimension` must be positive if `num_nodes > 0`.
    // Accessing `nodes.offset((node_idx * dimension + dim_idx) as isize)` is safe under these conditions.
    // Accessing `point.offset(dim_idx as isize)` is safe.
    // Writing to `*predicate` is safe.
    if num_nodes == 0 {
        *predicate = 0u8; // false as c_bool -> 0u8
        return;
    }

    for dim_idx in 0..dimension { 
        let mut min_val_dim = *nodes.offset((0 * dimension + dim_idx) as isize); 
        let mut max_val_dim = *nodes.offset((0 * dimension + dim_idx) as isize);

        for node_idx in 1..num_nodes { 
            let val = *nodes.offset((node_idx * dimension + dim_idx) as isize);
            if val < min_val_dim {
                min_val_dim = val;
            }
            if val > max_val_dim {
                max_val_dim = val;
            }
        }

        let p_val = *point.offset(dim_idx as isize);
        if p_val < min_val_dim || p_val > max_val_dim {
            *predicate = 0u8; // false as c_bool -> 0u8
            return;
        }
    }

    *predicate = 1u8; // true as c_bool -> 1u8
}

#[no_mangle]
pub unsafe extern "C" fn BEZ_vector_close(
    num_values: c_int,
    vec1: *const c_double, 
    vec2: *const c_double, 
    eps: c_double,
) -> u8 { // c_bool -> u8
    // Safety:
    // - `vec1` must be a valid pointer to at least `num_values` `c_double` elements.
    // - `vec2` must be a valid pointer to at least `num_values` `c_double` elements.
    // - `num_values` must accurately reflect the number of elements to compare.
    // Accessing `vec1.offset(i as isize)` and `vec2.offset(i as isize)` is safe under these conditions.
    let mut s1_sq = 0.0;
    let mut s2_sq = 0.0;
    let mut diff_sq = 0.0;

    for i in 0..num_values {
        let v1_i = *vec1.offset(i as isize);
        let v2_i = *vec2.offset(i as isize);
        s1_sq += v1_i * v1_i;
        s2_sq += v2_i * v2_i;
        let diff = v1_i - v2_i;
        diff_sq += diff * diff;
    }

    let size1 = s1_sq.sqrt();
    let size2 = s2_sq.sqrt();

    if size1 == 0.0 {
        (size2 <= eps) as u8 // c_bool -> u8
    } else if size2 == 0.0 {
        (size1 <= eps) as u8 // c_bool -> u8
    } else {
        (diff_sq.sqrt() <= eps * size1.min(size2)) as u8 // c_bool -> u8
    }
}

#[no_mangle]
pub extern "C" fn BEZ_in_interval(
    value: c_double,
    start: c_double,
    end: c_double,
) -> u8 { // c_bool -> u8
    (start <= value && value <= end) as u8 // c_bool -> u8
}

/// # Safety
/// - `points_ptr` must be a valid pointer to at least `2 * num_points` `f64` elements.
/// - `num_points` must be non-negative and accurately reflect the number of (x,y) pairs.
/// - `match_idx_0based` must be a valid pointer to an `i32`.
unsafe fn min_index_rs(num_points: i32, points_ptr: *const f64, match_idx_0based: &mut i32) {
    if num_points == 0 {
        return;
    }
    *match_idx_0based = 0; 
    for i in 1..num_points {
        // Safety: Accessing points_ptr offsets is safe if points_ptr is valid for 2*num_points elements
        // and i is within 0..num_points.
        let x_i = *points_ptr.offset((2 * i) as isize);
        let x_match = *points_ptr.offset((2 * (*match_idx_0based)) as isize);

        if x_i < x_match {
            *match_idx_0based = i;
        } else if x_i == x_match {
            let y_i = *points_ptr.offset((2 * i + 1) as isize);
            let y_match = *points_ptr.offset((2 * (*match_idx_0based) + 1) as isize);
            if y_i < y_match {
                *match_idx_0based = i;
            }
        }
    }
}

/// # Safety
/// - `points_ptr` must be a valid pointer to at least `2 * num_points` `f64` elements, and must be mutable.
/// - `num_points` must be non-negative and accurately reflect the number of (x,y) pairs.
/// - `num_uniques_rs` must be a valid pointer to an `i32`.
/// This function calls `min_index_rs` which has its own safety requirements.
unsafe fn sort_in_place_rs(num_points: i32, points_ptr: *mut f64, num_uniques_rs: &mut i32) {
    if num_points == 0 {
        *num_uniques_rs = 0;
        return;
    }

    *num_uniques_rs = num_points;
    
    // Place the "smallest" point first.
    let mut match_0based = 0;
    min_index_rs(*num_uniques_rs, points_ptr, &mut match_0based);
    if match_0based != 0 {
        // Swap points_ptr[0] with points_ptr[match_0based]
        let x0 = *points_ptr.offset(0);
        let y0 = *points_ptr.offset(1);
        *points_ptr.offset(0) = *points_ptr.offset((2 * match_0based) as isize);
        *points_ptr.offset(1) = *points_ptr.offset((2 * match_0based + 1) as isize);
        *points_ptr.offset((2 * match_0based) as isize) = x0;
        *points_ptr.offset((2 * match_0based + 1) as isize) = y0;
    }

    let mut i_1based = 2; // Fortran loop `i = 2, num_uniques`
    while i_1based <= *num_uniques_rs {
        let i_0based = i_1based -1; // Current position to fill (0-indexed)

        // Find min in points_ptr[i_0based .. *num_uniques_rs-1]
        // Number of elements in slice: *num_uniques_rs - i_0based
        // Pointer to start of slice: points_ptr.offset((2 * i_0based) as isize)
        let mut match_relative_0based = 0;
        min_index_rs(
            *num_uniques_rs - i_0based,
            points_ptr.offset((2 * i_0based) as isize),
            &mut match_relative_0based,
        );
        let match_absolute_0based = match_relative_0based + i_0based;

        // points_ptr[match_absolute_0based] is the smallest in the remainder.
        // Compare with points_ptr[i_0based - 1] (the previous unique point)
        let prev_x = *points_ptr.offset((2 * (i_0based - 1)) as isize);
        let prev_y = *points_ptr.offset((2 * (i_0based - 1) + 1) as isize);
        let current_match_x = *points_ptr.offset((2 * match_absolute_0based) as isize);
        let current_match_y = *points_ptr.offset((2 * match_absolute_0based + 1) as isize);

        if current_match_x == prev_x && current_match_y == prev_y {
            // Duplicate of previous unique point found at points_ptr[match_absolute_0based].
            // Move it to the end of the unique part by swapping with points_ptr[*num_uniques_rs - 1].
            if match_absolute_0based != (*num_uniques_rs - 1) { // Avoid self-swap
                let end_x = *points_ptr.offset((2 * (*num_uniques_rs - 1)) as isize);
                let end_y = *points_ptr.offset((2 * (*num_uniques_rs - 1) + 1) as isize);
                
                *points_ptr.offset((2 * match_absolute_0based) as isize) = end_x;
                *points_ptr.offset((2 * match_absolute_0based + 1) as isize) = end_y;
                // The duplicate (current_match_x,y) is not explicitly placed at the end,
                // it's overwritten and that part of array will be ignored.
            }
            *num_uniques_rs -= 1;
            // i_1based does not increment in this path, loop re-evaluates with smaller num_uniques_rs
        } else {
            // Not a duplicate of the previous unique point.
            // If smallest (at match_absolute_0based) is not already at current position (i_0based), swap them.
            if match_absolute_0based != i_0based {
                let i_x = *points_ptr.offset((2 * i_0based) as isize);
                let i_y = *points_ptr.offset((2 * i_0based + 1) as isize);
                *points_ptr.offset((2 * i_0based) as isize) = current_match_x;
                *points_ptr.offset((2 * i_0based + 1) as isize) = current_match_y;
                *points_ptr.offset((2 * match_absolute_0based) as isize) = i_x;
                *points_ptr.offset((2 * match_absolute_0based + 1) as isize) = i_y;
            }
            i_1based += 1; // Increment only if not a duplicate
        }
    }
}


fn in_sorted_rs(num_values: i32, values_ptr: *const i32, value_to_find: i32) -> bool {
    if num_values == 0 {
        return false;
    }
    let mut left = 0; 
    let mut right = num_values - 1; 

    while left < right {
        let midpoint = left + (right - left) / 2; 
        let mid_val = unsafe { *values_ptr.offset(midpoint as isize) };
        if value_to_find == mid_val {
            return true;
        } else if value_to_find < mid_val {
            right = midpoint - 1;
        } else {
            left = midpoint + 1;
        }
    }
    
    if left <= right && left < num_values { // Ensure 'left' is a valid index
         unsafe { *values_ptr.offset(left as isize) == value_to_find }
    } else {
        false
    }
}


#[no_mangle]
pub unsafe extern "C" fn BEZ_simple_convex_hull(
    num_points: c_int,
    points: *const c_double, 
    polygon_size: *mut c_int,
    polygon: *mut c_double, 
) {
    // Safety:
    // - `points` must be valid for reading `2 * num_points` doubles if `num_points > 0`.
    // - `polygon_size` must be a valid pointer to `c_int`.
    // - `polygon` must be valid for writing `2 * num_points` doubles in the worst case (all points are unique and form the hull).
    // - `num_points` must be non-negative.
    // - `ptr::copy_nonoverlapping` calls require source and dest to be valid and non-overlapping for the given count.
    // - `sort_in_place_rs` has its own safety requirements for `uniques_vec.as_mut_ptr()`.
    // - `in_sorted_rs` has its own safety requirements.
    // - Offsets into `uniques_vec` and `polygon` must be within bounds.
    if num_points == 0 {
        *polygon_size = 0;
        return;
    }

    let mut uniques_vec: Vec<f64> = vec![0.0; (num_points * 2) as usize];
    ptr::copy_nonoverlapping(points, uniques_vec.as_mut_ptr(), uniques_vec.len());
    
    let mut num_uniques_rs: c_int = 0; // This will be out param for sort_in_place_rs
    sort_in_place_rs(num_points, uniques_vec.as_mut_ptr(), &mut num_uniques_rs);

    if num_uniques_rs == 0 {
        *polygon_size = 0;
        return;
    } else if num_uniques_rs < 3 {
        *polygon_size = num_uniques_rs;
        ptr::copy_nonoverlapping(
            uniques_vec.as_ptr(),
            polygon,
            (num_uniques_rs * 2) as usize,
        );
        return;
    }

    let mut lower_hull_indices: Vec<i32> = Vec::with_capacity(num_uniques_rs as usize);
    let mut upper_hull_indices: Vec<i32> = Vec::with_capacity(num_uniques_rs as usize);

    // Lower hull (using 0-based indices for uniques_vec)
    lower_hull_indices.push(0); // Index of first point in uniques_vec
    lower_hull_indices.push(1); // Index of second point
    for i_idx in 2..num_uniques_rs { // i_idx is the index in uniques_vec
        let p3_x = uniques_vec[(2 * i_idx) as usize];
        let p3_y = uniques_vec[(2 * i_idx + 1) as usize];
        
        while lower_hull_indices.len() >= 2 {
            let p2_unique_idx = lower_hull_indices[lower_hull_indices.len() - 1];
            let p1_unique_idx = lower_hull_indices[lower_hull_indices.len() - 2];

            let p1_x = uniques_vec[(2 * p1_unique_idx) as usize];
            let p1_y = uniques_vec[(2 * p1_unique_idx + 1) as usize];
            let p2_x = uniques_vec[(2 * p2_unique_idx) as usize];
            let p2_y = uniques_vec[(2 * p2_unique_idx + 1) as usize];

            let cp = (p2_x - p1_x) * (p3_y - p1_y) - (p2_y - p1_y) * (p3_x - p1_x);
            if cp <= 0.0 {
                lower_hull_indices.pop();
            } else {
                break;
            }
        }
        lower_hull_indices.push(i_idx);
    }

    // Upper hull
    // Fortran: upper(1) = num_uniques (1-based index into uniques array)
    // Rust: upper_hull_indices.push(num_uniques_rs - 1) (0-based index)
    upper_hull_indices.push(num_uniques_rs - 1);
    // Fortran: upper(2) = num_uniques - 1 (implicit, handled by loop start)
    // The loop structure in Fortran is:
    // do i = num_uniques - 1, 1, -1  (i is 1-based index into uniques)
    //   if (i > 1 .AND. in_sorted(num_lower, lower(:num_lower), i)) then cycle
    // My loop: for i_idx in (0..num_uniques_rs - 1).rev() which is num_uniques_rs-2 down to 0.
    // This means my p3 will be uniques_vec[i_idx].
    // Fortran's `num_uniques - 1` is my `num_uniques_rs - 2`.
    // Fortran's `1` is my `0`.

    // Initial fill for upper hull (at least two points needed for cross product logic)
    // The first point added is uniques(num_uniques-1).
    // The loop then considers uniques(num_uniques-2) down to uniques(0).
    // The Fortran logic adds num_uniques as the first point, then num_uniques-1 as the second.
    // The loop starts from i = num_uniques - 2 (1-based).
    // So the first p3 considered is uniques(num_uniques-2).

    // Let's trace Fortran: uniques = [u0,u1,u2,u3,u4], num_uniques=5
    // lower_hull_indices = [0,1,2,4] (example)
    // upper(1) = 5 (means u4) -> upper_hull_indices = [4]
    // Loop i = 4 down to 1 (1-based for `uniques` array)
    // i = 4 (u3): not in_sorted(lower_hull excluding ends). p3 = u3.
    //    num_upper=1. Loop for pop: while(1>1) false. num_upper=2. upper(2)=4 (u3). -> upper_hull = [4,3]
    // i = 3 (u2): in_sorted. skip.
    // i = 2 (u1): not in_sorted. p3 = u1.
    //    num_upper=2. p1=u4, p2=u3. cp = cross(u3-u4, u1-u4). if <=0 pop. Suppose it is. num_upper=1. upper_hull=[4].
    //    Loop again: while(1>1) false. num_upper=2. upper(2)=2 (u1). -> upper_hull = [4,1]
    // i = 1 (u0): not in_sorted (it's an end of lower_hull). p3 = u0.
    //    num_upper=2. p1=u4, p2=u1. cp = cross(u1-u4, u0-u4). if <=0 pop. Suppose not. break.
    //    num_upper=3. upper(3)=1 (u0). -> upper_hull = [4,1,0]

    // Rust equivalent:
    // upper_hull_indices = [num_uniques_rs-1] (e.g. [4])
    // Consider adding num_uniques_rs-2:
    // This point must be added before loop if loop starts from num_uniques_rs-3
    if num_uniques_rs > 1 { // Need at least 2 points for initial upper hull segment
         // Check if num_uniques_rs-2 is part of internal lower_hull_indices
        let second_last_idx = num_uniques_rs - 2;
        let mut skip = false;
        if second_last_idx != lower_hull_indices[0] && second_last_idx != lower_hull_indices[lower_hull_indices.len()-1] {
            if in_sorted_rs(lower_hull_indices.len() as i32, lower_hull_indices.as_ptr(), second_last_idx) {
                 skip = true;
            }
        }
        if !skip {
            upper_hull_indices.push(second_last_idx);
        }
    }


    // Loop from third-to-last unique point, down to the first unique point.
    // (i.e. num_uniques_rs-3 down to 0)
    if num_uniques_rs > 2 { // Only if there are more points to consider
        for i_idx in (0..num_uniques_rs - 2).rev() { 
            let mut skip = false;
            // Fortran: if (i > 1 .AND. in_sorted(num_lower, lower(:num_lower), i_1based))
            // i_1based is i_idx + 1. Fortran `i > 1` means `i_idx+1 > 1` -> `i_idx > 0`.
            // So skip if i_idx is not the first unique point (0) AND i_idx is not the last unique point (num_uniques_rs-1)
            // AND i_idx is in the main body of lower_hull_indices.
            // The Fortran condition `i > 1` means the current point `uniques(i)` is not `uniques(1)` (the first point).
            // The `in_sorted` check is on `lower` which contains 1-based indices of `uniques`.
            // `lower(1)` is `1` (index of uniques(1)). `lower(num_lower)` is `num_uniques` (index of uniques(num_uniques)).
            // So, if `uniques(i)` is `uniques(lower(k))` for `1 < k < num_lower`.
            if i_idx != lower_hull_indices[0] && i_idx != lower_hull_indices[lower_hull_indices.len()-1] {
                 if in_sorted_rs(lower_hull_indices.len() as i32, lower_hull_indices.as_ptr(), i_idx) { // Pass 0-based i_idx
                     skip = true;
                }
            }
            if skip { continue; }

            let p3_x = uniques_vec[(2 * i_idx) as usize];
            let p3_y = uniques_vec[(2 * i_idx + 1) as usize];

            while upper_hull_indices.len() >= 2 {
                let p2_unique_idx = upper_hull_indices[upper_hull_indices.len() - 1];
                let p1_unique_idx = upper_hull_indices[upper_hull_indices.len() - 2];

                let p1_x = uniques_vec[(2 * p1_unique_idx) as usize];
                let p1_y = uniques_vec[(2 * p1_unique_idx + 1) as usize];
                let p2_x = uniques_vec[(2 * p2_unique_idx) as usize];
                let p2_y = uniques_vec[(2 * p2_unique_idx + 1) as usize];
                
                let cp = (p2_x - p1_x) * (p3_y - p1_y) - (p2_y - p1_y) * (p3_x - p1_x);
                if cp <= 0.0 { 
                    upper_hull_indices.pop();
                } else {
                    break;
                }
            }
            upper_hull_indices.push(i_idx);
        }
    }
    
    let mut current_polygon_size = 0;
    for k_idx in 0..lower_hull_indices.len() - 1 {
        let unique_idx = lower_hull_indices[k_idx];
        *polygon.offset((2 * current_polygon_size) as isize) = uniques_vec[(2 * unique_idx) as usize];
        *polygon.offset((2 * current_polygon_size + 1) as isize) = uniques_vec[(2 * unique_idx + 1) as usize];
        current_polygon_size += 1;
    }
    // If upper_hull_indices became empty (e.g. colinear points removed one by one)
    if !upper_hull_indices.is_empty() {
        for k_idx in 0..upper_hull_indices.len() - 1 {
            let unique_idx = upper_hull_indices[k_idx];
            *polygon.offset((2 * current_polygon_size) as isize) = uniques_vec[(2 * unique_idx) as usize];
            *polygon.offset((2 * current_polygon_size + 1) as isize) = uniques_vec[(2 * unique_idx + 1) as usize];
            current_polygon_size += 1;
        }
    }
    *polygon_size = current_polygon_size as c_int;
}


unsafe fn is_separating_rs(
    edge_direction_ptr: *const f64, 
    polygon_size1: i32,
    polygon1_ptr: *const f64, 
    polygon_size2: i32,
    polygon2_ptr: *const f64, 
) -> bool {
    // Safety:
    // - `edge_direction_ptr` must be valid for reading 2 `f64`s.
    // - `polygon1_ptr` must be valid for reading `2 * polygon_size1` `f64`s if `polygon_size1 > 0`.
    // - `polygon2_ptr` must be valid for reading `2 * polygon_size2` `f64`s if `polygon_size2 > 0`.
    // - `polygon_size1` and `polygon_size2` must be non-negative.
    let edge_dx = *edge_direction_ptr.offset(0);
    let edge_dy = *edge_direction_ptr.offset(1);

    let norm_squared = edge_dx * edge_dx + edge_dy * edge_dy;
    if norm_squared == 0.0 { 
        return false; 
    }

    let mut min_param1: f64 = 0.0;
    let mut max_param1: f64 = 0.0;
    for i in 0..polygon_size1 {
        let vx = *polygon1_ptr.offset((2 * i) as isize);
        let vy = *polygon1_ptr.offset((2 * i + 1) as isize);
        let param = (edge_dx * vy - edge_dy * vx) / norm_squared; // Projection onto axis (edge_dy, -edge_dx)
        if i == 0 {
            min_param1 = param;
            max_param1 = param;
        } else {
            if param < min_param1 { min_param1 = param; }
            if param > max_param1 { max_param1 = param; }
        }
    }

    let mut min_param2: f64 = 0.0;
    let mut max_param2: f64 = 0.0;
    for i in 0..polygon_size2 {
        let vx = *polygon2_ptr.offset((2 * i) as isize);
        let vy = *polygon2_ptr.offset((2 * i + 1) as isize);
        let param = (edge_dx * vy - edge_dy * vx) / norm_squared;
        if i == 0 {
            min_param2 = param;
            max_param2 = param;
        } else {
            if param < min_param2 { min_param2 = param; }
            if param > max_param2 { max_param2 = param; }
        }
    }
    
    min_param1 > max_param2 || max_param1 < min_param2
}


#[no_mangle]
pub unsafe extern "C" fn BEZ_polygon_collide(
    polygon_size1: c_int,
    polygon1: *const c_double, 
    polygon_size2: c_int,
    polygon2: *const c_double, 
    collision: *mut u8, // c_bool -> u8
) {
    // Safety:
    // - `polygon1` must be valid for reading `2 * polygon_size1` doubles if `polygon_size1 > 0`.
    // - `polygon2` must be valid for reading `2 * polygon_size2` doubles if `polygon_size2 > 0`.
    // - `collision` must be a valid pointer to a `u8`.
    // - `polygon_size1` and `polygon_size2` must be non-negative.
    // - `is_separating_rs` has its own safety requirements.
    //   `edge_direction.as_ptr()` is safe as `edge_direction` is a local array.
    // Polygons need at least 1 point for BEZ_bbox, but SAT typically assumes >= 3 points for closed shapes.
    // Fortran code implies polygons are just sequences of points, edges are formed implicitly.
    // If size < 1, it's an error or no collision. If size 1 or 2, they are points/lines, not areas.
    // The Fortran code doesn't explicitly check for polygon_size < 3.
    // Let's assume valid inputs as per Fortran's implicit assumptions for SAT.
    if polygon_size1 < 1 || polygon_size2 < 1 { 
        *collision = 0u8; // false as c_bool -> 0u8
        return;
    }
    
    // Check edges of polygon1
    // Fortran: edge_direction = polygon1(:, 1) - polygon1(:, polygon_size1)
    // Then: do i = 2, polygon_size1; edge_direction = polygon1(:, i) - polygon1(:, i - 1)
    for i_plus1 in 1..=polygon_size1 { // Iterate 1 to polygon_size1 (Fortran style)
        let p1_idx_0based = (i_plus1 -1) as isize; // Current vertex (0-indexed)
        let p2_idx_0based = if i_plus1 == 1 { (polygon_size1 - 1) as isize} else { (i_plus1 - 2) as isize}; // Previous vertex

        let p1_x = *polygon1.offset(2 * p1_idx_0based);
        let p1_y = *polygon1.offset(2 * p1_idx_0based + 1);
        let p2_x = *polygon1.offset(2 * p2_idx_0based);
        let p2_y = *polygon1.offset(2 * p2_idx_0based + 1);
        
        let edge_direction: [f64; 2] = [p1_x - p2_x, p1_y - p2_y];

        if is_separating_rs(edge_direction.as_ptr(), polygon_size1, polygon1, polygon_size2, polygon2) {
            *collision = 0u8; // false as c_bool -> 0u8
            return;
        }
    }

    // Check edges of polygon2
    for i_plus1 in 1..=polygon_size2 {
        let p1_idx_0based = (i_plus1 -1) as isize;
        let p2_idx_0based = if i_plus1 == 1 { (polygon_size2 - 1) as isize} else { (i_plus1 - 2) as isize};

        let p1_x = *polygon2.offset(2 * p1_idx_0based);
        let p1_y = *polygon2.offset(2 * p1_idx_0based + 1);
        let p2_x = *polygon2.offset(2 * p2_idx_0based);
        let p2_y = *polygon2.offset(2 * p2_idx_0based + 1);

        let edge_direction: [f64; 2] = [p1_x - p2_x, p1_y - p2_y];
        
        if is_separating_rs(edge_direction.as_ptr(), polygon_size1, polygon1, polygon_size2, polygon2) {
            *collision = 0u8; // false as c_bool -> 0u8
            return;
        }
    }

    *collision = 1u8; // true as c_bool -> 1u8
}

#[allow(clippy::many_single_char_names)]
unsafe fn solve2x2_rs(
    lhs_ptr: *const f64, 
    rhs_ptr: *const f64, 
    singular: &mut bool, // Stays as Rust bool for internal function
    x_val: &mut f64,
    y_val: &mut f64,
) {
    // Safety:
    // - `lhs_ptr` must be a valid pointer to at least 4 `f64` elements (for a 2x2 matrix).
    //   Assumed column-major: [a, c, b, d] -> lhs[0]=a, lhs[1]=c, lhs[2]=b, lhs[3]=d
    // - `rhs_ptr` must be a valid pointer to at least 2 `f64` elements.
    // - `singular`, `x_val`, `y_val` must be valid mutable references.
    let a = *lhs_ptr.offset(0); 
    let c = *lhs_ptr.offset(1); 
    let b = *lhs_ptr.offset(2); 
    let d = *lhs_ptr.offset(3); 

    let e = *rhs_ptr.offset(0); 
    let f = *rhs_ptr.offset(1); 

    if c.abs() > a.abs() {
        if c == 0.0 { // Should not happen if c.abs() > a.abs() unless a is also 0 and non-finite.
                      // Fortran does not check this inner condition. It relies on denominator check.
            *singular = true;
            return;
        }
        let ratio = a / c;
        let denominator = b - ratio * d;
        if denominator == 0.0 {
            *singular = true;
            return;
        }
        *y_val = (e - ratio * f) / denominator;
        *x_val = (f - d * *y_val) / c; 
        *singular = false;
    } else { // a.abs() >= c.abs()
        if a == 0.0 { // This means c must also be 0.0 if a.abs() >= c.abs().
            *singular = true;
            return;
        }
        let ratio = c / a;
        let denominator = d - ratio * b;
        if denominator == 0.0 {
            *singular = true;
            return;
        }
        *y_val = (f - ratio * e) / denominator;
        *x_val = (e - b * *y_val) / a; 
        *singular = false;
    }
}
