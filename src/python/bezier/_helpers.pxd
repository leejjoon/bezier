# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Cython wrapper for Rust FFI for helpers."""

ctypedef unsigned char uint8_t

cdef extern from "../../rust/target/include/bezier_rust_ffi.h":
    void BEZ_cross_product(
        const double* vec0, const double* vec1, double* result)
    void BEZ_bbox(
        int num_nodes, const double* nodes, double* left,
        double* right, double* bottom, double* top)
    void BEZ_wiggle_interval(
        double value, double* result, uint8_t* success)
    void BEZ_contains_nd(
        int num_nodes, int dimension,
        const double* nodes, const double* point, uint8_t* predicate)
    uint8_t BEZ_vector_close(
        int num_values, const double* vec1, const double* vec2,
        double eps)
    uint8_t BEZ_in_interval(
        double value, double start, double end_) # 'end' is a keyword in Cython, Rust uses 'end_' or similar if C maps 'end'
    void BEZ_simple_convex_hull(
        int num_points, const double* points, int* polygon_size, double* polygon)
    void BEZ_polygon_collide(
        int polygon_size1, const double* polygon1,
        int polygon_size2, const double* polygon2, uint8_t* collision)
