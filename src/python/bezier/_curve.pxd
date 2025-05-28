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

"""Cython wrapper for Rust FFI for curve module."""

ctypedef unsigned char uint8_t

cdef extern from "../../rust/target/include/bezier_rust_ffi.h":
    void BEZ_evaluate_curve_barycentric(
        int num_nodes, int dimension,
        const double* nodes, int num_vals, const double* lambda1,
        const double* lambda2, double* evaluated)
    void BEZ_evaluate_multi(
        int num_nodes, int dimension,
        const double* nodes, int num_vals, const double* s_vals,
        double* evaluated)
    void BEZ_specialize_curve(
        int num_nodes, int dimension,
        const double* nodes, double start_s, double end_s, # Changed start and end to value
        double* new_nodes)
    void BEZ_evaluate_hodograph(
        double s, int num_nodes, # s by value
        int dimension, const double* nodes, double* hodograph)
    void BEZ_subdivide_nodes_curve(
        int num_nodes, int dimension,
        const double* nodes, double* left_nodes, double* right_nodes)
    void BEZ_newton_refine_curve(
        int num_nodes, int dimension,
        const double* nodes, const double* point, double s, # s by value
        double* updated_s)
    void BEZ_locate_point_curve(
        int num_nodes, int dimension,
        const double* nodes, const double* point, double* s_approx)
    void BEZ_elevate_nodes_curve(
        int num_nodes, int dimension,
        const double* nodes, double* elevated)
    void BEZ_get_curvature(
        int num_nodes, const double* nodes, # dimension is implicitly 2
        const double* tangent_vec, double s, double* curvature) # s by value
    void BEZ_reduce_pseudo_inverse(
        int num_nodes, int dimension,
        const double* nodes, double* reduced, uint8_t* not_implemented)
    void BEZ_full_reduce(
        int num_nodes, int dimension,
        const double* nodes, int* num_reduced_nodes, double* reduced,
        uint8_t* not_implemented)
    void BEZ_compute_length(
        int num_nodes, int dimension,
        const double* nodes, double* length, int* error_val)
