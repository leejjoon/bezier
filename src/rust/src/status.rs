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

use libc;

// Values of Status enum:
// SUCCESS: Procedure exited with no error.
pub const STATUS_SUCCESS: libc::c_int = 0;
// BAD_MULTIPLICITY: An iterative method (e.g. Newton's method) failed to
//                   converge because it encountered a solution with an
//                   unsupported multiplicity (i.e. it was a triple root or
//                   higher). The multiplicity is detected by the observed
//                   rate of convergence. (Used by
//                   ``curve_intersection.full_newton_nonzero()``.)
pub const STATUS_BAD_MULTIPLICITY: libc::c_int = 1;
// NO_CONVERGE: An iterative method failed to converge in the maximum allowed
//              number of iterations. (Used by
//              ``curve_intersection.all_intersections()``.)
pub const STATUS_NO_CONVERGE: libc::c_int = 2;
// INSUFFICIENT_SPACE: Intended to be used by ABI versions of procedures. This
//                     will be used when the caller has not allocated enough
//                     space for the output values. (On the Fortran side, an
//                     ``allocatable`` array can be resized, but in the ABI,
//                     we can only use what the caller has passed.)
pub const STATUS_INSUFFICIENT_SPACE: libc::c_int = 3;
// SAME_CURVATURE: Classification of a curve-curve intersection failed due to
//                 two curves having identical at a tangent intersection.
pub const STATUS_SAME_CURVATURE: libc::c_int = 4;
// BAD_INTERIOR: Caused by a failure during the process of triangle-triangle
//               intersection. Occurs when the corners and edge-edge
//               intersections can't be converted into the curved polygon(s)
//               that make up the intersection of the two triangles.
pub const STATUS_BAD_INTERIOR: libc::c_int = 5;
// EDGE_END: A triangle-triangle intersection point occurs at the **end** of
//           an edge (only intersections at the beginning of an edge should
//           be used).
pub const STATUS_EDGE_END: libc::c_int = 6;
// SINGULAR: Signifies that an attempt was made to solve a linear system
//           that was singular.
pub const STATUS_SINGULAR: libc::c_int = 7;
// UNKNOWN: Signifies a block of code reached an "impossible" state. Either
//          the code has violated some mathematical invariant or the author
//          of that block of code misunderstood the possible states of the
//          system.
pub const STATUS_UNKNOWN: libc::c_int = 999;
