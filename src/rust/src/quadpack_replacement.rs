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

#![allow(dead_code)] // Will be used by other modules later

use libc::c_int;
use quadrature::{integrate, Output}; // Corrected imports, removed ErrorType
use crate::status;

// Default limit from quadrature::double_exponential, approximately 350.
// dqagse's limit is for subintervals, this crate's is for function evaluations.
const QUADRATURE_CRATE_EVAL_LIMIT_APPROX: u64 = 350;

#[allow(clippy::too_many_arguments)]
pub fn rust_dqagse<F>(
    func: F,
    a: f64,
    b: f64,
    epsabs: f64,
    epsrel: f64,
    limit: i32, // Corresponds to `limit` in dqagse, for max subintervals.
                // Not directly used by `quadrature` crate in the same way.
    // Outputs:
    result: &mut f64,
    abserr: &mut f64,
    // neval: &mut i32, // `quadrature::Output` has `iterations` (u64)
    ier: &mut c_int,
) where
    F: Fn(f64) -> f64,
{
    if limit <= 0 {
        *ier = status::STATUS_BAD_INTERIOR; // Or a more specific "invalid_input" status
                                        // DQAGSE uses IER=6 for invalid input.
                                        // STATUS_BAD_INTERIOR is not a perfect match.
                                        // Let's use a generic error for now or add one.
                                        // For now, reusing STATUS_UNKNOWN for bad input limit.
        *result = 0.0;
        *abserr = 0.0;
        // *neval = 0;
        *ier = status::STATUS_UNKNOWN; // Placeholder for invalid input
        return;
    }

    // The `quadrature` crate does not directly use `limit` in the way DQAGSE does
    // (number of subintervals). It has an internal evaluation limit (~350 evaluations).
    // We can check if the returned `iterations` exceeds what might be implied by `limit`,
    // but it's not a 1-to-1 mapping.
    // For now, `limit` parameter is noted but not directly passed or strictly enforced
    // in the same way as DQAGSE.

    // `quadrature::integrate` v0.1.2 only takes absolute error.
    // We will pass `epsabs` and then check `epsrel` against the result.
    let output: Output = integrate(func, a, b, epsabs);

    *result = output.integral; // Field name is `integral` not `estimate` in v0.1.2 Output
    *abserr = output.error_estimate; // Field name is `error_estimate` not `error`

    // *neval = output.num_evals as i32; // Field name is `num_evals`
    // The `iterations` field used before was from a misremembered/different crate's API.
    // quadrature::Output for 0.1.2: pub struct Output { pub integral: f64, pub error_estimate: f64, pub num_evals: u32 }
    // It does not have `error_type` or `iterations` as previously assumed.
    // This simplifies error checking significantly. The crate gives one error estimate.

    // Determine status (ier)
    // The Output struct for quadrature 0.1.2 is:
    // pub struct Output {
    //     pub integral: f64,
    //     pub error_estimate: f64,
    //     pub num_evals: u32,
    // }
    // There is no `error_type`. We must check both conditions.
    // Variables achieved_abs_err, achieved_rel_err, and success were removed as they were unused.
    // The logic is directly in success_check.
    
    // If the result is zero, the relative error check `(*abserr / *result).abs() <= epsrel` is problematic.
    // If *result is exactly 0.0, then success depends only on `*abserr <= epsabs`.
    // If *result is very small (but not zero), `epsrel * |*result|` could be tiny, making the condition hard to meet.
    // Let's refine success condition to match typical library behavior:
    // The error test is satisfied if EITHER abserr <= epsabs OR relerr <= epsrel.
    // However, QUADPACK's `dqagse` (which this aims to mimic behaviorally) uses:
    //   ERREST <= MAX(EPSABS, EPSREL*ABS(RESULT))
    // So, if *abserr <= epsabs, that's good enough.
    // If not, but *abserr <= epsrel * abs(*result), that's also good enough (provided *result isn't zero).
    
    // Refined success condition based on typical interpretation of epsabs and epsrel:
    // The integration is considered successful if the estimated absolute error (`*abserr`)
    // is less than or equal to `epsabs` OR it's less than or equal to `epsrel` times
    // the absolute value of the integral result (`*result`).
    // If `*result` is zero, then only the `epsabs` condition can apply.
    
    let success_check = if *result == 0.0 {
        *abserr <= epsabs
    } else {
        *abserr <= epsabs || *abserr <= epsrel * (*result).abs()
    };


    if success_check {
        *ier = status::STATUS_SUCCESS;
    } else {
        // If not successful, it's likely due to requested tolerance not being met.
        // This could be because the internal evaluation limit was hit.
        // `output.num_function_evaluations` can be checked against QUADRATURE_CRATE_EVAL_LIMIT_APPROX
        if output.num_function_evaluations as u64 >= QUADRATURE_CRATE_EVAL_LIMIT_APPROX {
            *ier = status::STATUS_NO_CONVERGE; // Max evaluations likely hit
        } else {
            // Tolerance not met for other reasons
            *ier = status::STATUS_NO_CONVERGE; 
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_integrate_x_squared() {
        let f = |x: f64| x * x;
        let a = 0.0;
        let b = 1.0;
        let epsabs = 1e-10;
        let epsrel = 1e-10;
        let limit = 50; // Arbitrary, as it's not directly used by `quadrature` in same way

        let mut result = 0.0;
        let mut abserr = 0.0;
        let mut ier = 0;

        rust_dqagse(f, a, b, epsabs, epsrel, limit, &mut result, &mut abserr, &mut ier);

        assert_eq!(ier, status::STATUS_SUCCESS);
        assert!((result - (1.0/3.0)).abs() < epsabs);
        assert!(abserr < epsabs);
    }

    #[test]
    fn test_integrate_sin_x() {
        let f = |x: f64| x.sin();
        let a = 0.0;
        let b = PI;
        let epsabs = 1e-10;
        let epsrel = 1e-10;
        let limit = 50;

        let mut result = 0.0;
        let mut abserr = 0.0;
        let mut ier = 0;

        rust_dqagse(f, a, b, epsabs, epsrel, limit, &mut result, &mut abserr, &mut ier);
        
        // Expected: -cos(PI) - (-cos(0)) = -(-1) - (-1) = 1 - (-1) = 2
        assert_eq!(ier, status::STATUS_SUCCESS);
        assert!((result - 2.0).abs() < epsabs);
        assert!(abserr < epsabs);
    }

    #[test]
    fn test_bad_limit() {
        let f = |x: f64| x * x;
        let a = 0.0;
        let b = 1.0;
        let epsabs = 1e-10;
        let epsrel = 1e-10;
        let limit = 0; // Invalid limit

        let mut result = 0.0;
        let mut abserr = 0.0;
        let mut ier = 0;

        rust_dqagse(f, a, b, epsabs, epsrel, limit, &mut result, &mut abserr, &mut ier);
        
        assert_eq!(ier, status::STATUS_UNKNOWN); // Or a more specific invalid input status
    }

    // Test that might hit evaluation limit or fail to converge (if possible to construct)
    // The `quadrature` crate has a hardcoded limit of ~350 evaluations.
    // A function that is highly oscillatory or has singularities might be hard for it.
    // Example: 1/sqrt(x) from 0 to 1. Integral is 2. `quadrature` might struggle at x=0.
    #[test]
    fn test_difficult_integral_sqrt_inv() {
        let f = |x: f64| if x == 0.0 { 0.0 } else { 1.0 / x.sqrt() }; // Handle singularity for direct eval
        let a = 0.0;
        let b = 1.0;
        let epsabs = 1e-6; // Looser tolerance
        let epsrel = 1e-6;
        let limit = 50;

        let mut result = 0.0;
        let mut abserr = 0.0;
        let mut ier = 0;

        rust_dqagse(f, a, b, epsabs, epsrel, limit, &mut result, &mut abserr, &mut ier);
        
        // `quadrature` might handle this due to double exponential method's strengths
        // near singularities if they are at the interval ends.
        // Expected result is 2.0
        if ier == status::STATUS_SUCCESS {
            assert!((result - 2.0).abs() < epsabs.max(epsrel * result.abs()));
            assert!(abserr < epsabs.max(epsrel * result.abs()));
        } else {
            // If it didn't succeed, check it's a convergence error
            assert_eq!(ier, status::STATUS_NO_CONVERGE);
        }
        // This test is more about seeing how the wrapper handles results from `quadrature`
        // when the integral is challenging.
    }
}
