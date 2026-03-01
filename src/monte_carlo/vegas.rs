use rand::Rng;

/// Result of a Vegas integration.
#[derive(Debug, Clone)]
pub struct VegasResult {
    pub integral: f64,
    pub error: f64,
    pub chi_squared_per_dof: f64,
    pub iterations_completed: usize,
}

/// Vegas adaptive Monte Carlo integrator.
///
/// Integrates a function f: [0,1]^d -> f64 using importance sampling
/// with an adaptive grid. The grid is refined after each iteration to
/// concentrate sampling points where the integrand is largest.
pub struct Vegas {
    n_dim: usize,
    n_bins: usize,
    n_iterations: usize,
    n_calls_per_iter: usize,
    alpha: f64,
    grid: Vec<Vec<f64>>,
}

impl Vegas {
    pub fn new(n_dim: usize) -> Self {
        let n_bins = 50;
        let grid = (0..n_dim)
            .map(|_| (0..=n_bins).map(|i| i as f64 / n_bins as f64).collect())
            .collect();

        Self {
            n_dim,
            n_bins,
            n_iterations: 10,
            n_calls_per_iter: 10000,
            alpha: 1.5,
            grid,
        }
    }

    pub fn with_calls(mut self, n: usize) -> Self {
        self.n_calls_per_iter = n;
        self
    }

    pub fn with_iterations(mut self, n: usize) -> Self {
        self.n_iterations = n;
        self
    }

    pub fn with_bins(mut self, n: usize) -> Self {
        self.n_bins = n;
        self.grid = (0..self.n_dim)
            .map(|_| (0..=n).map(|i| i as f64 / n as f64).collect())
            .collect();
        self
    }

    /// Integrate the function f over [0,1]^d.
    pub fn integrate<F>(&mut self, f: F, rng: &mut impl Rng) -> VegasResult
    where
        F: Fn(&[f64]) -> f64,
    {
        let mut results: Vec<(f64, f64)> = Vec::new(); // (integral, variance)
        let mut x = vec![0.0; self.n_dim];

        for _iter in 0..self.n_iterations {
            let mut sum = 0.0;
            let mut sum2 = 0.0;
            let mut bin_weights: Vec<Vec<f64>> =
                (0..self.n_dim).map(|_| vec![0.0; self.n_bins]).collect();

            for _ in 0..self.n_calls_per_iter {
                let mut jacobian = 1.0;

                for d in 0..self.n_dim {
                    // Sample a random bin
                    let r: f64 = rng.random_range(0.0..1.0);
                    let bin_f = r * self.n_bins as f64;
                    let bin = (bin_f as usize).min(self.n_bins - 1);
                    let bin_frac = bin_f - bin as f64;

                    // Map to x in [0, 1] using the grid
                    let lo = self.grid[d][bin];
                    let hi = self.grid[d][bin + 1];
                    let dx = hi - lo;
                    x[d] = lo + bin_frac * dx;
                    jacobian *= dx * self.n_bins as f64;
                }

                let fval = f(&x) * jacobian;
                sum += fval;
                sum2 += fval * fval;

                // Accumulate bin weights for grid adaptation
                for d in 0..self.n_dim {
                    let r: f64 = x[d];
                    // Find which bin this x falls in
                    let mut bin = 0;
                    for b in 0..self.n_bins {
                        if r >= self.grid[d][b] && r < self.grid[d][b + 1] {
                            bin = b;
                            break;
                        }
                    }
                    bin_weights[d][bin] += fval * fval;
                }
            }

            let n = self.n_calls_per_iter as f64;
            let mean = sum / n;
            let variance = (sum2 / n - mean * mean) / (n - 1.0);

            results.push((mean, variance));

            // Adapt the grid
            self.adapt_grid(&bin_weights);
        }

        // Combine results from all iterations
        let weighted_sum;
        let mut weight_total = 0.0;
        let mut chi2 = 0.0;

        // First pass: compute weighted average
        let mut combined_integral = 0.0;
        let mut combined_weight = 0.0;
        for &(mean, var) in &results {
            if var > 0.0 {
                let w = 1.0 / var;
                combined_integral += mean * w;
                combined_weight += w;
            }
        }

        if combined_weight > 0.0 {
            weighted_sum = combined_integral / combined_weight;
            weight_total = combined_weight;

            // Chi-squared
            for &(mean, var) in &results {
                if var > 0.0 {
                    let diff = mean - weighted_sum;
                    chi2 += diff * diff / var;
                }
            }
        } else {
            // Fallback: simple average
            weighted_sum = results.iter().map(|(m, _)| m).sum::<f64>() / results.len() as f64;
        }

        let error = if weight_total > 0.0 {
            (1.0 / weight_total).sqrt()
        } else {
            f64::INFINITY
        };

        let dof = (results.len() as f64 - 1.0).max(1.0);

        VegasResult {
            integral: weighted_sum,
            error,
            chi_squared_per_dof: chi2 / dof,
            iterations_completed: self.n_iterations,
        }
    }

    fn adapt_grid(&mut self, bin_weights: &[Vec<f64>]) {
        for d in 0..self.n_dim {
            // Smooth the bin weights
            let mut smoothed = vec![0.0; self.n_bins];
            for i in 0..self.n_bins {
                let prev = if i > 0 { bin_weights[d][i - 1] } else { 0.0 };
                let curr = bin_weights[d][i];
                let next = if i + 1 < self.n_bins {
                    bin_weights[d][i + 1]
                } else {
                    0.0
                };
                smoothed[i] = (prev + 6.0 * curr + next) / 8.0;
            }

            // Compute cumulative distribution
            let total: f64 = smoothed.iter().sum();
            if total <= 0.0 {
                continue;
            }

            // Apply damping
            let mut importance = vec![0.0; self.n_bins];
            for i in 0..self.n_bins {
                let r = smoothed[i] / total;
                if r > 0.0 {
                    importance[i] = ((r - 1.0) / r.ln()).powf(self.alpha);
                }
            }

            let imp_total: f64 = importance.iter().sum();
            if imp_total <= 0.0 {
                continue;
            }

            // Redistribute grid points according to importance
            let mut new_grid = vec![0.0; self.n_bins + 1];
            new_grid[0] = 0.0;
            new_grid[self.n_bins] = 1.0;

            let mut cumul = 0.0;
            let mut old_bin = 0;
            let target_per_bin = imp_total / self.n_bins as f64;

            for new_bin in 1..self.n_bins {
                let target = new_bin as f64 * target_per_bin;

                while old_bin < self.n_bins && cumul + importance[old_bin] < target {
                    cumul += importance[old_bin];
                    old_bin += 1;
                }

                if old_bin >= self.n_bins {
                    new_grid[new_bin] = 1.0;
                    continue;
                }

                let remaining = target - cumul;
                let frac = if importance[old_bin] > 0.0 {
                    remaining / importance[old_bin]
                } else {
                    0.0
                };

                let lo = self.grid[d][old_bin];
                let hi = self.grid[d][old_bin + 1];
                new_grid[new_bin] = lo + frac * (hi - lo);
            }

            self.grid[d] = new_grid;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vegas_constant() {
        let mut vegas = Vegas::new(1).with_calls(5000).with_iterations(5);
        let mut rng = rand::rng();
        let result = vegas.integrate(|_| 3.14, &mut rng);
        assert!(
            (result.integral - 3.14).abs() < 0.01,
            "integral = {}",
            result.integral
        );
    }

    #[test]
    fn test_vegas_linear() {
        // integral of x from 0 to 1 = 0.5
        let mut vegas = Vegas::new(1).with_calls(10000).with_iterations(5);
        let mut rng = rand::rng();
        let result = vegas.integrate(|x| x[0], &mut rng);
        assert!(
            (result.integral - 0.5).abs() < 0.01,
            "integral = {}",
            result.integral
        );
    }

    #[test]
    fn test_vegas_gaussian_peak() {
        // Integral of exp(-10*(x-0.5)^2) from 0 to 1
        // This has a peak at x=0.5, testing grid adaptation
        let mut vegas = Vegas::new(1).with_calls(20000).with_iterations(10);
        let mut rng = rand::rng();
        let result = vegas.integrate(
            |x| (-10.0 * (x[0] - 0.5).powi(2)).exp(),
            &mut rng,
        );
        // Analytic: sqrt(pi/10) * erf(sqrt(10)/2) ~ 0.5551
        assert!(
            (result.integral - 0.5551).abs() < 0.02,
            "integral = {} +/- {}",
            result.integral,
            result.error
        );
    }

    #[test]
    fn test_vegas_2d() {
        // integral of x*y over [0,1]^2 = 0.25
        let mut vegas = Vegas::new(2).with_calls(10000).with_iterations(5);
        let mut rng = rand::rng();
        let result = vegas.integrate(|x| x[0] * x[1], &mut rng);
        assert!(
            (result.integral - 0.25).abs() < 0.01,
            "integral = {}",
            result.integral
        );
    }
}
