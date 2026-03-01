/// Gauss-Legendre quadrature nodes and weights for numerical integration.
/// Suitable for smooth integrands like 2->2 cross sections over cos(theta).

/// Integrate f(x) from a to b using Gauss-Legendre quadrature with n points.
pub fn gauss_legendre(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let (nodes, weights) = gauss_legendre_nodes_weights(n);

    let mid = (b + a) / 2.0;
    let half_range = (b - a) / 2.0;

    let mut result = 0.0;
    for i in 0..n {
        let x = mid + half_range * nodes[i];
        result += weights[i] * f(x);
    }
    result * half_range
}

/// Generate Gauss-Legendre nodes and weights for n points on [-1, 1].
fn gauss_legendre_nodes_weights(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];

    let m = (n + 1) / 2;

    for i in 0..m {
        // Initial guess for roots of P_n(x)
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();

        loop {
            // Evaluate P_n(z) and P_{n-1}(z) using recurrence
            let mut p0 = 1.0;
            let mut p1 = z;
            for j in 2..=n {
                let pn = ((2 * j - 1) as f64 * z * p1 - (j - 1) as f64 * p0) / j as f64;
                p0 = p1;
                p1 = pn;
            }
            // p1 is now P_n(z)
            // Derivative: P_n'(z) = n * (z * P_n(z) - P_{n-1}(z)) / (z^2 - 1)
            let dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
            let dz = p1 / dp;
            z -= dz;
            if dz.abs() < 1e-15 {
                break;
            }
        }

        nodes[i] = -z;
        nodes[n - 1 - i] = z;

        // Weight
        let mut p0 = 1.0;
        let mut p1 = z;
        for j in 2..=n {
            let pn = ((2 * j - 1) as f64 * z * p1 - (j - 1) as f64 * p0) / j as f64;
            p0 = p1;
            p1 = pn;
        }
        let dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
        let w = 2.0 / ((1.0 - z * z) * dp * dp);
        weights[i] = w;
        weights[n - 1 - i] = w;
    }

    (nodes, weights)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integrate_constant() {
        let result = gauss_legendre(|_| 1.0, 0.0, 1.0, 4);
        assert!((result - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_integrate_linear() {
        let result = gauss_legendre(|x| x, 0.0, 1.0, 4);
        assert!((result - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_integrate_polynomial() {
        // x^4 from 0 to 1 = 0.2. Gauss-Legendre with n >= 3 should be exact.
        let result = gauss_legendre(|x| x.powi(4), 0.0, 1.0, 10);
        assert!((result - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_integrate_trig() {
        // sin(x) from 0 to pi = 2.0
        let result = gauss_legendre(|x| x.sin(), 0.0, std::f64::consts::PI, 20);
        assert!((result - 2.0).abs() < 1e-10);
    }
}
