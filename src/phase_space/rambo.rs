use std::f64::consts::PI;

use rand::Rng;

use crate::kinematics::FourMomentum;

/// RAMBO: Random Momenta Beautifully Organized.
///
/// Generates n isotropic massless four-momenta, then boosts and rescales
/// for the desired sqrt(s) and particle masses.
///
/// Returns (momenta, weight) where weight is the phase space weight
/// for importance sampling.
///
/// Reference: S.D. Ellis, R. Kleiss, W.J. Stirling, Comp. Phys. Comm. 40 (1986) 359.
pub fn generate_rambo(
    rng: &mut impl Rng,
    sqrt_s: f64,
    masses: &[f64],
) -> Option<(Vec<FourMomentum>, f64)> {
    let n = masses.len();
    if n < 2 {
        return None;
    }

    let mass_sum: f64 = masses.iter().sum();
    if sqrt_s <= mass_sum {
        return None;
    }

    // Step 1: Generate n isotropic massless four-momenta
    let mut q: Vec<FourMomentum> = Vec::with_capacity(n);
    for _ in 0..n {
        let cos_theta: f64 = rng.random_range(-1.0..1.0);
        let phi: f64 = rng.random_range(0.0..(2.0 * PI));
        let r1: f64 = rng.random_range(0.0_f64..1.0);
        let r2: f64 = rng.random_range(0.0_f64..1.0);
        let energy = -(r1 * r2).ln();

        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        q.push(FourMomentum::new(
            energy,
            energy * sin_theta * phi.cos(),
            energy * sin_theta * phi.sin(),
            energy * cos_theta,
        ));
    }

    // Step 2: Compute total momentum Q
    let total: FourMomentum = q.iter().copied().fold(FourMomentum::zero(), |a, b| a + b);
    let total_mass = total.invariant_mass();

    // Step 3: Boost all momenta to CM frame and rescale to desired sqrt(s)
    let beta = total.beta();
    let scale = sqrt_s / total_mass;

    let mut p: Vec<FourMomentum> = Vec::with_capacity(n);
    for qi in &q {
        let boosted = qi.boost(&beta);
        p.push(FourMomentum::new(
            boosted.e * scale,
            boosted.px * scale,
            boosted.py * scale,
            boosted.pz * scale,
        ));
    }

    // Step 4: If all massless, we're done. Otherwise, rescale for masses.
    let all_massless = masses.iter().all(|&m| m == 0.0);
    let weight;

    if all_massless {
        weight = massless_weight(n, sqrt_s);
    } else {
        // Massive RAMBO: iterative rescaling
        let result = massive_rescale(&mut p, sqrt_s, masses);
        if result.is_none() {
            return None;
        }
        weight = result.unwrap();
    }

    Some((p, weight))
}

/// Massless n-body phase space weight (unnormalized).
fn massless_weight(n: usize, sqrt_s: f64) -> f64 {
    let n = n as f64;
    let s = sqrt_s * sqrt_s;
    // Phase space volume factor
    (PI / 2.0).powf(n - 1.0) * s.powf(n - 2.0)
        / (factorial((n as u32) - 1) * factorial((n as u32) - 2)) as f64
        * (2.0 * PI).powf(4.0 - 3.0 * n)
}

/// Rescale massless momenta to massive ones using the iterative method.
/// Returns the weight correction factor.
fn massive_rescale(
    momenta: &mut Vec<FourMomentum>,
    sqrt_s: f64,
    masses: &[f64],
) -> Option<f64> {
    let n = momenta.len();

    // Find the rescaling parameter xi such that
    // sum_i sqrt(|xi*p_i|^2 + m_i^2) = sqrt(s)
    // Use Newton's method
    let mut xi = 1.0;
    for _ in 0..100 {
        let mut f = -sqrt_s;
        let mut df = 0.0;
        for i in 0..n {
            let p_mag_sq = momenta[i].px * momenta[i].px
                + momenta[i].py * momenta[i].py
                + momenta[i].pz * momenta[i].pz;
            let ei = (xi * xi * p_mag_sq + masses[i] * masses[i]).sqrt();
            f += ei;
            if ei > 1e-30 {
                df += xi * p_mag_sq / ei;
            }
        }
        if df.abs() < 1e-30 {
            return None;
        }
        let delta = -f / df;
        xi += delta;
        if delta.abs() < 1e-12 * xi.abs() {
            break;
        }
    }

    // Apply the rescaling
    let mut sum_p_over_e = 0.0;
    let mut product_p_over_e = 1.0;

    for i in 0..n {
        let px = xi * momenta[i].px;
        let py = xi * momenta[i].py;
        let pz = xi * momenta[i].pz;
        let p_mag = (px * px + py * py + pz * pz).sqrt();
        let ei = (p_mag * p_mag + masses[i] * masses[i]).sqrt();
        momenta[i] = FourMomentum::new(ei, px, py, pz);

        if ei > 1e-30 {
            let ratio = p_mag / ei;
            sum_p_over_e += p_mag * p_mag / ei;
            product_p_over_e *= ratio;
        }
    }

    // Weight correction for massive phase space
    let massless_w = massless_weight(n, sqrt_s);
    let weight = massless_w * product_p_over_e * (sqrt_s / sum_p_over_e).powi(n as i32 - 1)
        * xi.powi(2 * n as i32 - 3);

    Some(weight)
}

fn factorial(n: u32) -> u64 {
    (1..=n as u64).product()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rambo_massless_2body() {
        let mut rng = rand::rng();
        let sqrt_s = 100.0;
        let masses = [0.0, 0.0];

        let (momenta, weight) = generate_rambo(&mut rng, sqrt_s, &masses).unwrap();
        assert_eq!(momenta.len(), 2);
        assert!(weight > 0.0);

        // Check energy-momentum conservation
        let total = momenta[0] + momenta[1];
        assert!((total.e - sqrt_s).abs() < 1e-8);
        assert!(total.px.abs() < 1e-8);
        assert!(total.py.abs() < 1e-8);
        assert!(total.pz.abs() < 1e-8);
    }

    #[test]
    fn test_rambo_massive_2body() {
        let mut rng = rand::rng();
        let sqrt_s = 200.0;
        let masses = [80.377, 80.377]; // W+W-

        let (momenta, weight) = generate_rambo(&mut rng, sqrt_s, &masses).unwrap();
        assert_eq!(momenta.len(), 2);
        assert!(weight > 0.0);

        let total = momenta[0] + momenta[1];
        assert!((total.e - sqrt_s).abs() < 1e-6, "E = {}", total.e);
        assert!(total.px.abs() < 1e-6);
        assert!(total.py.abs() < 1e-6);
        assert!(total.pz.abs() < 1e-6);

        // Check on-shell
        for (i, p) in momenta.iter().enumerate() {
            assert!(
                (p.invariant_mass() - masses[i]).abs() < 1e-4,
                "mass[{}] = {}, expected {}",
                i,
                p.invariant_mass(),
                masses[i]
            );
        }
    }

    #[test]
    fn test_rambo_below_threshold() {
        let mut rng = rand::rng();
        assert!(generate_rambo(&mut rng, 100.0, &[80.0, 80.0]).is_none());
    }

    #[test]
    fn test_rambo_3body() {
        let mut rng = rand::rng();
        let sqrt_s = 500.0;
        let masses = [0.0, 0.0, 0.0];

        let (momenta, weight) = generate_rambo(&mut rng, sqrt_s, &masses).unwrap();
        assert_eq!(momenta.len(), 3);
        assert!(weight > 0.0);

        let total = momenta.iter().copied().fold(FourMomentum::zero(), |a, b| a + b);
        assert!((total.e - sqrt_s).abs() < 1e-6);
        assert!(total.px.abs() < 1e-6);
        assert!(total.py.abs() < 1e-6);
        assert!(total.pz.abs() < 1e-6);
    }
}
