//! Validation against B. Bajc and G. Senjanovic, "Seesaw at LHC",
//! JHEP 0708 (2007) 014 [hep-ph/0612029].
//!
//! Non-supersymmetric SU(5) with an adjoint fermion 24_F. Below M_GUT the
//! running is affected by the fermion and scalar weak triplets (1,3,0), the
//! fermion and scalar colour octets (8,1,0), the fermionic leptoquark
//! (3,2,-5/6) + conjugate at m_(3,2) <~ M_GUT^2 / M_Pl, and the scalar colour
//! triplet (3,1,-1/3) at m_T. Their one-loop results, eqs. (15) and (16):
//!
//!   exp[30 pi (1/a1 - 1/a2)] = (M_G/M_Z)^84 ((m3F^4 m3B)/M_Z^5)^5 (M_G/m32)^20 (M_G/m_T)
//!   exp[20 pi (1/a1 - 1/a3)] = (M_G/M_Z)^86 ((m8F^4 m8B)/M_Z^5)^5 (M_G/m32)^20 (M_G/m_T)^-1
//!
//! Benchmark: m3F = m3B = M_Z and m_T = M_G give M_GUT ~ 10^15.5 GeV, and the
//! octet masses come out in the range 10^5 - 10^8 GeV.

use reality::couplings::constants::M_Z;
use reality::couplings::unification::inverse_couplings_mz;
use reality::model::field::{Field, Kind, Su3};
use reality::model::predict::predict;
use reality::model::running::Running;
use reality::model::sm::standard_model;
use reality::model::spec::TopYukawa;
use reality::model::two_loop::TwoLoop;
use std::f64::consts::PI;

/// Planck mass (GeV).
const M_PLANCK: f64 = 1.22e19;

fn fields(m3: f64, m8f: f64, m8b: f64, m32: f64, m_t: f64) -> Vec<Field> {
    let mut f = standard_model();
    f.extend([
        Field::new("fermion triplet", Kind::Weyl, Su3::Singlet, 3, 0.0, 1).with_mass(m3),
        Field::new("scalar triplet", Kind::RealScalar, Su3::Singlet, 3, 0.0, 1).with_mass(m3),
        Field::new("fermion octet", Kind::Weyl, Su3::Octet, 1, 0.0, 1).with_mass(m8f),
        Field::new("scalar octet", Kind::RealScalar, Su3::Octet, 1, 0.0, 1).with_mass(m8b),
        Field::new("fermion leptoquark", Kind::Dirac, Su3::Triplet, 2, -5.0 / 6.0, 1).with_mass(m32),
        Field::new("scalar colour triplet", Kind::ComplexScalar, Su3::Triplet, 1, -1.0 / 3.0, 1).with_mass(m_t),
    ]);
    f
}

/// Self-consistent one-loop M_GUT with m_(3,2) = M_GUT^2 / M_Pl and m_T = M_GUT.
/// The octets do not affect where alpha_1 meets alpha_2.
fn one_loop_m_gut(m3: f64) -> f64 {
    let inv = inverse_couplings_mz();
    let mut m_g = 10f64.powf(15.5);
    for _ in 0..50 {
        let running = Running::new(&fields(m3, 1e7, 1e7, m_g * m_g / M_PLANCK, m_g));
        m_g = running.crossing(inv, 0, 1, [0.0; 3]).unwrap();
    }
    m_g
}

#[test]
fn test_eq15_gut_scale() {
    // Their eq. (15) with m3F = m3B = M_Z, m_T = M_G, m32 = M_G^2 / M_Pl reduces to
    // 64 x + 20 ln(M_Pl / M_Z) = 30 pi (1/a1 - 1/a2), with x = ln(M_G / M_Z).
    let inv = inverse_couplings_mz();
    let x = (30.0 * PI * (inv[0] - inv[1]) - 20.0 * (M_PLANCK / M_Z).ln()) / 64.0;
    let analytic = M_Z * x.exp();

    let engine = one_loop_m_gut(M_Z);
    assert!(
        (engine / analytic).ln().abs() < 1e-6,
        "engine M_GUT = {:e}, eq. (15) = {:e}",
        engine,
        analytic
    );
    // The paper quotes M_GUT ~ 10^15.5 GeV.
    let log = engine.log10();
    assert!(log > 15.3 && log < 15.6, "log10 M_GUT = {}", log);
}

#[test]
fn test_eq16_octet_mass_and_alpha_s() {
    // Their benchmark uses 1/a3(M_Z) = 8.55. Solve eq. (16) for a common octet
    // mass m8 = m8F = m8B, then check the engine predicts 1/a3 = 8.55 back.
    let inv = inverse_couplings_mz();
    let target_inv_a3 = 8.55;
    let m_g = one_loop_m_gut(M_Z);
    let m32 = m_g * m_g / M_PLANCK;
    // 20 pi (1/a1 - 1/a3) = 86 x + 25 ln(m8/M_Z) + 20 ln(M_G/m32) - ln(M_G/m_T), m_T = M_G.
    let x = (m_g / M_Z).ln();
    let y = (20.0 * PI * (inv[0] - target_inv_a3) - 86.0 * x - 20.0 * (m_g / m32).ln()) / 25.0;
    let m8 = M_Z * y.exp();
    assert!(m8 > 1e5 && m8 < 1e8, "octet mass {:e} outside the paper's 1e5 - 1e8 GeV", m8);

    let running = Running::new(&fields(M_Z, m8, m8, m32, m_g));
    let p = predict(&running, (0, 1), 2, [0.0; 3]).unwrap();
    assert!(
        (p.inverse_predicted - target_inv_a3).abs() < 1e-6,
        "engine 1/a3 = {}, paper 8.55",
        p.inverse_predicted
    );
    assert!((p.m_gut / m_g - 1.0).abs() < 1e-9);
}

#[test]
fn test_octet_range_from_scalar_octet_mass() {
    // The paper's 1e5 - 1e8 GeV range for the fermion octet: a heavy scalar octet
    // (at M_GUT) needs the lightest fermion octet, a degenerate one the heaviest.
    let inv = inverse_couplings_mz();
    let m_g = one_loop_m_gut(M_Z);
    let m32 = m_g * m_g / M_PLANCK;
    let x = (m_g / M_Z).ln();
    let rhs = 20.0 * PI * (inv[0] - 8.55) - 86.0 * x - 20.0 * (m_g / m32).ln();
    // m8B = M_G: 20 ln(m8F/M_Z) + 5 x = rhs.
    let light = M_Z * ((rhs - 5.0 * x) / 20.0).exp();
    // m8B = m8F: 25 ln(m8/M_Z) = rhs.
    let heavy = M_Z * (rhs / 25.0).exp();
    assert!(light > 1e5 && light < 1e6, "{:e}", light);
    assert!(heavy > 1e7 && heavy < 1e8, "{:e}", heavy);
}

#[test]
fn test_two_loop_triplet_mass() {
    // "The two loop effects relax this somewhat and for the above example of the
    // GUT scale the triplet mass increases to about 500 GeV." Find the triplet
    // mass that gives the one-loop M_GUT at two loops.
    let inv = inverse_couplings_mz();
    let target = one_loop_m_gut(M_Z).ln();
    let m8 = 1e7;
    let crossing_at = |m3: f64| {
        let mut m_g = target.exp();
        for _ in 0..30 {
            let two = TwoLoop::new(&fields(m3, m8, m8, m_g * m_g / M_PLANCK, m_g), TopYukawa::Auto);
            let p = TwoLoop::crossing(&two.trajectory(inv), 0, 1, [0.0; 3]).unwrap();
            m_g = p.ln_q.exp();
        }
        m_g.ln()
    };
    // Bisect in ln(m3): a heavier triplet lowers M_GUT.
    let (mut lo, mut hi) = (M_Z.ln(), 1e5f64.ln());
    assert!(crossing_at(lo.exp()) > target && crossing_at(hi.exp()) < target);
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if crossing_at(mid.exp()) > target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let m3 = (0.5 * (lo + hi)).exp();
    eprintln!("two-loop triplet mass for the one-loop M_GUT: {:.0} GeV (paper: about 500 GeV)", m3);
    assert!(m3 > 200.0 && m3 < 1500.0, "two-loop triplet mass {:.0} GeV", m3);
}
