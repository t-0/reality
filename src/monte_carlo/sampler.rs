use std::f64::consts::PI;

use rand::Rng;

use super::event::{Event, EventSample, Particle, ParticleStatus};
use crate::cross_section::total::total_cross_section_2to2;
use crate::kinematics::invariants::Mandelstam;
use crate::matrix_elements::MatrixElement;
use crate::phase_space::two_body;

/// Accept/reject Monte Carlo event generator for 2->2 processes.
pub struct EventGenerator<'a, M: MatrixElement> {
    matrix_element: &'a M,
    sqrt_s: f64,
    max_weight: f64,
}

impl<'a, M: MatrixElement> EventGenerator<'a, M> {
    pub fn new(me: &'a M, sqrt_s: f64) -> Self {
        Self {
            matrix_element: me,
            sqrt_s,
            max_weight: 0.0,
        }
    }

    /// Warmup phase: sample phase space to discover the maximum weight.
    pub fn warmup(&mut self, rng: &mut impl Rng, n_points: usize) {
        let s = self.sqrt_s * self.sqrt_s;
        let (m3, m4) = self.matrix_element.final_masses();

        for _ in 0..n_points {
            let cos_theta: f64 = rng.random_range(-1.0..1.0);
            let mandelstam =
                Mandelstam::from_s_costheta(s, cos_theta, [0.0, 0.0], [m3, m4]);
            let m2 = self.matrix_element.amplitude_squared(s, mandelstam.t);
            if m2 > self.max_weight {
                self.max_weight = m2;
            }
        }
        // Add 10% safety margin
        self.max_weight *= 1.1;
    }

    /// Generate a single unweighted event using accept/reject.
    pub fn generate_event(&self, rng: &mut impl Rng, event_number: u64) -> Event {
        let s = self.sqrt_s * self.sqrt_s;
        let (m3, m4) = self.matrix_element.final_masses();
        let initial = self.matrix_element.initial_state();
        let final_state = self.matrix_element.final_state();

        loop {
            let cos_theta: f64 = rng.random_range(-1.0..1.0);
            let phi: f64 = rng.random_range(0.0..(2.0 * PI));

            let mandelstam =
                Mandelstam::from_s_costheta(s, cos_theta, [0.0, 0.0], [m3, m4]);
            let m2 = self.matrix_element.amplitude_squared(s, mandelstam.t);

            // Accept/reject
            let r: f64 = rng.random_range(0.0..self.max_weight);
            if r < m2 {
                // Generate kinematics
                let (p3, p4) =
                    two_body::two_body_at_angle(self.sqrt_s, m3, m4, cos_theta, phi)
                        .unwrap();

                let half_e = self.sqrt_s / 2.0;
                return Event {
                    event_number,
                    weight: 1.0,
                    sqrt_s: self.sqrt_s,
                    particles: vec![
                        Particle {
                            pdg_id: initial[0].0,
                            status: ParticleStatus::Incoming,
                            momentum: [half_e, 0.0, 0.0, half_e],
                            mass: 0.0,
                        },
                        Particle {
                            pdg_id: initial[1].0,
                            status: ParticleStatus::Incoming,
                            momentum: [half_e, 0.0, 0.0, -half_e],
                            mass: 0.0,
                        },
                        Particle {
                            pdg_id: final_state[0].0,
                            status: ParticleStatus::Outgoing,
                            momentum: p3.as_array(),
                            mass: m3,
                        },
                        Particle {
                            pdg_id: final_state[1].0,
                            status: ParticleStatus::Outgoing,
                            momentum: p4.as_array(),
                            mass: m4,
                        },
                    ],
                };
            }
        }
    }

    /// Generate n unweighted events.
    pub fn generate_events(&self, rng: &mut impl Rng, n: usize) -> Vec<Event> {
        (0..n)
            .map(|i| self.generate_event(rng, i as u64))
            .collect()
    }

    /// Generate a complete EventSample with cross section calculation.
    pub fn generate_sample(&self, rng: &mut impl Rng, n_events: usize) -> EventSample {
        let sigma = total_cross_section_2to2(self.matrix_element, self.sqrt_s);
        let events = self.generate_events(rng, n_events);

        EventSample {
            process: self.matrix_element.name().to_string(),
            sqrt_s: self.sqrt_s,
            cross_section_nb: sigma,
            cross_section_error_nb: sigma * 0.01, // rough estimate
            n_events,
            events,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix_elements::ee_mumu::EpEmToMuPMuM;

    #[test]
    fn test_generate_events() {
        let me = EpEmToMuPMuM::qed_only();
        let mut gen = EventGenerator::new(&me, 91.2);
        let mut rng = rand::rng();

        gen.warmup(&mut rng, 10000);
        let events = gen.generate_events(&mut rng, 10);

        assert_eq!(events.len(), 10);
        for event in &events {
            assert_eq!(event.particles.len(), 4);
            let residual = event.check_momentum_conservation();
            assert!(
                residual < 1e-6,
                "Momentum not conserved: residual = {}",
                residual
            );
        }
    }

    #[test]
    fn test_generate_sample() {
        let me = EpEmToMuPMuM::qed_only();
        let mut gen = EventGenerator::new(&me, 200.0);
        let mut rng = rand::rng();

        gen.warmup(&mut rng, 5000);
        let sample = gen.generate_sample(&mut rng, 5);

        assert_eq!(sample.n_events, 5);
        assert!(sample.cross_section_nb > 0.0);
        assert_eq!(sample.events.len(), 5);
    }
}
