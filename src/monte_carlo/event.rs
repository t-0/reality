use serde::{Deserialize, Serialize};

/// Status of a particle in an event record.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ParticleStatus {
    Incoming,
    Outgoing,
    Intermediate,
}

/// A particle in the event record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Particle {
    pub pdg_id: i32,
    pub status: ParticleStatus,
    pub momentum: [f64; 4], // [E, px, py, pz]
    pub mass: f64,
}

/// A complete Monte Carlo event.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Event {
    pub event_number: u64,
    pub weight: f64,
    pub sqrt_s: f64,
    pub particles: Vec<Particle>,
}

impl Event {
    /// Check four-momentum conservation: sum(incoming) == sum(outgoing).
    pub fn check_momentum_conservation(&self) -> f64 {
        let mut sum = [0.0_f64; 4];
        for p in &self.particles {
            let sign = match p.status {
                ParticleStatus::Incoming => 1.0,
                ParticleStatus::Outgoing => -1.0,
                ParticleStatus::Intermediate => 0.0,
            };
            for i in 0..4 {
                sum[i] += sign * p.momentum[i];
            }
        }
        (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2] + sum[3] * sum[3]).sqrt()
    }
}

/// A sample of generated events with metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EventSample {
    pub process: String,
    pub sqrt_s: f64,
    pub cross_section_nb: f64,
    pub cross_section_error_nb: f64,
    pub n_events: usize,
    pub events: Vec<Event>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::couplings::constants::*;

    #[test]
    fn test_event_serialization() {
        let event = Event {
            event_number: 1,
            weight: 1.0,
            sqrt_s: 91.2,
            particles: vec![
                Particle {
                    pdg_id: -11,
                    status: ParticleStatus::Incoming,
                    momentum: [45.6, 0.0, 0.0, 45.6],
                    mass: M_ELECTRON,
                },
                Particle {
                    pdg_id: 11,
                    status: ParticleStatus::Incoming,
                    momentum: [45.6, 0.0, 0.0, -45.6],
                    mass: M_ELECTRON,
                },
                Particle {
                    pdg_id: -13,
                    status: ParticleStatus::Outgoing,
                    momentum: [45.6, 30.0, 20.0, 10.0],
                    mass: M_MUON,
                },
                Particle {
                    pdg_id: 13,
                    status: ParticleStatus::Outgoing,
                    momentum: [45.6, -30.0, -20.0, -10.0],
                    mass: M_MUON,
                },
            ],
        };

        let json = serde_json::to_string(&event).unwrap();
        let recovered: Event = serde_json::from_str(&json).unwrap();
        assert_eq!(recovered.event_number, 1);
        assert_eq!(recovered.particles.len(), 4);
    }

    #[test]
    fn test_momentum_conservation_check() {
        let event = Event {
            event_number: 1,
            weight: 1.0,
            sqrt_s: 10.0,
            particles: vec![
                Particle {
                    pdg_id: -11,
                    status: ParticleStatus::Incoming,
                    momentum: [5.0, 0.0, 0.0, 5.0],
                    mass: 0.0,
                },
                Particle {
                    pdg_id: 11,
                    status: ParticleStatus::Incoming,
                    momentum: [5.0, 0.0, 0.0, -5.0],
                    mass: 0.0,
                },
                Particle {
                    pdg_id: -13,
                    status: ParticleStatus::Outgoing,
                    momentum: [5.0, 3.0, 4.0, 0.0],
                    mass: 0.0,
                },
                Particle {
                    pdg_id: 13,
                    status: ParticleStatus::Outgoing,
                    momentum: [5.0, -3.0, -4.0, 0.0],
                    mass: 0.0,
                },
            ],
        };

        let residual = event.check_momentum_conservation();
        assert!(residual < 1e-10, "residual = {}", residual);
    }
}
