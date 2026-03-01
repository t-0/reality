use crate::particles::PdgId;

/// A single decay channel with its partial width and branching ratio.
#[derive(Debug, Clone)]
pub struct DecayChannelEntry {
    pub final_state: Vec<PdgId>,
    pub label: String,
    pub partial_width: f64,
    pub branching_ratio: f64,
}

/// A complete decay table for a particle.
#[derive(Debug, Clone)]
pub struct DecayTable {
    pub parent: PdgId,
    pub parent_name: String,
    pub channels: Vec<DecayChannelEntry>,
    pub total_width: f64,
}

impl DecayTable {
    /// Create a decay table from a list of channels (partial widths).
    /// Branching ratios are computed automatically.
    pub fn from_channels(
        parent: PdgId,
        parent_name: &str,
        channels: Vec<(Vec<PdgId>, String, f64)>,
    ) -> Self {
        let total_width: f64 = channels.iter().map(|(_, _, w)| w).sum();

        let entries = channels
            .into_iter()
            .map(|(fs, label, pw)| DecayChannelEntry {
                final_state: fs,
                label,
                partial_width: pw,
                branching_ratio: if total_width > 0.0 {
                    pw / total_width
                } else {
                    0.0
                },
            })
            .collect();

        Self {
            parent,
            parent_name: parent_name.to_string(),
            channels: entries,
            total_width,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_branching_ratios_sum_to_one() {
        let table = DecayTable::from_channels(
            PdgId::Z_BOSON,
            "Z",
            vec![
                (vec![], "a".into(), 1.0),
                (vec![], "b".into(), 2.0),
                (vec![], "c".into(), 3.0),
            ],
        );

        let br_sum: f64 = table.channels.iter().map(|c| c.branching_ratio).sum();
        assert!((br_sum - 1.0).abs() < 1e-10);
        assert!((table.total_width - 6.0).abs() < 1e-10);
    }
}
