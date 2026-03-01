use num_complex::Complex64;

/// CKM quark mixing matrix.
///
/// Matrix elements V_{ij} where i = (u, c, t) and j = (d, s, b).
/// Uses PDG 2024 central values.
pub struct CkmMatrix {
    /// elements[i][j]: i = up-type generation (0=u, 1=c, 2=t), j = down-type (0=d, 1=s, 2=b)
    elements: [[Complex64; 3]; 3],
}

/// Generation index for CKM matrix.
#[derive(Debug, Clone, Copy)]
pub enum Generation {
    First = 0,
    Second = 1,
    Third = 2,
}

impl CkmMatrix {
    /// Standard CKM matrix with PDG 2024 central values.
    /// Using the Wolfenstein parameterization magnitudes.
    pub fn standard() -> Self {
        // Magnitudes (PDG 2024):
        // |V_ud| = 0.97435  |V_us| = 0.22500  |V_ub| = 0.00369
        // |V_cd| = 0.22486  |V_cs| = 0.97349  |V_cb| = 0.04182
        // |V_td| = 0.00857  |V_ts| = 0.04110  |V_tb| = 0.999118
        //
        // For simplicity, we use real values (CP phase ignored for rate calculations).
        Self {
            elements: [
                [
                    Complex64::new(0.97435, 0.0),
                    Complex64::new(0.22500, 0.0),
                    Complex64::new(0.00369, 0.0),
                ],
                [
                    Complex64::new(0.22486, 0.0),
                    Complex64::new(0.97349, 0.0),
                    Complex64::new(0.04182, 0.0),
                ],
                [
                    Complex64::new(0.00857, 0.0),
                    Complex64::new(0.04110, 0.0),
                    Complex64::new(0.999118, 0.0),
                ],
            ],
        }
    }

    /// Get the CKM matrix element V_{ij}.
    pub fn element(&self, up_gen: Generation, down_gen: Generation) -> Complex64 {
        self.elements[up_gen as usize][down_gen as usize]
    }

    /// Get |V_{ij}|^2.
    pub fn abs_squared(&self, up_gen: Generation, down_gen: Generation) -> f64 {
        self.element(up_gen, down_gen).norm_sqr()
    }
}

/// Global CKM matrix instance.
pub fn ckm() -> CkmMatrix {
    CkmMatrix::standard()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unitarity_rows() {
        let v = CkmMatrix::standard();
        // Each row should have sum |V_ij|^2 ~ 1
        for i in [Generation::First, Generation::Second, Generation::Third] {
            let sum: f64 = [Generation::First, Generation::Second, Generation::Third]
                .iter()
                .map(|&j| v.abs_squared(i, j))
                .sum();
            assert!(
                (sum - 1.0).abs() < 0.01,
                "Row {:?} unitarity: sum = {}",
                i as u8,
                sum
            );
        }
    }

    #[test]
    fn test_unitarity_columns() {
        let v = CkmMatrix::standard();
        for j in [Generation::First, Generation::Second, Generation::Third] {
            let sum: f64 = [Generation::First, Generation::Second, Generation::Third]
                .iter()
                .map(|&i| v.abs_squared(i, j))
                .sum();
            assert!(
                (sum - 1.0).abs() < 0.01,
                "Column {:?} unitarity: sum = {}",
                j as u8,
                sum
            );
        }
    }

    #[test]
    fn test_vud() {
        let v = CkmMatrix::standard();
        let vud = v.element(Generation::First, Generation::First).norm();
        assert!((vud - 0.97435).abs() < 0.001);
    }
}
