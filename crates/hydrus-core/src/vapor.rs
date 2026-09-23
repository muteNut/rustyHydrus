//! Coupled water vapor flow physics (port of ConVapor, VaporContent, and xLatent).
//!
//! Evaluates isothermal vapor hydraulic conductivity (ConVh), thermal vapor
//! conductivity (ConVT), thermal liquid conductivity (ConLT), volumetric vapor
//! content (ThetaV), and latent heat of vaporization (Latent).

pub const DIFF0: f64 = 2.12e-5; // Diffusivity of water vapor in air [m2/s]
pub const G: f64 = 9.81; // Gravitational acceleration [m/s2]
pub const X_MOL: f64 = 0.018015; // Molecular weight of water [kg/mol]
pub const R_GAS: f64 = 8.314; // Universal gas constant [J/mol/K]
pub const FC: f64 = 0.02; // Mass fraction of clay in soil
pub const GAMMA0: f64 = 71.89; // Reference surface tension [g/s2]
pub const GWT: f64 = 7.0; // Gain factor [-]

/// Volumetric latent heat of vaporization of water [J/m3]
pub fn latent_heat_volumetric(temp_c: f64) -> f64 {
    let row = (1.0 - 7.37e-6 * (temp_c - 4.0).powi(2) + 3.79e-8 * (temp_c - 4.0).powi(3)) * 1000.0;
    let x_lw = 2.501e6 - 2369.2 * temp_c; // [J/kg]
    row * x_lw
}

/// Saturated vapor density [kg/m3]
pub fn saturated_vapor_density(t_kelv: f64) -> f64 {
    0.001 * (31.3716 - 6014.79 / t_kelv - 0.00792495 * t_kelv).exp() / t_kelv
}

/// Relative humidity from pressure head h [m] and temperature T [K]
pub fn relative_humidity(h_m: f64, t_kelv: f64) -> f64 {
    if h_m >= 0.0 {
        1.0
    } else {
        (h_m * X_MOL * G / R_GAS / t_kelv).exp()
    }
}

/// ConVapor: Compute ConLT, ConVT, and ConVh across the profile
///
/// * `h_new`: Pressure head [L]
/// * `temp`: Temperature [°C]
/// * `con`: Liquid hydraulic conductivity [L/T]
/// * `theta`: Volumetric water content [-]
/// * `ths`: Saturated water content [-]
/// * `x_conv`: Length conversion to meters (e.g. 100 for cm)
/// * `t_conv`: Time conversion to seconds (e.g. 1/86400 for days)
/// * `l_vapor`: Flag indicating if vapor transport is active
/// * `i_enhanc`: 1 = Campbell enhancement factor, 0 = 1.0
pub fn con_vapor(
    n: usize,
    mat: &[usize],
    h_new: &[f64],
    temp: &[f64],
    con: &[f64],
    theta: &[f64],
    ths: &[f64],
    con_lt: &mut [f64],
    con_vt: &mut [f64],
    con_vh: &mut [f64],
    x_conv: f64,
    t_conv: f64,
    l_vapor: bool,
    i_enhanc: i32,
) {
    for i in 0..n {
        let h = h_new[i] / x_conv; // [m]
        let con_lh = con[i] / x_conv * t_conv; // [m/s]
        let t = temp[i];
        let m = mat[i];
        let theta_s = ths[m];

        // Surface tension and derivative [g/s2]
        let d_gamma = -0.1425 - 0.000479 * t;
        con_lt[i] = 0.0;
        if h < 0.0 {
            con_lt[i] = con_lh * h * GWT * d_gamma / GAMMA0;
        }

        if l_vapor {
            let t_kelv = t + 273.15;
            let diff_t = DIFF0 * (t_kelv / 273.15).powi(2);
            let theta_a = (theta_s - theta[i]).max(0.0);
            let tau = if theta_a > 0.0 {
                theta_a.powf(7.0 / 3.0) / (theta_s * theta_s)
            } else {
                0.0
            };
            let diff = tau * theta_a * diff_t;
            let row = (1.0 - 7.37e-6 * (t - 4.0).powi(2) + 3.79e-8 * (t - 4.0).powi(3)) * 1000.0;
            let rovs = saturated_vapor_density(t_kelv);
            let hr = relative_humidity(h, t_kelv);

            // Isothermal vapor conductivity [m/s]
            con_vh[i] = diff / row * rovs * X_MOL * G / R_GAS / t_kelv * hr;

            // Thermal vapor conductivity [m2/s/K]
            let t_kelv1 = t_kelv + 1.0;
            let rovs1 = saturated_vapor_density(t_kelv1);
            let drovs = rovs1 - rovs;
            let mut eta = 1.0;
            if i_enhanc == 1 && theta_s > 0.0 {
                let sr = theta[i] / theta_s;
                eta = 9.5 + 3.0 * sr - 8.5 * (-((1.0 + 2.6 / FC.sqrt()) * sr).powi(4)).exp();
            }
            con_vt[i] = diff / row * eta * hr * drovs;

            // Convert to HYDRUS project units
            con_vh[i] = con_vh[i] * x_conv / t_conv;
            con_vt[i] = con_vt[i] * x_conv * x_conv / t_conv;
        }
        con_lt[i] = con_lt[i] * x_conv * x_conv / t_conv;
    }
}

/// VaporContent: Calculate equivalent volumetric vapor content and adjust capacity
pub fn vapor_content(
    n: usize,
    mat: &[usize],
    theta: &[f64],
    theta_v: &mut [f64],
    temp: &[f64],
    h_new: &[f64],
    ths: &[f64],
    cap: &mut [f64],
    x_conv: f64,
) {
    for i in 0..n {
        let h = h_new[i] / x_conv;
        let m = mat[i];
        let t_kelv = temp[i] + 273.15;
        let rovs = saturated_vapor_density(t_kelv);
        let hr = relative_humidity(h, t_kelv);
        let rov = rovs * hr;
        let row = (1.0 - 7.37e-6 * (temp[i] - 4.0).powi(2) + 3.79e-8 * (temp[i] - 4.0).powi(3)) * 1000.0;

        theta_v[i] = rov * (ths[m] - theta[i]).max(0.0) / row;
        cap[i] = (1.0 - rov / row) * cap[i] + theta_v[i] * X_MOL * G / R_GAS / t_kelv * x_conv;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thermodynamic_properties_at_20c() {
        let t = 20.0;
        let t_kelv = t + 273.15;

        // 1. Saturated vapor density at 20°C: ~0.0173 kg/m3 (17.3 g/m3)
        let rovs = saturated_vapor_density(t_kelv);
        assert!((rovs - 0.0173).abs() < 0.001);

        // 2. Relative humidity:
        // Saturated (h = 0) -> Hr = 1.0
        assert_eq!(relative_humidity(0.0, t_kelv), 1.0);

        // Capillary suction h = -100 m (-10,000 cm) -> exp(-100 * 0.018015 * 9.81 / (8.314 * 293.15))
        // exponent = -0.007251 -> Hr = 0.99277
        let hr_100m = relative_humidity(-100.0, t_kelv);
        let expected_hr = (-100.0 * X_MOL * G / R_GAS / t_kelv).exp();
        assert!((hr_100m - expected_hr).abs() < 1e-12);

        // High suction h = -10,000 m (-1,000,000 cm, near air-dry / residual)
        let hr_dry = relative_humidity(-10_000.0, t_kelv);
        let expected_dry = (-10_000.0 * X_MOL * G / R_GAS / t_kelv).exp();
        assert!((hr_dry - expected_dry).abs() < 1e-12);
        assert!(hr_dry > 0.0 && hr_dry < 1.0);

        // 3. Volumetric latent heat: row * (2.501e6 - 2369.2 * T)
        let lat = latent_heat_volumetric(t);
        let row_expected = (1.0 - 7.37e-6 * (t - 4.0).powi(2) + 3.79e-8 * (t - 4.0).powi(3)) * 1000.0;
        let expected_lat = row_expected * (2.501e6 - 2369.2 * t);
        assert!((lat - expected_lat).abs() < 1e-6);
    }

    #[test]
    fn test_convapor_dry_vs_saturated_soil() {
        let n = 2;
        let mat = vec![0, 0];
        let ths = vec![0.40];
        let x_conv = 100.0; // cm
        let t_conv = 1.0 / 86400.0; // days

        // Node 0: Dry soil (h = -500 cm, theta = 0.10)
        // Node 1: Saturated soil (h = 0 cm, theta = 0.40)
        let h_new = vec![-500.0, 0.0];
        let temp = vec![25.0, 20.0];
        let con = vec![1e-5, 25.0]; // cm/day
        let theta = vec![0.10, 0.40];

        let mut con_lt = vec![0.0; n];
        let mut con_vt = vec![0.0; n];
        let mut con_vh = vec![0.0; n];

        con_vapor(
            n,
            &mat,
            &h_new,
            &temp,
            &con,
            &theta,
            &ths,
            &mut con_lt,
            &mut con_vt,
            &mut con_vh,
            x_conv,
            t_conv,
            true,
            1, // Campbell enhancement
        );

        // Node 0 (dry): air-filled porosity is high (0.30) -> non-zero vapor conductivities
        assert!(con_vh[0] > 0.0, "Isothermal vapor conductivity should be > 0 in dry soil");
        assert!(con_vt[0] > 0.0, "Thermal vapor conductivity should be > 0 in dry soil");
        assert!(con_lt[0] > 0.0, "Thermal liquid conductivity should be > 0 under suction");

        // Node 1 (saturated): theta == ths -> air-filled porosity is 0 -> vapor flow shuts off
        assert_eq!(con_vh[1], 0.0, "Vapor conductivity must be 0 at complete saturation");
        assert_eq!(con_vt[1], 0.0, "Thermal vapor conductivity must be 0 at complete saturation");
        assert_eq!(con_lt[1], 0.0, "Thermal liquid conductivity is 0 when h = 0");
    }

    #[test]
    fn test_vapor_content_and_capacity_augmentation() {
        let n = 1;
        let mat = vec![0];
        let ths = vec![0.45];
        let x_conv = 100.0; // cm
        let theta = vec![0.15];
        let temp = vec![20.0];
        let h_new = vec![-1000.0]; // dry condition
        let mut theta_v = vec![0.0; n];
        let mut cap = vec![0.001; n]; // baseline liquid differential capacity dTheta/dh

        vapor_content(
            n,
            &mat,
            &theta,
            &mut theta_v,
            &temp,
            &h_new,
            &ths,
            &mut cap,
            x_conv,
        );

        // ThetaV = (rov / row) * (ths - theta)
        // With ths - theta = 0.30 and rov/row ~ 1.7e-5, ThetaV should be ~5e-6
        assert!(theta_v[0] > 0.0 && theta_v[0] < 1e-4);

        // Augmented capacity should remain strictly positive and finite
        assert!(cap[0] > 0.0 && cap[0].is_finite());
    }
}