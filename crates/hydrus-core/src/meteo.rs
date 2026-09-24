//! Meteorological ET (port of the non-daily, non-energy-balance parts of
//! TIME.FOR: SetMeteo, RadGlobal, Cloudiness, CropRes, RadLongNet, Aero).

use crate::model::{MeteoRecord, MeteoSettings};

fn ea(t: f64) -> f64 {
    0.6108 * ((17.27 * t) / (t + 237.3)).exp()
}

/// RadGlobal: potential (extraterrestrial) radiation.
fn rad_global(lat: f64, day_no: f64) -> (f64, f64, f64, f64) {
    let pi = std::f64::consts::PI;
    let sc = 118.08;
    let lat1 = (lat - lat.trunc()) * (5.0 / 3.0) + lat.trunc();
    let lat2 = lat1 * pi / 180.0;
    let sol_declin = (2.0 * pi / 365.0 * day_no - 1.39).sin() * 0.4093;
    let omega = (-sol_declin.tan() * lat2.tan()).clamp(-1.0, 1.0).acos();
    let xx = sol_declin.sin() * lat2.sin();
    let yy = sol_declin.cos() * lat2.cos();
    let dr = 1.0 + 0.033 * ((2.0 * pi / 365.0) * day_no).cos();
    let ra = sc / pi * dr * (omega * xx + omega.sin() * yy);
    (ra, omega, xx, yy)
}

/// Cloudiness: cloudiness factor from sunshine hours / cloudiness / transmission coeff. / solar radiation.
fn cloudiness(
    i_sun_sh: i32,
    sun_hours: f64,
    omega: f64,
    long_wave_b: f64,
    long_wave_a: f64,
    rad: f64,
    rad_cs: f64, // clear-sky radiation: Ra*(as + bs)
    ac: f64,
    bc: f64,
) -> f64 {
    let pi = std::f64::consts::PI;
    match i_sun_sh {
        0 => {
            let nn = 24.0 / pi * omega;
            let n_n = (sun_hours / nn).min(1.0);
            long_wave_b + long_wave_a * n_n
        }
        1 => sun_hours,
        2 => {
            let cover = (2.330 - 3.330 * sun_hours).clamp(0.0001, 1.0);
            long_wave_b + long_wave_a * (1.0 - cover)
        }
        _ => {
            // iSunSh = 3: cloudiness from measured solar radiation (Eq. 57), no clamp
            if rad_cs > 0.0 {
                ac * rad / rad_cs + bc
            } else {
                long_wave_b
            }
        }
    }
}

/// CropRes: crop canopy resistance and soil-cover fraction from LAI / crop height.
/// Returns (rc, scf).
fn crop_res(i_lai: i32, lai_in: f64, r_extinct: f64, i_crop: i32, crop_height: f64) -> (f64, f64) {
    let l_crop = i_crop != 0 && crop_height > 0.0;
    if !l_crop {
        return (0.0, 0.0);
    }
    let lai = match i_lai {
        1 => 0.24 * crop_height,                // clipped grass
        2 => 1.5 * crop_height.ln() + 5.5,      // alfalfa
        3 => {
            if lai_in < 1.0 {
                -((1.0 - lai_in).ln()) / r_extinct.max(0.1)
            } else {
                10.0
            }
        }
        _ => lai_in,
    };
    let mut rc = 0.0;
    let mut scf = 0.0;
    if lai > 0.0 {
        rc = 200.0 / lai;
        scf = (1.0 - (-r_extinct.max(0.1) * lai).exp()).max(0.0);
    }
    (rc, scf)
}

/// RadLongNet: net longwave radiation.
fn rad_long_net(t_max: f64, t_min: f64, a1: f64, b1: f64, ea_dew: f64, cloud_f: f64) -> f64 {
    let sigma = 0.00000000245 * ((t_max + 273.16).powi(4) + (t_min + 273.16).powi(4));
    let emissivity = a1 + b1 * ea_dew.max(0.0).sqrt();
    sigma * cloud_f * emissivity
}

/// Aero: aerodynamic and radiation terms of the Penman-Monteith equation.
#[allow(clippy::too_many_arguments)]
fn aero(
    rn: f64,
    rc: f64,
    crop_height: f64,
    wind_height: f64,
    temp_height: f64,
    wind_ms: f64,
    altitude: f64,
    t_max: f64,
    t_min: f64,
    ea_max: f64,
    ea_min: f64,
    ea_dew: f64,
) -> Option<(f64, f64)> {
    let ch = crop_height.max(0.1);
    let dl0 = 0.667 * ch;
    if dl0 >= wind_height || dl0 >= temp_height {
        return None;
    }
    let aer_dyn_res = ((wind_height - dl0) / (0.123 * ch)).ln() * ((temp_height - dl0) / (0.0123 * ch)).ln() / 0.41f64.powi(2);
    let raa = if wind_ms > 0.0 { aer_dyn_res / wind_ms } else { 0.0 };
    let aero_t_cff = 0.622 * 3.486 * 86400.0 / aer_dyn_res / 1.01;
    let p_atm = 101.3 * ((293.0 - 0.0065 * altitude) / 293.0).powf(5.253);
    let t_aver = (t_max + t_min) / 2.0;
    let lambda = 2.501 - 0.002361 * t_aver;
    let gamma = 0.0016286 * p_atm / lambda;
    let gamma1 = if raa > 0.0 { gamma * (1.0 + rc / raa) } else { gamma };
    let dlt = 2049.0 * ea_max / (t_max + 237.3).powi(2) + 2049.0 * ea_min / (t_min + 237.3).powi(2);
    let dl_dl = dlt / (dlt + gamma1);
    let gm_dl = gamma / (dlt + gamma1);
    let ea_mean = (ea_max + ea_min) / 2.0;
    let aero_term = gm_dl * aero_t_cff / (t_aver + 273.0) * wind_ms * (ea_mean - ea_dew);
    let rad_term = if lambda > 0.0 { dl_dl * rn / lambda } else { 0.0 };
    Some((aero_term, rad_term))
}

/// SetMeteo: potential evaporation and transpiration [mm/d] for one meteorological record.
pub fn potential_et(mp: &MeteoSettings, rec: &MeteoRecord, t_conv: f64) -> (f64, f64) {
    let tt_conv = 24.0 * 3600.0 * t_conv;
    let day_no = (rec.t / tt_conv).rem_euclid(365.0);
    let t_aver = (rec.t_max + rec.t_min) / 2.0;

    // Use dynamic daily crop properties when available (iCrop == 3)
    let crop_height = rec.crop_height.unwrap_or(mp.crop_height);
    let albedo = rec.albedo.unwrap_or(mp.albedo);
    let lai = rec.lai.unwrap_or(mp.lai);

    let (rc, scf) = crop_res(mp.i_lai, lai, mp.r_extinct, mp.i_crop, crop_height);

    if mp.hargreaves {
        let (ra, _, _, _) = rad_global(mp.latitude, day_no);
        let etcomb = 0.0023 * 0.408 * ra * (t_aver + 17.8) * (rec.t_max - rec.t_min).max(0.0).sqrt();
        return (etcomb * (1.0 - scf), etcomb * scf);
    }

    let ea_max = ea(rec.t_max);
    let ea_min = ea(rec.t_min);
    let ea_dew = if mp.i_rel_hum == 0 {
        rec.rh_mean / (50.0 / ea_min + 50.0 / ea_max)
    } else {
        rec.rh_mean
    };

    let rn = if mp.i_radiation != 2 {
        let (ra, omega, xx, yy) = rad_global(mp.latitude, day_no);

        // Detect sub-daily / hourly interval (iMetHour == 1 in TIME.FOR)
        let is_hourly = mp.records.len() > 1 && (mp.records[1].t - mp.records[0].t) <= (0.9999 * tt_conv);

        let cloud_f = if mp.i_radiation == 1 && mp.i_sun_sh == 3 && is_hourly {
            // Sub-daily clear-sky radiation diurnal variation (TIME.FOR lines 485-507)
            let pi = std::f64::consts::PI;
            let rad_cs = ra * (mp.short_wave_a + mp.short_wave_b * 1.0);

            let mut sum_sine = 0.0;
            for k in 1..=24 {
                let sine1 = xx + yy * (2.0 * pi / 24.0 * (k as f64 - 12.0)).cos();
                sum_sine += sine1.max(0.0) / 24.0;
            }

            let hour_no = 24.0 * day_no.fract();
            let sine = xx + yy * (2.0 * pi / 24.0 * (hour_no - 12.0)).cos();
            let rad_csh = if sum_sine > 0.0 {
                (sine * rad_cs / sum_sine).max(0.0)
            } else {
                0.0
            };

            if rad_csh <= 0.0001 {
                // Nighttime default
                mp.cloud_fact_ac * 0.6 + mp.cloud_fact_bc
            } else if rec.rad >= rad_csh {
                mp.cloud_fact_ac * 1.0 + mp.cloud_fact_bc
            } else {
                (mp.cloud_fact_ac * rec.rad / rad_csh + mp.cloud_fact_bc).max(0.01)
            }
        } else {
			cloudiness(
				mp.i_sun_sh,
				rec.sun_hours,
				omega,
				mp.long_wave_b,
				mp.long_wave_a,
				rec.rad,
				ra * (mp.short_wave_a + mp.short_wave_b),
				mp.cloud_fact_ac,
				mp.cloud_fact_bc,
			)
        };

        let rad = if mp.i_radiation == 0 {
            let rel_sun = ((cloud_f - mp.long_wave_b) / mp.long_wave_a.max(1e-12)).clamp(0.0, 1.0);
            ra * (mp.short_wave_a + mp.short_wave_b * rel_sun)
        } else {
            rec.rad
        };

        let rns = (1.0 - albedo) * rad;
        let rnl = rad_long_net(rec.t_max, rec.t_min, mp.long_wave_a1, mp.long_wave_b1, ea_dew, cloud_f);
        rns - rnl
    } else {
        rec.rad
    };

    let wind_ms = rec.wind_kmd / 86.4;
    let Some((aero_term, rad_term)) = aero(
        rn,
        rc,
        crop_height,
        mp.wind_height,
        mp.temp_height,
        wind_ms,
        mp.altitude,
        rec.t_max,
        rec.t_min,
        ea_max,
        ea_min,
        ea_dew,
    ) else {
        return (0.0, 0.0);
    };

    let mut etcomb = (rad_term + aero_term).max(0.0);
    let row = (1.0 - 7.37e-6 * (t_aver - 4.0).powi(2) + 3.79e-8 * (t_aver - 4.0).powi(3)) * 1000.0;
    etcomb = etcomb / row * 1000.0;
    (etcomb * (1.0 - scf), etcomb * scf)
}

/// Interception state container holding carry-over water between time steps.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct InterceptionState {
    pub excess: f64,
}

/// Canopy interception model (port of Intercept and IntercepS from TIME.FOR).
///
/// * `prec_in`: Precipitation intensity [L/T]
/// * `trans_p`: Potential transpiration [L/T]
/// * `lai`: Leaf area index [-]
/// * `a_interc`: Interception constant [L] or [L/T] (default ~0.25 mm)
/// * `scf`: Soil cover fraction [-]
/// * `state`: Carry-over intercepted water
///
/// Returns: `(net_prec, net_trans_p, actual_interc)` in [L/T].
pub fn calculate_interception(
    prec_in: f64,
    trans_p: f64,
    lai: f64,
    a_interc: f64,
    scf: f64,
    state: &mut InterceptionState,
) -> (f64, f64, f64) {
    if lai <= 0.0 || a_interc <= 0.0 {
        return (prec_in, trans_p, 0.0);
    }

    let mut new_interc = 0.0;
    if prec_in > 0.0 {
        let denom = 1.0 + (scf * prec_in) / (a_interc * lai);
        new_interc = (a_interc * lai * (1.0 - 1.0 / denom)).min(prec_in);
    }

    // Capacity limit
    let max_interc = a_interc * lai;
    if new_interc + state.excess > max_interc {
        new_interc = (max_interc - state.excess).max(0.0);
    }

    let net_prec = (prec_in - new_interc).max(0.0);
    let total_interc = state.excess + new_interc;
    state.excess = 0.0;

    // Transpiration consumes intercepted water first
    if trans_p < total_interc {
        state.excess = total_interc - trans_p;
    }
    let net_trans = (trans_p - total_interc).max(0.0);

    (net_prec, net_trans, total_interc)
}

/// Snow accumulation and melt tracking (port of subroutine Snow from TIME.FOR).
///
/// * `prec`: Incoming precipitation [L/T]
/// * `temp`: Air temperature [°C]
/// * `dt`: Time step length [T]
/// * `snow_mf`: Snowmelt degree-day factor [L/(T·°C)]
/// * `snow_layer`: Current snowpack equivalent depth [L]
/// * `r_evap`: Potential soil evaporation [L/T]
/// * `x_conv`: Length conversion to meters
///
/// Returns: `(net_prec, updated_snow_layer, net_evap, min_step_needed)`
pub fn calculate_snow(
    mut prec: f64,
    temp: f64,
    dt: f64,
    snow_mf: f64,
    mut snow_layer: f64,
    mut r_evap: f64,
    x_conv: f64,
    c_top: &mut [f64],
    c_t: &[f64],
) -> (f64, f64, f64, bool) {
    let prec_old = prec;
    let r_evap_old = r_evap;
    let safe_dt = dt.max(1e-8);

    // Rain vs snow fraction Q
    let q = if snow_layer < 0.001 * x_conv {
        if temp < -2.0 {
            1.0
        } else if temp < 2.0 {
            1.0 - ((temp + 2.0) / 4.0)
        } else {
            0.0
        }
    } else {
        1.0
    };

    let r_top = prec * (1.0 - q);
    let snow_f = prec * q;

    let mut snow_melt = if temp > 0.0 && snow_layer > 0.0 {
        temp * snow_mf * dt
    } else {
        0.0
    };

    snow_layer = snow_layer + snow_f * dt - snow_melt;

    if snow_layer < 0.0 {
        snow_melt += snow_layer;
        snow_layer = 0.0;
    } else if snow_layer > 0.0 && r_evap > 0.0 {
        if r_evap * dt < snow_layer {
            snow_layer -= r_evap * dt;
            r_evap = 0.0;
        } else {
            r_evap = (r_evap * dt - snow_layer) / safe_dt;
            snow_layer = 0.0;
        }
    }

    prec = r_top + snow_melt / safe_dt;

    // Check if time step needs reduction due to abrupt rain/melt flux changes
    let min_step = (prec_old - prec).abs() > prec.abs() * 0.2 && prec > 0.0;

    // Solute concentration mixing in snowpack
    if snow_layer > 0.001 * x_conv {
        let denom = snow_layer + dt * (prec_old - r_evap_old);
        if denom > 0.0 {
            for (ct, &cin) in c_top.iter_mut().zip(c_t.iter()) {
                *ct = (snow_layer * (*ct) + dt * prec_old * cin) / denom;
            }
        }
    }

    (prec, snow_layer, r_evap, min_step)
}

/// Atmospheric stability iteration for aerodynamic resistance to heat and vapor flow (AeroRes).
pub fn aero_res(
    temp_height: f64,
    wind_height: f64,
    wind_ms: f64,
    t_kelv_s: f64,
    t_kelv_a: f64,
) -> f64 {
    let r_k = 0.41;
    let pi = std::f64::consts::PI;
    let g = 9.81;
    let ca = 1200.0; // Volumetric heat capacity of air [J/m3/K]

    let zm = 0.001; // Roughness parameter for momentum [m]
    let zh = zm;
    let dl = 0.0;
    let t_height = temp_height / 100.0; // cm -> m
    let w_height = wind_height / 100.0;

    if wind_ms <= 0.0 {
        let diff0 = 2.12e-5;
        let diff_t = diff0 * (t_kelv_a / 273.15).powi(2);
        return t_height / diff_t;
    }

    if (t_kelv_s - t_kelv_a).abs() < 0.01 {
        // Neutral atmosphere
        return (1.0 / wind_ms / r_k / r_k)
            * ((t_height - dl) / zh).ln()
            * ((w_height - dl) / zm).ln();
    }

    let mut psim = 0.0;
    let mut psih = 0.0;
    let mut r_v = 0.0;
    let mut rv_min = 0.0;
    let mut rv_max = 0.0;

    for i in 1..=6 {
        let uu = wind_ms * r_k / (((w_height - dl + zm) / zm).ln() + psim);
        r_v = 1.0 / uu / r_k * (((t_height - dl + zh) / zh).ln() + psih);

        if i == 1 {
            rv_min = 0.1 * r_v;
            rv_max = 10.0 * r_v;
        } else if r_v > rv_max {
            return rv_max;
        } else if r_v < rv_min {
            return rv_min;
        }

        let denom = ca * (t_kelv_s - t_kelv_a) / r_v;
        let mo = if denom.abs() > 1e-12 {
            -ca * t_kelv_a * uu.powi(3) / r_k / g / denom
        } else {
            1e10
        };

        let zeta = (t_height - dl) / mo;
        if zeta < 0.0 {
            // Unstable
            let xx = (1.0 - 16.0 * zeta).max(0.0).powf(0.25);
            psih = -2.0 * ((1.0 + xx * xx) / 2.0).ln();
            psim = -2.0 * ((1.0 + xx) / 2.0).ln() - ((1.0 + xx * xx) / 2.0).ln() + 2.0 * xx.atan() - pi / 2.0;
        } else if zeta > 0.0 {
            // Stable
            if zeta < 1.0 {
                psih = 5.0 * zeta;
                psim = psih;
            } else {
                psih = 5.0;
                psim = psih;
            }
        }
    }
    r_v
}

/// Surface energy balance (port of subroutine Evapor from TIME.FOR).
/// Returns `(r_top_evap, heat_flux_w_m2, sens_flux, evap_kg_m2_s, r_v, r_s)`
pub fn surface_energy_balance(
    temp_s: f64,
    temp_a: f64,
    rh_mean: f64,
    rad: f64,
    h_top: f64,
    theta_top: f64,
    wind_ms: f64,
    wind_height: f64,
    temp_height: f64,
    x_conv: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let g = 9.81;
    let x_mol = 0.018015;
    let r_gas = 8.314;
    let ca = 1200.0;

    let t_kelv_s = temp_s + 273.15;
    let t_kelv_a = temp_a + 273.15;

    // Saturated and actual atmospheric vapor density [kg/m3]
    let rovs_a = 0.001 * (31.3716 - 6014.79 / t_kelv_a - 0.00792495 * t_kelv_a).exp() / t_kelv_a;
    let roa = (rh_mean / 100.0) * rovs_a;

    // Aerodynamic resistance r_v
    let r_v = aero_res(temp_height, wind_height, wind_ms, t_kelv_s, t_kelv_a);
    let r_h = r_v;

    // Soil surface resistance r_s (van de Griend and Owe, 1994)
    let r_s = if theta_top < 0.15 {
        10.0 * (35.63 * (0.15 - theta_top)).exp()
    } else {
        10.0
    };

    // Sensible heat flux [W/m2]
    let sens_flux = ca * (temp_s - temp_a) / r_h;

    // Soil surface vapor density with water potential effect
    let h_m = h_top / x_conv;
    let hr = (h_m * x_mol * g / r_gas / t_kelv_s).exp().clamp(0.0001, 1.0);
    let rovs_s = 0.001 * (31.3716 - 6014.79 / t_kelv_s - 0.00792495 * t_kelv_s).exp() / t_kelv_s;
    let rov = rovs_s * hr;

    // Evaporation flux [kg/(m2·s)]
    let evap_kg_m2_s = ((rov - roa) / (r_v + r_s)).max(0.0);

    // Latent heat of vaporization [J/kg]
    let lat = (2.501 - 0.002361 * temp_a) * 1e6;
    let latent_heat_flux = lat * evap_kg_m2_s; // [W/m2]

    // Net radiation converted from MJ/m2/d to W/m2 (1e6 / 86400)
    let rn_w_m2 = rad * (1e6 / 86400.0);

    // Soil heat flux G [W/m2]
    let heat_flux_w_m2 = rn_w_m2 - sens_flux - latent_heat_flux;

    // Water density [kg/m3]
    let row = (1.0 - 7.37e-6 * (temp_s - 4.0).powi(2) + 3.79e-8 * (temp_s - 4.0).powi(3)) * 1000.0;
    // Evaporation velocity [m/s]
    let evap_m_s = evap_kg_m2_s / row;

    (evap_m_s, heat_flux_w_m2, sens_flux, evap_kg_m2_s, r_v, r_s)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_snow_partition_and_melting() {
        let x_conv = 100.0; // cm
        let dt = 1.0; // days
        let snow_mf = 0.43; // cm/(day·°C)
        let mut c_top = vec![0.0];
        let c_t = vec![1.0];

        // 1. Below -2°C: All snow, 1.0 sublimated from 10.0 incoming -> 9.0 remaining
        let (p_net, snow_pack, evap_net, _) =
            calculate_snow(10.0, -5.0, dt, snow_mf, 0.0, 1.0, x_conv, &mut c_top, &c_t);
        assert_eq!(p_net, 0.0);
        assert!((snow_pack - 9.0).abs() < 1e-6);
        assert_eq!(evap_net, 0.0);

        // 2. Above 0°C with existing snowpack: Melt degree-day
        let (p_net, snow_pack_2, _, _) =
            calculate_snow(0.0, 5.0, dt, snow_mf, snow_pack, 0.0, x_conv, &mut c_top, &c_t);
        let expected_melt = 5.0 * 0.43 * 1.0;
        assert!((p_net - expected_melt).abs() < 1e-6);
        assert!((snow_pack_2 - (9.0 - expected_melt)).abs() < 1e-6);

        // 3. Above +2°C with no snowpack: Pure rain
        let (p_net, snow_pack_3, evap_net_3, _) =
            calculate_snow(5.0, 10.0, dt, snow_mf, 0.0, 2.0, x_conv, &mut c_top, &c_t);
        assert_eq!(p_net, 5.0);
        assert_eq!(snow_pack_3, 0.0);
        assert_eq!(evap_net_3, 2.0);
    }

    #[test]
    fn test_interception_capacity_and_transpiration_reduction() {
        let mut state = InterceptionState::default();
        let lai = 3.0;
        let a_interc = 0.25;
        let scf = 0.8;
        let prec = 10.0;
        let trans_p = 2.0;

        let (net_prec, net_trans, interc) =
            calculate_interception(prec, trans_p, lai, a_interc, scf, &mut state);

        assert!(net_prec < prec);
        assert!(interc > 0.0);
        assert!(net_trans < trans_p);
        assert!(interc <= a_interc * lai);
    }
	#[test]
    fn test_snow_and_interception_integration_pipeline() {
        let x_conv = 100.0; // cm
        let _t_conv = 1.0 / 86400.0; // days
        let dt = 1.0;

        // 1. Setup a test atmospheric & meteo configuration
        let mut mp = MeteoSettings::default();
        mp.lai = 2.0;
        mp.crop_height = 20.0;
        mp.records = vec![MeteoRecord {
            t: 1.0,
            rad: 15.0,
            t_max: -2.0,
            t_min: -8.0,
            rh_mean: 60.0,
            wind_kmd: 100.0,
            sun_hours: 6.0,
			..Default::default()
        }];

        // Resolve temperature from meteo record: (-2.0 + -8.0) / 2.0 = -5.0°C
        let m_rec = &mp.records[0];
        let temp_air = (m_rec.t_max + m_rec.t_min) / 2.0;
        assert_eq!(temp_air, -5.0);

        let mut prec = 10.0;
        let mut r_soil = 2.0;
        let mut snow_layer = 0.0;
        let snow_mf = 0.43;
        let mut c_top = vec![0.0];
        let c_t = vec![0.0];

        // 2. Snow processing: at -5°C, precipitation turns entirely into snowpack
        let (p_eff, s_layer, evap_eff, _) = calculate_snow(
            prec,
            temp_air,
            dt,
            snow_mf,
            snow_layer,
            r_soil,
            x_conv,
            &mut c_top,
            &c_t,
        );
        prec = p_eff;
        snow_layer = s_layer;
        r_soil = evap_eff;

        // Rain reaching surface should be zero; snowpack absorbs all 10 incoming units (minus 2 sublimated)
        assert_eq!(prec, 0.0);
        assert_eq!(evap_eff, 0.0);
        assert!((snow_layer - 8.0).abs() < 1e-6);

        // 3. Interception processing: with zero throughfall, canopy adds no net precipitation
        let mut interc_state = InterceptionState::default();
        let mut r_root = 1.5;
        let scf = (1.0 - (-0.463f64 * mp.lai).exp()).max(0.0);
        let a_interc = 0.25;

        let (p_net, tr_net, _) = calculate_interception(
            prec,
            r_root,
            mp.lai,
            a_interc,
            scf,
            &mut interc_state,
        );
        prec = p_net;
        r_root = tr_net;

        assert_eq!(prec, 0.0);
        assert_eq!(r_root, 1.5);

        // 4. Net top boundary flux (r_soil - prec)
        let r_top = r_soil.abs() - prec.abs();
        assert_eq!(r_top, 0.0);
    }
	
	#[test]
    fn test_dynamic_crop_overrides_and_snow_interception() {
        let x_conv = 100.0; // cm
        let t_conv = 1.0 / 86400.0; // days
        let dt = 1.0;

        // Base settings with static fallbacks
        let mp = MeteoSettings {
            crop_height: 10.0,
            albedo: 0.25,
            lai: 1.0,
            i_crop: 3, // Daily dynamic values from METEO.IN
            i_lai: 3,
            ..Default::default()
        };

        // Record with dynamic daily overrides (iCrop == 3)
        let rec = MeteoRecord {
            t: 1.0,
            rad: 20.0,
            t_max: 22.0,
            t_min: 12.0,
            rh_mean: 45.0,
            wind_kmd: 120.0,
            sun_hours: 9.0,
            crop_height: Some(50.0), // overrides mp.crop_height
            albedo: Some(0.18),      // overrides mp.albedo
            lai: Some(3.5),          // overrides mp.lai
            x_root: Some(60.0),      // dynamic root depth
        };

        // 1. Verify dynamic potential ET calculation runs with overrides
        let (evap_p, trans_p) = potential_et(&mp, &rec, t_conv);
        assert!(evap_p > 0.0);
        assert!(trans_p > 0.0);

        // 2. Verify dynamic rooting depth extraction
        assert_eq!(rec.x_root, Some(60.0));

        // 3. Verify canopy interception using the dynamic LAI
        let mut interc_state = InterceptionState::default();
        let dynamic_lai = rec.lai.unwrap();
        let scf = (1.0 - (-0.463f64 * dynamic_lai).exp()).max(0.0);
        let prec = 15.0;

        let (net_prec, net_trans, interc) = calculate_interception(
            prec,
            trans_p,
            dynamic_lai,
            0.25,
            scf,
            &mut interc_state,
        );

        assert!(net_prec < prec);
        assert!(interc > 0.0);
        assert!(net_trans <= trans_p);

        // 4. Verify snow degree-day melt with solute mixing slices
        let mut c_top = vec![0.0];
        let c_t = vec![2.0];
        let (p_melt, snow_rem, _, min_step) = calculate_snow(
            0.0,
            4.0, // 4°C melt condition
            dt,
            0.43,
            5.0, // 5 cm snowpack
            0.5, // 0.5 cm evaporation
            x_conv,
            &mut c_top,
            &c_t,
        );

        let expected_melt = 4.0 * 0.43 * dt;
		let expected_snow_rem = 5.0 - expected_melt - 0.5 * dt;
        assert!((p_melt - expected_melt).abs() < 1e-6);
		assert!((snow_rem - expected_snow_rem).abs() < 1e-6);
        assert!(min_step);
    }
	
	#[test]
    fn test_surface_energy_balance_and_aero_res() {
        let rv = aero_res(150.0, 200.0, 2.5, 20.0 + 273.15, 18.0 + 273.15);
        assert!(rv > 0.0 && rv < 500.0);

        let (evap_m_s, g_flux, sens_flux, evap_kg, _, _) = surface_energy_balance(
            22.0,   // Soil temp [C]
            20.0,   // Air temp [C]
            50.0,   // RH [%]
            15.0,   // Net radiation [MJ/m2/d]
            -100.0, // hTop [cm]
            0.20,   // ThetaTop
            2.0,    // Wind [m/s]
            200.0,
            150.0,
            100.0,  // xConv
        );

        assert!(evap_m_s > 0.0);
        assert!(evap_kg > 0.0);
        assert!(sens_flux > 0.0); // Soil is warmer than air -> sensible heat flux upwards
        // Net radiation = G + H + Latent
        let rn_w = 15.0 * (1e6 / 86400.0);
        let latent = evap_kg * (2.501 - 0.002361 * 20.0) * 1e6;
        assert!((rn_w - (g_flux + sens_flux + latent)).abs() < 1e-4);
    }
	
	#[test]
    fn test_hourly_radiation_imethour_diurnal_curve() {
        let t_conv = 1.0 / 86400.0; // days
        let mp = MeteoSettings {
            latitude: 33.58,
            altitude: 306.0,
            short_wave_a: 0.25,
            short_wave_b: 0.50,
            i_radiation: 1, // Solar radiation measured
            i_sun_sh: 3,     // Cloudiness from solar radiation
            cloud_fact_ac: 1.35,
            cloud_fact_bc: -0.35,
            records: vec![
                MeteoRecord {
                    t: 328.042, // ~1:00 AM (nighttime)
                    rad: 0.0,
                    t_max: 15.2,
                    t_min: 15.2,
                    rh_mean: 33.0,
                    wind_kmd: 138.24,
                    sun_hours: 0.545,
                    ..Default::default()
                },
                MeteoRecord {
                    t: 328.083, // ~2:00 AM (sub-daily interval confirms is_hourly = true)
                    rad: 0.0,
                    t_max: 14.9,
                    t_min: 14.9,
                    rh_mean: 32.0,
                    wind_kmd: 103.68,
                    sun_hours: 0.545,
                    ..Default::default()
                },
                MeteoRecord {
                    t: 328.500, // 12:00 PM (solar noon)
                    rad: 43.2864,
                    t_max: 29.6,
                    t_min: 29.6,
                    rh_mean: 14.0,
                    wind_kmd: 198.72,
                    sun_hours: 0.545,
                    ..Default::default()
                },
            ],
            ..Default::default()
        };

        // 1. Nighttime record (rad = 0.0): should evaluate nighttime cloudiness and minimal/zero ET
        let (evap_night, trans_night) = potential_et(&mp, &mp.records[0], t_conv);
        assert!(evap_night >= 0.0);
        assert!(trans_night >= 0.0);

        // 2. Solar noon record (rad = 43.28 MJ/m2/d): should compute positive daytime ET
        let (evap_noon, trans_noon) = potential_et(&mp, &mp.records[2], t_conv);
        assert!(evap_noon > evap_night);
        assert!(trans_noon >= 0.0);
    }
	
	#[test]
    fn test_lai_partitioning_formula() {
        let r_pet: f64 = 0.5; // cm/day total potential evapotranspiration
        let lai: f64 = 2.5;
        let r_extinct: f64 = 0.463;

        // SC = 1 - exp(-k * LAI)
        let sc: f64 = 1.0 - (-r_extinct * lai).exp();
        let expected_r_root: f64 = r_pet * sc;
        let expected_r_soil: f64 = r_pet - expected_r_root;

        assert!((expected_r_root + expected_r_soil - r_pet).abs() < 1e-10);
        assert!(expected_r_root > 0.0 && expected_r_soil > 0.0);
        // At LAI 2.5 with k=0.463, roughly 68.5% goes to transpiration and 31.5% to soil evaporation
        assert!((sc - 0.685).abs() < 0.01);
    }

    #[test]
    fn test_cyclic_boundary_time_wrapping() {
        // Simulating a 24-hour cycle (t_period = 1.0 day) repeated over 3 days
        let t_period: f64 = 1.0;
        let t_init: f64 = 0.0;

        let eval_times: [f64; 3] = [0.25, 1.25, 2.25];
        let wrapped_times: Vec<f64> = eval_times
            .iter()
            .map(|&t: &f64| t_init + (t - t_init).rem_euclid(t_period))
            .collect();

        assert!((wrapped_times[0] - 0.25).abs() < 1e-12);
        assert!((wrapped_times[1] - 0.25).abs() < 1e-12);
        assert!((wrapped_times[2] - 0.25).abs() < 1e-12);
    }
	#[test]
	fn test_boundary_condition_cycling_execution() {
		use crate::model::*;
		use crate::sim::Simulation;

		let mut prj = Project::default();
		prj.time.t_init = 0.0;
		prj.time.t_max = 10.0;
		prj.time.dt = 0.5;
		prj.time.dt_min = 0.01;
		prj.time.dt_max = 0.5;

		prj.water.bc.atmospheric = true;
		prj.water.bc.top_time_variable = true;
		prj.atmosphere.bc_cycles = true;

		// Define a 2-day repeating cycle of precipitation
		prj.atmosphere.records = vec![
			AtmRecord {
				t: 1.0,
				prec: 2.0,
				evap: 0.1,
				transp: 0.0,
				h_crit_a: 100000.0,
				..Default::default()
			},
			AtmRecord {
				t: 2.0,
				prec: 0.0,
				evap: 0.5,
				transp: 0.0,
				h_crit_a: 100000.0,
				..Default::default()
			},
		];

		let mut sim = Simulation::new(prj).expect("Simulation should initialize");

		// Run across multiple cycles up to t = 10.0
		let mut steps = 0;
		while sim.t < 10.0 && steps < 200 {
			sim.step();
			steps += 1;
		}

		assert!(sim.t >= 10.0, "Simulation should reach t = 10.0 via repeating cycles");
		assert!(!sim.failed);
	}
}