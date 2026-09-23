//! Meteorological ET (port of the non-daily, non-energy-balance parts of
//! TIME.FOR: SetMeteo, RadGlobal, Cloudiness, CropRes, RadLongNet, Aero).

use crate::model::{MeteoRecord, MeteoSettings};

fn ea(t: f64) -> f64 {
    0.6108 * ((17.27 * t) / (t + 237.3)).exp()
}

/// RadGlobal: potential (extraterrestrial) radiation.
fn rad_global(lat: f64, day_no: f64) -> (f64, f64) {
    let pi = std::f64::consts::PI;
    let sc = 118.08;
    let lat1 = (lat - lat.trunc()) * (5.0 / 3.0) + lat.trunc();
    let lat2 = lat1 * pi / 180.0;
    let sol_declin = (2.0 * pi / 365.0 * day_no - 1.39).sin() * 0.4093;
    let omega = (-sol_declin.tan() * lat2.tan()).acos();
    let xx = sol_declin.sin() * lat2.sin();
    let yy = sol_declin.cos() * lat2.cos();
    let dr = 1.0 + 0.033 * ((2.0 * pi / 365.0) * day_no).cos();
    let ra = sc / pi * dr * (omega * xx + omega.sin() * yy);
    (ra, omega)
}

/// Cloudiness: cloudiness factor from sunshine hours / cloudiness / transmission coeff. / solar radiation.
fn cloudiness(
    i_sun_sh: i32,
    sun_hours: f64,
    omega: f64,
    long_wave_b: f64,
    long_wave_a: f64,
    rad: f64,
    ra: f64,
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
            if ra > 0.0 {
                (ac * (rad / ra) + bc).clamp(0.0, 1.0)
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
    let (rc, scf) = crop_res(mp.i_lai, mp.lai, mp.r_extinct, mp.i_crop, mp.crop_height);

    if mp.hargreaves {
        let (ra, _) = rad_global(mp.latitude, day_no);
        let etcomb = 0.0023 * 0.408 * ra * (t_aver + 17.8) * (rec.t_max - rec.t_min).max(0.0).sqrt();
        return (etcomb * (1.0 - scf), etcomb * scf);
    }

    let ea_max = ea(rec.t_max);
    let ea_min = ea(rec.t_min);
    let ea_dew = if mp.i_rel_hum == 0 { rec.rh_mean / (50.0 / ea_min + 50.0 / ea_max) } else { rec.rh_mean };

    let rn = if mp.i_radiation != 2 {
        let (ra, omega) = rad_global(mp.latitude, day_no);
        let cloud_f = cloudiness(
            mp.i_sun_sh,
            rec.sun_hours,
            omega,
            mp.long_wave_b,
            mp.long_wave_a,
            rec.rad,
            ra,
            mp.cloud_fact_ac,
            mp.cloud_fact_bc,
        );
        let rad = if mp.i_radiation == 0 {
            let rel_sun = ((cloud_f - mp.long_wave_b) / mp.long_wave_a.max(1e-12)).clamp(0.0, 1.0);
            ra * (mp.short_wave_a + mp.short_wave_b * rel_sun)
        } else {
            rec.rad
        };
        let rns = (1.0 - mp.albedo) * rad;
        let rnl = rad_long_net(rec.t_max, rec.t_min, mp.long_wave_a1, mp.long_wave_b1, ea_dew, cloud_f);
        rns - rnl
    } else {
        rec.rad
    };

    let wind_ms = rec.wind_kmd / 86.4;
    let Some((aero_term, rad_term)) = aero(
        rn,
        rc,
        mp.crop_height,
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