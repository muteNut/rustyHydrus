//! Time-variable boundary conditions and time-step control (port of parts of TIME.FOR).

use crate::error::HydrusError;
use crate::model::SoilModel;
use crate::output::ALevel;
use crate::sim::Simulation;

impl Simulation {
    /// Read the next atmospheric record (SetBC).
    pub fn set_bc(&mut self) -> Result<(), HydrusError> {
        let n = self.n;
        let kod_top_old = self.kod_top;
        let recs = &self.prj.atmosphere.records;
        let n_rec = recs.len();

        if n_rec == 0 {
            return Ok(());
        }

        // When boundary cycles are enabled and the table ends, wrap back to the beginning
        if self.prj.atmosphere.bc_cycles && self.atm_idx >= n_rec {
            let t_first = recs.first().map(|r| r.t).unwrap_or(0.0);
            let t_last = recs.last().map(|r| r.t).unwrap_or(1.0);
            let cycle_period = (t_last - t_first).max(1e-6);

            self.atm_idx %= n_rec;
            let num_cycles = ((self.t - self.t_init) / cycle_period).floor();
            self.t_atm_old = self.t;
            self.t_atm1 = t_first + (num_cycles + 1.0) * cycle_period;
        }

        if self.atm_idx >= n_rec {
            // 'end' record: run stops at the last atmospheric time
            self.t_max = self.t_atm1;
            return Ok(());
        }

        let rec = recs[self.atm_idx].clone();
        self.atm_idx += 1;

        if self.prj.atmosphere.bc_cycles {
            let t_first = recs.first().map(|r| r.t).unwrap_or(0.0);
            let t_last = recs.last().map(|r| r.t).unwrap_or(1.0);
            let cycle_period = (t_last - t_first).max(1e-6);
            let num_cycles = ((self.t - self.t_init) / cycle_period).floor();
            self.t_atm1 = rec.t + num_cycles * cycle_period;
            if self.t_atm1 <= self.t && self.atm_idx < n_rec {
                self.t_atm1 = recs[self.atm_idx].t + num_cycles * cycle_period;
            }
        } else {
            self.t_atm1 = rec.t;
        }

        self.prec = rec.prec;
        let mut r_r = rec.transp;
        self.r_soil = rec.evap;
        if self.prj.atmosphere.lai_partitioning && self.prj.atmosphere.meteo.is_none() {
            let r_lai = r_r;
            let r_pet = self.r_soil;
            r_r = 0.0;
            if r_lai > 0.0 {
                r_r = r_pet * (1.0 - (-self.prj.atmosphere.extinction.max(0.1) * r_lai).exp()).max(0.0);
            }
            self.r_soil = r_pet - r_r;
        }
        let hca = rec.h_crit_a;
        if self.prj.atmosphere.has_root_depth {
            self.x_root = rec.x_root;
        }
        self.atm_record_scalars(&rec);
		
        if self.top_inf {
            let r_top_old = self.r_top;
            self.h_crit_a = -hca.abs();
            if self.l_var_bc {
                self.r_top = self.prec;
                if (r_top_old - self.r_top).abs() > self.r_top.abs() * 0.2 && self.r_top < 0.0 {
                    self.min_step = true;
                }
                self.kod_top = self.r_soil as i32;
                self.r_soil = 0.0;
                if self.kod_top == -1 && kod_top_old == 1 && self.prec > 0.0 && self.h_new[n - 1] > 0.0 {
                    self.h_new[n - 1] = -0.01 * self.x_conv;
                }
            } else {
                // Advance meteo index to current simulation time
                if let Some(ref mp) = self.prj.atmosphere.meteo {
                    while self.meteo_idx + 1 < mp.records.len() && self.t >= mp.records[self.meteo_idx].t {
                        self.meteo_idx += 1;
                    }
                }

                // Determine air temperature for snow accumulation/melting
                let temp_air = if self.l_temp {
                    rec.t_top
                } else if let Some(ref mp) = self.prj.atmosphere.meteo {
                    if let Some(m_rec) = mp.records.get(self.meteo_idx) {
                        (m_rec.t_max + m_rec.t_min) / 2.0
                    } else {
                        rec.t_top
                    }
                } else {
                    rec.t_top
                };

                // 1. Snow processing
                if self.prj.atmosphere.snow {
                    let mut c_top_dummy = vec![0.0; self.n_species()];
                    let c_t_dummy = vec![0.0; self.n_species()];
                    let (p_eff, s_layer, evap_eff, min_st) = crate::meteo::calculate_snow(
                        self.prec,
                        temp_air,
                        self.t_atm1 - self.t_atm_old,
                        self.prj.atmosphere.snow_mf,
                        self.snow_layer,
                        self.r_soil,
                        self.x_conv,
                        &mut c_top_dummy,
                        &c_t_dummy,
                    );
                    self.prec = p_eff;
                    self.snow_layer = s_layer;
                    self.r_soil = evap_eff;
                    if min_st {
                        self.min_step = true;
                    }
                }

                // 2. Potential ET / Surface Energy Balance from meteorological records
                if let Some(ref mp) = self.prj.atmosphere.meteo {
                    if let Some(m_rec) = mp.records.get(self.meteo_idx) {
                        if mp.l_en_bal {
                            // Coupled surface energy balance (TIME.FOR subroutine Evapor)
							let temp_s = if self.l_temp { self.temp(n - 1) } else { temp_air };
							let h_top_val = self.h_new[n - 1];
							let theta_top = self.th_new[n - 1];
                            let wind_ms = m_rec.wind_kmd / 86.4;

                            let (evap_m_s, _heat_flux_w, _sens_flux, _evap_kg, _rv, _rs) =
                                crate::meteo::surface_energy_balance(
                                    temp_s,
                                    temp_air,
                                    m_rec.rh_mean,
                                    m_rec.rad,
                                    h_top_val,
                                    theta_top,
                                    wind_ms,
                                    mp.wind_height,
                                    mp.temp_height,
                                    self.x_conv,
                                );

                            // Convert evaporation velocity [m/s] to simulation length/time units [L/T]
                            let evap_rate = evap_m_s * self.x_conv / self.t_conv;
                            self.r_soil = evap_rate;
                            self.r_root = 0.0;

                            if self.l_temp && self.prj.heat.k_top == -1 {
                                self.prj.heat.t_top = temp_s;
                            }
                        } else {
                            let (evap_p, trans_p) = crate::meteo::potential_et(mp, m_rec, self.t_conv);
                            let r_conv = 0.001 * self.x_conv;
                            let tt_conv = 24.0 * 3600.0 * self.t_conv;
                            self.r_soil = evap_p * r_conv / tt_conv;
                            self.r_root = trans_p * r_conv / tt_conv;
                        }

                        // Dynamic daily rooting depth (iCrop == 3)
                        if let Some(xr) = m_rec.x_root {
                            self.x_root = xr;
                        }
                    }
                }

                // 3. Canopy Interception
                if self.prj.atmosphere.interception {
                    let lai = if let Some(ref mp) = self.prj.atmosphere.meteo {
                        if let Some(m_rec) = mp.records.get(self.meteo_idx) {
                            m_rec.lai.unwrap_or(mp.lai)
                        } else {
                            mp.lai
                        }
                    } else {
                        1.0
                    };
                    let scf = (1.0 - (-self.prj.atmosphere.extinction.max(0.1) * lai).exp()).max(0.0);
                    let (p_net, tr_net, _) = crate::meteo::calculate_interception(
                        self.prec,
                        self.r_root,
                        lai,
                        self.prj.atmosphere.interception_a,
                        scf,
                        &mut self.interc_state,
                    );
                    self.prec = p_net;
                    self.r_root = tr_net;
                }

                // 4. Net surface flux
                self.r_top = self.r_soil.abs() - self.prec.abs();

                if (r_top_old - self.r_top).abs() > self.r_top.abs() * 0.2 && self.r_top < 0.0 {
                    self.min_step = true;
                }
                if self.r_top > 0.0 && r_top_old < 0.0 && !self.w_layer {
                    let mut x_limit = 0.0;
                    if self.model == SoilModel::VGAirEntry {
                        x_limit = -0.03 * self.x_conv;
                    }
                    if self.kod_top == 4 || self.h_new[n - 1] > x_limit {
                        if self.model != SoilModel::VGAirEntry {
                            x_limit = -0.01 * self.x_conv;
                        }
                        self.h_new[n - 1] = x_limit;
                        self.kod_top = -4;
                    }
                }
            }
            if self.kod_top == 3 || self.l_var_bc {
                self.h_top = rec.h_top;
            }
            if self.prj.atmosphere.meteo.is_none() {
                self.r_root = r_r.abs();
            }
        }
		
        if self.bot_inf {
            if (self.r_bot - rec.r_bot).abs() > self.r_bot.abs() * 0.2 {
                self.min_step = true;
            }
            self.r_bot = rec.r_bot;
            if (self.h_bot - rec.h_bot - self.gwl0l).abs() > self.h_bot.abs() * 0.2 {
                self.min_step = true;
            }
            self.h_bot = rec.h_bot + self.gwl0l;
        }
        Ok(())
    }

    pub fn tm_cont(&mut self, iter: usize, t_print: f64, dt_max_c: f64) {
        let dt_max_w = self.dt_max;
        let dt_max;
        if self.min_step {
            dt_max = dt_max_w.min(dt_max_c).min(self.dt_init).min(self.dt_opt);
            self.dt_opt = dt_max;
            self.min_step = false;
        } else {
            dt_max = dt_max_w.min(dt_max_c);
        }
        let t_fix = t_print.min(self.t_atm).min(self.t_max);
        if iter <= self.it_min && (t_fix - self.t) >= self.d_mul * self.dt_opt {
            self.dt_opt = dt_max.min(self.d_mul * self.dt_opt);
        }
        if iter >= self.it_max {
            self.dt_opt = self.dt_min.max(self.d_mul2 * self.dt_opt);
        }
        self.dt = self.dt_opt.min(t_fix - self.t);
        let mut i_step = 1i64;
        if self.dt > 0.0 {
            i_step = ((t_fix - self.t) / self.dt).round() as i64;
        }
        if i_step >= 1 && i_step <= 10 {
            self.dt = ((t_fix - self.t) / i_step as f64).min(dt_max);
        }
        if i_step == 1 {
            self.dt = t_fix - self.t;
            if self.dt - dt_max > self.dt_min {
                self.dt /= 2.0;
            }
        }
        if self.dt <= 0.0 {
            self.dt = self.dt_min / 3.0;
        }
    }

    fn daily_factor(&self, base: f64) -> f64 {
        let t_day = self.t / self.t_conv / 86400.0;
        let rem = t_day % 1.0;
        if rem <= 0.264 || rem >= 0.736 {
            0.24 * base
        } else {
            2.75 * base * (2.0 * std::f64::consts::PI * t_day - 6.0 * std::f64::consts::PI / 12.0).sin()
        }
    }

    pub fn daily_var_root_soil(&mut self) {
        self.r_root = self.daily_factor(self.r_root_d);
        self.r_soil = self.daily_factor(self.r_soil_d);
    }

    pub fn sin_prec(&mut self) {
        let dt = self.t_atm1 - self.t_atm_old;
        if self.prec_d > 0.0 && dt > 0.0 {
            self.prec = self.prec_d * (1.0 + (2.0 * std::f64::consts::PI * (self.t - self.t_atm_old) / dt - std::f64::consts::PI).cos());
        } else {
            self.prec = 0.0;
        }
    }

    pub fn a_level(&mut self) {
        let n = self.n;
        self.res.alevel.push(ALevel {
            t: self.t,
            cum: [self.cum_q[0], self.cum_q[1], self.cum_q[2], self.cum_q[3], self.cum_q[4]],
            h_top: self.h_new[n - 1],
            h_root: self.h_root,
            h_bot: self.h_new[0],
        });
    }
}
