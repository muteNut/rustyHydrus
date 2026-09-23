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
        if self.atm_idx >= recs.len() {
            // 'end' record: the run stops at the last atmospheric time
            self.t_max = self.t_atm1;
            return Ok(());
        }
        let rec = recs[self.atm_idx].clone();
        self.atm_idx += 1;
        self.t_atm1 = rec.t;
        self.prec = rec.prec;
        self.r_soil = rec.evap;
        let rr = rec.transp;
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
                self.r_top = self.r_soil.abs() - self.prec.abs();
				if let Some(ref mp) = self.prj.atmosphere.meteo {
					while self.meteo_idx + 1 < mp.records.len() && self.t >= mp.records[self.meteo_idx].t {
						self.meteo_idx += 1;
					}
					if let Some(rec) = mp.records.get(self.meteo_idx) {
						let (evap_p, trans_p) = crate::meteo::potential_et(mp, rec, self.t_conv);
						let r_conv = 0.001 * self.x_conv;
						let tt_conv = 24.0 * 3600.0 * self.t_conv;
						self.r_soil = evap_p * r_conv / tt_conv;
						self.r_root = trans_p * r_conv / tt_conv;
						self.r_top = self.r_soil.abs() - self.prec.abs();
					}
				}
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
				self.r_root = rr.abs();
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
