//! Heat transport (port of TEMPER.FOR, without vapor flow / energy balance).

use crate::error::HydrusError;
use crate::sim::Simulation;

pub struct HeatState {
    pub temp_o: Vec<f64>,
    pub temp_n: Vec<f64>,
    pub k_top: i32,
    pub k_bot: i32,
    pub t_top: f64,
    pub t_bot: f64,
    pub ampl: f64,
    pub heat_fl: f64,
}

pub fn thomas(b: &[f64], d: &mut [f64], e: &[f64], f: &mut [f64]) {
    // b: lower, d: diag, e: upper, f: rhs -> solution in f
    let n = d.len();
    for i in 1..n {
        let m = b[i] / d[i - 1];
        d[i] -= m * e[i - 1];
        f[i] -= m * f[i - 1];
    }
    f[n - 1] /= d[n - 1];
    for i in (0..n - 1).rev() {
        f[i] = (f[i] - e[i] * f[i + 1]) / d[i];
    }
}

impl Simulation {
    pub fn temp(&self, i: usize) -> f64 {
        self.heat.as_ref().map(|h| h.temp_n[i]).unwrap_or(0.0)
    }

    pub fn heat_init(&mut self) -> Result<(), HydrusError> {
        let ht = &self.prj.heat;
        let temp: Vec<f64> = self.prj.profile.nodes.iter().rev().map(|n| n.temp).collect();
        self.heat = Some(HeatState {
            temp_o: temp.clone(),
            temp_n: temp,
            k_top: ht.k_top,
            k_bot: ht.k_bot,
            t_top: ht.t_top,
            t_bot: ht.t_bot,
            ampl: ht.amplitude,
            heat_fl: 0.0,
        });
        Ok(())
    }

    fn temp_cap(&self, theta: &[f64], veloc: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = self.n;
        let ht = &self.prj.heat;
        let mut cap = vec![0.0; n];
        let mut cond = vec![0.0; n];
        for i in 0..n {
            let m = &ht.materials[self.mat[i]];
            let th = theta[i];
            cap[i] = m.cn * m.qn + m.co * m.qo + m.cw * th;
            cond[i] = (m.b1 + m.b2 * th + m.b3 * th.sqrt()).max(0.0);
            if ht.campbell {
                // Chung & Horton (Campbell) thermal conductivity; qn=solid fraction, b1..b3 = quartz, other minerals, clay
                let (qs, qz, qm, qc) = (m.qn, m.b1, m.b2, m.b3);
                let aa = (0.57 + 1.73 * qz + 0.93 * qm) / (1.0 - 0.74 * qz - 0.49 * qm) - 2.8 * qs * (1.0 - qs);
                let bb = 2.8 * qs;
                let xc = qc.max(0.005);
                let cc = 1.0 + 2.6 / xc.sqrt();
                let dd = 0.03 + 0.7 * qs * qs;
                let ee = 4.0;
                let lamb = aa + bb * th - (aa - dd) * (-(cc * th).powf(ee)).exp();
                cond[i] = (lamb * self.x_conv / self.t_conv.powi(3)).max(0.0);
            }
            cond[i] += m.cw * m.disp * veloc[i].abs();
        }
        (cap, cond)
    }

    /// Solve heat transport for one time step.
    /// Solve heat transport for one time step (mit Dampftransport/TempAdj).
    pub fn temper(&mut self) {
        let n = self.n;
        let dt = self.dt;
        let ht = self.prj.heat.clone();
        let cw = ht.materials[0].cw;
        let pi = std::f64::consts::PI;
        let (t, sink) = (self.t, self.sink.clone());
        let hs = self.heat.as_mut().unwrap();
        let mut t_top_a = hs.t_top;
        if ht.period > 0.0 {
            t_top_a = hs.t_top + hs.ampl * (2.0 * pi * t / ht.period - 7.0 * pi / 12.0).sin();
        }
        if hs.k_top < 0 {
            hs.heat_fl = cw * t_top_a * (self.v_new[n - 1] + self.v_old[n - 1]) / 2.0;
            // Dampffluss-Beitrag zur oberen Wärmestrom-Randbedingung
            if self.l_vapor {
                let cv = 1.8e6 / self.x_conv / self.t_conv.powi(2);
                let x_lat = crate::vapor::latent_heat_volumetric(hs.temp_n[n - 1]) / self.x_conv / self.t_conv.powi(2);
                let v_v_mean = (self.v_v_new[n - 1] + self.v_v_old[n - 1]) / 2.0;
                hs.heat_fl += cv * t_top_a * v_v_mean + x_lat * v_v_mean;
            }
        }
        let heat_fl = hs.heat_fl;
        let (k_top, k_bot, t_bot) = (hs.k_top, hs.k_bot, hs.t_bot);
        let temp_o = hs.temp_o.clone();
        let x = &self.x;
        let mut b = vec![0.0; n];
        let mut d = vec![0.0; n];
        let mut e = vec![0.0; n];
        let mut f = vec![0.0; n];

        // 1. Basis-Matrixaufbau für Level 1 (alt) und Level 2 (neu)
        for level in 1..=2 {
            let (cap, cond) = if level == 1 { self.temp_cap(&self.th_old, &self.v_old) } else { self.temp_cap(&self.th_new, &self.v_new) };
            let (v_new, v_old) = (&self.v_new, &self.v_old);
            let dx = x[1] - x[0];
            if k_bot > 0 {
                d[0] = 1.0;
                e[0] = 0.0;
                f[0] = t_bot;
            } else if k_bot < 0 {
                if level == 2 {
                    d[0] = dx / 2.0 / dt * cap[0] + (cond[0] + cond[1]) / dx / 4.0 + cw * (2.0 * v_new[0] + v_new[1]) / 12.0 + dx / 24.0 * cw * (3.0 * sink[0] + sink[1]);
                    e[0] = -(cond[0] + cond[1]) / 4.0 / dx + cw * (2.0 * v_new[1] + v_new[0]) / 12.0 + dx / 24.0 * cw * (sink[0] + sink[1]);
                } else {
                    f[0] = temp_o[0] * (dx / 2.0 / dt * cap[0] - (cond[0] + cond[1]) / dx / 4.0 - cw * (2.0 * v_old[0] + v_old[1]) / 12.0 - dx / 24.0 * cw * (3.0 * sink[0] + sink[1]))
                        + temp_o[1] * ((cond[0] + cond[1]) / 4.0 / dx - cw * (2.0 * v_old[1] + v_old[0]) / 12.0 - dx / 24.0 * cw * (sink[0] + sink[1]))
                        + t_bot * cw * (v_new[0] + v_old[0]) / 2.0;
                }
            } else {
                d[0] = -1.0;
                e[0] = 1.0;
                f[0] = 0.0;
            }
            for i in 1..n - 1 {
                let dxa = x[i] - x[i - 1];
                let dxb = x[i + 1] - x[i];
                let dx = (x[i + 1] - x[i - 1]) / 2.0;
                if level == 2 {
                    b[i] = -(cond[i] + cond[i - 1]) / 4.0 / dxa - cw * (v_new[i] + 2.0 * v_new[i - 1]) / 12.0 + dxa / 24.0 * cw * (sink[i - 1] + sink[i]);
                    d[i] = (cond[i - 1] + cond[i]) / 4.0 / dxa
                        + (cond[i] + cond[i + 1]) / 4.0 / dxb
                        + dx / dt * cap[i]
                        + cw * (v_new[i + 1] - v_new[i - 1]) / 12.0
                        + dxa / 24.0 * cw * (sink[i - 1] + 3.0 * sink[i])
                        + dxb / 24.0 * cw * (3.0 * sink[i] + sink[i + 1]);
                    e[i] = -(cond[i] + cond[i + 1]) / 4.0 / dxb + cw * (2.0 * v_new[i + 1] + v_new[i]) / 12.0 + dxb / 24.0 * cw * (sink[i + 1] + sink[i]);
                } else {
                    f[i] = temp_o[i - 1] * ((cond[i] + cond[i - 1]) / 4.0 / dxa + cw * (v_old[i] + 2.0 * v_old[i - 1]) / 12.0 - dxa / 24.0 * cw * (sink[i - 1] + sink[i]))
                        + temp_o[i]
                            * (-cw * (v_old[i + 1] - v_old[i - 1]) / 12.0 + dx / dt * cap[i]
                                - (cond[i + 1] + cond[i]) / 4.0 / dxb
                                - (cond[i] + cond[i - 1]) / 4.0 / dxa
                                - dxa / 24.0 * cw * (sink[i - 1] + 3.0 * sink[i])
                                - dxb / 24.0 * cw * (3.0 * sink[i] + sink[i + 1]))
                        + temp_o[i + 1] * ((cond[i + 1] + cond[i]) / 4.0 / dxb - cw * (2.0 * v_old[i + 1] + v_old[i]) / 12.0 - dxb / 24.0 * cw * (sink[i + 1] + sink[i]));
                }
            }
            let m = n - 1;
            if k_top > 0 {
                b[m] = 0.0;
                d[m] = 1.0;
                f[m] = t_top_a;
            } else if k_top < 0 {
                let dx = x[m] - x[m - 1];
                if level == 2 {
                    b[m] = -(cond[m] + cond[m - 1]) / 4.0 / dx - cw * (v_new[m] + 2.0 * v_new[m - 1]) / 12.0 + dx / 24.0 * cw * (sink[m - 1] + sink[m]);
                    d[m] = dx / 2.0 / dt * cap[m] + (cond[m - 1] + cond[m]) / 4.0 / dx - cw * (2.0 * v_new[m] + v_new[m - 1]) / 12.0 + dx / 24.0 * cw * (sink[m - 1] + 3.0 * sink[m]);
                } else {
                    f[m] = temp_o[m - 1] * ((cond[m] + cond[m - 1]) / 4.0 / dx + cw * (v_old[m] + 2.0 * v_old[m - 1]) / 12.0 - dx / 24.0 * cw * (sink[m - 1] + sink[m]))
                        + temp_o[m]
                            * (dx / 2.0 / dt * cap[m] - (cond[m - 1] + cond[m]) / 4.0 / dx + cw * (2.0 * v_old[m] + v_old[m - 1]) / 12.0
                                - dx / 24.0 * cw * (sink[m - 1] + 3.0 * sink[m]));
                    f[m] -= heat_fl;
                    f[m] -= heat_fl.min(0.0);
                }
            }
        }

        // ====================================================================
        // HIER WIRD EINGEFÜGT: TempAdj (Matrix-Korrektur für Wasserdampf)
        // ====================================================================
        if self.l_vapor {
            let cv = 1.8e6 / self.x_conv / self.t_conv.powi(2);
            let mut g0 = vec![0.0; n];
            for i_level in 1..=2 {
                for i in 0..n {
                    let v_v_grad = if i == 0 {
                        (self.v_v_old[1] - self.v_v_old[0]) / (x[1] - x[0])
                    } else if i == n - 1 {
                        (self.v_v_old[n - 1] - self.v_v_old[n - 2]) / (x[n - 1] - x[n - 2])
                    } else {
                        (self.v_v_old[i + 1] - self.v_v_old[i - 1]) / (x[i + 1] - x[i - 1]) * 2.0
                    };
                    let temp_eval = if i_level == 1 { self.heat.as_ref().unwrap().temp_o[i] } else { self.heat.as_ref().unwrap().temp_n[i] };
                    let lat = crate::vapor::latent_heat_volumetric(temp_eval) / self.x_conv / self.t_conv.powi(2);
                    let th_v_grad = (self.th_v_new[i] - self.th_v_old[i]) / dt;
                    g0[i] = -lat * (v_v_grad + th_v_grad);
                }

                // Unterer Rand
                let dx0 = x[1] - x[0];
                if k_bot < 0 {
                    if i_level == 1 {
                        f[0] += temp_o[0] * (-cv * (2.0 * self.v_v_old[0] + self.v_v_old[1]) / 12.0)
                            + temp_o[1] * (-cv * (2.0 * self.v_v_old[1] + self.v_v_old[0]) / 12.0)
                            + dx0 / 12.0 * (2.0 * g0[0] + g0[1])
                            + t_bot * cv * (self.v_v_new[0] + self.v_v_old[0]) / 2.0;
                    } else {
                        d[0] += cv * (2.0 * self.v_v_new[0] + self.v_v_new[1]) / 12.0;
                        e[0] += cv * (2.0 * self.v_v_new[1] + self.v_v_new[0]) / 12.0;
                        f[0] += dx0 / 12.0 * (2.0 * g0[0] + g0[1]);
                    }
                }

                // Innere Knoten
                for i in 1..n - 1 {
                    let dxa = x[i] - x[i - 1];
                    let dxb = x[i + 1] - x[i];
                    if i_level == 1 {
                        f[i] += temp_o[i - 1] * (cv * (self.v_v_old[i] + 2.0 * self.v_v_old[i - 1]) / 12.0)
                            + temp_o[i] * (-cv * (self.v_v_old[i + 1] - self.v_v_old[i - 1]) / 12.0)
                            + temp_o[i + 1] * (-cv * (2.0 * self.v_v_old[i + 1] + self.v_v_old[i]) / 12.0)
                            + dxa * (g0[i - 1] + 2.0 * g0[i]) / 12.0
                            + dxb * (2.0 * g0[i] + g0[i + 1]) / 12.0;
                    } else {
                        b[i] -= cv * (self.v_v_new[i] + 2.0 * self.v_v_new[i - 1]) / 12.0;
                        d[i] += cv * (self.v_v_new[i + 1] - self.v_v_new[i - 1]) / 12.0;
                        e[i] += cv * (2.0 * self.v_v_new[i + 1] + self.v_v_new[i]) / 12.0;
                        f[i] += dxa * (g0[i - 1] + 2.0 * g0[i]) / 12.0 + dxb * (2.0 * g0[i] + g0[i + 1]) / 12.0;
                    }
                }

                // Oberer Rand
                if k_top < 0 {
                    let dxm = x[n - 1] - x[n - 2];
                    let m = n - 1;
                    if i_level == 1 {
                        f[m] += temp_o[m - 1] * (cv * (self.v_v_old[m] + 2.0 * self.v_v_old[m - 1]) / 12.0)
                            + temp_o[m] * (cv * (2.0 * self.v_v_old[m] + self.v_v_old[m - 1]) / 12.0)
                            + dxm / 12.0 * (g0[m - 1] + 2.0 * g0[m]);
                        if k_top == -1 {
                            f[m] -= t_top_a * cv * (self.v_v_new[m] + self.v_v_old[m]) / 2.0;
                        }
                    } else {
                        b[m] -= cv * (self.v_v_new[m] + 2.0 * self.v_v_new[m - 1]) / 12.0;
                        d[m] -= cv * (2.0 * self.v_v_new[m] + self.v_v_new[m - 1]) / 12.0;
                        f[m] += dxm / 12.0 * (g0[m - 1] + 2.0 * g0[m]);
                    }
                }
            }
        }
        // ====================================================================

        // 3. Tridiagonales Gleichungssystem lösen
        let mut bb = b.clone();
        bb[0] = 0.0;
        let mut ee = e.clone();
        ee[n - 1] = 0.0;
        thomas(&bb, &mut d, &ee, &mut f);
        let hs = self.heat.as_mut().unwrap();
        hs.temp_n.copy_from_slice(&f);
    }
}
