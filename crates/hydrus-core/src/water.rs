//! Richards equation solver (port of WATFLOW.FOR).

use crate::material::*;
use crate::model::SoilModel;
use crate::sim::Simulation;

/// Tridiagonal system in the layout of the Fortran code.
struct Sys {
    p: Vec<f64>,
    r: Vec<f64>,
    s: Vec<f64>,
    pb: f64,
    rb: f64,
    sb: f64,
    pt: f64,
    rt: f64,
    st: f64,
}

pub fn fqh(gwl: f64, aqh: f64, bqh: f64) -> f64 {
    aqh * (bqh * gwl.abs()).exp()
}

impl crate::model::DrainSettings {
    /// Drainage flux, based on the SWAP model (FqDrain in Fortran).
    pub fn flux(&self, gwl: f64) -> f64 {
        let pi = 3.14159f64;
        let dh = gwl - self.z_bot;
        let sp = self.spacing;
        let mut zimp = 0.0;
        let mut dbot = 0.0;
        if self.position > 1 {
            zimp = self.base_gw.max(self.z_bot - 0.25 * sp);
            dbot = (self.z_bot - zimp).max(0.0);
        }
        if dh < 1e-10 {
            return 0.0;
        }
        let tot_res = match self.position {
            1 => sp * sp / (4.0 * self.kh_top * dh.abs()) + self.entrance_resistance,
            2 | 3 => {
                let x = 2.0 * pi * dbot / sp;
                let mut eqd;
                if x > 0.5 {
                    let mut fx = 0.0;
                    for i in [1.0f64, 3.0, 5.0] {
                        fx += 4.0 * (-2.0 * i * x).exp() / (i * (1.0 - (-2.0 * i * x).exp()));
                    }
                    eqd = pi * sp / 8.0 / ((sp / self.wet_perimeter).ln() + fx);
                } else if x < 1e-6 {
                    eqd = dbot;
                } else {
                    let fx = pi * pi / (4.0 * x) + (x / (2.0 * pi)).ln();
                    eqd = pi * sp / 8.0 / ((sp / self.wet_perimeter).ln() + fx);
                }
                if eqd > dbot {
                    eqd = dbot;
                }
                if self.position == 2 {
                    sp * sp / (8.0 * self.kh_top * eqd + 4.0 * self.kh_top * dh.abs()) + self.entrance_resistance
                } else {
                    sp * sp / (8.0 * self.kh_bot * eqd + 4.0 * self.kh_top * dh.abs()) + self.entrance_resistance
                }
            }
            4 => {
                let rver = (gwl - self.z_interface).max(0.0) / self.kv_top + (self.z_interface.min(gwl) - self.z_bot) / self.kv_bot;
                let rhor = sp * sp / (8.0 * self.kh_bot * dbot);
                let rrad = sp / (pi * (self.kh_bot * self.kv_bot).sqrt()) * (dbot / self.wet_perimeter).ln();
                rver + rhor + rrad + self.entrance_resistance
            }
            _ => {
                let rver = (gwl - self.z_bot) / self.kv_top;
                let rhor = sp * sp / (8.0 * self.kh_top * (self.z_bot - self.z_interface) + 8.0 * self.kh_bot * (self.z_interface - zimp));
                let rrad = sp / (pi * (self.kh_top * self.kv_top).sqrt()) * (self.geo_factor * (self.z_bot - self.z_interface) / self.wet_perimeter).ln();
                rver + rhor + rrad + self.entrance_resistance
            }
        };
        -dh / tot_res
    }
}

impl Simulation {
    fn lookup(&self, m: usize, him: f64) -> Option<(usize, f64)> {
        let t0 = &self.tabs[0];
        if self.l_table && him > t0.h[NTAB - 1] && him <= t0.h[0] {
            let alh1 = (-self.h_tab_first).log10();
            let dlh = ((-self.h_tab_last).log10() - alh1) / (NTAB as f64 - 1.0);
            let it = (((-him).log10() - alh1) / dlh).max(0.0) as usize;
            let it = it.min(NTAB - 2);
            let tb = &self.tabs[m];
            let dh = (him - tb.h[it]) / (tb.h[it + 1] - tb.h[it]);
            Some((it, dh))
        } else {
            None
        }
    }

    fn lookup_tabular(&self, m: usize, him: f64) -> (f64, f64, f64) {
        let tb = &self.tabs[m];
        let n_pts = tb.h.len();
        if n_pts == 0 {
            return (0.0, 0.0, 0.0);
        }
        if him >= tb.h[0] {
            return (tb.con[0], tb.cap[0], tb.the[0]);
        }
        if him <= tb.h[n_pts - 1] {
            return (tb.con[n_pts - 1], tb.cap[n_pts - 1], tb.the[n_pts - 1]);
        }

        // Bisection lookup in descending pressure head table
        let mut low = 0;
        let mut high = n_pts - 1;
        while high - low > 1 {
            let mid = (low + high) / 2;
            if tb.h[mid] <= him {
                high = mid;
            } else {
                low = mid;
            }
        }

        let dh = (him - tb.h[low]) / (tb.h[high] - tb.h[low]);
        let coni = tb.con[low] + (tb.con[high] - tb.con[low]) * dh;
        let capi = tb.cap[low] + (tb.cap[high] - tb.cap[low]) * dh;
        let thei = tb.the[low] + (tb.the[high] - tb.the[low]) * dh;
        (coni.max(1e-37), capi.max(0.0), thei.max(0.0))
    }
	
	/// Evaluate matrix hydraulic properties K_m(h_m) and C_m(h_m).
    pub fn set_mat_matrix(&mut self) {
        for i in 0..self.n {
            let m = self.mat[i];
            let par = &self.par_d[m];

            // If extra matrix params were not specified, fallback to primary soil retention params
            let (qr_m, qs_m, alpha_m, n_m, ks_m) = if par[7] > par[6] {
                (par[6], par[7], par[8], par[9], par[10])
            } else {
                (par[0], par[1], par[2], par[3], par[4])
            };

            let par_m: Par = [
                qr_m,
                qs_m,
                if alpha_m > 0.0 { alpha_m } else { par[2] },
                if n_m > 1.0 { n_m } else { par[3] },
                if ks_m > 0.0 { ks_m } else { par[4] },
                par[5],
                0.0, 0.0, 0.0, 0.0, 0.0,
            ];

            let hm = self.h_matrix_new[i].min(100.0);
            self.con_matrix[i] = fk(SoilModel::VanGenuchten, hm, &par_m).max(1e-20);
            self.cap_matrix[i] = fc(SoilModel::VanGenuchten, hm, &par_m).max(1e-12);
            let q_min = qr_m.min(qs_m);
            let q_max = qr_m.max(qs_m);
            self.th_matrix_new[i] = fq(SoilModel::VanGenuchten, hm, &par_m).clamp(q_min, q_max);
            if self.th_matrix_old[i] == 0.0 {
                self.th_matrix_old[i] = self.th_matrix_new[i];
            }
        }
    }

    /// Assemble tridiagonal system for the matrix domain.
    fn build_matrix_system(&mut self) -> Sys {
        let n = self.n;
        let dt = self.dt.max(1e-10);
        let grav = self.cos_alf;
        let w_m = (1.0 - self.w_fracture).clamp(0.001, 0.999);

        let mut p = vec![0.0; n];
        let mut r = vec![0.0; n];
        let mut s = vec![0.0; n];

        // Bottom node
        let dxb = (self.x[1] - self.x[0]).abs().max(1e-6);
        let dx = dxb / 2.0;
        let conb = (self.con_matrix[0] + self.con_matrix[1]) / 2.0;
        let b = conb * grav;
        s[0] = -conb / dxb;
        let f2 = (self.cap_matrix[0] * dx / dt).max(1e-12);
        let rb = conb / dxb + f2;
        let sb = -conb / dxb;
        let r_bot_m = if self.free_d { -conb * grav } else { 0.0 };

        let pb = b + f2 * self.h_matrix_new[0]
            - (self.th_matrix_new[0] - self.th_matrix_old[0]) * dx / dt
            + r_bot_m
            + (self.sink_im[0] / w_m) * dx;

        // Interior nodes
        for i in 1..n - 1 {
            let dxa = (self.x[i] - self.x[i - 1]).abs().max(1e-6);
            let dxb = (self.x[i + 1] - self.x[i]).abs().max(1e-6);
            let dx = (dxa + dxb) / 2.0;
            let cona = (self.con_matrix[i] + self.con_matrix[i - 1]) / 2.0;
            let conb = (self.con_matrix[i] + self.con_matrix[i + 1]) / 2.0;
            let b = (cona - conb) * grav;
            let a2 = cona / dxa + conb / dxb;
            let a3 = -conb / dxb;
            let f2 = (self.cap_matrix[i] * dx / dt).max(1e-12);
            r[i] = a2 + f2;
            p[i] = f2 * self.h_matrix_new[i]
                - (self.th_matrix_new[i] - self.th_matrix_old[i]) * dx / dt
                - b
                + (self.sink_im[i] / w_m) * dx;
            s[i] = a3;
        }

        // Top node
        let m = n - 1;
        let dxa = (self.x[m] - self.x[m - 1]).abs().max(1e-6);
        let dx = dxa / 2.0;
        let cona = (self.con_matrix[m] + self.con_matrix[m - 1]) / 2.0;
        let b = cona * grav;
        let f2 = (self.cap_matrix[m] * dx / dt).max(1e-12);
        let rt = cona / dxa + f2;
        let st = -cona / dxa;
        let pt = f2 * self.h_matrix_new[m]
            - (self.th_matrix_new[m] - self.th_matrix_old[m]) * dx / dt
            - b
            + (self.sink_im[m] / w_m) * dx;

        Sys { p, r, s, pb, rb, sb, pt, rt, st }
    }

    /// Solve the tridiagonal matrix system for the matrix domain.
    fn solve_matrix(&mut self, sys: &Sys) {
        let n = self.n;
        let rmin = 1e-100;
        let mut diag = sys.r.clone();
        let mut rhs = sys.p.clone();
        let upper = sys.s.clone();
        let mut lower = vec![0.0; n];
        for i in 1..n {
            lower[i] = sys.s[i - 1];
        }
        diag[0] = sys.rb;
        rhs[0] = sys.pb;
        diag[n - 1] = sys.rt;
        lower[n - 1] = sys.st;
        rhs[n - 1] = sys.pt;

        for i in 1..n {
            let mut d = diag[i - 1];
            if d.abs() < rmin {
                d = rmin.copysign(d);
            }
            let f = lower[i] / d;
            diag[i] -= f * upper[i - 1];
            rhs[i] -= f * rhs[i - 1];
        }
        let mut dn = diag[n - 1];
        if dn.abs() < rmin {
            dn = rmin.copysign(dn);
        }
        self.h_matrix_new[n - 1] = (rhs[n - 1] / dn).clamp(-1e5, 100.0);
        for i in (0..n - 1).rev() {
            let mut d = diag[i];
            if d.abs() < rmin {
                d = rmin.copysign(d);
            }
            self.h_matrix_new[i] = ((rhs[i] - upper[i] * self.h_matrix_new[i + 1]) / d).clamp(-1e5, 100.0);
        }
    }

    /// Hydraulic properties at every node (SetMat in Fortran).
    pub fn set_mat(&mut self, iter: usize) {
        let model = self.model;
        for i in 0..self.n {
            let m = self.mat[i];
            
            // In WATFLOW.FOR: Effective scaling factor A = ah(i) * (if kappa == 1 { ahW(m) } else { 1.0 })
            // However, when individual nodal scaling is active (ah != 1.0), ah already incorporates the scaling.
            let a_scale = if self.kappa[i] == 1 && (self.ah[i] - 1.0).abs() < 1e-5 {
                self.ah[i] * self.ah_w[m]
            } else {
                self.ah[i]
            };

            let hi1 = self.h_sat[m].min(self.h_temp[i] / a_scale);
            let hi2 = self.h_sat[m].min(self.h_new[i] / a_scale);
            let him = 0.1 * hi1 + 0.9 * hi2;
			
            if self.model == SoilModel::Tabular {
                let (coni, capi, thei) = self.lookup_tabular(m, him);
                self.con[i] = coni * self.ak[i] * self.ak_s[i];
                self.cap[i] = capi * self.ath[i] * self.ath_s[i];
                self.th_eq[i] = thei * self.ath[i] * self.ath_s[i];
                if iter == 0 {
                    self.con_o[i] = self.con[i];
                }
                continue;
            }

            let coni;
            if hi1 >= self.h_sat[m] && hi2 >= self.h_sat[m] {
                coni = self.con_sat[m];
            } else if let Some((it, dh)) = self.lookup(m, him) {
                let tb = &self.tabs[m];
                coni = tb.con[it] + (tb.con[it + 1] - tb.con[it]) * dh;
            } else {
                coni = fk(model, him, &self.par_d[m]);
            }

            let (capi, thei);
            if him >= self.h_sat[m] {
                capi = 0.0;
                thei = self.ths[m];
            } else if let Some((it, dh)) = self.lookup(m, him) {
                let tb = &self.tabs[m];
                capi = tb.cap[it] + (tb.cap[it + 1] - tb.cap[it]) * dh;
                thei = tb.the[it] + (tb.the[it + 1] - tb.the[it]) * dh;
            } else {
                capi = fc(model, him, &self.par_d[m]);
                thei = fq(model, him, &self.par_d[m]);
            }

            let (at, bt) = (1.0, 1.0);
            if self.kappa[i] == -1 {
                self.con[i] = coni * self.ak[i] * bt * self.ak_s[i];
                self.cap[i] = capi * self.ath[i] * self.ath_s[i] / self.ah[i] / at;
                self.th_eq[i] = self.thr[m] + (thei - self.thr[m]) * self.ath[i] * self.ath_s[i];
            } else {
                // When individual nodal scaling is active, Dxz is already scaled.
                let d_scale = if (self.ath[i] - 1.0).abs() < 1e-5 {
                    self.ath_w[m] * self.ath[i]
                } else {
                    self.ath[i]
                };
                let k_scale = if (self.ak[i] - 1.0).abs() < 1e-5 {
                    self.ak_w[m] * self.ak[i]
                } else {
                    self.ak[i]
                };

                self.con[i] = self.con_r[i] + coni * k_scale * bt * self.ak_s[i];
                self.cap[i] = capi * d_scale * self.ath_s[i] / a_scale / at;
                self.th_eq[i] = self.th_rr[i] + d_scale * self.ath_s[i] * (thei - self.thr[m]);
            }
            if iter == 0 {
                self.con_o[i] = self.con[i];
            }
        }
		if self.l_vapor || self.prj.water.l_w_dep {
            let temps: Vec<f64> = (0..self.n).map(|i| self.temp(i)).collect();
            crate::vapor::con_vapor(
                self.n,
                &self.mat,
                &self.h_new,
                &temps,
                &self.con,
                &self.th_eq,
                &self.ths,
                &mut self.con_lt,
                &mut self.con_vt,
                &mut self.con_vh,
                self.x_conv,
                self.t_conv,
                self.l_vapor,
                1,
            );
            if self.l_vapor {
                crate::vapor::vapor_content(
                    self.n,
                    &self.mat,
                    &self.th_eq,
                    &mut self.th_v_new,
                    &temps,
                    &self.h_new,
                    &self.ths,
                    &mut self.cap,
                    self.x_conv,
                );
            }
        }
    }

    /// Kool & Parker style reversal detection and scanning-curve parameters (Hyster in WATFLOW.FOR).
    fn hyster_kp(&mut self) {
        let model = self.model;
        for i in 0..self.n {
            self.kappa_o[i] = self.kappa[i];
            if (self.th_new[i] - self.th_old[i]) * self.kappa[i] as f64 >= -self.tol_th {
                continue;
            }
            self.kappa[i] = -self.kappa[i];
            let m = self.mat[i];
            let thr = self.par_d[m][0];
            let ths_d = self.par_d[m][1];
            let ths_w = self.par_w[m][1];
            let ks_d = self.par_d[m][4];
            let ks_w = self.par_w[m][4];
            let mut ths = ths_d;
            let mut ks = ks_d;
            if self.kappa[i] == 1 {
                if ths_w < 0.999 * ths_d {
                    let rr = 1.0 / (ths_d - ths_w) - 1.0 / (ths_d - thr);
                    ths = ths_d - (ths_d - self.th_old[i]) / (1.0 + rr * (ths_d - self.th_old[i]));
                }
                if ks_w < 0.999 * ks_d {
                    let rr = 1.0 / (ks_d - ks_w) - 1.0 / ks_d;
                    ks = ks_d - (ks_d - self.con_o[i]) / (1.0 + rr * (ks_d - self.con_o[i]));
                }
            }
            let hh = self.h_old[i] / self.ah[i];
            if self.kappa[i] == 1 {
                self.ath_s[i] = 1.0;
                let sew = if self.model == SoilModel::Tabular {
                    let (_, _, th_val) = self.lookup_tabular(m, hh);
                    if ths_w > thr { ((th_val - thr) / (ths_w - thr)).clamp(0.0, 1.0) } else { 1.0 }
                } else {
                    fs(model, hh, &self.par_w[m])
                };
                if sew < 0.999 {
                    self.ath_s[i] = (self.th_old[i] - ths) / (1.0 - sew) / (thr - ths_w);
                }
                self.th_rr[i] = ths - self.ath_s[i] * (ths_w - thr);
                self.ak_s[i] = 1.0;
                self.con_r[i] = 0.0;
                if self.i_hyst == 2 {
                    let kw = if self.model == SoilModel::Tabular {
                        let (k_val, _, _) = self.lookup_tabular(m, hh);
                        self.ak[i] * k_val
                    } else {
                        self.ak[i] * fk(model, hh, &self.par_w[m])
                    };
                    if kw < 0.999 * ks_w {
                        self.ak_s[i] = (self.con_o[i] - ks) / (kw - ks_w);
                    }
                    self.con_r[i] = ks - self.ak_s[i] * ks_w;
                }
            } else {
                let se_d = if self.model == SoilModel::Tabular {
                    let (_, _, th_val) = self.lookup_tabular(m, hh);
                    if ths_d > thr { ((th_val - thr) / (ths_d - thr)).clamp(1e-6, 1.0) } else { 1.0 }
                } else {
                    fs(model, hh, &self.par_d[m]).max(1e-6)
                };
                self.ath_s[i] = (self.th_old[i] - thr) / se_d / (ths_d - thr);
                self.th_rr[i] = thr;
                self.ak_s[i] = 1.0;
                self.con_r[i] = 0.0;
                if self.i_hyst == 2 {
                    let kd = if self.model == SoilModel::Tabular {
                        let (k_val, _, _) = self.lookup_tabular(m, hh);
                        k_val
                    } else {
                        fk(model, hh, &self.par_d[m])
                    };
                    if kd > 1e-30 {
                        self.ak_s[i] = self.con_o[i] / kd / self.ak[i];
                    }
                }
            }
		}
	}

    /// Velocity at node i using the Fortran "Veloc" formula.
    pub fn veloc(
        &self,
        h: &[f64],
        th_new: &[f64],
        th_old: &[f64],
        temp: &[f64],
    ) -> (Vec<f64>, Vec<f64>) {
        let n = self.n;
        let x = &self.x;
        let con = &self.con;
        let g = self.cos_alf;
        let dt = self.dt;
        let mut v = vec![0.0; n];
        let mut v_v = vec![0.0; n];

        // Top node (m = n - 1)
        let m = n - 1;
        let dxn = x[m] - x[m - 1];
        v[m] = -(con[m] + con[m - 1]) / 2.0 * ((h[m] - h[m - 1]) / dxn + g)
            - dxn / 2.0 * ((th_new[m] - th_old[m]) / dt + self.sink[m]);

        if self.prj.water.l_w_dep {
            v[m] -= (self.con_lt[m] + self.con_lt[m - 1]) / 2.0 * (temp[m] - temp[m - 1]) / dxn;
        }

        if self.l_vapor {
            v_v[m] = -(self.con_vh[m] + self.con_vh[m - 1]) / 2.0 * (h[m] - h[m - 1]) / dxn
                - (self.con_vt[m] + self.con_vt[m - 1]) / 2.0 * (temp[m] - temp[m - 1]) / dxn
                - dxn / 2.0 * (self.th_v_new[m] - self.th_v_old[m]) / dt;
        }

        // Interior nodes
        for i in 1..n - 1 {
            let dxa = x[i + 1] - x[i];
            let dxb = x[i] - x[i - 1];

            let va = -(con[i] + con[i + 1]) / 2.0 * ((h[i + 1] - h[i]) / dxa + g);
            let vb = -(con[i] + con[i - 1]) / 2.0 * ((h[i] - h[i - 1]) / dxb + g);
            let mut vi = (va * dxb + vb * dxa) / (dxa + dxb);

            if self.prj.water.l_w_dep {
                let v_ta = -(self.con_lt[i] + self.con_lt[i + 1]) / 2.0 * (temp[i + 1] - temp[i]) / dxa;
                let v_tb = -(self.con_lt[i] + self.con_lt[i - 1]) / 2.0 * (temp[i] - temp[i - 1]) / dxb;
                vi += (v_ta * dxb + v_tb * dxa) / (dxa + dxb);
            }
            v[i] = vi;

            if self.l_vapor {
                let mut v_va = -(self.con_vh[i] + self.con_vh[i + 1]) / 2.0 * (h[i + 1] - h[i]) / dxa;
                let mut v_vb = -(self.con_vh[i] + self.con_vh[i - 1]) / 2.0 * (h[i] - h[i - 1]) / dxb;
                v_va -= (self.con_vt[i] + self.con_vt[i + 1]) / 2.0 * (temp[i + 1] - temp[i]) / dxa;
                v_vb -= (self.con_vt[i] + self.con_vt[i - 1]) / 2.0 * (temp[i] - temp[i - 1]) / dxb;
                v_v[i] = (v_va * dxb + v_vb * dxa) / (dxa + dxb);
            }
        }

        // Bottom node (i = 0)
        let dx1 = x[1] - x[0];
        v[0] = -(con[0] + con[1]) / 2.0 * ((h[1] - h[0]) / dx1 + g)
            + dx1 / 2.0 * ((th_new[0] - th_old[0]) / dt + self.sink[0]);

        if self.prj.water.l_w_dep {
            v[0] -= (self.con_lt[0] + self.con_lt[1]) / 2.0 * (temp[1] - temp[0]) / dx1;
        }

        if self.l_vapor {
            v_v[0] = -(self.con_vh[0] + self.con_vh[1]) / 2.0 * (h[1] - h[0]) / dx1
                - (self.con_vt[0] + self.con_vt[1]) / 2.0 * (temp[1] - temp[0]) / dx1
                + dx1 / 2.0 * (self.th_v_new[0] - self.th_v_old[0]) / dt;
        }

        (v, v_v)
    }

    fn build_system(&mut self) -> Sys {
        let n = self.n;
        let dt = self.dt;
        let grav = self.cos_alf;

        for i in 0..n {
            self.th_new[i] = self.th_old[i] + (self.th_eq[i] - self.th_old[i]);
        }

        let mut p = vec![0.0; n];
        let mut r = vec![0.0; n];
        let mut s = vec![0.0; n];

        // --- Bottom BC ---
        let dxb = self.x[1] - self.x[0];
        let dx = dxb / 2.0;
        let mut conb = (self.con[0] + self.con[1]) / 2.0;
        let b = conb * grav;
		
        if self.l_vapor {
            conb += (self.con_vh[0] + self.con_vh[1]) / 2.0;
        }
        s[0] = -conb / dxb;
        if self.free_d {
            self.r_bot = -conb * grav;
        }
        let f2 = self.cap[0] * dx / dt;
        let rb = conb / dxb + f2;
        let sb = -conb / dxb;
        if self.gwl_f {
            self.r_bot = crate::water::fqh(self.h_new[0] - self.gwl0l, self.aqh, self.bqh);
        }
        if let Some(d) = &self.prj.water.bc.drains {
            self.r_bot = d.flux(self.x[0] + self.h_new[0]);
        }
		let w_f = self.w_fracture.max(0.001);
        let mut pb = b - self.sink[0] * dx + f2 * self.h_new[0] - (self.th_new[0] - self.th_old[0]) * dx / dt + self.r_bot;
		if self.i_dual_por > 0 {
			pb -= self.sink_im[0] * dx;
		}
		if self.l_dual_perm {
            pb -= (self.sink_im[0] / w_f) * dx;
        }
        if self.l_vapor || self.prj.water.l_w_dep {
            let mut con_tb = 0.0;
            if self.l_vapor {
                con_tb += (self.con_vt[0] + self.con_vt[1]) / 2.0;
            }
            if self.prj.water.l_w_dep {
                con_tb += (self.con_lt[0] + self.con_lt[1]) / 2.0;
            }
            pb += con_tb * (self.temp(1) - self.temp(0)) / dxb - (self.th_v_new[0] - self.th_v_old[0]) * dx / dt;
        }

        // --- Interior Nodes ---
        for i in 1..n - 1 {
            let dxa = self.x[i] - self.x[i - 1];
            let dxb = self.x[i + 1] - self.x[i];
            let dx = (dxa + dxb) / 2.0;
            let mut cona = (self.con[i] + self.con[i - 1]) / 2.0;
            let mut conb = (self.con[i] + self.con[i + 1]) / 2.0;
            let b = (cona - conb) * grav;
            if self.l_vapor {
                cona += (self.con_vh[i] + self.con_vh[i - 1]) / 2.0;
                conb += (self.con_vh[i] + self.con_vh[i + 1]) / 2.0;
            }
            let a2 = cona / dxa + conb / dxb;
            let a3 = -conb / dxb;
            let f2 = self.cap[i] * dx / dt;
            r[i] = a2 + f2;
            p[i] = f2 * self.h_new[i] - (self.th_new[i] - self.th_old[i]) * dx / dt - b - self.sink[i] * dx;
			if self.i_dual_por > 0 { p[i] -= self.sink_im[i] * dx; }
			if self.l_dual_perm   { p[i] -= (self.sink_im[i] / w_f) * dx; }
            s[i] = a3;
            if self.l_vapor || self.prj.water.l_w_dep {
                let mut con_ta = 0.0;
                let mut con_tb = 0.0;
                if self.l_vapor {
                    con_ta += (self.con_vt[i] + self.con_vt[i - 1]) / 2.0;
                    con_tb += (self.con_vt[i] + self.con_vt[i + 1]) / 2.0;
                }
                if self.prj.water.l_w_dep {
                    con_ta += (self.con_lt[i] + self.con_lt[i - 1]) / 2.0;
                    con_tb += (self.con_lt[i] + self.con_lt[i + 1]) / 2.0;
                }
                p[i] += con_tb * (self.temp(i + 1) - self.temp(i)) / dxb
                    - con_ta * (self.temp(i) - self.temp(i - 1)) / dxa
                    - (self.th_v_new[i] - self.th_v_old[i]) * dx / dt;
            }
        }

        // --- Top BC ---
        let m = n - 1;
        let dxa = self.x[m] - self.x[m - 1];
        let dx = dxa / 2.0;
        let mut cona = (self.con[m] + self.con[m - 1]) / 2.0;
        let b = cona * grav;
		
        if self.l_vapor {
            cona += (self.con_vh[m] + self.con_vh[m - 1]) / 2.0;
        }
        let f2 = self.cap[m] * dx / dt;
        let mut rt = cona / dxa + f2;
        let st = -cona / dxa;
        let mut pt = f2 * self.h_new[m] - (self.th_new[m] - self.th_old[m]) * dx / dt - self.sink[m] * dx - b;
		if self.i_dual_por > 0 {
            pt -= self.sink_im[m] * dx;
        }
		if self.l_dual_perm {
            pt -= (self.sink_im[m] / w_f) * dx;
        }
        if self.l_vapor || self.prj.water.l_w_dep {
            let mut con_ta = 0.0;
            if self.l_vapor {
                con_ta += (self.con_vt[m] + self.con_vt[m - 1]) / 2.0;
            }
            if self.prj.water.l_w_dep {
                con_ta += (self.con_lt[m] + self.con_lt[m - 1]) / 2.0;
            }
            pt -= con_ta * (self.temp(m) - self.temp(m - 1)) / dxa + (self.th_v_new[m] - self.th_v_old[m]) * dx / dt;
        }
        self.v_top = -st * self.h_new[m - 1] - rt * self.h_new[m] + pt;
        pt -= self.r_top;
        if self.w_layer {
            if self.h_new[m] > 0.0 {
                rt += 1.0 / dt;
            }
            pt += self.h_old[m].max(0.0) / dt;
        }

        Sys { p, r, s, pb, rb, sb, pt, rt, st }
    }

    fn solve(&mut self, sys: &Sys) {
        let n = self.n;
        let rmin = 1e-100;

        let mut p = sys.p.clone();
        let mut r = sys.r.clone();
        let s = &sys.s;

        let pb = sys.pb;
        let mut rb = sys.rb;
        let sb = sys.sb;

        let pt = sys.pt;
        let rt = sys.rt;
        let st = sys.st;

        // Forward
        if self.kod_bot >= 0 {
            p[1] -= s[0] * self.h_bot;
        } else {
            if rb.abs() < rmin {
                rb = rmin;
            }
            p[1] -= pb * s[0] / rb;
            r[1] -= sb * s[0] / rb;
        }

        for i in 2..n - 1 {
            let mut r_prev = r[i - 1];
            if r_prev.abs() < rmin {
                r_prev = rmin;
            }
            p[i] -= p[i - 1] * s[i - 1] / r_prev;
            r[i] -= s[i - 1] * s[i - 1] / r_prev;
        }

        if self.kod_top > 0 {
            p[n - 2] -= s[n - 2] * self.h_top;
        } else {
            let mut r_nm2 = r[n - 2];
            if r_nm2.abs() < rmin {
                r_nm2 = rmin;
            }
            p[n - 1] = pt - p[n - 2] * st / r_nm2;
            r[n - 1] = rt - s[n - 2] * st / r_nm2;
        }

        // Back
        let mut r_last = r[n - 2];
        if r_last.abs() < rmin {
            r_last = rmin;
        }

        if self.kod_top > 0 {
            self.h_new[n - 1] = self.h_top;
            self.h_new[n - 2] = p[n - 2] / r_last;
        } else {
            let mut r_nm1 = r[n - 1];
            if r_nm1.abs() < rmin {
                r_nm1 = rmin;
            }
            self.h_new[n - 1] = p[n - 1] / r_nm1;
            self.h_new[n - 2] = (p[n - 2] - s[n - 2] * self.h_new[n - 1]) / r_last;
        }

        for i in (1..n - 2).rev() {
            let mut r_i = r[i];
            if r_i.abs() < rmin {
                r_i = rmin;
            }
            self.h_new[i] = (p[i] - s[i] * self.h_new[i + 1]) / r_i;
        }

        if self.kod_bot >= 0 {
            self.h_new[0] = self.h_bot;
        } else {
            if rb.abs() < rmin {
                rb = rmin;
            }
            self.h_new[0] = (pb - sb * self.h_new[1]) / rb;
        }
    }

    /// Change boundary condition types depending on the solution (Shift in Fortran).
    fn shift(&mut self) {
        let n = self.n;
        let grav = self.cos_alf;
		
        if self.seep_f {
            let dx = self.x[1] - self.x[0];
            let v_bot = -(self.con[0] + self.con[1]) / 2.0 * ((self.h_new[1] - self.h_new[0]) / dx + grav)
                - dx / 2.0 * ((self.th_new[0] - self.th_old[0]) / self.dt + self.sink[0]);
            if self.kod_bot >= 0 {
                if v_bot > 0.0 {
                    self.kod_bot = -2;
                    self.r_bot = 0.0;
                }
            } else if self.h_new[0] >= self.h_seep {
                self.kod_bot = 2;
                self.h_bot = self.h_seep;
            }
        }
        if self.top_inf && (self.kod_top.abs() == 4 || (self.kod_top.abs() == 1 && self.r_top > 0.0)) {
            if self.kod_top > 0 {
                let m = n - 2;
                let dx = self.x[n - 1] - self.x[m];
                let mut v_top = -(self.con[n - 1] + self.con[m]) / 2.0 * ((self.h_new[n - 1] - self.h_new[m]) / dx + grav)
                    - (self.th_new[n - 1] - self.th_old[n - 1]) * dx / 2.0 / self.dt
                    - self.sink[n - 1] * dx / 2.0;

                if self.i_dual_por > 0 {
					v_top -= self.sink_im[n - 1] * dx / 2.0;
				}
				if self.prj.water.l_w_dep {
                    v_top -= (self.con_lt[n - 1] + self.con_lt[m]) / 2.0 * (self.temp(n - 1) - self.temp(m)) / dx;
                }
                if self.l_vapor {
                    v_top -= (self.con_vh[n - 1] + self.con_vh[m]) / 2.0 * (self.h_new[n - 1] - self.h_new[m]) / dx
                        + (self.con_vt[n - 1] + self.con_vt[m]) / 2.0 * (self.temp(n - 1) - self.temp(m)) / dx;
                }

                if v_top.abs() > self.r_top.abs() || v_top * self.r_top <= 0.0 {
                    if self.kod_top.abs() == 4 {
                        self.kod_top = -4;
                    }
                }
                if self.kod_top == 4 && self.h_new[n - 1] <= 0.99 * self.h_crit_a && self.r_top < 0.0 {
                    self.kod_top = -4;
                }
            } else {
                if !self.w_layer && self.h_new[n - 1] > 0.0 {
                    if self.kod_top.abs() == 4 {
                        self.kod_top = 4;
                    }
                    if self.kod_top.abs() == 1 {
                        self.kod_top = 1;
                    }
                    self.h_top = 0.0;
                }
                if self.h_new[n - 1] <= self.h_crit_a {
                    if self.kod_top.abs() == 4 {
                        self.kod_top = 4;
                    }
                    if self.kod_top.abs() == 1 {
                        self.kod_top = 1;
                    }
                    self.h_top = self.h_crit_a;
                }
            }
        }
    }

    /// Solve the Richards equation for one time step (WatFlow in Fortran).
    pub fn wat_flow(&mut self) {
        let n = self.n;
        let rmax = 1e10;
        self.it_cum_start();
        'outer: loop {
            self.iter_w = 0;
            self.convg = true;
			if self.i_dual_por > 0 {
                self.dual_por();
            }
            if self.l_dual_perm {
                self.dual_perm();
            }
            if self.w_layer && self.h_new[n - 1] > 0.0 && self.h_new[n - 1] < 0.00005 * self.x_conv && self.r_top >= 0.0 {
                let m = self.mat[n - 1];
                let hh = fh(self.model, 0.9999, &self.par_d[m]);
                self.h_new[n - 1] = hh;
                self.h_old[n - 1] = hh;
                self.h_temp[n - 1] = hh;
            }
            loop {
                if self.i_hyst == 3 {
                    self.lenhard_hyst(0, 2, crate::lenhard::ThetaTarget::Eq);
                } else {
                    self.set_mat(self.iter_w);
                    if self.iter_w == 2 && self.i_hyst > 0 {
                        self.hyster_kp();
                    }
                }
                let sys = self.build_system();
				self.shift();
				self.h_temp.copy_from_slice(&self.h_new);
                self.solve(&sys);
				if self.l_dual_perm {
                    self.set_mat_matrix();
                    let sys_m = self.build_matrix_system();
                    self.solve_matrix(&sys_m);
                    self.dual_perm();
                }
                for i in 0..n {
                    if self.h_new[i].abs() > rmax {
                        self.h_new[i] = rmax.copysign(self.h_new[i]);
                    }
                    if self.kod_top.abs() == 4 && self.h_new[i] < self.h_crit_a && i == n - 1 {
                        self.h_new[i] = self.h_crit_a;
                    }
                    if self.kod_top.abs() == 4 && self.h_new[i] < self.h_crit_a && (i as f64 + 1.0) > n as f64 * 9.0 / 10.0 && self.sink[i] <= 0.0 {
                        self.h_new[i] = self.h_crit_a;
                    }
                }
                self.iter_w += 1;
                self.it_cum += 1;
                // convergence test
                let mut it_crit = true;
                for i in 0..n {
                    let m = self.mat[i];
                    let mut eps_th = 0.0;
                    let mut eps_h = 0.0;
                    if self.h_temp[i] < self.h_sat[m] && self.h_new[i] < self.h_sat[m] {
                        let th = self.th_new[i]
                            + self.cap[i] * (self.h_new[i] - self.h_temp[i]) / (self.ths[m] - self.thr[m]) / self.ath[i];
                        eps_th = (self.th_new[i] - th).abs();
                    } else {
                        eps_h = (self.h_new[i] - self.h_temp[i]).abs();
                    }
                    if eps_th > self.tol_th || eps_h > self.tol_h || self.h_new[i].abs() > rmax * 0.999 {
                        it_crit = false;
                        if self.h_new[i].abs() > rmax * 0.999 {
                            self.iter_w = self.max_it;
                        }
                        break;
                    }
                }
                if !it_crit || self.iter_w <= 1 || (self.iter_w <= 2 && self.i_hyst > 0) {
                    if self.iter_w < self.max_it {
                        continue;
                    } else if self.dt <= self.dt_min {
                        self.convg = false;
                        return;
                    } else {
                        for i in 0..n {
                            if self.i_hyst > 0 {
                                self.kappa[i] = self.kappa_o[i];
                            }
                            self.h_new[i] = self.h_old[i];
                            self.h_temp[i] = self.h_old[i];
                            if let Some(h) = self.heat.as_mut() {
                                h.temp_n[i] = h.temp_o[i];
                            }
                        }
                        self.kod_top = self.k_top_old;
                        self.kod_bot = self.k_bot_old;
                        self.dt = (self.dt / 3.0).max(self.dt_min);
                        self.dt_opt = self.dt;
                        self.t = self.t_old + self.dt;
                        continue 'outer;
                    }
                }
                break;
            }
            for i in 0..n {
                self.th_new[i] += self.cap[i] * (self.h_new[i] - self.h_temp[i]);
            }
            if self.l_dual_perm {
                for i in 0..n {
                    let m = self.mat[i];
                    let par = &self.par_d[m];
                    let (qr_m, qs_m) = if par[7] > par[6] {
                        (par[6], par[7])
                    } else {
                        (par[0], par[1])
                    };
                    let q_min = qr_m.min(qs_m);
                    let q_max = qr_m.max(qs_m);
                    self.th_matrix_new[i] = (self.th_matrix_new[i] + self.cap_matrix[i] * (self.h_matrix_new[i] - self.h_matrix_old[i]))
                        .clamp(q_min, q_max);
                }
            }
            if self.w_layer && self.h_new[n - 1] > self.h_crit_s {
                self.kod_top = 4;
                self.h_top = self.h_crit_s;
            }
            if self.i_hyst == 3 {
                self.lenhard_hyst(0, 3, crate::lenhard::ThetaTarget::New);
            }
            return;
        }
    }
	
	/// Mass transfer between mobile and immobile water domains (DualPor in WATFLOW.FOR).
    pub fn dual_por(&mut self) {
        if self.i_dual_por == 0 {
            return;
        }
        let n = self.n;
        let dt = self.dt;
        self.w_transf = 0.0;

        for i in 0..n {
            let m = self.mat[i];
            let par = &self.par_d[m];
            let thr_m = self.thr[m];
            let ths_m = self.ths[m];
            let thr_im = par[6];
            let ths_im = par[7];

            let se_im = ((self.th_old_im[i] - thr_im) / (ths_im - thr_im)).clamp(0.0, 1.0);

            if self.i_dual_por == 1 {
                // Water Content driven exchange (Model 6)
                let se = ((self.th_old[i] - thr_m) / (ths_m - thr_m)).clamp(0.0, 1.0);
                let omega = par[8];
                self.sink_im[i] = omega * (se - se_im);

                let delta_th = (se - se_im) / (ths_m - thr_m + ths_im - thr_im)
                    * (ths_m - thr_m)
                    * (ths_im - thr_im);

                if self.sink_im[i] > 0.0 {
                    let mut tr_max_im = (ths_im - self.th_old_im[i]) / dt;
                    tr_max_im = tr_max_im.min(delta_th / dt);
                    if self.sink_im[i] > tr_max_im {
                        self.sink_im[i] = tr_max_im;
                    }
                } else if self.sink_im[i] < 0.0 {
                    let mut tr_max_im = -(ths_m - self.th_old[i]) / dt;
                    tr_max_im = tr_max_im.max(delta_th / dt);
                    if self.sink_im[i] < tr_max_im {
                        self.sink_im[i] = tr_max_im;
                    }
                }
            } else if self.i_dual_por == 2 {
                // Construct the matrix material parameters Par matching Fortran Par(1..6)
                let par_im: Par = [
                    par[6], // thr_im
                    par[7], // ths_im
                    par[8], // Alfa_im
                    par[9], // n_im
                    par[10], // Omega (used as Ks in the matrix retention/conductance)
                    par[5], // l (same tortuosity/connectivity as fracture)
                    0.0, 0.0, 0.0, 0.0, 0.0,
                ];

                // Immobile matrix head h_im from Se_im
                let h_im = fh(SoilModel::VanGenuchten, se_im, &par_im);

                // Immobile matrix conductivity and fracture conductivity evaluated with par_im
                let cond_m = fk(SoilModel::VanGenuchten, h_im, &par_im);
                let cond_f = fk(SoilModel::VanGenuchten, self.h_new[i], &par_im);

                self.sink_im[i] = 0.5 * (cond_m + cond_f) * (self.h_new[i] - h_im);

                // Near surface cutoff under strong evaporation
                if i == n - 1
                    && (self.h_crit_a - self.h_new[n - 1]).abs() < -0.001 * self.h_crit_a
                    && (self.h_crit_a - h_im).abs() < -0.01 * self.h_crit_a
                {
                    self.sink_im[i] = 0.0;
                }

                // Bounds checking against matrix/fracture capacities
                if self.sink_im[i] > 0.0 {
                    let tr_max_im = (ths_im - self.th_old_im[i]) / dt;
                    if self.sink_im[i] > tr_max_im {
                        self.sink_im[i] = tr_max_im;
                    }
                } else if self.sink_im[i] < 0.0 {
                    let tr_max_im = -(ths_m - self.th_old[i]) / dt;
                    if self.sink_im[i] < tr_max_im {
                        self.sink_im[i] = tr_max_im;
                    }
                }
            }

            self.th_new_im[i] = (self.th_old_im[i] + self.sink_im[i] * dt).clamp(thr_im, ths_im);

            if i >= 1 {
                self.w_transf += (self.sink_im[i - 1] + self.sink_im[i]) / 2.0 * (self.x[i] - self.x[i - 1]);
            }
        }
    }
	
	/// Calculate water transfer Gamma_w between fracture and matrix domains.
    pub fn dual_perm(&mut self) {
        if !self.l_dual_perm {
            return;
        }
        let n = self.n;
        self.w_transf = 0.0;

        for i in 0..n {
            let m = self.mat[i];
            let par = &self.par_d[m];
            let par_m: Par = [par[6], par[7], par[8], par[9], par[10], par[5], 0.0, 0.0, 0.0, 0.0, 0.0];

            let h_f = self.h_new[i];
            let h_m = self.h_matrix_new[i];

            // Effective interface conductivity K_a
            let k_m = fk(SoilModel::VanGenuchten, h_m, &par_m);
            let k_f_eval = fk(SoilModel::VanGenuchten, h_f, &par_m);
            let k_a = 0.5 * (k_m + k_f_eval);

            // Gamma_w = alpha_dw * K_a * (h_f - h_m)
            self.sink_im[i] = self.alpha_dw * k_a * (h_f - h_m);

            if i >= 1 {
                self.w_transf += (self.sink_im[i - 1] + self.sink_im[i]) / 2.0 * (self.x[i] - self.x[i - 1]);
            }
        }
    }
    fn it_cum_start(&mut self) {}
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Project;

    #[test]
    fn test_coupled_liquid_and_vapor_velocity() {
        let mut prj = Project::default();
        prj.processes.vapor = true;
        prj.processes.heat = true;
        
        let mut sim = Simulation::new(prj).expect("Simulation should initialize");
        
        // Setup simple 3-node column: bottom = -10, mid = -5, top = 0
        sim.n = 3;
        sim.x = vec![-10.0, -5.0, 0.0];
        sim.dt = 1.0;
        sim.cos_alf = 1.0;
        sim.sink = vec![0.0; 3];
        sim.l_vapor = true;

        let h = vec![-100.0, -90.0, -80.0];      // dh/dx = (-80 - -100)/10 = 2.0
        let th_new = vec![0.25; 3];
        let th_old = vec![0.25; 3];
        let temp = vec![15.0, 20.0, 25.0];        // dT/dx = (25 - 15)/10 = 1.0

        sim.con = vec![1.0; 3];                   // Liquid K = 1.0
        sim.con_vh = vec![0.1; 3];                // Isothermal vapor K_vh = 0.1
        sim.con_vt = vec![0.05; 3];               // Thermal vapor K_vT = 0.05
        sim.th_v_new = vec![0.001; 3];
        sim.th_v_old = vec![0.001; 3];

        let (v_liq, v_vap) = sim.veloc(&h, &th_new, &th_old, &temp);

        // 1. Liquid velocity at interior node: -K * (dh/dx + cos_alf)
        // dx = 5.0, dh/dx = 10 / 5 = 2.0. With gravity = 1.0 -> -1.0 * (2.0 + 1.0) = -3.0
        assert_eq!(v_liq.len(), 3);
        assert!((v_liq[1] - (-3.0)).abs() < 1e-6);

        // 2. Vapor velocity at interior node: -K_vh * (dh/dx) - K_vT * (dT/dx)
        // -0.1 * 2.0 - 0.05 * 1.0 = -0.20 - 0.05 = -0.25
        assert_eq!(v_vap.len(), 3);
        assert!((v_vap[1] - (-0.25)).abs() < 1e-6);

        // 3. Boundary nodes reflect matching gradient directions
        assert!((v_vap[0] - (-0.25)).abs() < 1e-6);
        assert!((v_vap[2] - (-0.25)).abs() < 1e-6);
    }
}