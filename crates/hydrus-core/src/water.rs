//! Richards equation solver (port of WATFLOW.FOR).

use crate::material::*;
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

    /// Hydraulic properties at every node (SetMat in Fortran).
    pub fn set_mat(&mut self, iter: usize) {
        let model = self.model;
        for i in 0..self.n {
            let m = self.mat[i];
            let (hi1, hi2);
            if self.kappa[i] == -1 {
                hi1 = self.h_sat[m].min(self.h_temp[i] / self.ah[i]);
                hi2 = self.h_sat[m].min(self.h_new[i] / self.ah[i]);
            } else {
                hi1 = self.h_sat[m].min(self.h_temp[i] / self.ah[i] / self.ah_w[m]);
                hi2 = self.h_sat[m].min(self.h_new[i] / self.ah[i] / self.ah_w[m]);
            }
            let him = 0.1 * hi1 + 0.9 * hi2;
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
                self.con[i] = self.con_r[i] + coni * self.ak[i] * bt * self.ak_s[i] * self.ak_w[m];
                self.cap[i] = capi * self.ath[i] * self.ath_s[i] * self.ath_w[m] / self.ah[i] / self.ah_w[m] / at;
                self.th_eq[i] = self.th_rr[i] + self.ath_w[m] * self.ath[i] * self.ath_s[i] * (thei - self.thr[m]);
            }
            if iter == 0 {
                self.con_o[i] = self.con[i];
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
                let sew = fs(model, hh, &self.par_w[m]);
                if sew < 0.999 {
                    self.ath_s[i] = (self.th_old[i] - ths) / (1.0 - sew) / (thr - ths_w);
                }
                self.th_rr[i] = ths - self.ath_s[i] * (ths_w - thr);
                self.ak_s[i] = 1.0;
                self.con_r[i] = 0.0;
                if self.i_hyst == 2 {
                    let kw = self.ak[i] * fk(model, hh, &self.par_w[m]);
                    if kw < 0.999 * ks_w {
                        self.ak_s[i] = (self.con_o[i] - ks) / (kw - ks_w);
                    }
                    self.con_r[i] = ks - self.ak_s[i] * ks_w;
                }
            } else {
                self.ath_s[i] = (self.th_old[i] - thr) / fs(model, hh, &self.par_d[m]) / (ths_d - thr);
                self.th_rr[i] = thr;
                self.ak_s[i] = 1.0;
                self.con_r[i] = 0.0;
                if self.i_hyst == 2 {
                    self.ak_s[i] = self.con_o[i] / fk(model, hh, &self.par_d[m]) / self.ak[i];
                }
            }
        }
    }

    /// Velocity at node i using the Fortran "Veloc" formula.
    pub fn veloc(&self, h: &[f64], th_new: &[f64], th_old: &[f64]) -> Vec<f64> {
        let n = self.n;
        let x = &self.x;
        let con = &self.con;
        let g = self.cos_alf;
        let mut v = vec![0.0; n];
        let m = n - 1;
        let dxn = x[m] - x[m - 1];
        v[m] = -(con[m] + con[m - 1]) / 2.0 * ((h[m] - h[m - 1]) / dxn + g)
            - dxn / 2.0 * ((th_new[m] - th_old[m]) / self.dt + self.sink[m]);
        for i in 1..n - 1 {
            let dxa = x[i + 1] - x[i];
            let dxb = x[i] - x[i - 1];
            let va = -(con[i] + con[i + 1]) / 2.0 * ((h[i + 1] - h[i]) / dxa + g);
            let vb = -(con[i] + con[i - 1]) / 2.0 * ((h[i] - h[i - 1]) / dxb + g);
            v[i] = (va * dxb + vb * dxa) / (dxa + dxb);
        }
        let dx1 = x[1] - x[0];
        v[0] = -(con[0] + con[1]) / 2.0 * ((h[1] - h[0]) / dx1 + g) + dx1 / 2.0 * ((th_new[0] - th_old[0]) / self.dt + self.sink[0]);
        v
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
        // bottom
        let dxb = self.x[1] - self.x[0];
        let dx = dxb / 2.0;
        let conb = (self.con[0] + self.con[1]) / 2.0;
        let b = conb * grav;
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
        let pb = b - self.sink[0] * dx + f2 * self.h_new[0] - (self.th_new[0] - self.th_old[0]) * dx / dt + self.r_bot;
        for i in 1..n - 1 {
            let dxa = self.x[i] - self.x[i - 1];
            let dxb = self.x[i + 1] - self.x[i];
            let dx = (dxa + dxb) / 2.0;
            let cona = (self.con[i] + self.con[i - 1]) / 2.0;
            let conb = (self.con[i] + self.con[i + 1]) / 2.0;
            let b = (cona - conb) * grav;
            let a2 = cona / dxa + conb / dxb;
            let a3 = -conb / dxb;
            let f2 = self.cap[i] * dx / dt;
            r[i] = a2 + f2;
            p[i] = f2 * self.h_new[i] - (self.th_new[i] - self.th_old[i]) * dx / dt - b - self.sink[i] * dx;
            s[i] = a3;
        }
        // top
        let m = n - 1;
        let dxa = self.x[m] - self.x[m - 1];
        let dx = dxa / 2.0;
        let cona = (self.con[m] + self.con[m - 1]) / 2.0;
        let b = cona * grav;
        let f2 = self.cap[m] * dx / dt;
        let mut rt = cona / dxa + f2;
        let st = -cona / dxa;
        let mut pt = f2 * self.h_new[m] - (self.th_new[m] - self.th_old[m]) * dx / dt - self.sink[m] * dx - b;
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
        // lower[i] = s[i-1]; diag = r; upper[i] = s[i]
        let mut diag = sys.r.clone();
        let mut rhs = sys.p.clone();
        let mut upper = sys.s.clone();
        let mut lower = vec![0.0; n];
        for i in 1..n {
            lower[i] = sys.s[i - 1];
        }
        // bottom row
        if self.kod_bot >= 0 {
            diag[0] = 1.0;
            upper[0] = 0.0;
            rhs[0] = self.h_bot;
        } else {
            diag[0] = sys.rb;
            upper[0] = sys.sb;
            rhs[0] = sys.pb;
        }
        // top row
        if self.kod_top > 0 {
            diag[n - 1] = 1.0;
            lower[n - 1] = 0.0;
            rhs[n - 1] = self.h_top;
        } else {
            diag[n - 1] = sys.rt;
            lower[n - 1] = sys.st;
            rhs[n - 1] = sys.pt;
        }
        // Thomas algorithm
        for i in 1..n {
            let mut d = diag[i - 1];
            if d.abs() < rmin {
                d = rmin;
            }
            let f = lower[i] / d;
            diag[i] -= f * upper[i - 1];
            rhs[i] -= f * rhs[i - 1];
        }
        let mut dn = diag[n - 1];
        if dn.abs() < rmin {
            dn = rmin;
        }
        self.h_new[n - 1] = rhs[n - 1] / dn;
        for i in (0..n - 1).rev() {
            let mut d = diag[i];
            if d.abs() < rmin {
                d = rmin;
            }
            self.h_new[i] = (rhs[i] - upper[i] * self.h_new[i + 1]) / d;
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
                let v_top = -(self.con[n - 1] + self.con[m]) / 2.0 * ((self.h_new[n - 1] - self.h_new[m]) / dx + grav)
                    - (self.th_new[n - 1] - self.th_old[n - 1]) * dx / 2.0 / self.dt
                    - self.sink[n - 1] * dx / 2.0;
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
            if self.w_layer && self.h_new[n - 1] > 0.0 && self.h_new[n - 1] < 0.00005 * self.x_conv && self.r_top >= 0.0 {
                let m = self.mat[n - 1];
                let hh = fh(self.model, 0.9999, &self.par_d[m]);
                self.h_new[n - 1] = hh;
                self.h_old[n - 1] = hh;
                self.h_temp[n - 1] = hh;
            }
            loop {
                self.set_mat(self.iter_w);
                if self.iter_w == 2 && self.i_hyst > 0 {
                    self.hyster_kp();
                }
                let sys = self.build_system();
                self.shift();
                self.h_temp.copy_from_slice(&self.h_new);
                self.solve(&sys);
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
            if self.w_layer && self.h_new[n - 1] > self.h_crit_s {
                self.kod_top = 4;
                self.h_top = self.h_crit_s;
            }
            return;
        }
    }

    fn it_cum_start(&mut self) {}
}
