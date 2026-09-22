//! Collection of output records (TLInf, NodOut, ObsNod, SubReg).

use crate::output::*;
use crate::sim::Simulation;

impl Simulation {
    /// Water-flow information for the current time level; also updates cumulative fluxes.
    pub fn tl_inf(&mut self, j_print: bool) {
        let n = self.n;
        let m = n - 1;
        let grav = self.cos_alf;
        let dxn = self.x[m] - self.x[m - 1];
        let dt = self.dt;
        let v_t = -(self.con[m] + self.con[m - 1]) / 2.0 * ((self.h_new[m] - self.h_new[m - 1]) / dxn + grav)
            - (self.th_new[m] - self.th_old[m]) * dxn / 2.0 / dt
            - self.sink[m] * dxn / 2.0;
        let dx1 = self.x[1] - self.x[0];
        let v_b = -(self.con[0] + self.con[1]) / 2.0 * ((self.h_new[1] - self.h_new[0]) / dx1 + grav)
            + (self.th_new[0] - self.th_old[0]) * dx1 / 2.0 / dt
            + self.sink[0] * dx1 / 2.0;
        self.v_top = v_t;
        self.v_bot = v_b;
        let mut run_off = 0.0;
        if (!self.w_layer || (self.w_layer && self.h_new[m] >= self.h_crit_s)) && self.r_top < 0.0 {
            run_off = (self.r_top - v_t).abs();
        }
        if run_off < 1e-5 {
            run_off = 0.0;
        }
        let (mut r_infil, mut r_evap) = (0.0, 0.0);
        let (prec, r_soil) = (self.prec, self.r_soil);
        if v_t < 0.0 && (prec > 0.0 || (self.w_layer && self.h_new[m] > 0.0)) {
            r_infil = -v_t + r_soil;
        }
        if v_t >= 0.0 && prec > 0.0 {
            r_infil = prec;
        }
        if v_t > 0.0 {
            r_evap = v_t + prec;
        }
        if v_t <= 0.0 && r_soil > 0.0 && prec > 0.0 {
            r_evap = r_soil;
        }
        if v_t < 0.0 && self.w_layer && self.h_new[m] > 0.0 {
            r_evap = r_soil;
        }
        if v_t < 0.0 && self.kod_top > 0 {
            r_infil = -v_t;
        }
        let q = &mut self.cum_q;
        q[0] += self.r_top * dt;
        q[1] += self.r_root * dt;
        q[2] += v_t * dt;
        q[3] += self.v_root * dt;
        q[4] += v_b * dt;
        q[5] += run_off * dt;
        q[6] += r_infil * dt;
        q[7] += r_evap * dt;
        self.w_cum_t += (v_b - v_t - self.v_root) * dt;
        self.w_cum_a += (v_b.abs() + v_t.abs() + self.v_root.abs()) * dt;
        let sol_extra = if self.l_chem { self.solute_cumulate(run_off) } else { vec![] };
        if !(j_print && (self.l_wat || self.t_level == 1 || self.l_end)) {
            return;
        }
        let mut volume = 0.0;
        for i in 0..n - 1 {
            volume += (self.x[i + 1] - self.x[i]) * (self.th_new[i] + self.th_new[i + 1]) / 2.0;
        }
        let rec = TLevel {
            tlevel: self.t_level,
            t: self.t,
            dt,
            iter_w: self.iter_w,
            iter_c: self.iter_c,
            it_cum: self.it_cum,
            kod_top: self.kod_top,
            kod_bot: self.kod_bot,
            converged: self.convg,
            r_top: self.r_top,
            r_root: self.r_root,
            v_top: v_t,
            v_root: self.v_root,
            v_bot: v_b,
            cum_r_top: self.cum_q[0],
            cum_r_root: self.cum_q[1],
            cum_v_top: self.cum_q[2],
            cum_v_root: self.cum_q[3],
            cum_v_bot: self.cum_q[4],
            h_top: self.h_new[m],
            h_root: self.h_root,
            h_bot: self.h_new[0],
            run_off,
            cum_run_off: self.cum_q[5],
            volume,
            cum_infil: self.cum_q[6],
            cum_evap: self.cum_q[7],
            precip: prec,
            peclet: self.peclet_courant().0,
            courant: self.peclet_courant().1,
            temp_top: self.temp(m),
            temp_bot: self.temp(0),
            solutes: sol_extra,
        };
        self.res.tlevel.push(rec);
    }

    /// Nodal profile output (NodOut). Fluxes use the NodOut formula.
    pub fn profile_out(&mut self, t: f64) {
        let n = self.n;
        let grav = self.cos_alf;
        let ns = self.n_species();
        let con_sn = self.con_sat[self.mat[n - 1]] * self.ak[n - 1];
        let mut nodes = Vec::with_capacity(n);
        for i in (0..n).rev() {
            let vi;
            if i == 0 {
                let dx = self.x[1] - self.x[0];
                vi = -(self.con[0] + self.con[1]) / 2.0 * ((self.h_new[1] - self.h_new[0]) / dx + grav);
            } else if i == n - 1 {
                let dx = self.x[n - 1] - self.x[n - 2];
                vi = -(self.con[n - 1] + self.con[n - 2]) / 2.0 * ((self.h_new[n - 1] - self.h_new[n - 2]) / dx + grav)
                    - (self.th_new[n - 1] - self.th_old[n - 1]) * dx / 2.0 / self.dt
                    - self.sink[n - 1] * dx / 2.0;
            } else {
                let dxa = self.x[i + 1] - self.x[i];
                let dxb = self.x[i] - self.x[i - 1];
                let va = -(self.con[i] + self.con[i + 1]) / 2.0 * ((self.h_new[i + 1] - self.h_new[i]) / dxa + grav);
                let vb = -(self.con[i] + self.con[i - 1]) / 2.0 * ((self.h_new[i] - self.h_new[i - 1]) / dxb + grav);
                vi = (va * dxa + vb * dxb) / (dxa + dxb);
            }
            nodes.push(NodeOut {
                node: n - i,
                depth: self.x[n - 1] - self.x[i],
                h: self.h_new[i],
                theta: self.th_new[i],
                k: self.con[i],
                c: self.cap[i],
                flux: vi,
                sink: self.sink[i],
                kappa: self.kappa[i],
                v_over_ks: vi / con_sn,
                temp: self.temp(i),
                conc: (0..ns).map(|j| self.conc(j, i)).collect(),
                sorb: (0..ns).map(|j| self.sorb(j, i)).collect(),
            });
        }
        self.res.profiles.push(ProfileOut { t, nodes });
    }

    pub fn obs_nod(&mut self) {
        if self.prj.profile.observation_nodes.is_empty() {
            return;
        }
        let n = self.n;
        let ns = self.n_species();
        let mut pts = vec![];
        for &nd in &self.prj.profile.observation_nodes {
            if nd < 1 || nd > n {
                continue;
            }
            let i = n - nd;
            pts.push(ObsPoint {
                node: nd,
                h: self.h_new[i],
                theta: self.th_new[i],
                temp: self.temp(i),
                flux: self.v_new.get(i).copied().unwrap_or(0.0),
                conc: (0..ns).map(|j| self.conc(j, i)).collect(),
            });
        }
        self.res.obs.push(ObsOut { t: self.t, points: pts });
    }

    /// Mass balance and sub-region summary (SubReg).
    pub fn sub_reg(&mut self, p_level: usize) {
        let n = self.n;
        let dt = self.dt;
        let n_lay = self.lay.iter().cloned().max().unwrap_or(1);
        let mut sub: Vec<SubRegion> = vec![SubRegion::default(); n_lay];
        let mut tot = SubRegion::default();
        let mut a_tot = 0.0;
        let mut delt_w = 0.0;
        let mut h_tot = 0.0;
        for i in (0..n - 1).rev() {
            let j = i + 1;
            let lay = self.lay[i] - 1;
            let dx = self.x[j] - self.x[i];
            sub[lay].area += dx;
            a_tot += dx;
            let he = (self.h_new[i] + self.h_new[j]) / 2.0;
            let vnew = dx * (self.th_new[i] + self.th_new[j]) / 2.0;
            let vold = dx * (self.th_old[i] + self.th_old[j]) / 2.0;
            tot.volume += vnew;
            tot.change += (vnew - vold) / dt;
            sub[lay].change += (vnew - vold) / dt;
            sub[lay].volume += vnew;
            h_tot += he * dx;
            sub[lay].h_mean += he * dx;
            if p_level == 0 {
                self.wat_in[i] = vnew;
            } else {
                delt_w += (self.wat_in[i] - vnew).abs();
            }
        }
        for s in sub.iter_mut() {
            if s.area > 0.0 {
                s.h_mean /= s.area;
            }
        }
        if a_tot > 0.0 {
            h_tot /= a_tot;
        }
        tot.area = a_tot;
        tot.h_mean = h_tot;
        let mut rec = BalanceOut { t: self.t, total: tot.clone(), sub, ..Default::default() };
        // fluxes at the boundary
        let grav = self.cos_alf;
        let dx1 = self.x[1] - self.x[0];
        rec.bot_flux = -(self.con[0] + self.con[1]) / 2.0 * ((self.h_new[1] - self.h_new[0]) / dx1 + grav);
        let dxn = self.x[n - 1] - self.x[n - 2];
        rec.top_flux = -(self.con[n - 1] + self.con[n - 2]) / 2.0 * ((self.h_new[n - 1] - self.h_new[n - 2]) / dxn + grav);
        if p_level == 0 {
            self.w_vol_i = tot.volume;
        } else {
            rec.wat_bal_t = tot.volume - self.w_vol_i - self.w_cum_t;
            let ww = delt_w.max(self.w_cum_a);
            if ww > 1e-25 {
                rec.wat_bal_r = rec.wat_bal_t.abs() / ww * 100.0;
            }
        }
        if self.l_chem {
            self.sub_reg_solute(p_level, &mut rec);
        }
        self.res.balance.push(rec);
    }
}
