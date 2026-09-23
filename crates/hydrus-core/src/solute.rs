//! Solute transport (port of SOLUTE.FOR): advection-dispersion with equilibrium and
//! non-equilibrium sorption (two-site / mobile-immobile), first-order decay chains,
//! zero-order production, gas-phase partitioning and passive root uptake.
//! Not ported: temperature / water-content dependent rates, virus/colloid filtration,
//! dual-porosity coupling, volatilisation boundary.

use crate::error::HydrusError;
use crate::heat::thomas;
use crate::model::*;
use crate::output::*;
use crate::sim::Simulation;

pub struct SoluteState {
    pub ns: usize,
    pub epsi: f64,
    pub upw: bool,
    pub art_d: bool,
    pub tort: bool,
    pub tort_model: i32,
    pub tol_a: f64,
    pub tol_r: f64,
    pub max_it_c: usize,
    pub pe_cr: f64,
    pub k_top: i32,
    pub k_bot: i32,
    pub t_pulse: f64,
	pub d_surf: f64,
    pub c_atm: f64,
	pub l_tdep: bool,
    pub t_dep: Vec<SpeciesTDep>,
    // parameters
    pub ro: Vec<f64>,
    pub disp_l: Vec<f64>,
    pub frac: Vec<f64>,
    pub th_im: Vec<f64>,
    pub diff_w: Vec<f64>,
    pub diff_g: Vec<f64>,
    pub pm: Vec<Vec<SpeciesMaterial>>,
    pub c_root_max: Vec<f64>,
    pub l_linear: Vec<bool>,
    pub l_equil: bool,
    pub l_mob_im: Vec<bool>,
    // state
    pub conc: Vec<Vec<f64>>,
    pub sorb: Vec<Vec<f64>>,
    pub c_top: Vec<f64>,
    pub c_bot: Vec<f64>,
    pub c_t: Vec<f64>,
    pub c_prev_o: Vec<f64>,
    pub c_new: Vec<f64>,
    pub c_temp: Vec<f64>,
    pub sorb_n: Vec<f64>,
	pub l_moist: bool,
	pub i_moist_dep: i32,
    pub w_dep: Vec<SpeciesWDep>,
	pub moist_tables: Vec<MoistDepTable>,
    pub th_ref: Vec<Vec<[f64; 9]>>,
	pub l_dual_neq: bool,
	pub l_nequil: bool,
    // work arrays
    pub th_n: Vec<f64>,
    pub th_o: Vec<f64>,
    pub v_o: Vec<f64>,
    pub v_n: Vec<f64>,
    pub disp: Vec<f64>,
    pub retard: Vec<f64>,
    pub g0: Vec<f64>,
    pub g1: Vec<f64>,
    pub q0: Vec<f64>,
    pub q1: Vec<f64>,
    pub wc: Vec<f64>,
    pub v_corr: Vec<f64>,
    pub s_sink: Vec<f64>,
    pub mb: Vec<f64>,
    pub md: Vec<f64>,
    pub me: Vec<f64>,
    pub mf: Vec<f64>,
    pub d_mob_i: f64,
    pub f1: f64,
    pub e1: f64,
    pub d1: f64,
    pub bn: f64,
    pub dn: f64,
    pub fn_: f64,
    // fluxes
    pub cv_top: Vec<f64>,
    pub cv_bot: Vec<f64>,
    pub cv_ch0: Vec<f64>,
    pub cv_ch1: Vec<f64>,
    pub cv_chr: Vec<f64>,
    pub cv_chim: Vec<f64>,
    pub cum_ch: Vec<[f64; 10]>,
    pub c_cum_t: Vec<f64>,
    pub c_cum_a: Vec<f64>,
    pub c_root: Vec<f64>,
    pub c_vol_i: Vec<f64>,
    pub sol_in: Vec<Vec<f64>>,
    pub peclet: f64,
    pub courant: f64,
	pub l_bact: bool,
	pub l_filtr: bool,
	pub sorb2: Vec<Vec<f64>>,
	pub sorb_n2: Vec<f64>,
	pub conc_m: Vec<Vec<f64>>,
    pub sorb_m: Vec<Vec<f64>>,	
}

impl Simulation {
    pub fn conc(&self, j: usize, i: usize) -> f64 {
        self.sol.as_ref().map(|s| s.conc.get(j).map(|c| c[i]).unwrap_or(0.0)).unwrap_or(0.0)
    }
    pub fn sorb(&self, j: usize, i: usize) -> f64 {
        self.sol.as_ref().map(|s| s.sorb.get(j).map(|c| c[i]).unwrap_or(0.0)).unwrap_or(0.0)
    }
    pub fn conc_matrix(&self, j: usize, i: usize) -> f64 {
        self.sol.as_ref().map(|s| s.conc_m.get(j).map(|c| c[i]).unwrap_or(0.0)).unwrap_or(0.0)
    }
    pub fn sorb_matrix(&self, j: usize, i: usize) -> f64 {
        self.sol.as_ref().map(|s| s.sorb_m.get(j).map(|c| c[i]).unwrap_or(0.0)).unwrap_or(0.0)
    }
    pub fn peclet_courant(&self) -> (f64, f64) {
        self.sol.as_ref().map(|s| (s.peclet, s.courant)).unwrap_or((0.0, 0.0))
    }
    pub fn set_c_root(&mut self, sums: &[f64], a: f64) {
        if let Some(s) = self.sol.as_mut() {
            for j in 0..s.ns.min(sums.len()) {
                s.c_root[j] = sums[j] / a;
            }
        }
    }

    pub fn atm_record_scalars(&mut self, rec: &AtmRecord) {
        if let Some(h) = self.heat.as_mut() {
            h.t_top = rec.t_top;
            h.t_bot = rec.t_bot;
            h.ampl = rec.ampl;
        }
        if let Some(s) = self.sol.as_mut() {
            for j in 0..s.ns {
                s.c_t[j] = rec.c_top.get(j).copied().unwrap_or(0.0);
                s.c_bot[j] = rec.c_bot.get(j).copied().unwrap_or(0.0);
            }
        }
    }

    pub fn set_chem_bc(&mut self) {
        let (prec, r_soil, kod_top) = (self.prec, self.r_soil, self.kod_top);
        let wl = self.w_layer && self.h_new[self.n - 1] > 0.0;
        if let Some(s) = self.sol.as_mut() {
            for j in 0..s.ns {
                if wl {
                    continue;
                }
                s.c_top[j] = s.c_t[j];
                if kod_top.abs() == 4 && s.k_top <= 0 {
                    if prec - r_soil > 0.0 {
                        s.c_top[j] = s.c_t[j] * prec / (prec - r_soil);
                    } else if r_soil > 0.0 {
                        s.c_top[j] = 0.0;
                    }
                }
            }
        }
    }

    /// Mass balance in the surface layer for concentration at the top.
    pub fn sol_wlayer_mix(&mut self) {
        let m = self.n - 1;
        if !(self.w_layer && self.h_new[m] > 0.0) {
            return;
        }
        let ht = self.h_new[m];
        let dt = self.dt;
        let prec = self.prec.max(0.0);
        let r_soil = self.r_soil.max(0.0);

        if let Some(s) = self.sol.as_mut() {
            for j in 0..s.ns {
                let c_rain = s.c_t[j];
                // Previous surface layer mass + precipitation mass
                let num = ht * s.c_top[j] + dt * prec * c_rain;
                let den = ht + dt * (prec - r_soil);
                if den > 1e-12 {
                    s.c_top[j] = (num / den).max(0.0);
                }
            }
        }
    }

    pub fn solute_init(&mut self) -> Result<(), HydrusError> {
        let n = self.n;
        let so = self.prj.solute.clone();
        let ns = so.species.len();
        let n_mat = self.n_mat;
        let mut l_equil = true;
        let mut l_mob = vec![false; n_mat];
        for m in 0..n_mat {
            let mm = &so.materials[m];
            if mm.frac < 1.0 || mm.th_immobile > 0.0 || so.l_bact || so.l_dual_neq || self.l_dual_perm {
                l_equil = false;
            }
            l_mob[m] = mm.th_immobile > 0.0;
            if !l_equil && mm.bulk_density == 0.0 {
                return Err(HydrusError::Invalid("Bulk density cannot be zero.".into()));
            }
        }
        let mut l_linear = vec![true; ns];
        for j in 0..ns {
            for m in 0..n_mat {
                let p = &so.species[j].per_material[m];
                if p.nu.abs() > 1e-12 || (p.beta - 1.0).abs() > 0.001 {
                    l_linear[j] = false;
                }
            }
        }
        let nodes: Vec<&ProfileNode> = self.prj.profile.nodes.iter().rev().collect();
        let mut conc = vec![vec![0.0; n]; ns];
        let mut sorb = vec![vec![0.0; n]; ns];
		let mut sorb2 = vec![vec![0.0; n]; ns];
        if so.l_bact || so.l_dual_neq {
            for j in 0..ns {
                for i in 0..n {
                    sorb2[j][i] = nodes[i].sorb.get(j).copied().unwrap_or(0.0);
                }
            }
        }
        for i in 0..n {
            for j in 0..ns {
                conc[j][i] = nodes[i].conc.get(j).copied().unwrap_or(0.0);
                sorb[j][i] = nodes[i].sorb.get(j).copied().unwrap_or(0.0);
            }
        }
        let ro: Vec<f64> = so.materials.iter().map(|m| m.bulk_density).collect();
        let frac: Vec<f64> = so.materials.iter().map(|m| m.frac).collect();
        // initial mass -> concentration
        if so.mass_init {
            for j in 0..ns {
                for i in 0..n {
                    if conc[j][i] > 0.0 {
                        let m = self.mat[i];
                        let p = &so.species[j].per_material[m];
                        let th = self.th_old[i];
                        let tha = self.ths[m] - th;
                        let mass = conc[j][i];

                        let ks_eff = if so.l_bact {
                            let mut k1 = 0.0;
                            let mut k2 = 0.0;
                            if p.r_kd1 > 0.0 && ro[m] > 0.0 {
                                k1 = th * p.r_ka1 / ro[m] / p.r_kd1;
                            }
                            if p.r_kd2 > 0.0 && ro[m] > 0.0 {
                                k2 = th * p.r_ka2 / ro[m] / p.r_kd2;
                            }
                            k1 + k2
                        } else {
                            p.ks
                        };

                        conc[j][i] = if l_linear[j] {
                            mass / (th + ro[m] * ks_eff + tha * p.henry)
                        } else {
                            let f = |c: f64| th * c + ro[m] * ks_eff * c.powf(p.beta) / (1.0 + p.nu * c.powf(p.beta)) + tha * p.henry * c;
                            let (mut lo, mut hi) = (1e-20f64, 1e20f64);
                            for _ in 0..300 {
                                let mid = (lo * hi).sqrt();
                                if f(mid) > mass {
                                    hi = mid
                                } else {
                                    lo = mid
                                }
                            }
                            (lo * hi).sqrt()
                        };
                    }
                }
            }
        }
        if !l_equil && (so.equil_init || so.mass_init) {
            for j in 0..ns {
                for i in 0..n {
                    let m = self.mat[i];
                    let p = &so.species[j].per_material[m];
                    let th = self.th_old[i];
                    let cc = conc[j][i];

                    if so.l_bact {
                        let ks1 = if p.r_kd1 > 0.0 && ro[m] > 0.0 { th * p.r_ka1 / ro[m] / p.r_kd1 } else { 0.0 };
                        let ks2 = if p.r_kd2 > 0.0 && ro[m] > 0.0 { th * p.r_ka2 / ro[m] / p.r_kd2 } else { 0.0 };
                        sorb[j][i] = ks1 * cc;
                        sorb2[j][i] = ks2 * cc;
                    } else if so.l_dual_neq {
                        sorb[j][i] = cc; // Physical immobile liquid in equilibrium with mobile liquid
                        let sc = if !l_linear[j] && cc > 0.0 { cc.powf(p.beta - 1.0) / (1.0 + p.nu * cc.powf(p.beta)) } else { 1.0 };
                        sorb2[j][i] = (1.0 - frac[m]) * sc * p.ks * cc; // Kinetic sorption site in mobile zone
                    } else if l_mob[m] {
                        sorb[j][i] = cc;
                    } else {
                        let sc = if !l_linear[j] && cc > 0.0 { cc.powf(p.beta - 1.0) / (1.0 + p.nu * cc.powf(p.beta)) } else { 1.0 };
                        sorb[j][i] = (1.0 - frac[m]) * sc * p.ks * cc;
                    }
                }
            }
        }
        let mut c_top = vec![0.0; ns];
        let mut c_bot = vec![0.0; ns];
        for j in 0..ns {
            c_top[j] = so.species[j].c_top;
            c_bot[j] = so.species[j].c_bot;
        }
        let z = |k: usize| vec![0.0; k];
        let conc_m = conc.clone();
        let sorb_m = sorb.clone();
		let mut th_ref = vec![vec![[0.0; 9]; n_mat]; ns];
        if so.l_moist {
            for j in 0..ns {
                if j < so.w_dep.len() {
                    let wdep = &so.w_dep[j];
                    for m in 0..n_mat {
                        for k in 0..9 {
                            let href = wdep.h_ref[k];
                            // Compute theta_ref using retention curve FQ for material m
                            th_ref[j][m][k] = self.theta_of(m, href);
                        }
                    }
                }
            }
        }
        self.sol = Some(SoluteState {
            ns,
            epsi: so.epsi,
            upw: so.upstream_weighting,
            art_d: so.artificial_dispersion,
            tort: so.tortuosity,
            tort_model: so.tort_model,
            tol_a: so.tol_abs,
            tol_r: so.tol_rel,
            max_it_c: so.max_iter.max(1),
            pe_cr: so.pe_cr.max(0.1),
            k_top: so.k_top,
            k_bot: so.k_bot,
            t_pulse: so.t_pulse,
            d_surf: so.d_surf,
            c_atm: so.c_atm,
			l_tdep: so.l_tdep,
            t_dep: so.t_dep.clone(),
            ro,
            disp_l: so.materials.iter().map(|m| m.disp_l).collect(),
            frac,
            th_im: so.materials.iter().map(|m| m.th_immobile).collect(),
            diff_w: so.species.iter().map(|s| s.diff_w).collect(),
            diff_g: so.species.iter().map(|s| s.diff_g).collect(),
            pm: so.species.iter().map(|s| s.per_material.clone()).collect(),
            c_root_max: (0..ns).map(|j| self.prj.root.c_root_max.get(j).copied().unwrap_or(0.0)).collect(),
            l_linear,
            l_equil,
            l_mob_im: l_mob,
            conc,
            sorb,
            c_top,
            c_bot,
            c_t: z(ns),
            c_prev_o: z(n),
            c_new: z(n),
            c_temp: z(n),
			l_moist: so.l_moist,
			i_moist_dep: so.i_moist_dep,
			w_dep: so.w_dep.clone(),
			moist_tables: so.moist_tables.clone(),
			th_ref,
			l_dual_neq: so.l_dual_neq,
			l_nequil: so.l_nequil,
            sorb_n: z(n),
            th_n: z(n),
            th_o: z(n),
            v_o: z(n),
            v_n: z(n),
            disp: z(n),
            retard: z(n),
            g0: z(n),
            g1: z(n),
            q0: z(n),
            q1: z(n),
            wc: z(n),
            v_corr: z(n),
            s_sink: z(n),
            mb: z(n),
            md: z(n),
            me: z(n),
            mf: z(n),
            d_mob_i: 1.0,
            f1: 0.0,
            e1: 0.0,
            d1: 0.0,
            bn: 0.0,
            dn: 0.0,
            fn_: 0.0,
            cv_top: z(ns),
            cv_bot: z(ns),
            cv_ch0: z(ns),
            cv_ch1: z(ns),
            cv_chr: z(ns),
            cv_chim: z(ns),
            cum_ch: vec![[0.0; 10]; ns],
            c_cum_t: z(ns),
            c_cum_a: z(ns),
            c_root: z(ns),
            c_vol_i: z(ns),
            sol_in: vec![vec![0.0; n]; ns],
            peclet: 0.0,
            courant: 0.0,
            l_bact: so.l_bact,
            l_filtr: so.l_filtr,
            sorb2,
            sorb_n2: z(n),
			conc_m,
            sorb_m,
        });
        Ok(())
    }

    pub fn solute(&mut self) -> Result<(), HydrusError> {
        let mut st = self.sol.take().expect("solute state");
        let r = self.solute_inner(&mut st);
        self.sol = Some(st);
        r
    }

    fn solute_inner(&mut self, st: &mut SoluteState) -> Result<(), HydrusError> {
        let n = self.n;
        let ns = st.ns;
        'restart: loop {
            let dt = self.dt;
            let epsi = st.epsi;
            let alf = 1.0 - epsi;
            st.peclet = 0.0;
            st.courant = 0.0;
            self.dt_max_c = 1e30;
            self.iter_c = 1;
            st.th_n.copy_from_slice(&self.th_new);
            st.th_o.copy_from_slice(&self.th_old);
            st.v_o.copy_from_slice(&self.v_old);
            st.v_n.copy_from_slice(&self.v_new);
			if st.l_bact {
                for i in 0..n {
                    let m = self.mat[i];
                    let th_c_param = st.th_im[m];
                    if th_c_param > 0.0 && m < self.par_d.len() {
                        let p = &self.par_d[m];
                        let (th_c_n, v_c_n) = exclusion(st.th_n[i], st.v_n[i], th_c_param, p[0], p[1], p[2], p[3], p[5]);
                        let (th_c_o, v_c_o) = exclusion(st.th_o[i], st.v_o[i], th_c_param, p[0], p[1], p[2], p[3], p[5]);
                        st.th_n[i] = th_c_n;
                        st.v_n[i] = v_c_n;
                        st.th_o[i] = th_c_o;
                        st.v_o[i] = v_c_o;
                    }
                }
            }
            for js in 0..ns {
                let mut iter = 0;
				st.c_prev_o.copy_from_slice(&st.conc[js]);
                st.cv_top[js] = 0.0;
                st.cv_bot[js] = 0.0;
                st.cv_ch0[js] = 0.0;
                st.cv_ch1[js] = 0.0;
                st.cv_chr[js] = 0.0;
                st.cv_chim[js] = 0.0;
                if self.t - st.t_pulse > self.dt_min && !self.atm_bc {
                    st.c_top[js] = 0.0;
                    st.c_bot[js] = 0.0;
                }
                if st.k_bot < 0 {
                    if self.l_vapor && self.r_bot.abs() < 1e-20 {
                        st.cv_bot[js] = 0.0;
                    } else if st.v_o[0] >= 0.0 {
                        st.cv_bot[js] = alf * st.c_bot[js] * st.v_o[0];
                    } else {
                        st.cv_bot[js] = alf * st.conc[js][0] * st.v_o[0];
                    }
                } else if st.k_bot == 0 {
                    st.cv_bot[js] = alf * st.conc[js][0] * st.v_o[0];
                }
                if st.k_top == -2 {
                    let m = n - 1;
                    let mat_top = self.mat[m];
                    let henry = st.pm[js][mat_top].henry;
                    let (g_vol, g_atm) = if st.d_surf > 1e-10 {
                        (st.diff_g[js] * henry / st.d_surf, st.diff_g[js] / st.d_surf)
                    } else {
                        (0.0, 0.0)
                    };
                    st.cv_top[js] += alf * g_vol * st.conc[js][m] - g_atm * st.c_atm;
                } else if st.k_top < 0 && self.t_level != 1 && st.v_o[n - 1] < 0.0 {
                    st.cv_top[js] = alf * st.c_top[js] * st.v_o[n - 1];
                }
                if !st.l_linear[js] {
                    st.c_new.copy_from_slice(&st.conc[js]);
                    if !st.l_equil {
                        st.sorb_n.copy_from_slice(&st.sorb[js]);
                    }
                }
                if self.sink_f {
                    for i in 0..n {
                        let c_val = st.conc[js][i].max(0.0);
                        // Passive uptake: S_w * min(c, cRootMax)
                        let mut sink_val = if self.beta[i] > 0.0 {
                            self.sink[i] * c_val.min(st.c_root_max[js])
                        } else {
                            0.0
                        };

                        // Active uptake (Michaelis-Menten kinetics)
                        if self.prj.root.l_act_rsu && js < self.prj.root.omega_act.len() {
                            let v_max = self.prj.root.omega_act[js];
                            let km = self.prj.root.r_km[js];
                            let c_min = self.prj.root.c_min[js];
                            if c_val > c_min && (km + c_val - c_min) > 1e-12 {
                                let r_act = self.beta[i] * v_max * (c_val - c_min) / (km + c_val - c_min);
                                sink_val += r_act.max(0.0);
                            }
                        }
                        st.s_sink[i] = sink_val;
                    }
                } else {
                    for v in st.s_sink.iter_mut() {
                        *v = 0.0;
                    }
                }
                loop {
                    iter += 1;
                    if !st.l_linear[js] {
                        st.c_temp.copy_from_slice(&st.c_new);
                    }
                    for level in 1..=2 {
                        let (pe, co, dtm) = self.sol_coeff(st, js, level, iter, dt);
                        st.peclet = st.peclet.max(pe);
                        st.courant = st.courant.max(co);
                        self.dt_max_c = self.dt_max_c.min(dtm);
                        self.sol_matset(st, js, level, dt);
                        for i in 0..n {
                            if level == 1 {
                                st.v_o[i] += st.v_corr[i];
                            } else {
                                st.v_n[i] += st.v_corr[i];
                            }
                        }
                        if level == 1 && iter == 1 {
                            self.sol_masstran(st, js, alf);
                        }
                    }
                    // solve
                    let mut d = st.md.clone();
                    let mut f = st.mf.clone();
                    let mut b = st.mb.clone();
                    let mut e = st.me.clone();
                    b[0] = 0.0;
                    e[n - 1] = 0.0;
                    thomas(&b, &mut d, &e, &mut f);
                    let mut lconv = true;
                    for i in 0..n {
                        if (ns > 1 && iter == 1) || false {
                            st.c_prev_o[i] = st.conc[js][i];
                        }
                        if st.l_linear[js] {
                            st.conc[js][i] = f[i].max(0.0);
                            if st.conc[js][i] < 1e-30 && st.conc[js][i] > 0.0 {
                                st.conc[js][i] = 0.0;
                            }
                        } else {
                            st.c_new[i] = f[i];
                            if st.c_new[i] < 1e-30 {
                                st.c_new[i] = 0.0;
                            }
                            if (st.c_new[i] - st.c_temp[i]).abs() > st.tol_a + st.tol_r * st.conc[js][i] {
                                lconv = false;
                            }
                        }
                    }
                    if !st.l_linear[js] {
                        if !lconv {
                            if iter < st.max_it_c {
                                continue;
                            } else if dt > self.dt_min && !self.l_wat {
                                let dt_old = self.dt;
                                self.dt = (self.dt / 3.0).max(self.dt_min);
                                self.dt_opt = self.dt;
                                self.t = self.t - dt_old + self.dt;
                                continue 'restart;
                            } else {
                                return Err(HydrusError::Numerical(format!("Solute transport did not converge at t = {}", self.t)));
                            }
                        }
                        for i in 0..n {
                            st.conc[js][i] = st.c_new[i];
                            if !st.l_equil {
                                st.sorb[js][i] = st.sorb_n[i];
                            }
                        }
                    }
                    break;
                }
                if !st.l_equil && st.l_linear[js] {
                    self.sol_sorbconc(st, js, dt);
                }
                self.sol_masstran(st, js, epsi);
                // fluxes across the boundaries
                if st.k_top == -2 {
                    let m = n - 1;
                    let mat_top = self.mat[m];
                    let henry = st.pm[js][mat_top].henry;
                    let (g_vol, g_atm) = if st.d_surf > 1e-10 {
                        (st.diff_g[js] * henry / st.d_surf, st.diff_g[js] / st.d_surf)
                    } else {
                        (0.0, 0.0)
                    };
                    st.cv_top[js] += epsi * g_vol * st.conc[js][m] - g_atm * st.c_atm;
                } else if st.k_top < 0 {
                    if self.t_level != 1 {
                        if st.v_n[n - 1] < 0.0 {
                            st.cv_top[js] += epsi * st.v_n[n - 1] * st.c_top[js];
                        }
                    } else if st.v_n[n - 1] < 0.0 {
                        st.cv_top[js] += st.v_n[n - 1] * st.c_top[js];
                    }
                } else {
                    st.cv_top[js] = st.fn_ - st.bn * st.conc[js][n - 2] - st.dn * st.conc[js][n - 1];
                }
                if st.k_bot < 0 {
                    if self.l_vapor && self.r_bot.abs() < 1e-20 {
                        st.cv_bot[js] = 0.0;
                    } else if st.v_n[0] >= 0.0 {
                        st.cv_bot[js] += epsi * st.c_bot[js] * st.v_n[0];
                    } else {
                        st.cv_bot[js] += epsi * st.conc[js][0] * st.v_n[0];
                    }
                } else if st.k_bot == 0 {
                    st.cv_bot[js] += epsi * st.conc[js][0] * st.v_n[0];
                } else {
                    st.cv_bot[js] = st.d1 * st.conc[js][0] + st.e1 * st.conc[js][1] - st.f1;
                }
                self.iter_c = self.iter_c.max(iter);
                if st.cv_top[js].abs() < 1e-30 {
                    st.cv_top[js] = 0.0;
                }
                if st.cv_bot[js].abs() < 1e-30 {
                    st.cv_bot[js] = 0.0;
                }
            }
            break;
        }
        Ok(())
    }

    /// Dispersion coefficients, retardation, source/decay coefficients (Coeff + Disper + NEquil + PeCour).
    fn sol_coeff(&self, st: &mut SoluteState, js: usize, level: usize, iter: usize, dt: f64) -> (f64, f64, f64) {
        let n = self.n;
        let last = level == 2;
        let epsi = st.epsi;
        let mut peclet = 0.0f64;
        let mut courant = 0.0f64;
        let mut dt_max_c = 1e30f64;
        for i in (0..n).rev() {
            let j = i + 1;
            let m = self.mat[i];
            let parent = if js > 0 { Some(st.pm[js - 1][m].clone()) } else { None };
            let mob = st.l_mob_im[m];
            let th_imob = st.th_im[m];
            let (th_arr, v_arr) = if last { (&st.th_n, &st.v_n) } else { (&st.th_o, &st.v_o) };
            let mut th_w = th_arr[i];
            let th_g = (self.ths[m] - th_w).max(0.0);
            if mob {
                th_w = (th_w - th_imob).max(0.001);
            }
            let v = v_arr[i];
            let (mut vj, mut thj) = (0.0, 0.0);
            if i != n - 1 {
                vj = v_arr[j];
                thj = th_arr[j];
                if mob {
                    thj = (thj - th_imob).max(0.001);
                }
            }
            let ro = st.ro[m];
            let frac = st.frac[m];
			let p_base = st.pm[js][m].clone();
            let mut p = p_base.clone();
			let mut diff_w = st.diff_w[js];
			let mut diff_g = st.diff_g[js];

            if st.l_tdep && js < st.t_dep.len() {
                let td = &st.t_dep[js];
                let temp_val = if last {
                    self.heat.as_ref().map(|h| h.temp_n[i]).unwrap_or(20.0)
                } else {
                    self.heat.as_ref().map(|h| h.temp_o[i]).unwrap_or(20.0)
                };
                let tr = 293.15;
                let r_gas = 8.314;
                let tt = (temp_val + 273.15 - tr) / (r_gas * (temp_val + 273.15) * tr);

                let scale = |base: f64, ea: f64| if ea.abs() > 1e-12 { base * (ea * tt).exp() } else { base };

                diff_w = scale(diff_w, td.diff_w);
				diff_g = scale(diff_g, td.diff_g);
				p.ks = scale(p_base.ks, td.ks);
                p.nu = scale(p_base.nu, td.nu);
                p.henry = scale(p_base.henry, td.henry);
                p.mu_w = scale(p_base.mu_w, td.mu_w);
                p.mu_s = scale(p_base.mu_s, td.mu_s);
                p.mu_g = scale(p_base.mu_g, td.mu_g);
                p.gam_w = scale(p_base.gam_w, td.gam_w);
                p.gam_s = scale(p_base.gam_s, td.gam_s);
                p.gam_g = scale(p_base.gam_g, td.gam_g);
                p.mu0_w = scale(p_base.mu0_w, td.mu0_w);
                p.mu0_s = scale(p_base.mu0_s, td.mu0_s);
                p.mu0_g = scale(p_base.mu0_g, td.mu0_g);
                p.omega = scale(p_base.omega, td.omega);
                if st.l_bact {
                    p.s_max2 = scale(p_base.s_max2, td.mu0_s);
                    p.r_ka2 = scale(p_base.r_ka2, td.mu0_g);
                    p.r_kd2 = scale(p_base.r_kd2, td.omega);
                    p.s_max1 = scale(p_base.s_max1, td.ks);
                    p.r_ka1 = scale(p_base.r_ka1, td.nu);
                    p.r_kd1 = scale(p_base.r_kd1, td.henry);
                }
            }
            if st.l_moist {
                let scale_m = |rate: f64, k: usize| -> f64 {
                    if rate.abs() < 1e-30 {
                        return rate;
                    }
                    if st.i_moist_dep == 2 {
                        // Tabular interpolation (MoistDep.in)
                        if js < st.moist_tables.len() {
                            let factor = interp_moist_factor(&st.moist_tables[js], k, th_w);
                            rate * factor
                        } else {
                            rate
                        }
                    } else {
                        // Walker power law (i_moist_dep == 1)
                        if js < st.w_dep.len() {
                            let b = st.w_dep[js].exp_b[k];
                            let th_ref = st.th_ref[js][m][k];
                            if b.abs() > 1e-12 && th_ref > 1e-6 && th_w > 1e-6 {
                                rate * (th_w / th_ref).powf(b)
                            } else {
                                rate
                            }
                        } else {
                            rate
                        }
                    }
                };

                p.mu_w = scale_m(p.mu_w, 0);
                p.mu_s = scale_m(p.mu_s, 1);
                p.mu_g = scale_m(p.mu_g, 2);
                p.gam_w = scale_m(p.gam_w, 3);
                p.gam_s = scale_m(p.gam_s, 4);
                p.gam_g = scale_m(p.gam_g, 5);
                p.mu0_w = scale_m(p.mu0_w, 6);
                p.mu0_s = scale_m(p.mu0_s, 7);
                p.mu0_g = scale_m(p.mu0_g, 8);
            }

            let (xks, xnu, fexp, henry) = (p.ks, p.nu, p.beta, p.henry);
            let (gam_l, gam_s, gam_g, gam_l1, gam_s1, gam_g1) = (p.mu_w, p.mu_s, p.mu_g, p.gam_w, p.gam_s, p.gam_g);
            let (xmu_l, xmu_s, xmu_g, omega) = (p.mu0_w, p.mu0_s, p.mu0_g, p.omega);
            let (mut dsconc, mut dconc, mut sconc, mut sconc_o) = (1.0, 1.0, 1.0, 1.0);
            let (mut sconc_s, mut dsconc_s, mut sconc_os) = (1.0, 1.0, 1.0);
            let mut sconc_p = 1.0;
            let mut sconc_ps = 1.0;
            let mut c_mid = 0.0;
            if !st.l_linear[js] {
                let c0 = st.conc[js][i];
                c_mid = (c0 + st.c_new[i]) / 2.0;
                let cc = if last { st.c_new[i] } else { c0 };
                if cc > 0.0 {
                    dsconc = fexp * cc.powf(fexp - 1.0) / (1.0 + xnu * cc.powf(fexp)).powi(2);
                    sconc = cc.powf(fexp - 1.0) / (1.0 + xnu * cc.powf(fexp));
                }
                if c_mid > 0.0 {
                    dconc = fexp * c_mid.powf(fexp - 1.0) / (1.0 + xnu * c_mid.powf(fexp)).powi(2);
                }
                if last && !st.l_equil && c0 > 0.0 {
                    sconc_o = c0.powf(fexp - 1.0) / (1.0 + xnu * c0.powf(fexp));
                }
                if mob {
                    let s0 = st.sorb[js][i];
                    let ss = if last { st.sorb_n[i] } else { s0 };
                    if ss > 0.0 {
                        dsconc_s = fexp * ss.powf(fexp - 1.0) / (1.0 + xnu * ss.powf(fexp)).powi(2);
                        sconc_s = ss.powf(fexp - 1.0) / (1.0 + xnu * ss.powf(fexp));
                    }
                    if last && !st.l_equil && s0 > 0.0 {
                        sconc_os = s0.powf(fexp - 1.0) / (1.0 + xnu * s0.powf(fexp));
                    }
                }
            }
            if js > 0 && !st.l_linear[js - 1] {
                let pp = parent.as_ref().unwrap();
                let cprev = if last { st.conc[js - 1][i] } else { st.c_prev_o[i] };
                if cprev > 0.0 {
                    sconc_p = cprev.powf(pp.beta - 1.0) / (1.0 + pp.nu * cprev.powf(pp.beta));
                }
                let sp = st.sorb[js - 1][i];
                if sp > 0.0 {
                    sconc_ps = sp.powf(pp.beta - 1.0) / (1.0 + pp.nu * sp.powf(pp.beta));
                }
            }
            let _ = sconc_ps;
            let retard = (ro * frac * xks * dconc + th_g * henry) / th_w + 1.0;
            st.retard[i] = retard;
            // ---- dispersion
            let (mut tau_w, mut tau_g) = (1.0, 1.0);
            if st.tort {
                let mut ths = self.ths[m];
                if mob {
                    ths = (ths - th_imob).max(0.001);
                }
                if st.tort_model == 0 {
                    tau_w = th_w.powf(7.0 / 3.0) / ths.powi(2);
                    tau_g = th_g.powf(7.0 / 3.0) / ths.powi(2);
                } else {
                    tau_w = 0.66 * (th_w / ths).powf(8.0 / 3.0);
                    tau_g = th_g.powf(1.5) / ths;
                }
            }
            let mut disp = st.disp_l[m] * v.abs() / th_w + diff_w * tau_w + th_g / th_w * diff_g * henry * tau_g;
            if !st.art_d && !st.upw {
                let mut fi = 0.0;
                if c_mid > 0.0 {
                    fi = 6.0 * th_w * ro * xks * c_mid.powf(fexp - 1.0) * (fexp / (1.0 + xnu * c_mid.powf(fexp)).powi(2) - 1.0 / (1.0 + xnu * c_mid.powf(fexp)));
                }
                let dpom = (dt / (6.0 * th_w * (th_w + ro * frac * xks * dsconc + th_g * henry) + fi)).max(0.0);
                if !last {
                    disp += v * v * dpom;
                } else {
                    disp = (disp - v * v * dpom).max(disp / 2.0);
                }
            }
            if st.art_d {
                let mut dd = 0.0;
                if st.pe_cr != 0.0 && v.abs() > 1e-15 {
                    dd = v * v * dt / th_w / th_w / retard / st.pe_cr;
                }
                if dd > disp {
                    disp = dd;
                }
            }
            st.disp[i] = disp;
			// ---- Colloid / bacteria kinetics preparation
            let mut r_ka1 = p.r_ka1;
            let mut r_ka2 = p.r_ka2;
            let mut psi1 = 1.0;
            let mut psi2 = 1.0;
            if st.l_bact {
                if st.l_filtr {
                    let (ka1, ka2) = deposit(
                        p.d_c,
                        p.d_p,
                        p.r_ka1,
                        p.r_ka2,
                        th_w,
                        v,
                        self.temp(i),
                        self.x_conv,
                        self.t_conv,
                    );
                    r_ka1 = ka1;
                    r_ka2 = ka2;
                }
                let ss1 = if last { st.sorb_n[i] } else { st.sorb[js][i] };
                let ss2 = if last { st.sorb_n2[i] } else { st.sorb2[js][i] };
                let z_dist = self.depth_of(i);
                if p.i_psi1 > 0 {
                    psi1 = blocking(p.i_psi1, p.s_max1, z_dist, ss1, p.d_c, p.s_max2);
                }
                if p.i_psi2 > 0 {
                    psi2 = blocking(p.i_psi2, p.s_max2, z_dist, ss2, p.d_c, p.s_max2);
                }
            }

            // ---- non-equilibrium sorbed phase
            let mut ssorb = st.sorb[js][i];
            let mut ssorb2 = st.sorb2[js][i];
            if !st.l_equil {
                let gamma_w = if self.i_dual_por > 0 { self.sink_im[i] } else { 0.0 };
                let tr_adv_m = if gamma_w > 0.0 { gamma_w } else { 0.0 };
                let tr_adv_im = if gamma_w < 0.0 { -gamma_w } else { 0.0 };
                let omega_eff = omega + tr_adv_m;

                let a_mob = th_imob + (1.0 - frac) * ro * xks * dsconc_s;
                let b_mob = th_imob * (gam_l + gam_l1) + (1.0 - frac) * ro * (gam_s + gam_s1) * xks * sconc_s + tr_adv_im;
                let d_mob = 2.0 * a_mob + dt * (omega_eff + b_mob);
                st.d_mob_i = d_mob;

                let c_parent_im = if st.l_nequil && js > 0 {
                    let pp = &st.pm[js - 1][m];
                    let r_w = if pp.gam_w > 0.0 { pp.gam_w } else { pp.mu_w };
                    let r_s = if pp.gam_s > 0.0 { pp.gam_s } else { pp.mu_s };
                    st.sorb[js - 1][i] * (th_imob * r_w + (1.0 - frac) * ro * r_s * pp.ks * sconc_ps)
                } else {
                    0.0
                };

                let c_parent_kin = if st.l_nequil && js > 0 {
                    let pp = &st.pm[js - 1][m];
                    let r_s = if pp.gam_s > 0.0 { pp.gam_s } else { pp.mu_s };
                    r_s * ro * (if st.l_dual_neq { st.sorb2[js - 1][i] } else { st.sorb[js - 1][i] })
                } else {
                    0.0
                };

                if last && iter == 1 {
                    if st.l_bact {
                        if st.l_linear[js] {
                            st.sorb[js][i] = ((2.0 - dt * (p.r_kd1 + gam_s + gam_s1)) * st.sorb[js][i]
                                + dt * r_ka1 * th_w * st.conc[js][i] / ro)
                                / (2.0 + dt * (p.r_kd1 + gam_s + gam_s1));
                            ssorb = st.sorb[js][i];
                            st.sorb2[js][i] = ((2.0 - dt * (p.r_kd2 + gam_s + gam_s1)) * st.sorb2[js][i]
                                + dt * r_ka2 * th_w * st.conc[js][i] / ro)
                                / (2.0 + dt * (p.r_kd2 + gam_s + gam_s1));
                            ssorb2 = st.sorb2[js][i];
                        } else {
                            let cc = st.c_new[i];
                            st.sorb_n[i] = st.sorb[js][i]
                                + dt * (epsi * (r_ka1 * th_w / ro * psi1 * cc - (p.r_kd1 + gam_s + gam_s1) * st.sorb_n[i])
                                    + (1.0 - epsi) * (r_ka1 * th_w / ro * psi1 * st.conc[js][i] - (p.r_kd1 + gam_s + gam_s1) * st.sorb[js][i]));
                            ssorb = st.sorb_n[i];
                            st.sorb_n2[i] = st.sorb2[js][i]
                                + dt * (epsi * (r_ka2 * th_w / ro * psi2 * cc - (p.r_kd2 + gam_s + gam_s1) * st.sorb_n2[i])
                                    + (1.0 - epsi) * (r_ka2 * th_w / ro * psi2 * st.conc[js][i] - (p.r_kd2 + gam_s + gam_s1) * st.sorb2[js][i]));
                            ssorb2 = st.sorb_n2[i];
                        }
                    } else if st.l_dual_neq {
                        let e_mob = th_imob * xmu_l + (1.0 - frac) * ro * xmu_s + c_parent_im;
                        if st.l_linear[js] {
                            let g_mob0 = (2.0 * a_mob - dt * (omega_eff + b_mob)) / d_mob;
                            st.sorb[js][i] = st.sorb[js][i] * g_mob0 + dt * (omega_eff * st.conc[js][i] + 2.0 * e_mob) / d_mob;
                            ssorb = st.sorb[js][i];

                            let s_old2 = st.sorb2[js][i];
                            st.sorb2[js][i] = ((2.0 - (omega + gam_s + gam_s1) * dt) * s_old2
                                + dt * (1.0 - frac) * omega * xks * st.conc[js][i]
                                + dt * (1.0 - frac) * (2.0 * xmu_s)
                                + 2.0 * dt * c_parent_kin)
                                / (2.0 + dt * (omega + gam_s + gam_s1));
                            ssorb2 = st.sorb2[js][i];
                        } else {
                            let cc = st.c_new[i];
                            let s_old = st.sorb[js][i];
                            st.sorb_n[i] = s_old
                                + dt / a_mob
                                    * (epsi * (omega_eff * cc - (omega_eff + b_mob) * st.sorb_n[i] + e_mob)
                                        + (1.0 - epsi) * (omega_eff * st.conc[js][i] - (omega_eff + b_mob) * s_old + e_mob));
                            ssorb = st.sorb_n[i];

                            let s_old2 = st.sorb2[js][i];
                            st.sorb_n2[i] = s_old2
                                + dt * (epsi * (omega * ((1.0 - frac) * sconc * xks * cc - st.sorb_n2[i]) - (gam_s + gam_s1) * st.sorb_n2[i] + (1.0 - frac) * xmu_s + c_parent_kin)
                                    + (1.0 - epsi) * (omega * ((1.0 - frac) * sconc_o * xks * st.conc[js][i] - s_old2) - (gam_s + gam_s1) * s_old2 + (1.0 - frac) * xmu_s + c_parent_kin));
                            ssorb2 = st.sorb_n2[i];
                        }
                    } else if mob {
                        let e_mob = th_imob * xmu_l + (1.0 - frac) * ro * xmu_s + c_parent_im;
                        if st.l_linear[js] {
                            let g_mob0 = (2.0 * a_mob - dt * (omega_eff + b_mob)) / d_mob;
                            st.sorb[js][i] = st.sorb[js][i] * g_mob0 + dt * (omega_eff * st.conc[js][i] + 2.0 * e_mob) / d_mob;
                            ssorb = st.sorb[js][i];
                        } else {
                            let cc = st.c_new[i];
                            let s_old = st.sorb[js][i];
                            st.sorb_n[i] = s_old
                                + dt / a_mob
                                    * (epsi * (omega_eff * cc - (omega_eff + b_mob) * st.sorb_n[i] + e_mob)
                                        + (1.0 - epsi) * (omega_eff * st.conc[js][i] - (omega_eff + b_mob) * s_old + e_mob));
                            ssorb = st.sorb_n[i];
                        }
                    } else if st.l_linear[js] {
                        let s_old = st.sorb[js][i];
                        st.sorb[js][i] = ((2.0 - (omega + gam_s + gam_s1) * dt) * s_old
                            + dt * (1.0 - frac) * omega * xks * st.conc[js][i]
                            + dt * (1.0 - frac) * (2.0 * xmu_s)
                            + 2.0 * dt * c_parent_kin)
                            / (2.0 + dt * (omega + gam_s + gam_s1));
                        ssorb = st.sorb[js][i];
                    } else {
                        let s_old = st.sorb[js][i];
                        let cc = st.c_new[i];
                        st.sorb_n[i] = s_old
                            + dt * (epsi * (omega * ((1.0 - frac) * sconc * xks * cc - st.sorb_n[i]) - (gam_s + gam_s1) * st.sorb_n[i] + (1.0 - frac) * xmu_s + c_parent_kin)
                                + (1.0 - epsi) * (omega * ((1.0 - frac) * sconc_o * xks * st.conc[js][i] - s_old) - (gam_s + gam_s1) * s_old + (1.0 - frac) * xmu_s + c_parent_kin));
                        ssorb = st.sorb_n[i];
                    }
                } else if !st.l_linear[js] {
                    ssorb = st.sorb_n[i];
                    ssorb2 = st.sorb_n2[i];
                }
            }
            let _ = sconc_os;

            // ---- zero-order coefficient
            let gamma_w = if self.i_dual_por > 0 || self.l_dual_perm { self.sink_im[i] } else { 0.0 };
            let tr_adv_im = if gamma_w < 0.0 { -gamma_w } else { 0.0 };
            let omega_im_eff = omega + tr_adv_im;

            let mut g0 = xmu_l * th_w + frac * ro * xmu_s + th_g * xmu_g - st.s_sink[i];
            let mut q0 = xmu_l * th_w + ro * xmu_s + th_g * xmu_g;
            if self.l_dual_perm {
                let w_f = self.w_fracture.max(0.001);
                let c_matrix = st.conc_m[js][i];
                let alpha_eff = (omega + self.alpha_dw) / w_f;
                let gamma_adv = if gamma_w < 0.0 { -gamma_w * c_matrix / w_f } else { 0.0 };
                g0 += alpha_eff * c_matrix + gamma_adv;
            } else if !st.l_equil {
                if st.l_bact {
                    let ss1 = if last { st.sorb_n[i] } else { st.sorb[js][i] };
                    let ss2 = if last { st.sorb_n2[i] } else { st.sorb2[js][i] };
                    g0 += p.r_kd1 * ro * ss1 + p.r_kd2 * ro * ss2;
                    q0 += p.r_kd1 * ro * ss1 + p.r_kd2 * ro * ss2;
                } else if st.l_dual_neq {
                    g0 += omega_im_eff * ssorb + omega * ro * ssorb2;
                    q0 += omega_im_eff * ssorb + omega * ro * ssorb2;
                } else if mob {
                    g0 += omega_im_eff * ssorb;
                } else {
                    g0 += omega * ro * ssorb;
                }
            }
            if let Some(pp) = &parent {
                let cprev = if last { st.conc[js - 1][i] } else { st.c_prev_o[i] };
                let r_w = if pp.gam_w > 0.0 { pp.gam_w } else { pp.mu_w };
                let r_s = if pp.gam_s > 0.0 { pp.gam_s } else { pp.mu_s };
                let r_g = if pp.gam_g > 0.0 { pp.gam_g } else { pp.mu_g };

                let mut cg = cprev * (r_w * th_w + ro * frac * pp.ks * r_s * sconc_p + th_g * pp.henry * r_g);
                let mut cg1 = cg;
                if !st.l_equil {
                    if mob || self.i_dual_por > 0 {
                        let aa = st.sorb[js - 1][i] * (th_imob * r_w + (1.0 - frac) * ro * r_s * pp.ks * sconc_ps);
                        if !st.l_nequil {
                            cg += aa;
                        }
                        cg1 += aa;
                        if st.l_dual_neq {
                            let aa2 = r_s * ro * st.sorb2[js - 1][i];
                            if !st.l_nequil {
                                cg += aa2;
                            }
                            cg1 += aa2;
                        }
                    } else if !st.l_bact {
                        let aa = r_s * ro * st.sorb[js - 1][i];
                        if !st.l_nequil {
                            cg += aa;
                        }
                        cg1 += aa;
                    }
                }
                g0 += cg;
                q0 += cg1;
            }
            st.g0[i] = g0;
            st.q0[i] = q0;

            // ---- first-order coefficient
            let tr_adv_m = if gamma_w > 0.0 { gamma_w } else { 0.0 };
            let omega_m_eff = omega + tr_adv_m;

            let mut g1 = -(gam_l + gam_l1) * th_w - (gam_s + gam_s1) * ro * frac * xks * sconc - (gam_g + gam_g1) * th_g * henry;
            if self.l_dual_perm {
                let w_f = self.w_fracture.max(0.001);
                let alpha_eff = (omega + self.alpha_dw) / w_f;
                let gamma_adv = if gamma_w > 0.0 { gamma_w / w_f } else { 0.0 };
                g1 -= alpha_eff + gamma_adv;
            } else if !st.l_equil {
                if st.l_bact {
                    g1 -= th_w * (r_ka1 * psi1 + r_ka2 * psi2);
                    if last && st.l_linear[js] {
                        g1 += dt * th_w * (p.r_kd1 * r_ka1 / (2.0 + dt * (p.r_kd1 + gam_s)) + p.r_kd2 * r_ka2 / (2.0 + dt * (p.r_kd2 + gam_s)));
                    }
                } else if mob {
                    g1 -= omega_m_eff;
                    if last && st.l_linear[js] {
                        g1 += omega_im_eff * dt * omega_m_eff / st.d_mob_i;
                    }
				} else if st.l_dual_neq {
                    g1 -= omega_m_eff + omega * ro * (1.0 - frac) * sconc * xks;
                    if last && st.l_linear[js] {
                        g1 += omega_im_eff * dt * omega_m_eff / st.d_mob_i
                            + omega * ro * (dt * omega * (1.0 - frac) * xks / (2.0 + dt * (omega + gam_s + gam_s1)));
                    }
                } else {
                    g1 -= omega * ro * (1.0 - frac) * sconc * xks;
                    if last && st.l_linear[js] {
                        g1 += omega * ro * (dt * omega * (1.0 - frac) * xks / (2.0 + dt * (omega + gam_s + gam_s1)));
                    }
                }
            }
            st.g1[i] = g1;
			
            // ---- velocity correction (gas-phase diffusion, varying Henry constant with temperature)
            let temp_at = |idx: usize| -> f64 {
                if last {
                    self.heat.as_ref().map(|h| h.temp_n[idx]).unwrap_or(20.0)
                } else {
                    self.heat.as_ref().map(|h| h.temp_o[idx]).unwrap_or(20.0)
                }
            };
            let henry_at = |idx: usize| -> f64 {
                let m_idx = self.mat[idx];
                let base_h = st.pm[js][m_idx].henry;
                if st.l_tdep && js < st.t_dep.len() {
                    let ea_h = st.t_dep[js].henry;
                    if ea_h.abs() > 1e-12 {
                        let t_k = temp_at(idx) + 273.15;
                        let tr = 293.15;
                        let r_gas = 8.314;
                        let tt = (t_k - tr) / (r_gas * t_k * tr);
                        base_h * (ea_h * tt).exp()
                    } else {
                        base_h
                    }
                } else {
                    base_h
                }
            };

            let der_k;
            if i == 0 {
                der_k = (henry_at(j) - henry) / (self.x[1] - self.x[0]);
            } else if i == n - 1 {
                der_k = (henry - henry_at(i - 1)) / (self.x[n - 1] - self.x[n - 2]);
            } else {
                der_k = (henry_at(j) - henry_at(i - 1)) / ((self.x[j] - self.x[i - 1]) / 2.0);
            }
            let vcorr = th_g * diff_g * tau_g * der_k;
            st.v_corr[i] = vcorr;
            if level == 1 {
                st.v_o[i] -= vcorr;
            } else {
                st.v_n[i] -= vcorr;
            }
            // ---- Peclet / Courant / upwind weights
            if i != n - 1 {
                let dx = self.x[j] - self.x[i];
                let mut vv = 0.0;
                let mut vv1 = 0.0;
                if th_w > 1e-6 && thj > 1e-6 {
                    vv = (v.abs() / th_w + vj.abs() / thj) / 2.0;
                    vv1 = (v / th_w + vj / thj) / 2.0;
                }
                let dd = (st.disp[i] + st.disp[j]) / 2.0;
                if last {
                    let mut pec = 99999.0;
                    let mut dtmax = 1e30;
                    let vmax = (v.abs() + vj.abs()) / (th_w + thj);
                    let rmin = st.retard[i].min(st.retard[j]);
                    if dd > 0.0 {
                        pec = vv.abs() * dx / dd;
                    }
                    let cour = vmax * dt / dx / rmin;
                    peclet = peclet.max(pec);
                    courant = courant.max(cour);
                    let mut cour1 = 1.0;
                    if !st.upw && !st.art_d && pec != 99999.0 {
                        cour1 = (1.0f64).min(st.pe_cr / pec.max(0.5));
                    }
                    if epsi < 1.0 && vmax > 1e-20 {
                        dtmax = cour1 * dx * rmin / vmax;
                    }
                    dt_max_c = dt_max_c.min(dtmax);
                } else if st.upw && iter == 1 {
                    let mut pe2 = 11.0;
                    if dd > 0.0 {
                        pe2 = dx * vv1 / dd / 2.0;
                    }
                    st.wc[i] = if vv.abs() < 1e-30 {
                        0.0
                    } else if pe2.abs() > 10.0 {
                        if vv1 > 0.0 { 1.0 } else { -1.0 }
                    } else {
                        (1.0 / pe2.tanh() - 1.0 / pe2).clamp(-1.0, 1.0)
                    };
                }
            }
        }
        (peclet, courant, dt_max_c)
    }

    /// Assemble the tridiagonal system (MatSet).
    fn sol_matset(&self, st: &mut SoluteState, js: usize, level: usize, dt: f64) {
        let n = self.n;
        let epsi = st.epsi;
        let alf = 1.0 - epsi;
        let mut tho = st.th_o.clone();
        let mut thn = st.th_n.clone();
        for i in 0..n {
            let m = self.mat[i];
            if st.l_mob_im[m] {
                tho[i] = (tho[i] - st.th_im[m]).max(0.001);
                thn[i] = (thn[i] - st.th_im[m]).max(0.001);
            }
        }
        let (vo, vn, wc, disp, retard, g0, g1) = (&st.v_o, &st.v_n, &st.wc, &st.disp, &st.retard, &st.g0, &st.g1);
        let c = &st.conc[js];
        let x = &self.x;
        let (mut b, mut d, mut e, mut f) = (st.mb.clone(), st.md.clone(), st.me.clone(), st.mf.clone());
        let mut b1 = x[1] - x[0];
        let (mut f1, mut e1, mut d1) = (st.f1, st.e1, st.d1);
        let (mut bn, mut dn, mut fnn) = (st.bn, st.dn, st.fn_);
        if level == 1 {
            f1 = c[0] * (b1 / 2.0 / dt * tho[0] * retard[0]
                + alf * (-(tho[0] * disp[0] + tho[1] * disp[1]) / b1 / 2.0 - ((2.0 + 3.0 * wc[0]) * vo[0] + vo[1]) / 6.0 + b1 / 12.0 * (3.0 * g1[0] + g1[1])))
                + c[1] * alf * ((tho[0] * disp[0] + tho[1] * disp[1]) / b1 / 2.0 - (vo[0] + (2.0 - 3.0 * wc[0]) * vo[1]) / 6.0 + b1 / 12.0 * (g1[0] + g1[1]))
                + alf * b1 / 6.0 * (2.0 * g0[0] + g0[1]);
            if st.k_bot == -1 {
                f[0] = f1 + alf * st.c_bot[js] * vo[0];
            }
        } else {
            e1 = epsi * (-(thn[0] * disp[0] + thn[1] * disp[1]) / b1 / 2.0 + (vn[0] + (2.0 - 3.0 * wc[0]) * vn[1]) / 6.0 - b1 / 12.0 * (g1[0] + g1[1]));
            d1 = b1 / 2.0 / dt * thn[0] * retard[0]
                + epsi * ((thn[0] * disp[0] + thn[1] * disp[1]) / b1 / 2.0 + ((2.0 + 3.0 * wc[0]) * vn[0] + vn[1]) / 6.0 - b1 / 12.0 * (3.0 * g1[0] + g1[1]));
            let f2 = epsi * b1 / 6.0 * (2.0 * g0[0] + g0[1]);
            f1 += f2;
            if st.k_bot == 1 {
                d[0] = 1.0;
                e[0] = 0.0;
                f[0] = st.c_bot[js];
            }
            if st.k_bot == -1 {
                if vn[0] > 0.0 {
                    e[0] = e1;
                    d[0] = d1;
                    f[0] = f[0] + f2 + epsi * st.c_bot[js] * vn[0];
                } else {
                    d[0] = -1.0;
                    e[0] = 1.0;
                    f[0] = 0.0;
                }
            }
            if st.k_bot == 0 {
                d[0] = -1.0;
                e[0] = 1.0;
                f[0] = 0.0;
            }
        }
        for i in 1..n - 1 {
            let a1 = b1;
            b1 = x[i + 1] - x[i];
            let dx = (x[i + 1] - x[i - 1]) / 2.0;
            if level == 1 {
                f[i] = c[i - 1]
                    * alf
                    * ((tho[i - 1] * disp[i - 1] + tho[i] * disp[i]) / a1 / 2.0 + ((2.0 + 3.0 * wc[i - 1]) * vo[i - 1] + vo[i]) / 6.0 + a1 / 12.0 * (g1[i - 1] + g1[i]))
                    + c[i]
                        * (dx / dt * tho[i] * retard[i]
                            + alf
                                * (-(tho[i - 1] * disp[i - 1] + tho[i] * disp[i]) / a1 / 2.0 - (tho[i + 1] * disp[i + 1] + tho[i] * disp[i]) / b1 / 2.0
                                    - (vo[i + 1] + 3.0 * (wc[i - 1] + wc[i]) * vo[i] - vo[i - 1]) / 6.0
                                    + (a1 * (g1[i - 1] + 3.0 * g1[i]) + b1 * (3.0 * g1[i] + g1[i + 1])) / 12.0))
                    + c[i + 1]
                        * alf
                        * ((tho[i + 1] * disp[i + 1] + tho[i] * disp[i]) / b1 / 2.0 - (vo[i] + (2.0 - 3.0 * wc[i]) * vo[i + 1]) / 6.0 + b1 / 12.0 * (g1[i] + g1[i + 1]))
                    + alf * (a1 * (g0[i - 1] + 2.0 * g0[i]) + b1 * (2.0 * g0[i] + g0[i + 1])) / 6.0;
            } else {
                b[i] = epsi * (-(thn[i - 1] * disp[i - 1] + thn[i] * disp[i]) / a1 / 2.0 - ((2.0 + 3.0 * wc[i - 1]) * vn[i - 1] + vn[i]) / 6.0 - a1 / 12.0 * (g1[i - 1] + g1[i]));
                d[i] = dx / dt * thn[i] * retard[i]
                    + epsi
                        * ((thn[i - 1] * disp[i - 1] + thn[i] * disp[i]) / a1 / 2.0 + (thn[i + 1] * disp[i + 1] + thn[i] * disp[i]) / b1 / 2.0
                            + (vn[i + 1] + 3.0 * (wc[i - 1] + wc[i]) * vn[i] - vn[i - 1]) / 6.0
                            - (a1 * (g1[i - 1] + 3.0 * g1[i]) + b1 * (3.0 * g1[i] + g1[i + 1])) / 12.0);
                e[i] = epsi * (-(thn[i + 1] * disp[i + 1] + thn[i] * disp[i]) / b1 / 2.0 + (vn[i] + (2.0 - 3.0 * wc[i]) * vn[i + 1]) / 6.0 - b1 / 12.0 * (g1[i] + g1[i + 1]));
                f[i] += epsi * (a1 * (g0[i - 1] + 2.0 * g0[i]) + b1 * (2.0 * g0[i] + g0[i + 1])) / 6.0;
            }
        }
        let m = n - 1;
        let mat_top = self.mat[m];
        let henry_top = st.pm[js][mat_top].henry;
        let (g_vol, g_atm) = if st.k_top == -2 && st.d_surf > 1e-10 {
            (st.diff_g[js] * henry_top / st.d_surf, st.diff_g[js] / st.d_surf)
        } else {
            (0.0, 0.0)
        };

        if level == 1 {
            fnn = c[m - 1]
                * alf
                * ((tho[m - 1] * disp[m - 1] + tho[m] * disp[m]) / b1 / 2.0 + ((2.0 + 3.0 * wc[m - 1]) * vo[m - 1] + vo[m]) / 6.0 + b1 / 12.0 * (g1[m - 1] + g1[m]))
                + c[m]
                    * (b1 / 2.0 / dt * tho[m] * retard[m]
                        + alf
                            * (-(tho[m - 1] * disp[m - 1] + tho[m] * disp[m]) / b1 / 2.0 + (vo[m - 1] + (2.0 - 3.0 * wc[m - 1]) * vo[m]) / 6.0
                                + b1 / 12.0 * (g1[m - 1] + 3.0 * g1[m])))
                + alf * b1 / 6.0 * (g0[m - 1] + 2.0 * g0[m]);
            if st.k_top <= 0 {
                f[m] = fnn;
                if vo[m] < 0.0 {
                    f[m] -= alf * vo[m] * st.c_top[js];
                }
                if st.k_top == -2 {
                    f[m] += -alf * g_vol * c[m] + g_atm * st.c_atm;
                }
            }
        } else {
            bn = epsi * (-(thn[m - 1] * disp[m - 1] + thn[m] * disp[m]) / b1 / 2.0 - ((2.0 + 3.0 * wc[m - 1]) * vn[m - 1] + vn[m]) / 6.0 - b1 / 12.0 * (g1[m - 1] + g1[m]));
            dn = b1 / 2.0 / dt * thn[m] * retard[m]
                + epsi * ((thn[m - 1] * disp[m - 1] + thn[m] * disp[m]) / b1 / 2.0 - (vn[m - 1] + (2.0 - 3.0 * wc[m - 1]) * vn[m]) / 6.0 - b1 / 12.0 * (g1[m - 1] + 3.0 * g1[m]));
            let fe = epsi * b1 / 6.0 * (g0[m - 1] + 2.0 * g0[m]);
            fnn += fe;
            if st.k_top > 0 {
                b[m] = 0.0;
                d[m] = 1.0;
                f[m] = st.c_top[js];
            } else {
                b[m] = bn;
                d[m] = dn;
                f[m] += fe;
                if vn[m] < 0.0 {
                    f[m] -= epsi * vn[m] * st.c_top[js];
                }
                if st.k_top == -2 {
                    d[m] += epsi * g_vol;
                }
            }
        }
        st.mb = b;
        st.md = d;
        st.me = e;
        st.mf = f;
        st.f1 = f1;
        st.e1 = e1;
        st.d1 = d1;
        st.bn = bn;
        st.dn = dn;
        st.fn_ = fnn;
    }

    /// Mass-transfer fluxes (MassTran). `epsi` is the time weight of this call.
    fn sol_masstran(&self, st: &mut SoluteState, js: usize, epsi: f64) {
        let n = self.n;
        for i in 0..n - 1 {
            let j = i + 1;
            let dx = self.x[j] - self.x[i];
            st.cv_ch0[js] += epsi * dx * (st.q0[i] + st.q0[j]) / 2.0;
            st.cv_ch1[js] += epsi * dx * (st.q1[i] + st.q1[j]) / 2.0;
            st.cv_chr[js] += epsi * dx * (st.s_sink[i] + st.s_sink[j]) / 2.0;

            if self.l_dual_perm {
                let alpha_eff_i = self.sol.as_ref().map(|s| s.pm[js][self.mat[i]].omega).unwrap_or(0.0) + self.alpha_dw;
                let alpha_eff_j = self.sol.as_ref().map(|s| s.pm[js][self.mat[j]].omega).unwrap_or(0.0) + self.alpha_dw;
                let transf_i = alpha_eff_i * (st.conc[js][i] - st.conc_m[js][i]);
                let transf_j = alpha_eff_j * (st.conc[js][j] - st.conc_m[js][j]);
                st.cv_chim[js] += epsi * dx / 2.0 * (transf_i + transf_j);
            } else if !st.l_equil {
                let (mi, mj) = (self.mat[i], self.mat[j]);
                let (pi_, pj) = (&st.pm[js][mi], &st.pm[js][mj]);
                if st.l_dual_neq {
                    // 1. Mobile-immobile transfer
                    let transf_m_im = pi_.omega * (st.conc[js][i] - st.sorb[js][i]) + pj.omega * (st.conc[js][j] - st.sorb[js][j]);
                    // 2. Mobile kinetic sorption transfer
                    let se_i = (1.0 - st.frac[mi]) * pi_.ks * st.conc[js][i];
                    let se_j = (1.0 - st.frac[mj]) * pj.ks * st.conc[js][j];
                    let transf_kin = st.ro[mi] * pi_.omega * (se_i - st.sorb2[js][i]) + st.ro[mj] * pj.omega * (se_j - st.sorb2[js][j]);
                    st.cv_chim[js] += epsi * dx / 2.0 * (transf_m_im + transf_kin);
                } else if st.l_mob_im[mi] {
                    st.cv_chim[js] += epsi * dx / 2.0 * (pi_.omega * (st.conc[js][i] - st.sorb[js][i]) + pj.omega * (st.conc[js][j] - st.sorb[js][j]));
                } else {
                    let se = |p: &SpeciesMaterial, fr: f64, cc: f64| if cc > 0.0 { (1.0 - fr) * p.ks * cc.powf(p.beta) / (1.0 + p.nu * cc.powf(p.beta)) } else { 0.0 };
                    let sei = se(pi_, st.frac[mi], st.conc[js][i]);
                    let sej = se(pj, st.frac[mj], st.conc[js][j]);
                    st.cv_chim[js] += epsi * dx / 2.0
                        * (st.ro[mi] * pi_.omega * (sei - st.sorb[js][i]) + st.ro[mj] * pj.omega * (sej - st.sorb[js][j]));
                }
            }
        }
    }

    /// Sorbed concentration at the end of the step for linear non-equilibrium transport.
    fn sol_sorbconc(&self, st: &mut SoluteState, js: usize, dt: f64) {
        for i in 0..self.n {
            let m = self.mat[i];
			
            let p_base = &st.pm[js][m];
            let mut p = p_base.clone();
            if st.l_tdep && js < st.t_dep.len() {
                let td = &st.t_dep[js];
                let temp_val = self.heat.as_ref().map(|h| h.temp_n[i]).unwrap_or(20.0);
                let tr = 293.15;
                let r_gas = 8.314;
                let tt = (temp_val + 273.15 - tr) / (r_gas * (temp_val + 273.15) * tr);
                let scale = |base: f64, ea: f64| if ea.abs() > 1e-12 { base * (ea * tt).exp() } else { base };

                p.mu_s = scale(p_base.mu_s, td.mu_s);
                p.gam_s = scale(p_base.gam_s, td.gam_s);
                p.omega = scale(p_base.omega, td.omega);
                p.ks = scale(p_base.ks, td.ks);
                if st.l_bact {
                    p.r_ka1 = scale(p_base.r_ka1, td.nu);
                    p.r_kd1 = scale(p_base.r_kd1, td.henry);
                    p.r_ka2 = scale(p_base.r_ka2, td.mu0_g);
                    p.r_kd2 = scale(p_base.r_kd2, td.omega);
                }
            }
            let th_w = st.th_n[i];
            let theta = (th_w - st.th_im[m]).max(0.001);
            if st.l_moist {
                let scale_m = |rate: f64, k: usize| -> f64 {
                    if rate.abs() < 1e-30 {
                        return rate;
                    }
                    if st.i_moist_dep == 2 {
                        if js < st.moist_tables.len() {
                            let factor = interp_moist_factor(&st.moist_tables[js], k, theta);
                            rate * factor
                        } else {
                            rate
                        }
                    } else if js < st.w_dep.len() {
                        let b = st.w_dep[js].exp_b[k];
                        let th_ref = st.th_ref[js][m][k];
                        if b.abs() > 1e-12 && th_ref > 1e-6 && theta > 1e-6 {
                            rate * (theta / th_ref).powf(b)
                        } else {
                            rate
                        }
                    } else {
                        rate
                    }
                };
                p.mu_s = scale_m(p.mu_s, 1);
                p.gam_s = scale_m(p.gam_s, 4);
            }

            let (frac, ro) = (st.frac[m], st.ro[m]);
            let (gam_s, gam_s1, omega, xks) = (p.mu_s, p.gam_s, p.omega, p.ks);

            if st.l_bact {
                let mut r_ka1 = p.r_ka1;
                let mut r_ka2 = p.r_ka2;
                if st.l_filtr {
                    let (ka1, ka2) = deposit(p.d_c, p.d_p, p.r_ka1, p.r_ka2, theta, st.v_n[i], self.temp(i), self.x_conv, self.t_conv);
                    r_ka1 = ka1;
                    r_ka2 = ka2;
                }
                st.sorb[js][i] += dt * r_ka1 * theta * st.conc[js][i] / ro / (2.0 + dt * (p.r_kd1 + gam_s + gam_s1));
                st.sorb2[js][i] += dt * r_ka2 * theta * st.conc[js][i] / ro / (2.0 + dt * (p.r_kd2 + gam_s + gam_s1));
            } else if st.l_mob_im[m] || self.i_dual_por > 0 {
                let gamma_w = if self.i_dual_por > 0 { self.sink_im[i] } else { 0.0 };
                let tr_adv_m = if gamma_w > 0.0 { gamma_w } else { 0.0 };
                let tr_adv_im = if gamma_w < 0.0 { -gamma_w } else { 0.0 };
                let omega_eff = omega + tr_adv_m;

                let th_im = st.th_im[m];
                let a_mob = th_im + (1.0 - frac) * ro * xks;
                let b_mob = th_im * (p.mu_w + p.gam_w) + (gam_s + gam_s1) * ro * (1.0 - frac) * xks + tr_adv_im;
                let d_mob = 2.0 * a_mob + dt * (omega_eff + b_mob);

                st.sorb[js][i] += dt * omega_eff * st.conc[js][i] / d_mob;
            } else if st.l_dual_neq {
                let gamma_w = if self.i_dual_por > 0 { self.sink_im[i] } else { 0.0 };
                let tr_adv_m = if gamma_w > 0.0 { gamma_w } else { 0.0 };
                let tr_adv_im = if gamma_w < 0.0 { -gamma_w } else { 0.0 };
                let omega_eff = omega + tr_adv_m;

                let th_im = st.th_im[m];
                let a_mob = th_im + (1.0 - frac) * ro * xks;
                let b_mob = th_im * (p.mu_w + p.gam_w) + (gam_s + gam_s1) * ro * (1.0 - frac) * xks + tr_adv_im;
                let d_mob = 2.0 * a_mob + dt * (omega_eff + b_mob);
                
                // 1. Update physical immobile domain
                st.sorb[js][i] += dt * omega_eff * st.conc[js][i] / d_mob;
                // 2. Update chemical kinetic mobile domain
                st.sorb2[js][i] += dt * omega * (1.0 - frac) * xks * st.conc[js][i] / (2.0 + dt * (omega + gam_s + gam_s1));
            } else {
                st.sorb[js][i] += dt * omega * (1.0 - frac) * xks * st.conc[js][i] / (2.0 + dt * (omega + gam_s + gam_s1));
            }

            if self.l_dual_perm {
                let c_f_new = st.conc[js][i];
                let c_f_old = st.c_prev_o[i];
                let c_f_mid = (c_f_old + c_f_new) / 2.0;
                let c_m = st.conc_m[js][i];

                let alpha = omega + self.alpha_dw;
                let gamma_w = self.sink_im[i];
                let c_star = if gamma_w >= 0.0 { c_f_mid } else { c_m };
                let gamma_s = gamma_w * c_star + alpha * (c_f_mid - c_m);

                let w_m = (1.0 - self.w_fracture).max(0.001);
                let th_matrix = self.th_matrix_new[i].max(0.001);

                // Exact match to explicit c_m sink in fracture FE equation: delta_M_m = dt * gamma_s
                let delta_c_m = (dt * gamma_s) / (w_m * th_matrix);
                st.conc_m[js][i] = (c_m + delta_c_m).max(0.0);
            }
        }
    }

    /// Cumulative boundary fluxes (TLInf, solute part).
    pub fn solute_cumulate(&mut self, run_off: f64) -> Vec<SoluteTLevel> {
        let dt = self.dt;
        let n = self.n;
        let st = self.sol.as_mut().unwrap();
        let mut out = vec![];
        for j in 0..st.ns {
            st.cum_ch[j][0] -= st.cv_top[j] * dt;
            st.cum_ch[j][1] += st.cv_bot[j] * dt;
            st.cum_ch[j][2] += st.cv_ch0[j] * dt;
            st.cum_ch[j][3] += st.cv_ch1[j] * dt;
            st.cum_ch[j][4] += st.cv_chr[j] * dt;
            st.cum_ch[j][5] += st.cv_chim[j] * dt;
            st.c_cum_t[j] += (st.cv_top[j] - st.cv_bot[j] - st.cv_ch0[j] - st.cv_ch1[j] + st.cv_chr[j]) * dt;
            st.c_cum_a[j] += (st.cv_bot[j].abs() + st.cv_top[j].abs() + st.cv_ch0[j].abs() + st.cv_ch1[j].abs() + st.cv_chr[j].abs()) * dt;
            st.cum_ch[j][9] += run_off * st.c_top[j] * dt;
            out.push(SoluteTLevel {
                cv_top: -st.cv_top[j],
                cv_bot: st.cv_bot[j],
                cum_top: st.cum_ch[j][0],
                cum_bot: st.cum_ch[j][1],
                cum_ch0: st.cum_ch[j][2],
                cum_ch1: st.cum_ch[j][3],
                c_top: st.conc[j][n - 1],
                c_root: st.c_root[j],
                c_bot: st.conc[j][0],
                cv_root: st.cv_chr[j],
                cum_root: st.cum_ch[j][4],
                cum_neq: st.cum_ch[j][5],
            });
        }
        out
    }

    /// Solute mass in sub-regions and solute mass balance (SubReg, solute part).
    pub fn sub_reg_solute(&mut self, p_level: usize, rec: &mut BalanceOut) {
        let n = self.n;
        let st = self.sol.as_mut().unwrap();
        let n_lay = rec.sub.len();
        for s in rec.sub.iter_mut() {
            s.c_vol = vec![0.0; st.ns];
            s.c_mean = vec![0.0; st.ns];
        }
        rec.total.c_vol = vec![0.0; st.ns];
        rec.total.c_mean = vec![0.0; st.ns];
        rec.sol_bal_t = vec![0.0; st.ns];
        rec.sol_bal_r = vec![0.0; st.ns];

        let w_f = if self.l_dual_perm { self.w_fracture } else { 1.0 };
        let w_m = (1.0 - w_f).max(0.0);

        for js in 0..st.ns {
            let mut vol = 0.0;
            let mut vol_im = 0.0;
            let mut delt_c = 0.0;
            let mut c_tot = 0.0;
            let mut a_tot = 0.0;
            let mut sub_v = vec![0.0; n_lay];
            let mut sub_c = vec![0.0; n_lay];
            let mut sub_a = vec![0.0; n_lay];

            for i in (0..n - 1).rev() {
                let j = i + 1;
                let (mi, mj) = (self.mat[i], self.mat[j]);
                let lay = self.lay[i] - 1;
                let dx = self.x[j] - self.x[i];
                let (pi_, pj) = (&st.pm[js][mi], &st.pm[js][mj]);

                let side_f = |m: usize, p: &SpeciesMaterial, th: f64, cc: f64| {
                    let mut thw = th;
                    let thg = (self.ths[m] - thw).max(0.0);
                    if st.l_mob_im[m] {
                        thw = (thw - st.th_im[m]).max(0.001);
                    }
                    let c1 = if !st.l_linear[js] && cc > 0.0 {
                        cc.powf(p.beta - 1.0) / (1.0 + p.nu * cc.powf(p.beta))
                    } else {
                        1.0
                    };
                    cc * (w_f * thw + w_f * st.frac[m] * st.ro[m] * p.ks * c1 + w_f * thg * p.henry)
                };

                let cnew_f = dx / 2.0 * (side_f(mi, pi_, self.th_new[i], st.conc[js][i]) + side_f(mj, pj, self.th_new[j], st.conc[js][j]));
                vol += cnew_f;
                sub_v[lay] += cnew_f;

                // Matrix continuum mass under dual-permeability
                if self.l_dual_perm {
					let side_m = |m: usize, p: &SpeciesMaterial, th_m: f64, cc_m: f64, s_m: f64| {
						w_m * th_m * cc_m + w_m * st.ro[m] * (s_m + (1.0 - st.frac[m]) * p.ks * cc_m)
					};
					let mass_m = dx / 2.0 * (side_m(mi, pi_, self.th_matrix_new[i], st.conc_m[js][i], st.sorb_m[js][i]) + side_m(mj, pj, self.th_matrix_new[j], st.conc_m[js][j], st.sorb_m[js][j]));
					vol += mass_m;
					sub_v[lay] += mass_m;
				}

                if !st.l_equil {
                    let sdside = |idx: usize, m: usize, p: &SpeciesMaterial, s1: f64, s2: f64| {
                        if st.l_dual_neq {
                            let th_im_dyn = if self.i_dual_por > 0 { self.th_new_im[idx] } else { st.th_im[m] };
                            let immob_mass = s1 * th_im_dyn;
                            let kin_mass = st.ro[m] * s2;
                            immob_mass + kin_mass
                        } else if st.l_mob_im[m] || self.i_dual_por > 0 {
                            let th_im_dyn = if self.i_dual_por > 0 { self.th_new_im[idx] } else { st.th_im[m] };
                            let s1_sc = if !st.l_linear[js] && s1 > 0.0 {
                                s1.powf(p.beta - 1.0) / (1.0 + p.nu * s1.powf(p.beta))
                            } else {
                                1.0
                            };
                            s1 * (th_im_dyn + (1.0 - st.frac[m]) * st.ro[m] * p.ks * s1_sc)
                        } else {
                            st.ro[m] * s1
                        }
                    };
                    vol_im += dx / 2.0 * (
                        sdside(i, mi, pi_, st.sorb[js][i], st.sorb2[js][i])
                        + sdside(j, mj, pj, st.sorb[js][j], st.sorb2[js][j])
                    );
                }

                let ce = (st.conc[js][i] + st.conc[js][j]) / 2.0;
                c_tot += ce * dx;
                a_tot += dx;
                sub_c[lay] += ce * dx;
                sub_a[lay] += dx;

                if p_level == 0 {
                    st.sol_in[js][i] = cnew_f;
                } else {
                    delt_c += (st.sol_in[js][i] - cnew_f).abs();
                }
            }

            for l in 0..n_lay {
                rec.sub[l].c_vol[js] = sub_v[l];
                rec.sub[l].c_mean[js] = if sub_a[l] > 0.0 { sub_c[l] / sub_a[l] } else { 0.0 };
            }
            rec.total.c_vol[js] = vol + vol_im;
            rec.total.c_mean[js] = if a_tot > 0.0 { c_tot / a_tot } else { 0.0 };

            if p_level == 0 {
                st.c_vol_i[js] = vol + vol_im;
            } else {
                let bal_t = vol + vol_im - st.c_vol_i[js] + st.c_cum_t[js];
                rec.sol_bal_t[js] = bal_t;
                let cc = delt_c.max(st.c_cum_a[js]);
                if cc > 1e-25 {
                    rec.sol_bal_r[js] = bal_t.abs() / cc * 100.0;
                }
            }
        }
    }
	/// Calculate flux concentration for a species (port of FluxConc in SOLUTE.FOR).
    pub fn flux_conc(&self, js: usize) -> Vec<f64> {
        let n = self.n;
        let mut conc_f = vec![0.0; n];
        let Some(st) = self.sol.as_ref() else {
            return conc_f;
        };

        for i in 0..n {
            let m = self.mat[i];
            let mut th_w = st.th_n[i];
            if st.l_mob_im[m] {
                th_w = (th_w - st.th_im[m]).max(0.001);
            }

            let qw = st.v_n[i];
            let c_grad = if i == 0 {
                (st.conc[js][1] - st.conc[js][0]) / (self.x[1] - self.x[0])
            } else if i == n - 1 {
                (st.conc[js][n - 1] - st.conc[js][n - 2]) / (self.x[n - 1] - self.x[n - 2])
            } else {
                (st.conc[js][i + 1] - st.conc[js][i - 1]) / (self.x[i + 1] - self.x[i - 1])
            };

            let mut cf = st.conc[js][i];
            if qw.abs() > 1e-20 {
                cf = st.conc[js][i] - (st.disp[i] * th_w / qw) * c_grad;
            }
            conc_f[i] = cf.max(0.0);
        }

        conc_f
    }
}

/// Calculates blocking factor psi for attachment/straining (port of Blocking in SOLUTE.FOR).
pub fn blocking(i_psi: i32, s_max: f64, x_depth: f64, ss: f64, d_c: f64, s_max2: f64) -> f64 {
	let mut psi = 1.0;
	match i_psi {
		1 => {
			// Langmuirian blocking
			if s_max > 0.0 {
				psi = (1.0 - ss / s_max).max(0.0);
			}
		}
		2 => {
			// Ripening
			if s_max > 0.0 {
				psi = (ss.powf(s_max)).max(psi);
			}
		}
		3 => {
			// Random Sequential Adsorption (RSA)
			if s_max > 0.0 {
				let b_inf = 1.0 / s_max;
				let s_inf = 0.546;
				let m_inf = d_c;
				let const_val = s_inf * b_inf * ss;
				if ss <= 0.8 * s_max {
					psi = 1.0 - 4.0 * const_val + 3.08 * const_val.powi(2) + 1.4069 * const_val.powi(3);
				} else {
					psi = (1.0 - b_inf * ss).powi(3) / (2.0 * m_inf.powi(2) * b_inf.powi(3));
				}
				psi = psi.max(0.0);
			}
		}
		4 => {
			// Depth-dependent straining
			if s_max > 0.0 && d_c > 0.0 {
				psi = ((x_depth.abs() + d_c) / d_c).powf(-s_max);
			}
		}
		5 => {
			// Combined straining and Langmuirian blocking
			if s_max > 0.0 && d_c > 0.0 {
				psi = ((x_depth.abs() + d_c) / d_c).powf(-s_max);
			}
			if s_max2 > 0.0 {
				psi *= (1.0 - ss / s_max2).max(0.0);
			}
		}
		_ => {}
	}
	psi
}

/// Single-collector clean-bed filtration deposition rate (port of Deposit in SOLUTE.FOR).
pub fn deposit(
	d_c1: f64,
	d_p1: f64,
	alfa1: f64,
	alfa2: f64,
	theta: f64,
	q: f64,
	temp_c: f64,
	x_conv: f64,
	t_conv: f64,
) -> (f64, f64) {
	let d_c = d_c1 / x_conv;
	let d_p = d_p1 / x_conv;
	if d_p <= 0.0 && d_c <= 0.0 {
		return (0.0, 0.0);
	}
	let pi = std::f64::consts::PI;
	let mu = 0.00093; // Pa·s
	let bk = 1.38048e-23; // J/K
	let h_hamaker = 1e-20; // J
	let g = 9.81; // m/s^2
	let ro_p = 1080.0; // kg/m^3
	let ro_f = 998.0; // kg/m^3

	let veloc = (q / x_conv * t_conv).abs();
	let p_veloc = if theta > 1e-6 { veloc / theta } else { 0.0 };

	let eta = if veloc > 0.0 {
		let gamma = (1.0 - theta).max(0.0).powf(1.0 / 3.0);
		let as_corr = 2.0 * (1.0 - gamma.powi(5)) / (2.0 - 3.0 * gamma + 3.0 * gamma.powi(5) - 2.0 * gamma.powi(6));

		let n_pe = 3.0 * pi * mu * d_p * d_c * veloc / (bk * (temp_c + 273.15));
		let e_diff = 4.0 * as_corr.powf(1.0 / 3.0) * n_pe.powf(-2.0 / 3.0);

		let n_lo = 4.0 * h_hamaker / (9.0 * pi * mu * d_p.powi(2) * veloc);
		let n_r = d_p / d_c;
		let e_inter = as_corr * n_lo.powf(1.0 / 8.0) * n_r.powf(15.0 / 8.0);

		let n_g = g * (ro_p - ro_f) * d_p.powi(2) / (18.0 * mu * veloc);
		let e_grav = 0.00338 * as_corr * n_g.powf(1.2) * n_r.powf(-0.4);

		e_diff + e_inter + e_grav
	} else {
		0.0
	};

	let ka1 = 3.0 * (1.0 - theta) / (2.0 * d_c) * eta * alfa1 * p_veloc / t_conv;
	let ka2 = 3.0 * (1.0 - theta) / (2.0 * d_c) * eta * alfa2 * p_veloc / t_conv;
	(ka1, ka2)
}
/// Calculates accessible water content and accelerated colloid velocity via Burdine model
/// (port of subroutine Exclusion in SOLUTE.FOR).
pub fn exclusion(
    theta: f64,
    v: f64,
    th_c_param: f64,
    qr: f64,
    qs: f64,
    _alpha: f64,
    n: f64,
    l_param: f64,
) -> (f64, f64) {
    if th_c_param.abs() < 1e-20 || qs <= 0.0 {
        return (theta, v);
    }

    let swr = (qr / qs).clamp(0.0, 1.0);
    let sw_c = (th_c_param / qs).clamp(0.0, 1.0);

    // Accessible water content to colloid
    let th_c = (theta - qs * sw_c).max(0.001);

    let sw = (theta / qs).clamp(0.0, 1.0);
    let sw_eff = ((sw - swr) / (1.0 - swr)).clamp(0.0, 1.0);
    let sw_c_eff = ((sw_c - swr) / (1.0 - swr)).clamp(0.0, 1.0);

    let vgn1 = n;
    let vgm1 = 1.0 - 1.0 / vgn1;
    let vgn2 = vgn1 + 1.0;
    let vgm2 = 1.0 - 2.0 / vgn2;

    // Colloid relative permeability according to Burdine model
    let kr_c = if sw_eff > sw_c_eff && vgm2 > 0.0 {
        let term1 = (1.0 - sw_c_eff.powf(1.0 / vgm2)).max(0.0).powf(vgm2);
        let term2 = (1.0 - sw_eff.powf(1.0 / vgm2)).max(0.0).powf(vgm2);
        sw_eff.powi(2) * (term1 - term2)
    } else {
        0.0
    };

    // Mualem water relative permeability
    let vel_c = if sw_eff > 0.0 && vgm1 > 0.0 {
        let term_w = 1.0 - (1.0 - sw_eff.powf(1.0 / vgm1)).max(0.0).powf(vgm1);
        let krw = sw_eff.powf(l_param) * term_w.powi(2);
        if krw > 1e-30 {
            v * kr_c / krw
        } else {
            v
        }
    } else {
        v
    };

    (th_c, vel_c)
}

/// Interpolate moisture reduction factor from a user table (MoistDep.in in HYDRUS-1D).
pub fn interp_moist_factor(tab: &MoistDepTable, par_idx: usize, theta: f64) -> f64 {
    let n = tab.theta.len();
    if n == 0 {
        return 1.0;
    }
    let fac = &tab.factors[par_idx];
    if fac.is_empty() {
        return 1.0;
    }
    if theta <= tab.theta[0] {
        return fac[0];
    }
    if theta >= tab.theta[n - 1] {
        return fac[n - 1];
    }
    let mut low = 0;
    let mut high = n - 1;
    while high - low > 1 {
        let mid = (low + high) / 2;
        if tab.theta[mid] <= theta {
            low = mid;
        } else {
            high = mid;
        }
    }
    let dth = tab.theta[high] - tab.theta[low];
    if dth.abs() < 1e-12 {
        return fac[low];
    }
    let ratio = (theta - tab.theta[low]) / dth;
    fac[low] + ratio * (fac[high] - fac[low])
}