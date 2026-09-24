//! Root water uptake and root growth (port of SINK.FOR).

use crate::material::fq;
use crate::model::*;
use crate::sim::Simulation;

fn falfa_feddes(t_pot: f64, h: f64, p0: f64, p1: f64, p2h: f64, p2l: f64, p3: f64, r2h: f64, r2l: f64) -> f64 {
    let mut p2 = 0.0;
    if t_pot < r2l {
        p2 = p2l;
    }
    if t_pot > r2h {
        p2 = p2h;
    }
    if t_pot >= r2l && t_pot <= r2h {
        p2 = p2h + (r2h - t_pot) / (r2h - r2l) * (p2l - p2h);
    }
    let mut a = 0.0;
    if h > p3 && h < p2 {
        a = (h - p3) / (p2 - p3);
    }
    if h >= p2 && h <= p1 {
        a = 1.0;
    }
    if h > p1 && h < p0 && p0 - p1 > 0.0 {
        a = (h - p0) / (p1 - p0);
    }
    if h >= p2 && p1 == 0.0 && p0 == 0.0 {
        a = 1.0;
    }
    a
}

fn fsalfa(s_shaped: bool, c_red: f64, c50: f64, p3c: f64) -> f64 {
    if s_shaped {
        if c50.abs() > 0.0 {
            1.0 / (1.0 + (c_red / c50).powf(p3c))
        } else {
            0.0
        }
    } else if c_red <= c50 {
        1.0
    } else {
        (1.0 - (c_red - c50) * p3c * 0.01).max(0.0)
    }
}

impl Simulation {
    /// Root water uptake sink term (SetSnk).
    pub fn set_snk(&mut self) {
        let n = self.n;
        let t_pot = self.r_root;
        let n_step = if self.omega_c < 1.0 { 2 } else { 1 };
        let mut omega = 0.0;
        self.v_root = 0.0;
        self.h_root = 0.0;
        let mut a_root = 0.0;
        let mut c_root_sum = vec![0.0; self.n_species()];
        let root = self.prj.root.clone();
        let ns = self.n_species();

        // Node 0 is skipped matching SINK.FOR loop do i = 2, N
        self.sink[0] = 0.0;

        for i_step in 1..=n_step {
            for i in 1..n {
                let dxm = if i == n - 1 {
                    (self.x[i] - self.x[i - 1]) / 2.0
                } else {
                    (self.x[i + 1] - self.x[i - 1]) / 2.0
                };

                if self.beta[i] > 0.0 {
                    let m = self.mat[i];
                    let mut h_red = self.h_new[i];
                    let mut s_alfa = 1.0;
                    if let (Some(ss), true) = (&root.solute_stress, self.l_chem) {
                        let mut c_red = 0.0;
                        for j in 0..ns {
                            c_red += ss.a_osm.get(j).copied().unwrap_or(0.0) * self.conc(j, i);
                        }
                        if ss.additive {
                            h_red += c_red;
                        } else {
                            s_alfa = fsalfa(ss.s_shaped, c_red, ss.c50, ss.p3c);
                        }
                    }
                    let (alfa, p_min) = match &root.stress {
                        RootStress::Feddes { p0, p2h, p2l, p3, r2h, r2l, p_optm } => {
                            let p1 = p_optm.get(m).copied().unwrap_or(-25.0);
                            (falfa_feddes(t_pot, h_red, *p0, p1, *p2h, *p2l, *p3, *r2h, *r2l), *p3)
                        }
                        RootStress::SShaped { p50, exponent } => (1.0 / (1.0 + (h_red / p50).powf(*exponent)), 10.0 * p50),
                    };

                    if i_step != n_step {
                        omega += alfa * s_alfa * self.beta[i] * dxm;
                        continue;
                    }

                    let mut compen = 1.0;
                    if omega < self.omega_c && omega > 0.0 {
                        compen = self.omega_c;
                    } else if omega >= self.omega_c {
                        compen = omega;
                    }

                    self.sink[i] = alfa * s_alfa * self.beta[i] * t_pot / compen;
                    if self.th_new[i] - 0.00025 < self.par_d[m][0] {
                        self.sink[i] = 0.0;
                    }
                    let th_limit = fq(self.model, p_min, &self.par_d[m]);
                    self.sink[i] = self.sink[i].min((0.5 * (self.th_new[i] - th_limit) / self.dt).max(0.0));
                    self.v_root += self.sink[i] * dxm;
                    self.h_root += self.h_new[i] * dxm;
                    for (j, c) in c_root_sum.iter_mut().enumerate() {
                        *c += self.conc(j, i) * dxm;
                    }
                    a_root += dxm;
                } else {
                    self.sink[i] = 0.0;
                }

                if self.beta[i] < 0.0 {
                    // Eddy Woehling's modification: source term at the bottom
                    let m = self.mat[i];
                    self.sink[i] = self.beta[i] * self.r_bot;
                    self.sink[i] = self.sink[i].max(0.5 * (self.th_new[i] - self.par_d[m][1]) / self.dt);
                }
            }
        }

        if a_root > 0.001 {
            self.h_root /= a_root;
            self.set_c_root(&c_root_sum, a_root);
        }
    }

    /// Root-zone distribution beta(x) from the root depth (SetRG).
    pub fn set_rg(&mut self) {
        let n = self.n;
        let mut x_r = self.x_root;
        match &self.prj.root.growth {
            RootGrowthMode::FromAtmosphere => {}
            RootGrowthMode::Table(tab) => {
                if !tab.is_empty() {
                    let t = self.t;
                    x_r = if t <= tab[0].0 {
                        tab[0].1
                    } else if t >= tab[tab.len() - 1].0 {
                        tab[tab.len() - 1].1
                    } else {
                        let mut v = tab[0].1;
                        for w in tab.windows(2) {
                            if t > w[0].0 && t <= w[1].0 {
                                v = w[0].1 + (w[1].1 - w[0].1) * (t - w[0].0) / (w[1].0 - w[0].0);
                            }
                        }
                        v
                    };
                    self.x_root = x_r;
                }
            }
            RootGrowthMode::Logistic { t_min, t_med, t_harv, x_min, x_med, x_max, period, .. } => {
                let (t_min, t_med, t_harv, mut x_min, x_med, x_max, period) = (*t_min, *t_med, *t_harv, *x_min, *x_med, *x_max, *period);
                let t_root = if period > 0.0 { self.t % period } else { self.t };
                if t_root < t_min || t_root > t_harv {
                    for b in self.beta.iter_mut() {
                        *b = 0.0;
                    }
                    return;
                }
                if x_min <= 0.001 {
                    x_min = 0.001;
                }
                let rtm = t_med - t_min;
                let rgr = -(1.0 / rtm) * ((0.0001f64).max(x_min * (x_max - x_med)) / (x_med * (x_max - x_min))).ln();
                self.rg = rgr;
                let tt = t_root - t_min;
                x_r = (x_max * x_min) / (x_min + (x_max - x_min) * (-rgr * tt).exp());
                self.x_root = x_r;
            }
        }
        let xtop = self.x[n - 1];
        let mut sbeta = 0.0;
        for i in 1..n - 1 {
            self.beta[i] = if self.x[i] < xtop - x_r {
                0.0
            } else if self.x[i] < xtop - 0.2 * x_r {
                2.08333 / x_r * (1.0 - (xtop - self.x[i]) / x_r)
            } else {
                1.66667 / x_r
            };
            sbeta += self.beta[i] * (self.x[i + 1] - self.x[i - 1]) / 2.0;
        }
        if sbeta < 0.0001 {
            self.beta[n - 2] = 1.0 / ((self.x[n - 1] - self.x[n - 3]) / 2.0);
        } else {
            for i in 1..n - 1 {
                self.beta[i] /= sbeta;
            }
        }
    }
}
