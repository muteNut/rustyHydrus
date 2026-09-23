//! Soil hydraulic property functions (port of MATERIAL.FOR).
//!
//! `par` uses the HYDRUS layout (0-based here):
//! `[Qr, Qs, Alfa, n, Ks, l, P7, P8, P9, P10, P11]`
//! * modified VG: P7..P10 = Qm, Qa, Qk, Kk
//! * VG with air entry: P7 = Qm (computed)
//! * Durner: P7..P9 = w2, alpha2, n2

use crate::model::{SoilMaterial, SoilModel};
use serde::{Deserialize, Serialize};

pub type Par = [f64; 11];

const EPS1: f64 = 0.999999999999999;

pub fn par_of(m: &SoilMaterial, model: SoilModel, x_conv: f64) -> Par {
    let mut p = [0.0; 11];
    p[0] = m.qr;
    p[1] = m.qs;
    p[2] = m.alpha;
    p[3] = m.n;
    p[4] = m.ks;
    p[5] = m.l;
    match model {
        SoilModel::ModifiedVG => {
            p[6] = m.extra[0].max(m.qs);
            p[7] = m.extra[1].min(m.qr);
            p[8] = m.extra[2].min(m.qs);
            p[9] = m.extra[3].min(m.ks);
        }
        SoilModel::VGAirEntry => {
            let h_entry = 0.02 * x_conv;
            p[6] = m.qr + (m.qs - m.qr) * (1.0 + (m.alpha * h_entry).powf(m.n)).powf(1.0 - 1.0 / m.n);
        }
        SoilModel::Durner => {
            p[6] = m.extra[0];
            p[7] = m.extra[1];
            p[8] = m.extra[2];
        }
		SoilModel::DualPorosityW => {
            p[6] = m.extra[0]; // thr_im
            p[7] = m.extra[1]; // ths_im
            p[8] = m.extra[2]; // Omega
        }
        SoilModel::DualPorosityH => {
            p[6] = m.extra[0]; // thr_im
            p[7] = m.extra[1]; // ths_im
            p[8] = m.extra[2]; // Alfa_im
            p[9] = m.extra[3]; // n_im
            p[10] = m.extra[4]; // Omega
        }
		SoilModel::DualPermeability => {
            p[6] = m.extra[0]; // thr_m
            p[7] = m.extra[1]; // ths_m
            p[8] = m.extra[2]; // alpha_m
            p[9] = m.extra[3]; // n_m
            p[10] = m.extra[4]; // ks_m
        }
        _ => {}
    }
    p
}

fn vg_terms(model: SoilModel, p: &Par) -> (f64, f64, f64, f64, f64, f64) {
    // returns Qm, Qa, Qk, Kk, Qs, Ks
    let (qr, qs, ks) = (p[0], p[1], p[4].max(1e-37));
    let (mut qm, mut qa, mut qk, mut kk) = (qs, qr, qs, ks);
    if model == SoilModel::ModifiedVG {
        qm = p[6];
        qa = p[7];
        qk = p[8];
        kk = p[9];
    }
    if model == SoilModel::VGAirEntry {
        qm = p[6];
    }
    (qm, qa, qk, kk, qs, ks)
}

fn qnorm(x: f64) -> f64 {
    let z = (x / 2f64.sqrt()).abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let erfc = t
        * (-z * z - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
            .exp();
    let erfc = if x < 0.0 { 2.0 - erfc } else { erfc };
    erfc / 2.0
}

fn hmin(alfa: f64, n: f64) -> f64 {
    -(1e300f64.powf(1.0 / n)) / alfa.max(1.0)
}

/// Hydraulic conductivity K(h)
pub fn fk(model: SoilModel, h: f64, p: &Par) -> f64 {
    let alfa = p[2];
    let n = p[3];
    let bpar = p[5];
    match model {
        SoilModel::VanGenuchten
		| SoilModel::ModifiedVG
		| SoilModel::VGAirEntry
		| SoilModel::DualPorosityW
		| SoilModel::DualPorosityH
		| SoilModel::DualPermeability
		| SoilModel::Tabular => {
            let (qm, qa, qk, kk, qs, ks) = vg_terms(model, p);
            let ppar = 2.0;
            let m = 1.0 - 1.0 / n;
            let hh = h.max(hmin(alfa, n));
            let qees = ((qs - qa) / (qm - qa)).min(EPS1);
            let qeek = ((qk - qa) / (qm - qa)).min(qees);
            let hs = -1.0 / alfa * (qees.powf(-1.0 / m) - 1.0).powf(1.0 / n);
            let hk = -1.0 / alfa * (qeek.powf(-1.0 / m) - 1.0).powf(1.0 / n);
            let mut fkv = 0.0;
            if h < hk {
                let qee = (1.0 + (-alfa * hh).powf(n)).powf(-m);
                let qe = (qm - qa) / (qs - qa) * qee;
                let qek = (qm - qa) / (qs - qa) * qeek;
                let mut ffq = 1.0 - (1.0 - qee.powf(1.0 / m)).powf(m);
                let ffqk = 1.0 - (1.0 - qeek.powf(1.0 / m)).powf(m);
                if ffq <= 0.0 {
                    ffq = m * qee.powf(1.0 / m);
                }
                let mut kr = (qe / qek).powf(bpar) * (ffq / ffqk).powf(ppar) * kk / ks;
                if model == SoilModel::VanGenuchten {
                    kr = qe.powf(bpar) * ffq.powf(ppar);
                }
                fkv = (ks * kr).max(1e-37);
            }
            if h >= hk && h < hs {
                let kr = (1.0 - kk / ks) / (hs - hk) * (h - hs) + 1.0;
                fkv = ks * kr;
            }
            if h >= hs {
                fkv = ks;
            }
            fkv
        }
        SoilModel::BrooksCorey => {
            let ks = p[4].max(1e-37);
            let lambda = 2.0;
            let hs = -1.0 / alfa;
            if h < hs {
                let kr = 1.0 / (-alfa * h).powf(n * (bpar + lambda) + 2.0);
                (ks * kr).max(1e-37)
            } else {
                ks
            }
        }
        SoilModel::Kosugi => {
            let ks = p[4].max(1e-37);
            if h < 0.0 {
                let qee = qnorm((-h / alfa).ln() / n);
                let t = qnorm((-h / alfa).ln() / n + n);
                let kr = qee.powf(bpar) * t * t;
                (ks * kr).max(1e-37)
            } else {
                ks
            }
        }
        SoilModel::Durner => {
            let ks = p[4].max(1e-37);
            let (w2, alfa2, n2) = (p[6], p[7], p[8]);
            let m = 1.0 - 1.0 / n;
            let m2 = 1.0 - 1.0 / n2;
            let w1 = 1.0 - w2;
            if h >= 0.0 {
                return ks;
            }
            let sw1 = w1 * (1.0 + (-alfa * h).powf(n)).powf(-m);
            let sw2 = w2 * (1.0 + (-alfa2 * h).powf(n2)).powf(-m2);
            let qe = sw1 + sw2;
            let sv1 = (-alfa * h).powf(n - 1.0);
            let sv2 = (-alfa2 * h).powf(n2 - 1.0);
            let sk1 = w1 * alfa * (1.0 - sv1 * (1.0 + (-alfa * h).powf(n)).powf(-m));
            let sk2 = w2 * alfa2 * (1.0 - sv2 * (1.0 + (-alfa2 * h).powf(n2)).powf(-m2));
            let numer = sk1 + sk2;
            let denom = w1 * alfa + w2 * alfa2;
            let kr = if denom != 0.0 { qe.powf(bpar) * (numer / denom).powi(2) } else { 0.0 };
            (ks * kr).max(1e-37)
        }
    }
}

/// Soil water capacity C(h) = dθ/dh
pub fn fc(model: SoilModel, h: f64, p: &Par) -> f64 {
    let (qr, qs, alfa, n) = (p[0], p[1], p[2], p[3]);
    match model {
        SoilModel::VanGenuchten
		| SoilModel::ModifiedVG
		| SoilModel::VGAirEntry
		| SoilModel::DualPorosityW
		| SoilModel::DualPorosityH
		| SoilModel::DualPermeability
		| SoilModel::Tabular => {
            let (qm, qa, _qk, _kk, _qs, _ks) = vg_terms(model, p);
            let m = 1.0 - 1.0 / n;
            let hh = h.max(hmin(alfa, n));
            let qees = ((qs - qa) / (qm - qa)).min(EPS1);
            let hs = -1.0 / alfa * (qees.powf(-1.0 / m) - 1.0).powf(1.0 / n);
            if h < hs {
                let c1 = (1.0 + (-alfa * hh).powf(n)).powf(-m - 1.0);
                let c2 = (qm - qa) * m * n * alfa.powf(n) * (-hh).powf(n - 1.0) * c1;
                c2.max(1e-37)
            } else {
                0.0
            }
        }
        SoilModel::BrooksCorey => {
            let hs = -1.0 / alfa;
            if h < hs {
                ((qs - qr) * n * alfa.powf(-n) * (-h).powf(-n - 1.0)).max(1e-37)
            } else {
                0.0
            }
        }
        SoilModel::Kosugi => {
            if h < 0.0 {
                let t = (-(-h / alfa).ln().powi(2) / (2.0 * n * n)).exp();
                ((qs - qr) / (2.0 * 3.141592654f64).sqrt() / n / (-h) * t).max(1e-37)
            } else {
                0.0
            }
        }
        SoilModel::Durner => {
            let (w2, alfa2, n2) = (p[6], p[7], p[8]);
            let m = 1.0 - 1.0 / n;
            let m2 = 1.0 - 1.0 / n2;
            let w1 = 1.0 - w2;
            if h >= 0.0 {
                return 0.0;
            }
            let c1a = (1.0 + (-alfa * h).powf(n)).powf(-m - 1.0);
            let c1b = (1.0 + (-alfa2 * h).powf(n2)).powf(-m2 - 1.0);
            (qs - qr) * m * n * alfa.powf(n) * (-h).powf(n - 1.0) * c1a * w1
                + (qs - qr) * m2 * n2 * alfa2.powf(n2) * (-h).powf(n2 - 1.0) * c1b * w2
        }
    }
}

/// Water content θ(h)
pub fn fq(model: SoilModel, h: f64, p: &Par) -> f64 {
    let (qr, qs, alfa, n) = (p[0], p[1], p[2], p[3]);
    match model {
        SoilModel::VanGenuchten
		| SoilModel::ModifiedVG
		| SoilModel::VGAirEntry
		| SoilModel::DualPorosityW
		| SoilModel::DualPorosityH
		| SoilModel::DualPermeability
		| SoilModel::Tabular => {
            let (qm, qa, ..) = vg_terms(model, p);
            let m = 1.0 - 1.0 / n;
            let hh = h.max(hmin(alfa, n));
            let qees = ((qs - qa) / (qm - qa)).min(EPS1);
            let hs = -1.0 / alfa * (qees.powf(-1.0 / m) - 1.0).powf(1.0 / n);
            if h < hs {
                let qee = (1.0 + (-alfa * hh).powf(n)).powf(-m);
                (qa + (qm - qa) * qee).max(1e-37)
            } else {
                qs
            }
        }
        SoilModel::BrooksCorey => {
            let hs = -1.0 / alfa;
            if h < hs {
                (qr + (qs - qr) * (-alfa * h).powf(-n)).max(1e-37)
            } else {
                qs
            }
        }
        SoilModel::Kosugi => {
            if h < 0.0 {
                (qr + (qs - qr) * qnorm((-h / alfa).ln() / n)).max(1e-37)
            } else {
                qs
            }
        }
        SoilModel::Durner => {
            let (w2, alfa2, n2) = (p[6], p[7], p[8]);
            let m = 1.0 - 1.0 / n;
            let m2 = 1.0 - 1.0 / n2;
            let w1 = 1.0 - w2;
            if h >= 0.0 {
                return qs;
            }
            let sw1 = w1 * (1.0 + (-alfa * h).powf(n)).powf(-m);
            let sw2 = w2 * (1.0 + (-alfa2 * h).powf(n2)).powf(-m2);
            (qr + (qs - qr) * (sw1 + sw2)).max(1e-37)
        }
    }
}

/// Effective saturation S_e(h)
pub fn fs(model: SoilModel, h: f64, p: &Par) -> f64 {
    let (qr, qs, alfa, n) = (p[0], p[1], p[2], p[3]);
    match model {
        SoilModel::VanGenuchten
		| SoilModel::ModifiedVG
		| SoilModel::VGAirEntry
		| SoilModel::DualPorosityW
		| SoilModel::DualPorosityH
		| SoilModel::DualPermeability
		| SoilModel::Tabular => {
            let (qm, qa, ..) = vg_terms(model, p);
            let m = 1.0 - 1.0 / n;
            let qees = ((qs - qa) / (qm - qa)).min(EPS1);
            let hs = -1.0 / alfa * (qees.powf(-1.0 / m) - 1.0).powf(1.0 / n);
            if h < hs {
                let hh = h.max(hmin(alfa, n));
                let qee = (1.0 + (-alfa * hh).powf(n)).powf(-m);
                (qee * (qm - qa) / (qs - qa)).max(1e-37)
            } else {
                1.0
            }
        }
        SoilModel::BrooksCorey => {
            if h < -1.0 / alfa {
                (-alfa * h).powf(-n).max(1e-37)
            } else {
                1.0
            }
        }
        SoilModel::Kosugi => {
            if h < 0.0 {
                qnorm((-h / alfa).ln() / n).max(1e-37)
            } else {
                1.0
            }
        }
        SoilModel::Durner => {
            let _ = qr;
            (fq(model, h, p) - p[0]) / (p[1] - p[0])
        }
    }
}

/// Pressure head h(Se)
pub fn fh(model: SoilModel, qe: f64, p: &Par) -> f64 {
    let (_qr, qs, alfa, n) = (p[0], p[1], p[2], p[3]);
    match model {
        SoilModel::VanGenuchten
		| SoilModel::ModifiedVG
		| SoilModel::VGAirEntry
		| SoilModel::DualPorosityW
		| SoilModel::DualPorosityH
		| SoilModel::DualPermeability
		| SoilModel::Tabular => {
            let (qm, qa, ..) = vg_terms(model, p);
            let m = 1.0 - 1.0 / n;
            let hm = hmin(alfa, n);
            let qeem = (1.0 + (-alfa * hm).powf(n)).powf(-m);
            let qee = (qe * (qs - qa) / (qm - qa)).max(qeem).min(EPS1);
            (-1.0 / alfa * (qee.powf(-1.0 / m) - 1.0).powf(1.0 / n)).max(-1e37)
        }
        SoilModel::BrooksCorey => (-1.0 / alfa * qe.max(1e-10).powf(-1.0 / n)).max(-1e37),
        SoilModel::Kosugi => {
            if qe > 0.9999 {
                0.0
            } else if qe < 0.00001 {
                -1e8
            } else {
                let y = qe * 2.0;
                let pp = if y < 1.0 { (-(y / 2.0).ln()).sqrt() } else { (-(1.0 - y / 2.0).ln()).sqrt() };
                let mut x = pp
                    - (1.881796 + 0.9425908 * pp + 0.0546028 * pp.powi(3))
                        / (1.0 + 2.356868 * pp + 0.3087091 * pp * pp + 0.0937563 * pp.powi(3) + 0.021914 * pp.powi(4));
                if y >= 1.0 {
                    x = -x;
                }
                -alfa * (2f64.sqrt() * n * x).exp()
            }
        }
        SoilModel::Durner => {
            if qe > 0.9999 {
                0.0
            } else if qe < 0.00001 {
                -1e8
            } else {
                // bisection in log10(-h) on Se(h) - qe
                let f = |lh: f64| {
                    let h = -(10f64.powf(lh));
                    fs(model, h, p) - qe
                };
                let (mut a, mut b) = (-6.0, 6.0);
                let fa = f(a);
                if fa * f(b) > 0.0 {
                    return -1e8;
                }
                for _ in 0..200 {
                    let c = 0.5 * (a + b);
                    if f(c) * fa > 0.0 {
                        a = c;
                    } else {
                        b = c;
                    }
                }
                -(10f64.powf(0.5 * (a + b))).max(-1e37)
            }
        }
    }
}

pub const NTAB: usize = 100;

/// Pre-tabulated hydraulic properties for one material (GenMat in Fortran)
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MatTable {
    pub h: Vec<f64>,
    pub con: Vec<f64>,
    pub cap: Vec<f64>,
    pub the: Vec<f64>,
}

impl MatTable {
    pub fn build(model: SoilModel, p: &Par, h_tab1: f64, h_tab_n: f64) -> MatTable {
        let dlh = ((-h_tab_n).log10() - (-h_tab1).log10()) / (NTAB as f64 - 1.0);
        let mut t = MatTable { h: vec![], con: vec![], cap: vec![], the: vec![] };
        for i in 0..NTAB {
            let alh = (-h_tab1).log10() + i as f64 * dlh;
            let h = -(10f64.powf(alh));
            t.h.push(h);
            t.con.push(fk(model, h, p));
            t.cap.push(fc(model, h, p));
            t.the.push(fq(model, h, p));
        }
        t
    }
}
