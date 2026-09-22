//! Bob Lenhard air-entrapment hysteresis (port of HYSTER.FOR, iHyst = 3).
//!
//! This implements the reversal-point ("scanning curve") bookkeeping used by
//! HYSTER.FOR's Hyst/HysterIni/Path/HawPath/Drain/Dry/Wet/UpdateHyst routines,
//! ported as literally as practical (including the odd but intentional
//! `PHSW(n) = -hOld(n)` reset performed at the top of every call, which is
//! what tells HawPath whether the trial head is wetting or drying relative
//! to the last converged state).
//!
//! IPATH (number of reversal curves remembered) is hard-coded to 7, matching
//! `iPath=7` in HYSTER.FOR's `Hyst` subroutine.

use crate::sim::Simulation;

/// Number of reversal curves tracked per node (`iPath` in HYSTER.FOR).
const IPATH: usize = 7;

/// Persistent per-node hysteresis state (the `/glob/` common block in HYSTER.FOR).
/// Arrays are 1-based: index 0 is unused padding so that Fortran's `RHSW(N,k)`
/// maps directly to `rhsw[n][k]`.
#[derive(Clone)]
pub struct LenhardState {
    pub rhsw: Vec<[f64; IPATH + 1]>,
    pub rasw: Vec<[f64; IPATH + 1]>,
    pub sarw: Vec<f64>,
    pub rswaw: Vec<f64>,
    pub phsw: Vec<f64>,
    pub mpsw: Vec<i32>,
    pub jjh: Vec<i32>,
    pub ipsw: Vec<i32>,
}

impl LenhardState {
    pub fn new(n: usize) -> Self {
        LenhardState {
            rhsw: vec![[0.0; IPATH + 1]; n],
            rasw: vec![[0.0; IPATH + 1]; n],
            sarw: vec![0.0; n],
            rswaw: vec![0.0; n],
            phsw: vec![0.0; n],
            mpsw: vec![3; n],
            jjh: vec![0; n],
            ipsw: vec![1; n],
        }
    }
}

/// Where the computed water content for this call should be written
/// (mirrors which Fortran array — ThOld / ThEq / ThNew — is passed as
/// `Theta` into `Hyst` at the various call sites).
#[derive(Clone, Copy)]
pub enum ThetaTarget {
    Old,
    Eq,
    New,
}

/// Per-material properties used by one `Hyst` node evaluation
/// (the `/properties/` common block in HYSTER.FOR).
struct Props {
    alphad: f64,
    alphai: f64,
    xn: f64,
    xm: f64,
    xxm: f64,
    xnw: f64,
    xmw: f64,
    xxmw: f64,
    sm: f64,
    sarwi: f64,
}

/// Transient per-call scratch values (the `/local/` common block in HYSTER.FOR).
#[derive(Default, Clone, Copy)]
struct LocalVars {
    sw: f64,
    sraw: f64,
    #[allow(dead_code)]
    sat: f64,
    rsw: f64,
    raw: f64,
    permw: f64,
    mpsw1: i32,
    jjjh: i32,
    ipsw1: i32,
    ilsw: i32,
    esw: f64,
    esat: f64,
    asw: f64,
    dww: f64,
}

impl Simulation {
    fn ls(&self) -> &LenhardState {
        self.lenhard.as_ref().expect("lenhard state initialised")
    }
    fn ls_mut(&mut self) -> &mut LenhardState {
        self.lenhard.as_mut().expect("lenhard state initialised")
    }

    fn lenhard_init_state(&mut self) {
        if self.lenhard.is_none() {
            self.lenhard = Some(LenhardState::new(self.n));
        }
    }

    /// HysterIni: initialise reversal-point arrays for the given initial
    /// branch (`iHyst` 1 = main drainage, 2 = main imbibition, 3 = primary
    /// drainage with maximum entrapped air — the latter is not reachable
    /// from `IKappa` alone but is kept for completeness).
    fn lenhard_hyster_ini(&mut self) {
        let n = self.n;
        // IPATH = 7 is odd -> Fortran's "label 160" branch of HysterIni.
        for ij in (4..=IPATH - 1).step_by(2) {
            for ji in 0..n {
                self.ls_mut().rasw[ji][ij] = 1.0;
            }
        }
        for ij in (3..=IPATH).step_by(2) {
            for ji in 0..n {
                self.ls_mut().rasw[ji][ij] = 0.0;
            }
        }
        for i in 0..n {
            let ls = self.ls_mut();
            ls.rasw[i][1] = 1.0;
            ls.rhsw[i][1] = 0.0;
            ls.mpsw[i] = 3;
            ls.sarw[i] = 0.0;
            ls.rswaw[i] = 1.0;
        }
        let i_hyst = self.lenhard_i_hyst;
        for i in 0..n {
            let ls = self.ls_mut();
            match i_hyst {
                1 => {
                    ls.rhsw[i][2] = 0.0;
                    ls.rasw[i][2] = 1.0;
                    ls.jjh[i] = 0;
                    ls.phsw[i] = 0.0;
                    ls.rswaw[i] = 1.0;
                    ls.ipsw[i] = 1;
                }
                2 => {
                    ls.rhsw[i][2] = 1.0e3;
                    ls.rasw[i][2] = 0.0;
                    ls.phsw[i] = 1.0e3;
                    ls.jjh[i] = 1;
                    ls.rswaw[i] = 0.0;
                    ls.ipsw[i] = 2;
                }
                3 => {
                    ls.rhsw[i][2] = 1.0e5;
                    ls.rhsw[i][3] = 0.0;
                    ls.rasw[i][2] = 0.0;
                    ls.rasw[i][3] = 1.0;
                    ls.phsw[i] = 0.0;
                    ls.jjh[i] = 0;
                    ls.rswaw[i] = 0.0;
                    ls.ipsw[i] = 3;
                }
                _ => {}
            }
        }
    }

    /// Hyst: evaluate hysteretic K-S-P relations for every node.
    ///
    /// `i_kappa`: initial-branch flag (-1 drying, +1 wetting), only used
    /// when `i_kod == 1`. `i_kod`: 1 = initialisation (calls HysterIni and
    /// commits reversal points), 2 = trial evaluation inside a Picard
    /// iteration (no commit), 3 = final evaluation after convergence
    /// (commits reversal points). `target` selects which Theta array
    /// (Old/Eq/New) receives the computed water content, matching the three
    /// call sites in HYDRUS.FOR / WATFLOW.FOR.
    pub fn lenhard_hyst(&mut self, i_kappa: i32, i_kod: i32, target: ThetaTarget) {
        let n = self.n;
        self.lenhard_init_state();
        if i_kod == 1 {
            self.lenhard_i_hyst = if i_kappa == -1 { 1 } else { 2 };
            self.lenhard_hyster_ini();
        }
        for i in 0..n {
            let m = self.mat[i];
            let alphad = self.par_d[m][2];
            let alphai = self.par_w[m][2];
            let xn = self.par_d[m][3];
            let xm = 1.0 - 1.0 / xn;
            let xxm = 1.0 / xm;
            let xnw = self.par_w[m][3];
            let xmw = 1.0 - 1.0 / xnw;
            let xxmw = 1.0 / xmw;
            let sm = self.par_d[m][0] / self.par_d[m][1];
            let sarwi = (self.par_d[m][1] - self.par_w[m][1]) / (self.par_d[m][1] - self.par_d[m][0]);
            let p = Props { alphad, alphai, xn, xm, xxm, xnw, xmw, xxmw, sm, sarwi };

            let h_w = self.h_new[i];
            let h_a = 0.0f64;
            // "Initialize arrays for hysteretic saturations" in HYSTER.FOR:
            // PHSW(n) is reset from the last converged head on every call,
            // before HawPath uses it to decide drying vs. wetting.
            self.ls_mut().phsw[i] = -self.h_old[i];

            let lv = self.lenhard_path(i, h_w, h_a, &p);

            if i_kod == 1 || i_kod == 3 {
                self.lenhard_update(h_w, h_a, i, &p, lv);
            }

            let theta = lv.sw * self.par_d[m][1];
            match target {
                ThetaTarget::Old => self.th_old[i] = theta,
                ThetaTarget::Eq => self.th_eq[i] = theta,
                ThetaTarget::New => self.th_new[i] = theta,
            }
            self.con[i] = self.par_d[m][4] * lv.permw;
            self.cap[i] = lv.dww * self.par_d[m][1];
            self.kappa[i] = if self.ls().jjh[i] == 0 { -1 } else { 1 };
        }
    }

    /// Path: assign SARWI-derived quantities, seed SARW on a fresh RSWAW,
    /// then call HawPath.
    fn lenhard_path(&self, i: usize, h_w: f64, h_a: f64, p: &Props) -> LocalVars {
        let mut lv = LocalVars::default();
        if p.sarwi == 0.0 {
            lv.raw = 0.0;
            lv.sraw = 0.0;
        } else {
            lv.raw = 1.0 / p.sarwi - 1.0;
        }
        let haw = (h_a - h_w).max(0.0);
        lv.jjjh = self.ls().jjh[i];
        lv.ipsw1 = self.ls().ipsw[i];
        lv.mpsw1 = self.ls().mpsw[i];
        lv.rsw = self.ls().rswaw[i];
        // Note: SARW(n) is only re-seeded here (in-place); since we don't
        // commit it back to persistent state except via UpdateHyst, this
        // mirrors Fortran's transient read of SARW(N) for this call.
        lv.sraw = if lv.rsw == 0.0 { p.sarwi } else { self.ls().sarw[i] };
        self.lenhard_hawpath(i, haw, p, &mut lv);
        lv
    }

    /// HawPath: determine which branch (main drainage, wetting scanning,
    /// drying scanning) the trial head lies on and dispatch to Drain/Dry/Wet.
    fn lenhard_hawpath(&self, i: usize, haw: f64, p: &Props, lv: &mut LocalVars) {
        #[derive(Clone, Copy, PartialEq)]
        enum Entry {
            Top,
            Scan,
        }
        let mut entry = Entry::Top;
        loop {
            if entry == Entry::Top {
                let rhsw2 = self.ls().rhsw[i][2];
                if haw >= rhsw2 {
                    self.lenhard_drain(haw, p, lv);
                    lv.ilsw = 0;
                    lv.jjjh = 0;
                    lv.ipsw1 = 1;
                    lv.mpsw1 = 3;
                    return;
                }
            }
            let phsw_i = self.ls().phsw[i];
            let wetting_cond = ((haw < phsw_i && lv.jjjh == 0) || (haw <= phsw_i && lv.jjjh == 1)) && lv.mpsw1 == 3
                || lv.mpsw1 == 1
                || (lv.jjjh == 1 && haw == 0.0);
            if wetting_cond {
                let rhsw_ipath = self.ls().rhsw[i][IPATH];
                if lv.ipsw1 == IPATH as i32 && haw > rhsw_ipath {
                    lv.mpsw1 = 0;
                    lv.ipsw1 = IPATH as i32 - 1;
                    entry = Entry::Top;
                    continue;
                }
                if haw == 0.0 {
                    if lv.ipsw1 < IPATH as i32 {
                        lv.mpsw1 = 3;
                    }
                    lv.ipsw1 = 2;
                    lv.jjjh = 1;
                    self.lenhard_wet(i, haw, p, lv);
                    return;
                }
                if lv.jjjh == 0 && lv.mpsw1 == 3 {
                    lv.ipsw1 += 1;
                    if lv.ipsw1 >= IPATH as i32 {
                        lv.mpsw1 = 1;
                        lv.ipsw1 = IPATH as i32;
                    }
                }
                let ipsw2 = lv.ipsw1;
                let mut l = 0i32;
                while l <= ipsw2 - 1 {
                    let idx1 = (lv.ipsw1 - l) as usize;
                    let idx2 = (lv.ipsw1 - 1 - l) as usize;
                    let rh1 = self.ls().rhsw[i][idx1];
                    let rh2 = self.ls().rhsw[i][idx2];
                    if haw < rh1 && haw > rh2 {
                        lv.ipsw1 -= l;
                        if lv.ipsw1 < IPATH as i32 {
                            lv.mpsw1 = 3;
                        }
                        self.lenhard_wet(i, haw, p, lv);
                        lv.ilsw = l;
                        lv.jjjh = 1;
                        return;
                    }
                    l += 2;
                }
                // Passed through the wetting-scanning loop without a match
                // (mirrors HAWPATH's error branch) — leave lv as computed
                // so far; con/cap/theta for this node keep their previous
                // values for this call.
                return;
            } else {
                let drying_cond = ((haw > phsw_i && lv.jjjh == 1) || (haw >= phsw_i && lv.jjjh == 0)) && lv.mpsw1 == 3
                    || lv.mpsw1 == 0;
                if drying_cond {
                    if lv.ipsw1 == 1 {
                        self.lenhard_drain(haw, p, lv);
                        lv.ilsw = 0;
                        lv.jjjh = 0;
                        lv.ipsw1 = 1;
                        lv.mpsw1 = 3;
                        return;
                    }
                    let rhsw_ipath = self.ls().rhsw[i][IPATH];
                    if lv.ipsw1 == IPATH as i32 && haw < rhsw_ipath {
                        lv.mpsw1 = 1;
                        lv.ipsw1 = IPATH as i32 - 1;
                        entry = Entry::Scan;
                        continue;
                    }
                    if lv.jjjh == 1 && lv.mpsw1 == 3 {
                        lv.ipsw1 += 1;
                        if lv.ipsw1 >= IPATH as i32 {
                            lv.mpsw1 = 0;
                            lv.ipsw1 = IPATH as i32;
                        }
                    }
                    let ipsw2 = lv.ipsw1;
                    let mut l = 0i32;
                    while l <= ipsw2 - 2 {
                        let idx1 = (lv.ipsw1 - l) as usize;
                        let idx2 = (lv.ipsw1 - 1 - l) as usize;
                        let rh1 = self.ls().rhsw[i][idx1];
                        let rh2 = self.ls().rhsw[i][idx2];
                        if haw > rh1 && haw < rh2 {
                            lv.ipsw1 -= l;
                            if lv.ipsw1 < IPATH as i32 {
                                lv.mpsw1 = 3;
                            }
                            self.lenhard_dry(i, haw, p, lv);
                            lv.ilsw = l;
                            lv.jjjh = 0;
                            return;
                        }
                        l += 2;
                    }
                    return;
                } else {
                    // Neither wetting nor drying condition matched — this is
                    // HAWPATH's terminal error branch; leave lv unchanged.
                    return;
                }
            }
        }
    }

    /// DRAIN: main drainage k-S-P relations.
    fn lenhard_drain(&self, haw: f64, p: &Props, lv: &mut LocalVars) {
        let asw = (1.0 + (p.alphad * haw).powf(p.xn)).powf(-p.xm).min(1.0);
        lv.asw = asw;
        if asw < lv.rsw {
            lv.rsw = asw;
            lv.sraw = if p.sarwi != 0.0 { (1.0 - asw) / (1.0 + lv.raw * (1.0 - asw)) } else { 0.0 };
        }
        // iPath != 1 in the Lenhard module (iPath is fixed to 7), so ESAT is
        // always 0 here, matching the `ELSE` branch of DRAIN's IPATH check.
        lv.esat = 0.0;
        lv.esw = asw - lv.esat;
        lv.sw = lv.esw * (1.0 - p.sm) + p.sm;
        lv.sat = lv.esat * (1.0 - p.sm);
        let x1 = asw.powf(p.xxm);
        lv.dww = (1.0 - p.sm) * p.alphad * (p.xn - 1.0) * x1 * (1.0 - x1).powf(p.xm);
        let p0 = 1.0 - (1.0 - lv.esw.powf(p.xxm)).powf(p.xm);
        lv.permw = (lv.esw.max(0.0)).powf(0.5) * p0 * p0;
        lv.permw = lv.permw.clamp(0.0, 1.0);
    }

    /// DRY: hysteretic k-S-P relations on a drying scanning curve.
    fn lenhard_dry(&self, i: usize, haw: f64, p: &Props, lv: &mut LocalVars) {
        let ip = lv.ipsw1 as usize;
        let (rh_ip, rh_ip1, ra_ip, ra_ip1) = {
            let ls = self.ls();
            (ls.rhsw[i][ip], ls.rhsw[i][ip - 1], ls.rasw[i][ip], ls.rasw[i][ip - 1])
        };
        let ridswd = (1.0 + (p.alphad * rh_ip).powf(p.xn)).powf(-p.xm);
        let rdiswd = (1.0 + (p.alphad * rh_ip1).powf(p.xn)).powf(-p.xm);
        let swd = (1.0 + (p.alphad * haw).powf(p.xn)).powf(-p.xm);
        let mut asw = ((swd - rdiswd) * (ra_ip - ra_ip1)) / (ridswd - rdiswd) + ra_ip1;
        asw = asw.min(1.0);
        lv.asw = asw;
        lv.esat = lv.sraw * ((asw - lv.rsw) / (1.0 - lv.rsw));
        lv.esw = asw - lv.esat;
        lv.sw = lv.esw * (1.0 - p.sm) + p.sm;
        lv.sat = lv.esat * (1.0 - p.sm);
        let p0 = 1.0 - (1.0 - asw.powf(p.xxm)).powf(p.xm);
        let csatw = lv.sraw / (1.0 - lv.rsw);
        let mut p3 = (1.0 - lv.rsw.powf(p.xxm)).powf(p.xm) - (1.0 - asw.powf(p.xxm)).powf(p.xm);
        p3 *= csatw;
        lv.permw = ((lv.esw.max(0.0)).powf(0.5) * (p0 - p3) * (p0 - p3)).clamp(0.0, 1.0);
        let x1 = swd.powf(p.xxm);
        let x7 = (1.0 - csatw).max(0.0);
        let x8 = (ra_ip - ra_ip1) / (ridswd - rdiswd);
        lv.dww = (1.0 - p.sm) * p.alphad * (p.xn - 1.0) * x1 * x7 * x8 * (1.0 - x1).powf(p.xm);
    }

    /// WET: hysteretic k-S-P relations on a wetting scanning curve.
    fn lenhard_wet(&self, i: usize, haw: f64, p: &Props, lv: &mut LocalVars) {
        let ip = lv.ipsw1 as usize;
        let (rh_ip, rh_ip1, ra_ip, ra_ip1) = {
            let ls = self.ls();
            (ls.rhsw[i][ip], ls.rhsw[i][ip - 1], ls.rasw[i][ip], ls.rasw[i][ip - 1])
        };
        let ridswi = (1.0 + (p.alphai * rh_ip1).powf(p.xnw)).powf(-p.xmw);
        let rdiswi = (1.0 + (p.alphai * rh_ip).powf(p.xnw)).powf(-p.xmw);
        let swi = (1.0 + (p.alphai * haw).powf(p.xnw)).powf(-p.xmw);
        let mut asw = ((swi - ridswi) * (ra_ip - ra_ip1)) / (rdiswi - ridswi) + ra_ip1;
        asw = asw.min(1.0);
        lv.asw = asw;
        lv.esat = lv.sraw * ((asw - lv.rsw) / (1.0 - lv.rsw));
        lv.esw = asw - lv.esat;
        lv.sw = lv.esw * (1.0 - p.sm) + p.sm;
        lv.sat = lv.esat * (1.0 - p.sm);
        let p0 = 1.0 - (1.0 - asw.powf(p.xxm)).powf(p.xm);
        let csatw = lv.sraw / (1.0 - lv.rsw);
        let mut p3 = (1.0 - lv.rsw.powf(p.xxm)).powf(p.xm) - (1.0 - asw.powf(p.xxm)).powf(p.xm);
        p3 *= csatw;
        lv.permw = ((lv.esw.max(0.0)).powf(0.5) * (p0 - p3) * (p0 - p3)).clamp(0.0, 1.0);
        let x1 = swi.powf(p.xxmw);
        let x7 = (1.0 - csatw).max(0.0);
        let x8 = (ra_ip - ra_ip1) / (rdiswi - ridswi);
        lv.dww = (1.0 - p.sm) * p.alphai * (p.xnw - 1.0) * x1 * x8 * x7 * (1.0 - x1).powf(p.xmw);
    }

    /// UpdateHyst: after convergence (or on initialisation), commit the
    /// reversal points implied by this call's trial values.
    fn lenhard_update(&mut self, h_w: f64, h_a: f64, i: usize, p: &Props, lv0: LocalVars) {
        let mut lv = lv0;
        let haw = (h_a - h_w).max(0.0);

        if lv.asw < self.ls().rswaw[i] {
            self.ls_mut().rswaw[i] = lv.asw;
            if p.sarwi != 0.0 {
                let v = (1.0 - lv.esw) / (1.0 + lv.raw * (1.0 - lv.esw));
                self.ls_mut().sarw[i] = v;
            }
        }
        // iPath != 1 always here (fixed to 7), so we never take the early
        // "GOTO 200" return that the fluid-entrapment-only option would.

        if haw == 0.0 && lv.ipsw1 != 1 {
            lv.ipsw1 = 2;
            lv.jjjh = 1;
            let mut m = lv.ipsw1 + 3;
            while m <= IPATH as i32 {
                self.ls_mut().rasw[i][m as usize] = 0.0;
                m += 2;
            }
            let mut m = lv.ipsw1 + 2;
            while m <= IPATH as i32 {
                self.ls_mut().rasw[i][m as usize] = 1.0;
                m += 2;
            }
            let idx = (lv.ipsw1 + 1) as usize;
            let ls = self.ls_mut();
            ls.rhsw[i][idx] = haw;
            ls.rasw[i][idx] = lv.asw;
            ls.phsw[i] = haw;
        }

        if lv.ilsw > 0 && lv.jjjh == 0 {
            let mut m = lv.ipsw1 + 3;
            while m <= IPATH as i32 {
                self.ls_mut().rasw[i][m as usize] = 1.0;
                m += 2;
            }
            let mut m = lv.ipsw1 + 2;
            while m <= IPATH as i32 {
                self.ls_mut().rasw[i][m as usize] = 0.0;
                m += 2;
            }
            let idx = (lv.ipsw1 + 1) as usize;
            let ls = self.ls_mut();
            ls.rhsw[i][idx] = haw;
            ls.rasw[i][idx] = lv.asw;
            ls.phsw[i] = haw;
        } else if lv.ilsw > 0 && lv.jjjh == 1 {
            let mut m = lv.ipsw1 + 3;
            while m <= IPATH as i32 {
                self.ls_mut().rasw[i][m as usize] = 0.0;
                m += 2;
            }
            let mut m = lv.ipsw1 + 2;
            while m <= IPATH as i32 {
                self.ls_mut().rasw[i][m as usize] = 1.0;
                m += 2;
            }
            let idx = (lv.ipsw1 + 1) as usize;
            let ls = self.ls_mut();
            ls.rhsw[i][idx] = haw;
            ls.rasw[i][idx] = lv.asw;
            ls.phsw[i] = haw;
        }

        {
            let rhsw2 = self.ls().rhsw[i][2];
            let rasw3 = self.ls().rasw[i][3];
            if haw >= rhsw2 && rasw3 > 0.0 {
                lv.ipsw1 = 1;
                lv.jjjh = 0;
                let mut m = lv.ipsw1 + 2;
                while m <= IPATH as i32 {
                    self.ls_mut().rasw[i][m as usize] = 0.0;
                    m += 2;
                }
                let mut m = lv.ipsw1 + 3;
                while m <= IPATH as i32 {
                    self.ls_mut().rasw[i][m as usize] = 1.0;
                    m += 2;
                }
                let idx = (lv.ipsw1 + 1) as usize;
                let ls = self.ls_mut();
                ls.rhsw[i][idx] = haw;
                ls.rasw[i][idx] = lv.asw;
                ls.phsw[i] = haw;
            }
        }

        {
            let ls = self.ls_mut();
            ls.jjh[i] = lv.jjjh;
            ls.ipsw[i] = lv.ipsw1;
            ls.mpsw[i] = lv.mpsw1;
        }

        if lv.ipsw1 == IPATH as i32 {
            return;
        }

        let idx = (lv.ipsw1 + 1) as usize;
        if lv.jjjh == 1 && lv.asw > self.ls().rasw[i][idx] {
            let ls = self.ls_mut();
            ls.rhsw[i][idx] = haw;
            ls.rasw[i][idx] = lv.asw;
            ls.phsw[i] = haw;
        }
        if lv.jjjh == 0 && lv.asw < self.ls().rasw[i][idx] {
            let ls = self.ls_mut();
            ls.rhsw[i][idx] = haw;
            ls.rasw[i][idx] = lv.asw;
            ls.phsw[i] = haw;
        }
        if lv.ipsw1 == 1 {
            let rh2 = self.ls().rhsw[i][2];
            self.ls_mut().phsw[i] = rh2;
        }
    }
}