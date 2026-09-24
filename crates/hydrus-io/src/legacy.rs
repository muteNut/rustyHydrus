//! Reader for HYDRUS-1D (version 3/4) project directories.
//!
//! The reading order and the list-directed semantics follow INPUT.FOR
//! (BasInf, NodInf, MatIn, TmIn, RootIn, TempIn, ChemIn, SinkIn) and SetBC.

use hydrus_core::material::MatTable;
use hydrus_core::*;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};

type R<T> = Result<T, HydrusError>;

fn perr<T>(msg: impl Into<String>) -> R<T> {
    Err(HydrusError::Parse(msg.into()))
}

struct Rd {
    name: &'static str,
    lines: Vec<String>,
    pos: usize,
}

struct Toks(Vec<String>);

impl Toks {
    fn f(&self, i: usize) -> R<f64> {
        let s = self.0.get(i).ok_or_else(|| HydrusError::Parse(format!("missing value #{}", i + 1)))?;
        parse_f(s)
    }
    fn i(&self, i: usize) -> R<i32> {
        Ok(self.f(i)?.round() as i32)
    }
    fn b(&self, i: usize) -> R<bool> {
        let s = self.0.get(i).ok_or_else(|| HydrusError::Parse(format!("missing value #{}", i + 1)))?;
        let t = s.trim_matches('.').to_ascii_lowercase();
        match t.as_str() {
            "t" | "true" => Ok(true),
            "f" | "false" => Ok(false),
            _ => perr(format!("expected logical, got '{}'", s)),
        }
    }
    fn len(&self) -> usize {
        self.0.len()
    }
}

fn parse_f(s: &str) -> R<f64> {
    let t = s.replace(['d', 'D'], "e");
    t.parse::<f64>().map_err(|_| HydrusError::Parse(format!("cannot parse number '{}'", s)))
}

impl Rd {
    fn open(path: &Path, name: &'static str) -> R<Rd> {
        let bytes = std::fs::read(path).map_err(|e| HydrusError::Io(format!("{}: {}", path.display(), e)))?;
        let text = String::from_utf8_lossy(&bytes).to_string();
        Ok(Rd { name, lines: text.lines().map(|s| s.to_string()).collect(), pos: 0 })
    }
    fn version(&mut self) -> i32 {
        if let Some(l) = self.lines.get(0) {
            if l.starts_with("Pcp_File_Version=") {
                let v = l["Pcp_File_Version=".len()..].trim().parse().unwrap_or(0);
                self.pos = 1;
                return v;
            }
        }
        self.pos = 0;
        0
    }
    fn skip(&mut self) -> R<()> {
        if self.pos >= self.lines.len() {
            return perr(format!("{}: unexpected end of file", self.name));
        }
        self.pos += 1;
        Ok(())
    }
    fn line(&mut self) -> R<String> {
        if self.pos >= self.lines.len() {
            return perr(format!("{}: unexpected end of file", self.name));
        }
        self.pos += 1;
        Ok(self.lines[self.pos - 1].clone())
    }
    fn peek(&self) -> Option<&str> {
        self.lines.get(self.pos).map(|s| s.as_str())
    }
    /// List-directed read of `n` values (continues on following records if needed).
    fn read(&mut self, n: usize) -> R<Toks> {
        let mut out: Vec<String> = vec![];
        let mut first = true;
        while out.len() < n {
            let l = self.line()?;
            let mut got = 0;
            for tok in l.split(|c: char| c.is_whitespace() || c == ',').filter(|s| !s.is_empty()) {
                if tok.starts_with('/') {
                    break;
                }
                // repeat count r*value
                if let Some((r, v)) = tok.split_once('*') {
                    if let Ok(cnt) = r.parse::<usize>() {
                        for _ in 0..cnt {
                            out.push(v.to_string());
                            got += 1;
                        }
                        continue;
                    }
                }
                out.push(tok.to_string());
                got += 1;
            }
            if first && got == 0 && n > 0 {
                // blank line inside list-directed read: keep reading
            }
            first = false;
        }
        Ok(Toks(out))
    }
    /// Read as many as available on one record, at least `min`.
    fn read_opt(&mut self, min: usize) -> R<Toks> {
        let l = self.line()?;
        let v: Vec<String> = l.split(|c: char| c.is_whitespace() || c == ',').filter(|s| !s.is_empty()).map(|s| s.to_string()).collect();
        if v.len() < min {
            self.pos -= 1;
            return self.read(min);
        }
        Ok(Toks(v))
    }
}

/// Read a HYDRUS-1D project directory (containing SELECTOR.IN, PROFILE.DAT, ATMOSPH.IN).
pub fn read_legacy_project(dir: &Path) -> R<Project> {
    let find = |name: &str| -> R<std::path::PathBuf> {
        let rd = std::fs::read_dir(dir).map_err(|e| HydrusError::Io(format!("{}: {}", dir.display(), e)))?;
        for e in rd.flatten() {
            if e.file_name().to_string_lossy().eq_ignore_ascii_case(name) {
                return Ok(e.path());
            }
        }
        Err(HydrusError::Io(format!("{} not found in {}", name, dir.display())))
    };
    let mut prj = Project::default();
    let mut s = Rd::open(&find("selector.in")?, "Selector.in")?;
    let ver = s.version();

    // ---------------- BasInf
    s.skip()?;
    s.skip()?;
    let hed = s.line()?;
    prj.title = hed.trim().to_string();
    s.skip()?;
    let lunit = s.line()?.trim().to_ascii_lowercase();
    let tunit = s.line()?.trim().to_ascii_lowercase();
    let munit = s.line()?.trim().to_string();
    prj.units.length = match lunit.as_str() {
        "mm" => LengthUnit::Mm,
        "m" => LengthUnit::M,
        _ => LengthUnit::Cm,
    };
    prj.units.time = match tunit.as_str() {
        "sec" | "s" => TimeUnit::Sec,
        "min" => TimeUnit::Min,
        "hours" | "hour" | "h" => TimeUnit::Hours,
        "years" | "year" | "y" => TimeUnit::Years,
        _ => TimeUnit::Days,
    };
    prj.units.mass = munit;
    s.skip()?;
    let t = s.read(10)?;
    let (l_wat, l_chem, l_temp, sink_f, l_root, short_o) = (t.b(0)?, t.b(1)?, t.b(2)?, t.b(3)?, t.b(4)?, t.b(5)?);
    let (l_wdep, atm_bc, l_equil) = (t.b(6)?, t.b(8)?, t.b(9)?);
    let (mut l_snow, mut l_meteo, mut l_vapor, mut l_act_rsu, mut _l_flux) = (false, false, false, false, false);
    if ver == 3 {
        s.skip()?;
        let t = s.read(4)?;
        l_snow = t.b(0)?;
    } else if ver >= 4 {
        s.skip()?;
        let t = s.read_opt(6)?;
        l_snow = t.b(0).unwrap_or(false);
        if t.len() > 2 { l_meteo = t.b(2).unwrap_or(false); }
        if t.len() > 3 { l_vapor = t.b(3).unwrap_or(false); }
        if t.len() > 4 { l_act_rsu = t.b(4).unwrap_or(false); }
        if t.len() > 5 { _l_flux = t.b(5).unwrap_or(false); }
    }
    prj.atmosphere.snow = l_snow;
	prj.root.l_act_rsu = l_act_rsu;
	
    prj.processes = Processes {
        water_flow: l_wat,
        solute: l_chem,
        heat: l_temp,
        root_water_uptake: sink_f,
        root_growth: l_root,
        equilibrium_adsorption: l_equil,
        short_output: short_o,
		vapor: l_vapor,
    };
    s.skip()?;
    let t = s.read(3)?;
    let n_mat = t.i(0)? as usize;
    let _n_lay = t.i(1)?;
    prj.cos_alpha = t.f(2)?;
    if n_mat > 20 {
        return perr("NMat > 20");
    }
    s.skip()?;
    s.skip()?;
    let t = s.read(3)?;
    prj.water.max_iter = t.i(0)? as usize;
    prj.water.tol_th = t.f(1)?;
    prj.water.tol_h = t.f(2)?;
	prj.water.l_w_dep = l_wdep;
    s.skip()?;
    let t = s.read(4)?;
    let mut bc = WaterBc { atmospheric: atm_bc, ..Default::default() };
    bc.top_time_variable = t.b(0)?;
    bc.surface_layer = t.b(1)?;
    bc.kod_top = t.i(2)?;
    prj.water.init_in_water_content = t.b(3)?;
    s.skip()?;
    if ver <= 3 {
        let t = s.read_opt(5)?;
        bc.bot_time_variable = t.b(0)?;
        bc.gwl_flux = t.b(1)?;
        bc.free_drainage = t.b(2)?;
        bc.seepage_face = t.b(3)?;
        bc.kod_bot = t.i(4)?;
        let q = if t.len() > 5 { t.b(5)? } else { false };
        if q {
            bc.drains = Some(DrainSettings::default());
        }
    } else {
        let t = s.read(7)?;
        bc.bot_time_variable = t.b(0)?;
        bc.gwl_flux = t.b(1)?;
        bc.free_drainage = t.b(2)?;
        bc.seepage_face = t.b(3)?;
        bc.kod_bot = t.i(4)?;
        if t.b(5)? {
            bc.drains = Some(DrainSettings::default());
        }
        bc.h_seep = t.f(6)?;
    }
    let has_drain = bc.drains.is_some();
    if (!bc.top_time_variable && bc.kod_top == -1)
        || (!bc.bot_time_variable && bc.kod_bot == -1 && !bc.gwl_flux && !bc.free_drainage && !bc.seepage_face && !has_drain)
    {
        s.skip()?;
        let t = s.read(3)?;
        bc.r_top = t.f(0)?;
        bc.r_bot = t.f(1)?;
        bc.r_root = t.f(2)?;
    }
    if bc.gwl_flux {
        s.skip()?;
        let t = s.read(3)?;
        bc.gwl0l = t.f(0)?;
        bc.aqh = t.f(1)?;
        bc.bqh = t.f(2)?;
    }
    if let Some(d) = bc.drains.as_mut() {
        s.skip()?;
        d.position = s.read(1)?.i(0)?;
        s.skip()?;
        let t = s.read(3)?;
        d.z_bot = -t.f(0)?.abs();
        d.spacing = t.f(1)?;
        d.entrance_resistance = t.f(2)?;
        s.skip()?;
        match d.position {
            1 => {
                d.kh_top = s.read(1)?.f(0)?;
            }
            2 => {
                let t = s.read(3)?;
                d.base_gw = t.f(0)?;
                d.kh_top = t.f(1)?;
                d.wet_perimeter = t.f(2)?;
            }
            3 => {
                let t = s.read(4)?;
                d.base_gw = t.f(0)?;
                d.kh_top = t.f(1)?;
                d.kh_bot = t.f(2)?;
                d.wet_perimeter = t.f(3)?;
            }
            4 => {
                let t = s.read(6)?;
                d.base_gw = t.f(0)?;
                d.kv_top = t.f(1)?;
                d.kv_bot = t.f(2)?;
                d.kh_bot = t.f(3)?;
                d.wet_perimeter = t.f(4)?;
                d.z_interface = t.f(5)?;
            }
            _ => {
                let t = s.read(7)?;
                d.base_gw = t.f(0)?;
                d.kh_top = t.f(1)?;
                d.kv_top = t.f(2)?;
                d.kh_bot = t.f(3)?;
                d.wet_perimeter = t.f(4)?;
                d.z_interface = t.f(5)?;
                d.geo_factor = t.f(6)?;
            }
        }
        d.base_gw = -d.base_gw.abs();
        d.z_interface = -d.z_interface.abs();
    }
    prj.water.bc = bc;

    // ---------------- NodInf (Profile.dat) - needs NS and flags
    let (ns_profile, nobs) = read_profile(&find("profile.dat")?, &mut prj, l_chem, l_temp, l_equil)?;
    let _ = nobs;

    // ---------------- MatIn
    s.skip()?;
    let t = s.read(2)?;
    prj.water.h_tab1 = t.f(0)?;
    prj.water.h_tab_n = t.f(1)?;
    s.skip()?;
    let t = s.read(2)?;
    let i_model = t.i(0)?;
    let i_hyst = t.i(1)?;  

    prj.water.model = SoilModel::from_code(i_model)
        .ok_or_else(|| HydrusError::Unsupported(format!("hydraulic model {}", i_model)))?;

    // Check for external tabular model (Model 10)
    if prj.water.model == SoilModel::Tabular {
        let mater_path = find("mater.in")?;
        prj.water.tabs = Some(read_mater_in(&mater_path, n_mat)?);
    }

    prj.water.hysteresis = Hysteresis::from_code(i_hyst);
    if i_hyst > 0 {
        s.skip()?;
        prj.water.init_kappa = s.read(1)?.i(0)?;
    }
    s.skip()?;
    let npar = match prj.water.model {
        SoilModel::ModifiedVG => 10,
        SoilModel::Durner => 9,
        SoilModel::DualPorosityW => 9,
        SoilModel::DualPorosityH => 11,
		SoilModel::DualPermeability => 11,
        SoilModel::Tabular => 6, // Standard 6 parameters in Selector.in (Qr, Qs, Alfa, n, Ks, l)
        _ => 6,
    };

    prj.water.materials.clear();
    for m in 0..n_mat {
        let mut mat = SoilMaterial {
            name: format!("Material {}", m + 1),
            ..Default::default()
        };
        if i_hyst == 0 {
            let t = s.read(npar)?;
            mat.qr = t.f(0)?;
            mat.qs = t.f(1)?;
            mat.alpha = t.f(2)?;
            mat.n = t.f(3)?;
            mat.ks = t.f(4)?;
            mat.l = t.f(5)?;
            if npar > 6 {
                for k in 0..(npar - 6) {
                    mat.extra[k] = t.f(6 + k)?;
                }
            }
        } else {
            let t = s.read(10)?;
            mat.qr = t.f(0)?;
            mat.qs = t.f(1)?;
            mat.alpha = t.f(2)?;
            mat.n = t.f(3)?;
            mat.ks = t.f(4)?;
            mat.l = t.f(5)?;
            mat.qm = t.f(6)?;
            mat.qs_w = t.f(7)?;
            mat.alpha_w = t.f(8)?;
            mat.ks_w = t.f(9)?;
        }
        prj.water.materials.push(mat);
    }

    // ---------------- TmIn
    s.skip()?;
    s.skip()?;
    let t = s.read(8)?;
    let tm = &mut prj.time;
    tm.dt = t.f(0)?;
    tm.dt_min = t.f(1)?;
    tm.dt_max = t.f(2)?;
    tm.d_mul = t.f(3)?;
    tm.d_mul2 = t.f(4)?;
    tm.it_min = t.i(5)? as usize;
    tm.it_max = t.i(6)? as usize;
    let mpl = t.i(7)? as usize;
    s.skip()?;
    let t = s.read(2)?;
    tm.t_init = t.f(0)?;
    tm.t_max = t.f(1)?;
    if ver >= 3 {
        s.skip()?;
        let t = s.read_opt(3)?;
        tm.print_at_interval = t.b(0).unwrap_or(false);
        tm.print_step = t.i(1).unwrap_or(1).max(1) as usize;
        tm.print_interval = t.f(2).unwrap_or(0.0);
    }
    s.skip()?;
    let t = s.read(mpl)?;
    tm.print_times = (0..mpl).map(|i| t.f(i)).collect::<R<Vec<_>>>()?;

    let top_inf = prj.water.bc.top_time_variable;
    let bot_inf = prj.water.bc.bot_time_variable;
    let atm = prj.water.bc.atmospheric;
    let ns = if l_chem { ns_profile.max(1) } else { 0 };
    let mut rec_has_root_depth = false;
    let mut i_root_in = -1;

    // ---------------- Root growth
    if l_root {
        s.skip()?;
        if ver >= 4 {
            s.skip()?;
            i_root_in = s.read(1)?.i(0)?;
            if i_root_in == 1 {
                s.skip()?;
                let n = s.read(1)?.i(0)? as usize;
                s.skip()?;
                let mut tab = vec![];
                for _ in 0..n {
                    let t = s.read(2)?;
                    tab.push((t.f(0)?, t.f(1)?));
                }
                prj.root.growth = RootGrowthMode::Table(tab);
            }
        }
        if ver < 4 || i_root_in == 2 {
            i_root_in = 2;
            s.skip()?;
            let t = if ver < 4 { s.read(7)? } else { s.read(8)? };
            let irfak = t.i(0)?;
            let (t_min, mut t_med, t_harv) = (t.f(1)?, t.f(2)?, t.f(3)?);
            let (x_min, mut x_med, x_max) = (t.f(4)?, t.f(5)?, t.f(6)?);
            let period = if ver >= 4 { t.f(7)? } else { 1e30 };
            if irfak == 1 {
                t_med = (t_harv + t_min) / 2.0;
                x_med = (x_max + x_min) / 2.0;
            }
            prj.root.growth = RootGrowthMode::Logistic { t_min, t_med, t_harv, x_min, x_med, x_max, period, fit_from_half: irfak == 1 };
        }
        if i_root_in == 0 {
            prj.root.growth = RootGrowthMode::FromAtmosphere;
            rec_has_root_depth = true;
        }
    }

    // ---------------- Atmosph.in
    if top_inf || bot_inf || atm {
        read_atmosphere(&find("atmosph.in")?, &mut prj, l_temp, l_chem, ns, i_root_in, rec_has_root_depth)?;
    }
	
	// ---------------- Meteo.in
    if l_meteo {
		prj.water.bc.atmospheric = true;
        if let Ok(meteo_path) = find("meteo.in") {
            read_meteo(&meteo_path, &mut prj)?;
        }
    }

    // ---------------- Heat
    if l_temp {
        s.skip()?;
        s.skip()?;
        let mut mats = vec![];
        for _ in 0..n_mat {
            let t = s.read(9)?;
            mats.push(HeatMaterial { qn: t.f(0)?, qo: t.f(1)?, disp: t.f(2)?, b1: t.f(3)?, b2: t.f(4)?, b3: t.f(5)?, cn: t.f(6)?, co: t.f(7)?, cw: t.f(8)? });
        }
        s.skip()?;
        let mut ht = HeatSettings { materials: mats, ..Default::default() };
        if ver <= 2 {
            let t = s.read(2)?;
            ht.amplitude = t.f(0)?;
            ht.period = t.f(1)?;
        } else if ver == 3 {
            let t = s.read(4)?;
            ht.amplitude = t.f(0)?;
            ht.period = t.f(1)?;
            ht.campbell = t.i(2)? == 1;
        } else {
            let t = s.read(5)?;
            ht.amplitude = t.f(0)?;
            ht.period = t.f(1)?;
            ht.campbell = t.i(2)? == 1;
            if t.b(4)? {
                s.skip()?;
                s.skip()?;
            }
        }
        s.skip()?;
        let t = s.read(4)?;
        ht.k_top = t.i(0)?;
        ht.t_top = t.f(1)?;
        ht.k_bot = t.i(2)?;
        ht.t_bot = t.f(3)?;
        prj.heat = ht;
    }

    // ---------------- Solute
    if l_chem {
        read_chem(&mut s, &mut prj, ver, n_mat, l_equil)?;
        if prj.solute.i_moist_dep == 2 {
            if let Ok(moist_path) = find("moistdep.in") {
                let _ = crate::moist_dep::read_moist_dep(&moist_path, &mut prj);
            }
        }
    }

    // ---------------- Sink
    if sink_f {
        s.skip()?;
        s.skip()?;
        let ncr = prj.n_solutes();
        let t = if ver <= 2 { s.read(1 + ncr)? } else { s.read(1 + ncr + 1)? };
        let i_mo_sink = t.i(0)?;
        let mut c_root_max = vec![];
        for j in 0..ncr {
            c_root_max.push(t.f(1 + j)?);
        }
        if ver > 2 {
            prj.root.omega_c = t.f(1 + ncr)?;
        }
        prj.root.c_root_max = c_root_max;
        s.skip()?;
        if i_mo_sink == 0 {
            let t = s.read(6)?;
            let (p0, p2h, p2l, p3, r2h, r2l) = (-t.f(0)?.abs(), -t.f(1)?.abs(), -t.f(2)?.abs(), -t.f(3)?.abs(), t.f(4)?, t.f(5)?);
            s.skip()?;
            let t = s.read(n_mat)?;
            let p_optm = (0..n_mat).map(|i| t.f(i).map(|v| -v.abs())).collect::<R<Vec<_>>>()?;
            prj.root.stress = RootStress::Feddes { p0, p2h, p2l, p3, r2h, r2l, p_optm };
        } else {
            let t = s.read(2)?;
            prj.root.stress = RootStress::SShaped { p50: t.f(0)?, exponent: t.f(1)? };
        }
        if l_chem {
            s.skip()?;
            let l_sol_red = s.read(1)?.b(0)?;
            if l_sol_red {
                s.skip()?;
                let add = s.read(1)?.b(0)?;
                s.skip()?;
                if add {
                    let t = s.read(ncr)?;
                    let a_osm = (0..ncr).map(|j| t.f(j)).collect::<R<Vec<_>>>()?;
                    prj.root.solute_stress = Some(SoluteStress { additive: true, a_osm, c50: 0.0, p3c: 0.0, s_shaped: false });
                } else {
                    let t = s.read(2 + ncr + 1)?;
                    let c50 = t.f(0)?;
                    let p3c = t.f(1)?;
                    let a_osm = (0..ncr).map(|j| t.f(2 + j)).collect::<R<Vec<_>>>()?;
                    let ms = t.i(2 + ncr)?;
                    prj.root.solute_stress = Some(SoluteStress { additive: false, a_osm, c50, p3c, s_shaped: ms != 0 });
                }
            }
        }
		if ncr > 1 {
			prj.root.l_act_rsu = false; // Fortran: only for NS = 1
		}
		if prj.root.l_act_rsu && ncr == 1 {
			s.skip()?;
			let t = s.read(5)?; // OmegaS, SPot, rKM, cMin, lOmegaW
			prj.root.omega_act = vec![t.f(0)?];
			prj.root.s_pot = t.f(1)?;
			prj.root.r_km = vec![t.f(2)?];
			prj.root.c_min = vec![t.f(3)?];
			prj.root.l_omega_w = t.b(4)?;
		}
    }

    // observation nodes etc. done in read_profile. Layers:
    let errs = prj.validate();
    if !errs.is_empty() {
        return perr(errs.join("; "));
    }
    Ok(prj)
}

fn read_profile(path: &Path, prj: &mut Project, l_chem: bool, l_temp: bool, l_equil: bool) -> R<(usize, usize)> {
    let mut p = Rd::open(path, "Profile.dat")?;
    let _ver = p.version();
    let nhead = p.read(1)?.i(0)?;
    for _ in 0..nhead {
        p.skip()?;
    }
    let t = p.read(3)?;
    let num_np = t.i(0)? as usize;
    let ns = t.i(2)?.max(0) as usize;
    if num_np > 1001 {
        return perr("Profile.dat: more than 1001 nodes");
    }
    let ns_eff = if l_chem { ns } else { 0 };
    // read node lines, keeping the file order (top -> bottom) and interpolating skipped nodes
    #[derive(Clone)]
    struct Raw {
        x: f64,
        h: f64,
        m: usize,
        l: usize,
        b: f64,
        ax: f64,
        bx: f64,
        dx: f64,
        te: f64,
        c: Vec<f64>,
        s: Vec<f64>,
    }
    let mut nodes: Vec<Option<Raw>> = vec![None; num_np];
    let mut prev: Option<(usize, Raw)> = None;
    let mut cur = 0usize; // next expected 1-based index - 1
    while cur < num_np {
        let ncols = if !l_chem && !l_temp {
            8
        } else if !l_chem {
            9
        } else if l_equil {
            9 + ns_eff
        } else {
            9 + 2 * ns_eff
        };
        let t = p.read(ncols + 1)?;
        let idx = t.i(0)? as usize;
        let mut r = Raw {
            x: t.f(1)?,
            h: t.f(2)?,
            m: t.i(3)? as usize,
            l: t.i(4)? as usize,
            b: t.f(5)?,
            ax: t.f(6)?,
            bx: t.f(7)?,
            dx: t.f(8)?,
            te: 20.0,
            c: vec![0.0; ns_eff],
            s: vec![0.0; ns_eff],
        };
        if l_temp || l_chem {
            r.te = t.f(9)?;
        }
        if l_chem {
            for j in 0..ns_eff {
                r.c[j] = t.f(10 + j)?;
            }
            if !l_equil {
                for j in 0..ns_eff {
                    r.s[j] = t.f(10 + ns_eff + j)?;
                }
            }
        }
        if idx < 1 || idx > num_np || idx - 1 < cur {
            return perr(format!("Profile.dat: node {} out of order", idx));
        }
        if let Some((pi, pr)) = &prev {
            // linear interpolation for nodes cur..idx-1
            let span = (idx - 1 - *pi) as f64;
            for k in (*pi + 1)..(idx - 1) {
                let f = (k - *pi) as f64 / span;
                let lerp = |a: f64, b: f64| a + (b - a) * f;
                nodes[k] = Some(Raw {
                    x: lerp(pr.x, r.x),
                    h: lerp(pr.h, r.h),
                    m: pr.m,
                    l: pr.l,
                    b: lerp(pr.b, r.b),
                    ax: lerp(pr.ax, r.ax),
                    bx: lerp(pr.bx, r.bx),
                    dx: lerp(pr.dx, r.dx),
                    te: lerp(pr.te, r.te),
                    c: (0..ns_eff).map(|j| lerp(pr.c[j], r.c[j])).collect(),
                    s: (0..ns_eff).map(|j| lerp(pr.s[j], r.s[j])).collect(),
                });
            }
        }
        nodes[idx - 1] = Some(r.clone());
        prev = Some((idx - 1, r));
        cur = idx;
    }
    let nobs = p.read(1)?.i(0)?.max(0) as usize;
    let obs = if nobs > 0 {
        let t = p.read(nobs)?;
        (0..nobs).map(|i| t.i(i).map(|v| v as usize)).collect::<R<Vec<_>>>()?
    } else {
        vec![]
    };
    let x0 = nodes[0].as_ref().map(|r| r.x).unwrap_or(0.0);
    prj.profile.nodes = nodes
        .into_iter()
        .map(|r| {
            let r = r.unwrap();
            ProfileNode {
                depth: x0 - r.x,
                h: r.h,
                mat: r.m.max(1),
                layer: r.l.max(1),
                beta: r.b,
                ah: r.ax,
                ak: r.bx,
                ath: r.dx,
                temp: r.te,
                conc: r.c,
                sorb: r.s,
            }
        })
        .collect();
    prj.profile.observation_nodes = obs;
    Ok((ns, nobs))
}

fn read_atmosphere(path: &Path, prj: &mut Project, l_temp: bool, l_chem: bool, ns: usize, i_root_in: i32, has_root_depth: bool) -> R<()> {
    let mut a = Rd::open(path, "Atmosph.in")?;
    let ver = a.version();
    a.skip()?;
    a.skip()?;
    let _max_al = a.read(1)?;
    let at = &mut prj.atmosphere;
    if ver == 4 {
        a.skip()?;
        let t = a.read(5)?;
        at.daily_variation = t.b(0)?;
        at.sinusoidal_precip = t.b(1)?;
        at.lai_partitioning = t.b(2)?;
        at.bc_cycles = t.b(3)?;
        at.interception = t.b(4)?;

        if at.lai_partitioning {
            a.skip()?;
            at.extinction = a.read(1)?.f(0)?;
        }
        
        if at.interception {
            a.skip()?;
            at.interception_a = a.read(1)?.f(0)?;
        }
    }
    
    a.skip()?;
    at.h_crit_s = a.read(1)?.f(0)?;
    a.skip()?;
    at.has_root_depth = has_root_depth || i_root_in == 0;
    at.records.clear();
    loop {
        match a.peek() {
            None => break,
            Some(l) if l.trim_start().to_ascii_lowercase().starts_with("end") => break,
            _ => {}
        }
        let mut ncol = 8;
        if l_temp {
            ncol += 3;
        }
        if l_chem {
            ncol += 2 * ns;
        }
        if at.has_root_depth {
            ncol += 1;
        }
        let t = a.read(ncol)?;
        let mut r = AtmRecord {
            t: t.f(0)?,
            prec: t.f(1)?,
            evap: t.f(2)?,
            transp: t.f(3)?,
            h_crit_a: t.f(4)?,
            r_bot: t.f(5)?,
            h_bot: t.f(6)?,
            h_top: t.f(7)?,
            ..Default::default()
        };
        let mut k = 8;
        if l_temp {
            r.t_top = t.f(k)?;
            r.t_bot = t.f(k + 1)?;
            r.ampl = t.f(k + 2)?;
            k += 3;
        }
        if l_chem {
            for _ in 0..ns {
                r.c_top.push(t.f(k)?);
                r.c_bot.push(t.f(k + 1)?);
                k += 2;
            }
        }
        if at.has_root_depth {
            r.x_root = t.f(k)?;
        }
        at.records.push(r);
    }
    Ok(())
}

fn read_chem(s: &mut Rd, prj: &mut Project, ver: i32, n_mat: usize, l_equil_flag: bool) -> R<()> {
    let _ = l_equil_flag;
    s.skip()?;
    s.skip()?;
    let mut sol = SoluteSettings::default();
    let (ns, l_tdep, l_tort);
    if ver <= 2 {
        let t = s.read(10)?;
        sol.epsi = t.f(0)?;
        sol.upstream_weighting = t.b(1)?;
        sol.artificial_dispersion = t.b(2)?;
        l_tdep = t.b(3)?;
        sol.tol_abs = t.f(4)?;
        sol.tol_rel = t.f(5)?;
        sol.max_iter = t.i(6)? as usize;
        sol.pe_cr = t.f(7)?;
        ns = t.i(8)? as usize;
        l_tort = t.b(9)?;
    } else {
        let t = s.read(12)?;
        sol.epsi = t.f(0)?;
        sol.upstream_weighting = t.b(1)?;
        sol.artificial_dispersion = t.b(2)?;
        l_tdep = t.b(3)?;
        sol.tol_abs = t.f(4)?;
        sol.tol_rel = t.f(5)?;
        sol.max_iter = t.i(6)? as usize;
        sol.pe_cr = t.f(7)?;
        ns = t.i(8)? as usize;
        l_tort = t.b(9)?;
        sol.l_bact = t.i(10).unwrap_or(0) == 1;
        sol.l_filtr = t.b(11).unwrap_or(false);
    }
    sol.tortuosity = l_tort;
    sol.l_tdep = l_tdep;
    let (mut l_moist, mut l_dual_neq, mut l_mass_ini, mut l_eq_init, mut l_var) = (false, false, false, false, false);
    if ver >= 4 {
        s.skip()?;
        let t = s.read_opt(6)?;
        l_moist = t.b(1)?;
        l_dual_neq = t.b(2)?;
        l_mass_ini = t.b(3)?;
        l_eq_init = t.b(4)?;
        l_var = t.b(5)?;
        if t.len() > 6 {
            sol.l_nequil = t.b(6).unwrap_or(false);
        }
        if t.len() > 7 {
            sol.i_moist_dep = t.i(7).unwrap_or(if l_moist { 1 } else { 0 });
        }
        if t.len() > 8 {
            sol.i_conc_type = t.i(8).unwrap_or(1);
        }
    }
    sol.l_moist = l_moist;
    sol.mass_init = l_mass_ini;
    sol.equil_init = l_eq_init;
    sol.tort_model = if l_var { 1 } else { 0 };
	sol.l_dual_neq = l_dual_neq;
    if ns == 0 || ns > 11 {
        return perr("invalid number of solutes");
    }
    s.skip()?;
    sol.materials.clear();
    for _ in 0..n_mat {
        let t = s.read(4)?;
        sol.materials.push(SoluteMaterial { bulk_density: t.f(0)?, disp_l: t.f(1)?, frac: t.f(2)?, th_immobile: t.f(3)? });
    }
    sol.species.clear();
    for j in 0..ns {
        s.skip()?;
        let t = s.read(2)?;
        let mut sp = Species { name: format!("Solute {}", j + 1), diff_w: t.f(0)?, diff_g: t.f(1)?, per_material: vec![], c_top: 0.0, c_bot: 0.0 };
        s.skip()?;
        for _ in 0..n_mat {
            let t = s.read_opt(14)?;
            let mut sp_mat = SpeciesMaterial {
                ks: t.f(0)?,
                nu: t.f(1)?,
                beta: t.f(2)?,
                henry: t.f(3)?,
                mu_w: t.f(4)?,
                mu_s: t.f(5)?,
                mu_g: t.f(6)?,
                gam_w: t.f(7)?,
                gam_s: t.f(8)?,
                gam_g: t.f(9)?,
                mu0_w: t.f(10)?,
                mu0_s: t.f(11)?,
                mu0_g: t.f(12)?,
                omega: t.f(13)?,
                ..Default::default()
            };
            if sol.l_bact {
                sp_mat.s_max2 = sp_mat.mu0_s;
                sp_mat.r_ka2 = sp_mat.mu0_g;
                sp_mat.r_kd2 = sp_mat.omega;
                sp_mat.r_ka1 = sp_mat.ks;
                sp_mat.r_kd1 = sp_mat.nu;
            }
            if t.len() >= 24 {
                sp_mat.s_max1 = t.f(14)?;
                sp_mat.r_ka1 = t.f(15)?;
                sp_mat.r_kd1 = t.f(16)?;
                sp_mat.s_max2 = t.f(17)?;
                sp_mat.r_ka2 = t.f(18)?;
                sp_mat.r_kd2 = t.f(19)?;
                sp_mat.i_psi1 = t.i(20)?;
                sp_mat.i_psi2 = t.i(21)?;
                sp_mat.d_c = t.f(22)?;
                sp_mat.d_p = t.f(23)?;
            }
            sp.per_material.push(sp_mat);
        }
        sol.species.push(sp);
    }

    sol.t_dep.clear();
	sol.w_dep.clear();
   if sol.l_tdep {
        for jj in 0..ns {
            if jj == 0 {
                s.skip()?; // Header: Temperature dependence
            }
            s.skip()?; // Header: Dif.w. Dif.g.
            let t_diff = s.read(2)?;
            s.skip()?; // Header: Ks Nu Beta Henry ...
            let t_par = s.read(14)?;
            sol.t_dep.push(SpeciesTDep {
                diff_w: t_diff.f(0)?,
                diff_g: t_diff.f(1)?,
                ks: t_par.f(0)?,
                nu: t_par.f(1)?,
                beta: t_par.f(2)?,
                henry: t_par.f(3)?,
                mu_w: t_par.f(4)?,
                mu_s: t_par.f(5)?,
                mu_g: t_par.f(6)?,
                gam_w: t_par.f(7)?,
                gam_s: t_par.f(8)?,
                gam_g: t_par.f(9)?,
                mu0_w: t_par.f(10)?,
                mu0_s: t_par.f(11)?,
                mu0_g: t_par.f(12)?,
                omega: t_par.f(13)?,
            });
        }
    }
    if sol.l_moist {
        for jj in 0..ns {
            if jj == 0 {
                s.skip()?; // Header: Water content dependence
            }
            s.skip()?; // Header
            let _n_par2 = s.read(1)?.i(0)?;
            s.skip()?; // Header: Exponents B
            let t_exp = s.read(9)?;
            s.skip()?; // Header: Reference h
            let t_href = s.read(9)?;
            let mut wdep = SpeciesWDep::default();
            for k in 0..9 {
                wdep.exp_b[k] = t_exp.f(k)?;
                wdep.h_ref[k] = t_href.f(k)?;
            }
            sol.w_dep.push(wdep);
        }
    }
    s.skip()?;
    let t = s.read(2 + 2 * ns)?;
    sol.k_top = t.i(0)?;
    for j in 0..ns {
        sol.species[j].c_top = t.f(1 + j)?;
    }
    sol.k_bot = t.i(1 + ns)?;
    for j in 0..ns {
        sol.species[j].c_bot = t.f(2 + ns + j)?;
    }
    if sol.k_top == -2 {
        s.skip()?;
        let t = s.read(2)?;
        sol.d_surf = t.f(0)?;
        sol.c_atm = t.f(1)?;
    }
    s.skip()?;
    sol.t_pulse = s.read(1)?.f(0)?;
    prj.solute = sol;
    Ok(())
}
fn read_meteo(path: &Path, prj: &mut Project) -> R<()> {
    let mut m = Rd::open(path, "Meteo.in")?;
    let _ver = m.version();

    m.skip()?;
    m.skip()?;
    let t = m.read(3)?;
    let num_records = t.i(0)? as usize;
    let i_radiation = t.i(1)?;
    let hargreaves = t.b(2)?;

    // Line 4 & 5: Flags (lEnBal, lDaily, ...)
    m.skip()?;
    let flags = m.read(10)?;
	let l_en_bal = flags.b(0)?;
    let l_daily = flags.b(1)?;

    // Latitude, Altitude
    m.skip()?;
    let t = m.read(2)?;
    let latitude = t.f(0)?;
    let altitude = t.f(1)?;

    // ShortWaveRadA, ShortWaveRadB
    m.skip()?;
    let t = m.read(2)?;
    let short_wave_a = t.f(0)?;
    let short_wave_b = t.f(1)?;

    // LongWaveRadA, LongWaveRadB
    m.skip()?;
    let t = m.read(2)?;
    let long_wave_a = t.f(0)?;
    let long_wave_b = t.f(1)?;

    // LongWaveRadA1, LongWaveRadB1
    m.skip()?;
    let t = m.read(2)?;
    let long_wave_a1 = t.f(0)?;
    let long_wave_b1 = t.f(1)?;

    // WindHeight, TempHeight
    m.skip()?;
    let t = m.read(2)?;
    let wind_height = t.f(0)?;
    let temp_height = t.f(1)?;

    // iCrop, SunShine, RelativeHum
    m.skip()?;
    let t = m.read(3)?;
    let i_crop = t.i(0)?;
    let i_sun_sh = t.i(1)?;
    let i_rel_hum = t.i(2)?;

    // Albedo
    m.skip()?;
    let albedo = m.read(1)?.f(0)?;

    // Daily values header rows (2 lines)
    m.skip()?;
    m.skip()?;
	if let Some(next_line) = m.peek() {
        if next_line.trim().starts_with('[') {
            m.skip()?; // Skip units row if present
        }
    }

    let x_conv = prj.units.x_conv(); // For crop height unit conversion (to cm)
    let mut records = Vec::with_capacity(num_records);
    loop {
        match m.peek() {
            None => break,
            Some(l) if l.trim_start().to_ascii_lowercase().starts_with("end") => break,
            _ => {}
        }

        // iCrop = 3 (daily) has 11 columns; otherwise standard 7 columns
        let ncols = if i_crop == 3 { 11 } else { 7 };
        let t = match m.read_opt(ncols) {
            Ok(toks) => toks,
            Err(_) => break,
        };
        if t.len() < 7 {
            break;
        }

        let mut r = MeteoRecord {
            t: t.f(0)?,
            rad: t.f(1)?,
            t_max: t.f(2)?,
            t_min: t.f(3)?,
            rh_mean: t.f(4)?,
            wind_kmd: t.f(5)?,
            sun_hours: t.f(6)?,
            crop_height: None,
            albedo: None,
            lai: None,
            x_root: None,
        };

        if i_crop == 3 && t.len() >= 11 {
            // Conversion to cm as in TIME.FOR: CropHeight * 100. / xConv
            r.crop_height = Some(t.f(7)? * 100.0 / x_conv);
            r.albedo = Some(t.f(8)?);
            r.lai = Some(t.f(9)?);
            r.x_root = Some(t.f(10)?);
        }

        records.push(r);

        if records.len() >= num_records {
            break;
        }
    }

    let settings = MeteoSettings {
		latitude,
        altitude,
        short_wave_a,
        short_wave_b,
        long_wave_a,
        long_wave_b,
        long_wave_a1,
        long_wave_b1,
        wind_height,
        temp_height,
        i_radiation,
        i_sun_sh,
        i_rel_hum,
        hargreaves,
        l_en_bal,
        l_daily,
        i_crop,
        albedo,
        records,
        ..Default::default()
    };

    prj.atmosphere.meteo = Some(settings);
    Ok(())
}

/// Read an external user-defined material table file (Mater.in, Model 10).
pub fn read_mater_in(path: &Path, n_mat: usize) -> R<Vec<MatTable>> {
    let mut rd = Rd::open(path, "Mater.in")?;
    let mut tabs = Vec::with_capacity(n_mat);

    for m in 0..n_mat {
        let header = rd.line()?;
        if header.trim().is_empty() {
            return perr(format!("Mater.in: empty header for material {}", m + 1));
        }

        let n_points = rd.read(1)?.i(0)? as usize;
        if n_points < 2 {
            return perr(format!(
                "Mater.in: material {} requires at least 2 points, found {}",
                m + 1,
                n_points
            ));
        }

        let mut t = MatTable {
            h: Vec::with_capacity(n_points),
            the: Vec::with_capacity(n_points),
            con: Vec::with_capacity(n_points),
            cap: Vec::with_capacity(n_points),
        };

        for pt in 0..n_points {
            let row = rd.read(4).map_err(|e| {
                HydrusError::Parse(format!(
                    "Mater.in: failed reading row {} for material {}: {}",
                    pt + 1,
                    m + 1,
                    e
                ))
            })?;
            t.h.push(row.f(0)?);
            t.the.push(row.f(1)?);
            t.con.push(row.f(2)?);
            t.cap.push(row.f(3)?);
        }

        tabs.push(t);
    }
    Ok(tabs)
}

// ---------------------------------------------------------------------------
// Legacy Output File Readers (Nod_Inf.out)
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Default)]
pub struct LegacyProfileNode {
    pub node: usize,
    pub depth: f64,
    pub head: f64,
    pub moisture: f64,
    pub k: f64,
    pub c: f64,
    pub flux: f64,
    pub sink: f64,
    pub kappa: i32,
    pub temp: f64,
    pub conc: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct LegacyProfileTimeBlock {
    pub time: f64,
    pub nodes: Vec<LegacyProfileNode>,
}

/// Tokenizes a line from a legacy HYDRUS Fortran output file, safely separating 
/// contiguous scientific-notation tokens (e.g. "0.00-7.2970E+02" or "1.00E+01-2.00E+01").
pub fn sanitize_fortran_line(line: &str) -> String {
    let mut out = String::with_capacity(line.len() + 16);
    let bytes = line.as_bytes();
    for i in 0..bytes.len() {
        let b = bytes[i];
        if b == b'-' && i > 0 {
            let prev = bytes[i - 1];
            if prev != b' ' && prev != b'e' && prev != b'E' && prev != b'd' && prev != b'D' && prev != b',' {
                out.push(' ');
            }
        }
        out.push(b as char);
    }
    out
}

/// Reads a legacy Nod_Inf.out file, parsing all profile time blocks with complete column fidelity.
pub fn read_nod_inf(path: &Path) -> R<Vec<LegacyProfileTimeBlock>> {
    let file = File::open(path).map_err(|e| HydrusError::Io(format!("{}: {}", path.display(), e)))?;
    let reader = BufReader::new(file);

    let mut blocks: Vec<LegacyProfileTimeBlock> = Vec::new();
    let mut cur_time: Option<f64> = None;
    let mut cur_nodes: Vec<LegacyProfileNode> = Vec::new();

    let mut head_col = 2;
    let mut depth_col = 1;
    let mut node_col = 0;
    let mut th_col = 3;
    let mut k_col = 4;
    let mut c_col = 5;
    let mut flux_col = 6;
    let mut sink_col = 7;
    let mut kappa_col = 8;
    let mut temp_col = 10;

    for line_res in reader.lines() {
        let raw_line = line_res.map_err(|e| HydrusError::Io(e.to_string()))?;
        let trimmed = raw_line.trim();

        if trimmed.is_empty() || trimmed.starts_with("====") || trimmed.starts_with("----") {
            continue;
        }

        let lower = trimmed.to_ascii_lowercase();

        // New time block header
        if lower.starts_with("time:") {
            if let Some(t) = cur_time {
                if !cur_nodes.is_empty() {
                    blocks.push(LegacyProfileTimeBlock {
                        time: t,
                        nodes: std::mem::take(&mut cur_nodes),
                    });
                }
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            cur_time = parts.get(1).and_then(|s| s.parse::<f64>().ok());
            continue;
        }

        if lower.starts_with("end") {
            break;
        }

        let upper = trimmed.to_ascii_uppercase();
        // Dynamically detect column header layouts (even if preceded by comment asterisks)
        if upper.contains("NODE") && (upper.contains("DEPTH") || upper.contains("HEAD")) {
            let clean_line = trimmed.replace('*', " ");
            let headers: Vec<&str> = clean_line.split_whitespace().collect();
            for (idx, &h) in headers.iter().enumerate() {
                let clean = h.trim_matches(|c: char| !c.is_alphanumeric());
                if clean.eq_ignore_ascii_case("NODE") {
                    node_col = idx;
                } else if clean.eq_ignore_ascii_case("DEPTH") {
                    depth_col = idx;
                } else if clean.eq_ignore_ascii_case("HEAD") || clean.eq_ignore_ascii_case("H") {
                    head_col = idx;
                } else if clean.eq_ignore_ascii_case("MOISTURE") || clean.eq_ignore_ascii_case("THETA") {
                    th_col = idx;
                } else if clean.eq_ignore_ascii_case("K") {
                    k_col = idx;
                } else if clean.eq_ignore_ascii_case("C") {
                    c_col = idx;
                } else if clean.eq_ignore_ascii_case("FLUX") {
                    flux_col = idx;
                } else if clean.eq_ignore_ascii_case("SINK") {
                    sink_col = idx;
                } else if clean.eq_ignore_ascii_case("KAPPA") {
                    kappa_col = idx;
                } else if clean.eq_ignore_ascii_case("TEMP") {
                    temp_col = idx;
                }
            }
            continue;
        }

        if trimmed.starts_with('[') || trimmed.contains("[L]") || trimmed.contains("[T]") {
            continue;
        }

        // Clean out observation asterisks that prefix or suffix node numbers (e.g. "* 22" or "22*")
        let clean_row = trimmed.replace('*', " ");
        let sanitized = sanitize_fortran_line(&clean_row);
        let vals: Vec<f64> = sanitized
            .split_whitespace()
            .filter_map(|s| s.replace(['d', 'D'], "e").parse::<f64>().ok())
            .collect();

        if vals.len() >= 6 && vals.len() > head_col {
            let node = vals.get(node_col).copied().unwrap_or(0.0).round() as usize;
            if node == 0 {
                continue;
            }

            // Prevent trailing observation/summary tables from re-adding or overriding nodes in this block
            if cur_nodes.iter().any(|n| n.node == node) {
                continue;
            }

            let depth = vals.get(depth_col).copied().unwrap_or(0.0);
            let head = vals[head_col];
            let moisture = vals.get(th_col).copied().unwrap_or(0.0);
            let k = vals.get(k_col).copied().unwrap_or(0.0);
            let c = vals.get(c_col).copied().unwrap_or(0.0);
            let flux = vals.get(flux_col).copied().unwrap_or(0.0);
            let sink = vals.get(sink_col).copied().unwrap_or(0.0);
            let kappa = vals.get(kappa_col).copied().unwrap_or(-1.0).round() as i32;
            let temp = vals.get(temp_col).copied().unwrap_or(0.0);

            let conc = if vals.len() > 11 {
                vals[11..].to_vec()
            } else {
                Vec::new()
            };

            cur_nodes.push(LegacyProfileNode {
                node,
                depth,
                head,
                moisture,
                k,
                c,
                flux,
                sink,
                kappa,
                temp,
                conc,
            });
        }
    }

    if let Some(t) = cur_time {
        if !cur_nodes.is_empty() {
            blocks.push(LegacyProfileTimeBlock {
                time: t,
                nodes: cur_nodes,
            });
        }
    }

    Ok(blocks)
}