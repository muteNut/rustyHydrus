	//! Project data model. All quantities are expressed in the user's units
//! (length/time/mass units chosen in [`Units`]), exactly as in HYDRUS-1D.
//!
//! Node ordering in the *project* is top -> bottom (like Profile.dat and
//! the GUI). The solver internally reorders to bottom -> top.

use crate::material::MatTable;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum LengthUnit {
    Mm,
    Cm,
    M,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum TimeUnit {
    Sec,
    Min,
    Hours,
    Days,
    Years,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Units {
    pub length: LengthUnit,
    pub time: TimeUnit,
    pub mass: String,
}
impl Default for Units {
    fn default() -> Self {
        Units { length: LengthUnit::Cm, time: TimeUnit::Days, mass: "mmol".into() }
    }
}
impl Units {
    /// Conversion from metres to the length unit.
    pub fn x_conv(&self) -> f64 {
        match self.length {
            LengthUnit::M => 1.0,
            LengthUnit::Cm => 100.0,
            LengthUnit::Mm => 1000.0,
        }
    }
    /// Conversion from seconds to the time unit.
    pub fn t_conv(&self) -> f64 {
        match self.time {
            TimeUnit::Sec => 1.0,
            TimeUnit::Min => 1.0 / 60.0,
            TimeUnit::Hours => 1.0 / 3600.0,
            TimeUnit::Days => 1.0 / 86400.0,
            TimeUnit::Years => 1.0 / (86400.0 * 365.0),
        }
    }
    pub fn length_str(&self) -> &'static str {
        match self.length {
            LengthUnit::Mm => "mm",
            LengthUnit::Cm => "cm",
            LengthUnit::M => "m",
        }
    }
    pub fn time_str(&self) -> &'static str {
        match self.time {
            TimeUnit::Sec => "s",
            TimeUnit::Min => "min",
            TimeUnit::Hours => "h",
            TimeUnit::Days => "d",
            TimeUnit::Years => "yr",
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Processes {
    pub water_flow: bool,
    pub solute: bool,
    pub heat: bool,
    pub root_water_uptake: bool,
    pub root_growth: bool,
    /// Adsorption is an equilibrium process (lEquil). Recomputed from solute parameters.
    pub equilibrium_adsorption: bool,
    /// Report only at print times (lShort).
    pub short_output: bool,
    pub vapor: bool,
}

impl Default for Processes {
    fn default() -> Self {
        Processes {
            water_flow: true,
            solute: false,
            heat: false,
            root_water_uptake: false,
            root_growth: false,
            equilibrium_adsorption: true,
            short_output: false,
            vapor: false,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SoilModel {
    /// van Genuchten - Mualem (iModel 0)
    VanGenuchten,
    /// Modified van Genuchten, Vogel & Cislerova (1)
    ModifiedVG,
    /// Brooks & Corey (2)
    BrooksCorey,
    /// van Genuchten with air-entry value of 2 cm (3)
    VGAirEntry,
    /// Kosugi lognormal (4)
    Kosugi,
    /// Durner dual-porosity function (5)
    Durner,
    /// Dual-porosity model with water-content driven exchange (6)
    DualPorosityW,
    /// Dual-porosity model with pressure-head driven exchange (7)
    DualPorosityH,
	DualPermeability,
    /// External user tabular retention and conductivity curves (Mater.in) (10)
    Tabular,
}

impl SoilModel {
    pub fn code(self) -> i32 {
        match self {
            SoilModel::VanGenuchten => 0,
            SoilModel::ModifiedVG => 1,
            SoilModel::BrooksCorey => 2,
            SoilModel::VGAirEntry => 3,
            SoilModel::Kosugi => 4,
            SoilModel::Durner => 5,
            SoilModel::DualPorosityW => 6,
            SoilModel::DualPorosityH => 7,
			SoilModel::DualPermeability => 8,
            SoilModel::Tabular => 10,
        }
    }

    pub fn from_code(c: i32) -> Option<Self> {
        Some(match c {
            0 => SoilModel::VanGenuchten,
            1 => SoilModel::ModifiedVG,
            2 => SoilModel::BrooksCorey,
            3 => SoilModel::VGAirEntry,
            4 => SoilModel::Kosugi,
            5 => SoilModel::Durner,
            6 => SoilModel::DualPorosityW,
            7 => SoilModel::DualPorosityH,
			8 => SoilModel::DualPermeability,
            10 => SoilModel::Tabular,
            _ => return None,
        })
    }

    pub fn name(self) -> &'static str {
        match self {
            SoilModel::VanGenuchten => "van Genuchten - Mualem",
            SoilModel::ModifiedVG => "Modified van Genuchten (Vogel & Cislerova)",
            SoilModel::BrooksCorey => "Brooks & Corey",
            SoilModel::VGAirEntry => "van Genuchten with air-entry value (2 cm)",
            SoilModel::Kosugi => "Kosugi log-normal",
            SoilModel::Durner => "Durner dual-porosity",
            SoilModel::DualPorosityW => "Dual-porosity (water content driven)",
            SoilModel::DualPorosityH => "Dual-porosity (pressure head driven)",
            SoilModel::DualPermeability => "Dual-permeability",
            SoilModel::Tabular => "Tabular (Mater.in)",
        }
    }

    /// Labels of the extra parameters stored in `SoilMaterial::extra`.
    pub fn extra_labels(self) -> &'static [&'static str] {
        match self {
            SoilModel::ModifiedVG => &["θm", "θa", "θk", "Kk"],
            SoilModel::Durner => &["w2", "α2", "n2"],
            SoilModel::DualPorosityW => &["θr,im", "θs,im", "ω"],
            SoilModel::DualPorosityH => &["θr,im", "θs,im", "α_im", "n_im", "ω"],
            SoilModel::DualPermeability => &["θr,m", "θs,m", "α_m", "n_m", "Ks,m"],
            _ => &[],
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Hysteresis {
    None,
    /// Kool & Parker style scanning curves, retention curve only (iHyst 1)
    Retention,
    /// ... in retention curve and hydraulic conductivity (iHyst 2)
    RetentionAndConductivity,
    /// Lenhard et al. air-entrapment model (iHyst 3)
    Lenhard,
}
impl Hysteresis {
    pub fn code(self) -> i32 {
        match self {
            Hysteresis::None => 0,
            Hysteresis::Retention => 1,
            Hysteresis::RetentionAndConductivity => 2,
            Hysteresis::Lenhard => 3,
        }
    }
    pub fn from_code(c: i32) -> Self {
        match c {
            1 => Hysteresis::Retention,
            2 => Hysteresis::RetentionAndConductivity,
            3 => Hysteresis::Lenhard,
            _ => Hysteresis::None,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SoilMaterial {
    pub name: String,
    pub qr: f64,
    pub qs: f64,
    pub alpha: f64,
    pub n: f64,
    pub ks: f64,
    pub l: f64,
    /// Modified VG: [Qm, Qa, Qk, Kk]; Durner: [w2, alpha2, n2]; Dual-porosity: [thr_im, ths_im, omega, alpha_im, ks_im]
    pub extra: [f64; 5],
    pub qm: f64,
    pub qs_w: f64,
    pub alpha_w: f64,
    pub ks_w: f64,
}

impl Default for SoilMaterial {
    fn default() -> Self {
        SoilMaterial {
            name: "Loam".into(),
            qr: 0.078,
            qs: 0.43,
            alpha: 0.036,
            n: 1.56,
            ks: 24.96,
            l: 0.5,
            extra: [0.43, 0.078, 0.43, 24.96, 0.0],
            qm: 0.43,
            qs_w: 0.43,
            alpha_w: 0.072,
            ks_w: 24.96,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DrainSettings {
    /// 1..5, see FqDrain (SWAP-based drainage formulae)
    pub position: i32,
    pub z_bot: f64,
    pub spacing: f64,
    pub entrance_resistance: f64,
    pub base_gw: f64,
    pub kh_top: f64,
    pub kh_bot: f64,
    pub kv_top: f64,
    pub kv_bot: f64,
    pub wet_perimeter: f64,
    pub z_interface: f64,
    pub geo_factor: f64,
}
impl Default for DrainSettings {
    fn default() -> Self {
        DrainSettings {
            position: 1,
            z_bot: -100.0,
            spacing: 1000.0,
            entrance_resistance: 0.0,
            base_gw: -150.0,
            kh_top: 10.0,
            kh_bot: 10.0,
            kv_top: 10.0,
            kv_bot: 10.0,
            wet_perimeter: 10.0,
            z_interface: -50.0,
            geo_factor: 1.0,
        }
    }
}

/// Water-flow boundary conditions; the encoding follows HYDRUS' Selector.in.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WaterBc {
    /// Atmospheric boundary condition at the surface (AtmBC).
    pub atmospheric: bool,
    /// Time-variable top boundary condition read from the atmospheric file (TopInF).
    pub top_time_variable: bool,
    /// Allow a surface layer / ponding (WLayer).
    pub surface_layer: bool,
    /// 1 = constant pressure head, -1 = constant flux, 0 = variable pressure head / flux.
    pub kod_top: i32,
    pub bot_time_variable: bool,
    /// Boundary flux as a function of the groundwater level (qGWLF).
    pub gwl_flux: bool,
    pub free_drainage: bool,
    pub seepage_face: bool,
    /// 1 = constant head, -1 = constant flux, 0 = variable.
    pub kod_bot: i32,
    pub h_seep: f64,
    pub drains: Option<DrainSettings>,
    /// Constant boundary fluxes (used when not time-variable).
    pub r_top: f64,
    pub r_bot: f64,
    pub r_root: f64,
    pub gwl0l: f64,
    pub aqh: f64,
    pub bqh: f64,
}
impl Default for WaterBc {
    fn default() -> Self {
        WaterBc {
            atmospheric: false,
            top_time_variable: false,
            surface_layer: false,
            kod_top: -1,
            bot_time_variable: false,
            gwl_flux: false,
            free_drainage: true,
            seepage_face: false,
            kod_bot: -1,
            h_seep: 0.0,
            drains: None,
            r_top: 0.0,
            r_bot: 0.0,
            r_root: 0.0,
            gwl0l: 0.0,
            aqh: 0.0,
            bqh: 0.0,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WaterFlow {
    pub max_iter: usize,
    pub tol_th: f64,
    pub tol_h: f64,
    pub model: SoilModel,
    pub hysteresis: Hysteresis,
    /// Initial branch: -1 drying, +1 wetting (IKappa)
    pub init_kappa: i32,
    pub h_tab1: f64,
    pub h_tab_n: f64,
    pub materials: Vec<SoilMaterial>,
    /// Initial condition given as water content instead of pressure head (lInitW)
    pub init_in_water_content: bool,
    pub bc: WaterBc,
	pub l_w_dep: bool,
	pub tabs: Option<Vec<MatTable>>,
}
impl Default for WaterFlow {
    fn default() -> Self {
        WaterFlow {
            max_iter: 10,
            tol_th: 0.001,
            tol_h: 1.0,
            model: SoilModel::VanGenuchten,
            hysteresis: Hysteresis::None,
            init_kappa: -1,
            h_tab1: 1e-6,
            h_tab_n: 1e5,
            materials: vec![SoilMaterial::default()],
            init_in_water_content: false,
            bc: WaterBc::default(),
			l_w_dep: false,
			tabs: None,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TimeSettings {
    pub dt: f64,
    pub dt_min: f64,
    pub dt_max: f64,
    pub d_mul: f64,
    pub d_mul2: f64,
    pub it_min: usize,
    pub it_max: usize,
    pub t_init: f64,
    pub t_max: f64,
    pub print_times: Vec<f64>,
    /// Print every n-th time step (nPrStep)
    pub print_step: usize,
    /// Print at regular intervals (lPrintD / tPrintInt)
    pub print_at_interval: bool,
    pub print_interval: f64,
}
impl Default for TimeSettings {
    fn default() -> Self {
        TimeSettings {
            dt: 0.001,
            dt_min: 1e-5,
            dt_max: 5.0,
            d_mul: 1.3,
            d_mul2: 0.7,
            it_min: 3,
            it_max: 7,
            t_init: 0.0,
            t_max: 1.0,
            print_times: vec![0.5, 1.0],
            print_step: 1,
            print_at_interval: false,
            print_interval: 1.0,
        }
    }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AtmRecord {
    pub t: f64,
    pub prec: f64,
    /// Potential surface evaporation (or the KodTop code for variable BCs)
    pub evap: f64,
    /// Potential transpiration
    pub transp: f64,
    pub h_crit_a: f64,
    pub r_bot: f64,
    pub h_bot: f64,
    pub h_top: f64,
    pub t_top: f64,
    pub t_bot: f64,
    pub ampl: f64,
    pub c_top: Vec<f64>,
    pub c_bot: Vec<f64>,
    pub x_root: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Atmosphere {
    pub records: Vec<AtmRecord>,
    pub h_crit_s: f64,
    pub daily_variation: bool,
    pub sinusoidal_precip: bool,
    pub lai_partitioning: bool,
    pub extinction: f64,
    pub interception: bool,
    pub interception_a: f64,
    pub snow: bool,
    pub snow_mf: f64,
    pub has_root_depth: bool,
    pub bc_cycles: bool,
    pub meteo: Option<MeteoSettings>,
}
impl Default for Atmosphere {
    fn default() -> Self {
        Atmosphere {
            records: vec![],
            h_crit_s: 0.0,
            daily_variation: false,
            sinusoidal_precip: false,
            lai_partitioning: false,
            extinction: 0.463,
            interception: false,
            interception_a: 0.25,
            snow: false,
            snow_mf: 0.43,
            has_root_depth: false,
            bc_cycles: false,
            meteo: None,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum RootStress {
    /// Feddes et al. (1978) piecewise-linear
    Feddes { p0: f64, p2h: f64, p2l: f64, p3: f64, r2h: f64, r2l: f64, p_optm: Vec<f64> },
    /// van Genuchten (1987) S-shaped
    SShaped { p50: f64, exponent: f64 },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SoluteStress {
    /// Additive (true) or multiplicative (false) osmotic stress
    pub additive: bool,
    pub a_osm: Vec<f64>,
    pub c50: f64,
    pub p3c: f64,
    /// S-shaped (true) vs threshold-slope (false) for multiplicative
    pub s_shaped: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum RootGrowthMode {
    /// Rooting depth read from the atmospheric time series (iRootIn = 0)
    FromAtmosphere,
    /// Table of (time, root depth) (iRootIn = 1)
    Table(Vec<(f64, f64)>),
    /// Verhulst-Pearl logistic growth (iRootIn = 2)
    Logistic { t_min: f64, t_med: f64, t_harv: f64, x_min: f64, x_med: f64, x_max: f64, period: f64, fit_from_half: bool },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RootUptake {
    pub stress: RootStress,
    pub omega_c: f64,
    pub solute_stress: Option<SoluteStress>,
    pub growth: RootGrowthMode,
    pub c_root_max: Vec<f64>,
}
impl Default for RootUptake {
    fn default() -> Self {
        RootUptake {
            stress: RootStress::Feddes {
                p0: -10.0,
                p2h: -200.0,
                p2l: -800.0,
                p3: -8000.0,
                r2h: 0.5,
                r2l: 0.1,
                p_optm: vec![-25.0],
            },
            omega_c: 1.0,
            solute_stress: None,
            growth: RootGrowthMode::FromAtmosphere,
            c_root_max: vec![0.0],
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HeatMaterial {
    /// Volume fractions of solid (Qn) and organic (Qo) phases
    pub qn: f64,
    pub qo: f64,
    /// thermal dispersivity
    pub disp: f64,
    /// thermal-conductivity coefficients b1,b2,b3
    pub b1: f64,
    pub b2: f64,
    pub b3: f64,
    /// volumetric heat capacities of solid, organic, water
    pub cn: f64,
    pub co: f64,
    pub cw: f64,
}
impl Default for HeatMaterial {
    fn default() -> Self {
        // SI-derived defaults in J, m, s -> converted by the user; here for cm/day:
        HeatMaterial { qn: 0.6, qo: 0.0, disp: 5.0, b1: 1.42e17, b2: 7.5e17, b3: 2.99e18, cn: 1.43e14, co: 1.87e14, cw: 3.12e14 }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HeatSettings {
    pub materials: Vec<HeatMaterial>,
    pub amplitude: f64,
    pub period: f64,
    /// Use the Chung & Horton / Campbell thermal conductivity function
    pub campbell: bool,
    /// >0 Dirichlet temperature; <0 heat flux
    pub k_top: i32,
    pub t_top: f64,
    pub k_bot: i32,
    pub t_bot: f64,
}
impl Default for HeatSettings {
    fn default() -> Self {
        HeatSettings { materials: vec![HeatMaterial::default()], amplitude: 0.0, period: 1.0, campbell: false, k_top: 1, t_top: 20.0, k_bot: 0, t_bot: 20.0 }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SoluteMaterial {
    pub bulk_density: f64,
    pub disp_l: f64,
    /// fraction of sorption sites in contact with mobile water (equilibrium sites for two-site model)
    pub frac: f64,
    /// immobile water content (mobile-immobile model)
    pub th_immobile: f64,
}
impl Default for SoluteMaterial {
    fn default() -> Self {
        SoluteMaterial { bulk_density: 1.5, disp_l: 1.0, frac: 1.0, th_immobile: 0.0 }
    }
}

/// Per (species, material) reaction / sorption parameters (14 values in Selector.in).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SpeciesMaterial {
    pub ks: f64,
    pub nu: f64,
    pub beta: f64,
    pub henry: f64,
    pub mu_w: f64,
    pub mu_s: f64,
    pub mu_g: f64,
    pub gam_w: f64,
    pub gam_s: f64,
    pub gam_g: f64,
    pub mu0_w: f64,
    pub mu0_s: f64,
    pub mu0_g: f64,
    pub omega: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Species {
    pub name: String,
    pub diff_w: f64,
    pub diff_g: f64,
    pub per_material: Vec<SpeciesMaterial>,
    /// Top BC: >0 concentration BC (Dirichlet), <=0 flux BC (Cauchy)
    pub c_top: f64,
    pub c_bot: f64,
}
impl Default for Species {
    fn default() -> Self {
        Species {
            name: "Tracer".into(),
            diff_w: 0.0,
            diff_g: 0.0,
            per_material: vec![SpeciesMaterial { beta: 1.0, ..Default::default() }],
            c_top: 0.0,
            c_bot: 0.0,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SoluteSettings {
    /// Time weighting: 0 explicit, 0.5 Crank-Nicolson, 1 implicit
    pub epsi: f64,
    pub upstream_weighting: bool,
    pub artificial_dispersion: bool,
    pub tortuosity: bool,
    /// tortuosity model: 0 Millington & Quirk, 1 Moldrup
    pub tort_model: i32,
    pub tol_abs: f64,
    pub tol_rel: f64,
    pub max_iter: usize,
    pub pe_cr: f64,
    pub materials: Vec<SoluteMaterial>,
    pub species: Vec<Species>,
    /// Top BC: 1 Dirichlet, -1 Cauchy (concentration flux), -2 volatilization
    pub k_top: i32,
    pub k_bot: i32,
    pub t_pulse: f64,
    /// Non-equilibrium initial condition in equilibrium with the liquid phase
    pub equil_init: bool,
    /// Initial condition given as total mass
    pub mass_init: bool,
    pub d_surf: f64,
    pub c_atm: f64,
}
impl Default for SoluteSettings {
    fn default() -> Self {
        SoluteSettings {
            epsi: 0.5,
            upstream_weighting: false,
            artificial_dispersion: true,
            tortuosity: true,
            tort_model: 0,
            tol_abs: 0.0,
            tol_rel: 0.0,
            max_iter: 1,
            pe_cr: 2.0,
            materials: vec![SoluteMaterial::default()],
            species: vec![Species::default()],
            k_top: -1,
            k_bot: 0,
            t_pulse: 1.0,
            equil_init: false,
            mass_init: false,
            d_surf: 0.0,
            c_atm: 0.0,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ProfileNode {
    /// Depth below the surface (positive downward), same length unit as the project.
    pub depth: f64,
    /// Initial pressure head (or water content if `init_in_water_content`)
    pub h: f64,
    /// Material number (1-based)
    pub mat: usize,
    /// Layer / sub-region number (1-based)
    pub layer: usize,
    /// Root water uptake distribution (b(x))
    pub beta: f64,
    /// Scaling factors: pressure head, conductivity, water content
    pub ah: f64,
    pub ak: f64,
    pub ath: f64,
    pub temp: f64,
    pub conc: Vec<f64>,
    pub sorb: Vec<f64>,
}
impl ProfileNode {
    pub fn new(depth: f64, h: f64, mat: usize) -> Self {
        ProfileNode { depth, h, mat, layer: 1, beta: 0.0, ah: 1.0, ak: 1.0, ath: 1.0, temp: 20.0, conc: vec![], sorb: vec![] }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Profile {
    pub nodes: Vec<ProfileNode>,
    /// 1-based observation node numbers counted from the top
    pub observation_nodes: Vec<usize>,
}
impl Default for Profile {
    fn default() -> Self {
        let n = 101;
        let nodes = (0..n).map(|i| ProfileNode::new(i as f64, -100.0, 1)).collect();
        Profile { nodes, observation_nodes: vec![] }
    }
}
impl Profile {
    pub fn depth(&self) -> f64 {
        self.nodes.last().map(|n| n.depth).unwrap_or(0.0) - self.nodes.first().map(|n| n.depth).unwrap_or(0.0)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Project {
    pub title: String,
    pub units: Units,
    pub processes: Processes,
    /// cosine of the angle between flow direction and vertical axis (1 = vertical)
    pub cos_alpha: f64,
    pub water: WaterFlow,
    pub time: TimeSettings,
    pub atmosphere: Atmosphere,
    pub root: RootUptake,
    pub heat: HeatSettings,
    pub solute: SoluteSettings,
    pub profile: Profile,
}

impl Default for Project {
    fn default() -> Self {
        Project {
            title: "New project".into(),
            units: Units::default(),
            processes: Processes::default(),
            cos_alpha: 1.0,
            water: WaterFlow::default(),
            time: TimeSettings::default(),
            atmosphere: Atmosphere::default(),
            root: RootUptake::default(),
            heat: HeatSettings::default(),
            solute: SoluteSettings::default(),
            profile: Profile::default(),
        }
    }
}

impl Project {
    pub fn n_solutes(&self) -> usize {
        if self.processes.solute {
            self.solute.species.len()
        } else {
            0
        }
    }

    /// Regenerate node depths with `n` equidistant nodes over `depth`
    /// (keeping the properties of the closest existing node).
    pub fn regrid_uniform(&mut self, depth: f64, n: usize) {
        let old = self.profile.nodes.clone();
        let d0 = old.first().map(|x| x.depth).unwrap_or(0.0);
        let mut nodes = Vec::with_capacity(n);
        for i in 0..n {
            let z = d0 + depth * i as f64 / (n as f64 - 1.0);
            let mut best = old[0].clone();
            let mut bd = f64::MAX;
            for o in &old {
                let d = (o.depth - z).abs();
                if d < bd {
                    bd = d;
                    best = o.clone();
                }
            }
            best.depth = z;
            nodes.push(best);
        }
        self.profile.nodes = nodes;
        let nn = n;
        self.profile.observation_nodes.retain(|&k| k >= 1 && k <= nn);
    }

    /// Basic validation; returns a list of human-readable problems.
    pub fn validate(&self) -> Vec<String> {
        let mut e = vec![];
        let n = self.profile.nodes.len();
		if let Some(meteo) = &self.atmosphere.meteo {
			if meteo.records.is_empty() {
				e.push("Meteorological ET is active, but no meteo records were provided.".into());
			}
			let dl0 = 0.667 * meteo.crop_height.max(0.1);
			if meteo.wind_height <= dl0 || meteo.temp_height <= dl0 {
				e.push(format!(
				"Measurement heights (wind: {} cm, temp: {} cm) must be greater than displacement height ({:.2} cm).",
				meteo.wind_height, meteo.temp_height, dl0
				));
			}
			if meteo.i_radiation == 0 && meteo.long_wave_a.abs() < 1e-6 {
				e.push("Cloudiness parameter long_wave_a (a1) cannot be zero when potential radiation is selected.".into());
			}
			if meteo.i_sun_sh == 3 && meteo.cloud_fact_ac.abs() < 1e-6 && meteo.cloud_fact_bc.abs() < 1e-6 {
				e.push("Cloudiness parameters (ac, bc) must be set when estimating from solar radiation.".into());
			}
		}
        if n < 3 {
            e.push("The profile needs at least 3 nodes.".into());
        }
        if n > 1001 {
            e.push("At most 1001 nodes are supported.".into());
        }
        if self.water.materials.is_empty() {
            e.push("At least one soil material is required.".into());
        }
        for (i, nd) in self.profile.nodes.iter().enumerate() {
            if nd.mat == 0 || nd.mat > self.water.materials.len() {
                e.push(format!("Node {} refers to undefined material {}.", i + 1, nd.mat));
                break;
            }
        }
        for w in self.profile.nodes.windows(2) {
            if w[1].depth <= w[0].depth {
                e.push("Node depths must be strictly increasing (top to bottom).".into());
                break;
            }
        }
        for (i, m) in self.water.materials.iter().enumerate() {
            if m.qs <= m.qr {
                e.push(format!("Material {}: θs must be larger than θr.", i + 1));
            }
            if m.alpha <= 0.0 || m.n <= 1.0 && self.water.model != SoilModel::BrooksCorey && self.water.model != SoilModel::Kosugi {
                e.push(format!("Material {}: α must be > 0 and n > 1.", i + 1));
            }
            if m.ks <= 0.0 {
                e.push(format!("Material {}: Ks must be > 0.", i + 1));
            }
        }
        if self.time.t_max <= self.time.t_init {
            e.push("Final time must be larger than the initial time.".into());
        }
        if self.time.dt <= 0.0 || self.time.dt_min <= 0.0 || self.time.dt_max < self.time.dt_min {
            e.push("Check the time-step settings (dt, dtMin, dtMax).".into());
        }
        let bc = &self.water.bc;
        if (bc.top_time_variable || bc.bot_time_variable || bc.atmospheric) && self.atmosphere.records.is_empty() {
            e.push("Time-variable boundary conditions require at least one atmospheric record.".into());
        }
        if self.processes.solute {
            if self.solute.species.is_empty() {
                e.push("Solute transport is on but no solute is defined.".into());
            }
            if self.solute.materials.len() < self.water.materials.len() {
                e.push("Solute parameters are missing for some materials.".into());
            }
        }
        if self.processes.heat && self.heat.materials.len() < self.water.materials.len() {
            e.push("Heat parameters are missing for some materials.".into());
        }
        e
    }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct MeteoRecord {
    pub t: f64,
    pub rad: f64,
    pub t_max: f64,
    pub t_min: f64,
    pub rh_mean: f64,
    pub wind_kmd: f64,
	pub sun_hours: f64,
    pub crop_height: Option<f64>,
    pub albedo: Option<f64>,
    pub lai: Option<f64>,
    pub x_root: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MeteoSettings {
    pub latitude: f64,
    pub altitude: f64,
    pub short_wave_a: f64,
    pub short_wave_b: f64,
    pub long_wave_a: f64,
    pub long_wave_b: f64,
    pub long_wave_a1: f64,
    pub long_wave_b1: f64,
    pub cloud_fact_ac: f64,
    pub cloud_fact_bc: f64,
    pub wind_height: f64,
    pub temp_height: f64,
    pub i_radiation: i32,
    pub i_sun_sh: i32,
    pub i_rel_hum: i32,
    pub hargreaves: bool,
	pub l_en_bal: bool,
    pub l_daily: bool,
    pub i_crop: i32,
    pub crop_height: f64,
    pub albedo: f64,
    pub i_lai: i32,
    pub lai: f64,
    pub r_extinct: f64,
    pub records: Vec<MeteoRecord>,
}

impl Default for MeteoSettings {
    fn default() -> Self {
        MeteoSettings {
            latitude: 40.0,
            altitude: 110.0,
            short_wave_a: 0.25,
            short_wave_b: 0.50,
            long_wave_a: 0.90,
            long_wave_b: 0.10,
            long_wave_a1: 0.34,
            long_wave_b1: -0.139,
            cloud_fact_ac: 1.35,
            cloud_fact_bc: -0.35,
            wind_height: 200.0,
            temp_height: 200.0,
            i_radiation: 1, // Solar Radiation default as shown in dialog
            i_sun_sh: 3,     // Solar Radiation cloudiness default as shown in dialog
            i_rel_hum: 0,
            hargreaves: false,
			l_en_bal: false,
			l_daily: false,
            i_crop: 1,
            crop_height: 0.0,
            albedo: 0.23,
            i_lai: 1,       // From Crop Height, Clipped Grass
            lai: 0.0,
            r_extinct: 0.463,
            records: vec![],
        }
    }
}
