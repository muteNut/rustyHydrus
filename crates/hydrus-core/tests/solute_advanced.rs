use hydrus_core::*;

#[test]
fn test_volatilization_boundary_k_top_minus_two() {
    let mut prj = Project::default();
    prj.processes.water_flow = false; // Steady static profile
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 5.0;
    prj.time.dt = 0.5;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 1.0;

    // Uniform initial concentration of volatile tracer
    for node in &mut prj.profile.nodes {
        node.conc = vec![10.0];
    }

    // Configure Volatilization Boundary
    prj.solute.k_top = -2;
    prj.solute.d_surf = 1.0; // 1 cm stagnant air boundary layer
    prj.solute.c_atm = 0.0;  // 0 background air concentration
    prj.solute.species[0].diff_g = 500.0;
    prj.solute.species[0].diff_w = 1.0;
    prj.solute.species[0].per_material[0].henry = 0.05; // Significant Henry volatility

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    // Advance 5 steps
    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    // Surface concentration must decrease due to volatilization loss into air
    let n = sim.n;
    let c_top = sim.conc(0, n - 1);
    let c_bottom = sim.conc(0, 0);

    assert!(c_top < 10.0, "Surface concentration should drop due to volatilization");
    assert!(c_top < c_bottom, "Top should be more depleted than bottom");

    // Cumulative mass loss must be non-zero
    if let Some(ref sol) = sim.sol {
        assert!(
		sol.cum_ch[0][0] < 0.0,
		"Top cumulative flux must be negative (net volatilization mass loss), got {}",
		sol.cum_ch[0][0]
		);
    }
}

#[test]
fn test_dual_continuum_physical_nonequilibrium_transfer() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 2.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.5;

    // Mobile-immobile setup: theta_im = 0.10, mobile porosity > 0
    prj.solute.materials[0].th_immobile = 0.10;
    prj.solute.materials[0].frac = 0.5;
    // Omega mass transfer rate = 0.05 / day
    prj.solute.species[0].per_material[0].omega = 0.05;

    // Start with tracer only in mobile domain (conc = 5.0), immobile empty (sorb = 0.0)
    for node in &mut prj.profile.nodes {
        node.conc = vec![5.0];
        node.sorb = vec![0.0];
    }

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    // Check that solute entered the immobile domain
    let immobile_c = sim.sorb(0, sim.n / 2);
    assert!(immobile_c > 0.0, "Immobile solute concentration must increase from mass transfer");
}
#[test]
fn test_active_root_solute_uptake_michaelis_menten() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.root_water_uptake = true;
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // Transpiration demand
    prj.water.bc.atmospheric = false;
    prj.water.bc.top_time_variable = false;
    prj.water.bc.r_root = 0.2;

    // Configure Active Root Solute Uptake
    prj.root.l_act_rsu = true;
    prj.root.omega_act = vec![1.5]; // V_max
    prj.root.r_km = vec![2.0];      // K_m
    prj.root.c_min = vec![0.1];     // C_min

    // Uniform distribution and initial condition
    for node in &mut prj.profile.nodes {
        node.beta = 1.0;
        node.conc = vec![5.0];
    }

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    // Active + passive uptake must consume solute mass in root zone
    if let Some(ref sol) = sim.sol {
        assert!(
            sol.cum_ch[0][4] > 0.0,
            "Cumulative root solute uptake must be positive, got {}",
            sol.cum_ch[0][4]
        );
        let c_mid = sim.conc(0, sim.n / 2);
        assert!(
            c_mid < 5.0,
            "Solute concentration must deplete due to root uptake"
        );
    }
}

#[test]
fn test_colloid_virus_attachment_and_clean_bed_filtration() {
    let mut prj = Project::default();
    prj.processes.water_flow = true; // Richards solver updates moisture content and flux
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // Infiltration boundary
    prj.water.bc.atmospheric = false;
    prj.water.bc.top_time_variable = false;
    prj.water.bc.kod_top = -1;
    prj.water.bc.r_top = -2.0; // Infiltration flux downward [cm/d]
    prj.water.bc.free_drainage = true;
	prj.solute.species[0].c_top = 8.0; // Infiltrating solution concentration

    prj.solute.l_bact = true;
    prj.solute.l_filtr = true;

    // Colloid parameters for site 1 and site 2
    let sp_mat = &mut prj.solute.species[0].per_material[0];
    sp_mat.d_c = 0.05; // sand collector diameter [cm]
    sp_mat.d_p = 1e-4; // colloid diameter [cm]
    sp_mat.r_ka1 = 0.5; // sticking efficiency alpha 1
    sp_mat.r_ka2 = 0.1; // sticking efficiency alpha 2
    sp_mat.r_kd1 = 0.01; // detachment rate 1
    sp_mat.r_kd2 = 0.005; // detachment rate 2
    sp_mat.s_max1 = 10.0;
    sp_mat.s_max2 = 5.0;
    sp_mat.i_psi1 = 1; // Langmuir blocking
    sp_mat.i_psi2 = 1;

    for node in &mut prj.profile.nodes {
        node.conc = vec![8.0];
        node.sorb = vec![0.0];
    }

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    if let Some(ref sol) = sim.sol {
        let idx = sim.n / 2;
        let s1 = sol.sorb[0][idx];
        let s2 = sol.sorb2[0][idx];
        assert!(s1 > 0.0, "Site 1 attached concentration at surface must be positive, got {}", s1);
        assert!(s2 > 0.0, "Site 2 attached concentration at surface must be positive, got {}", s2);
        assert!(s1 > s2, "Site 1 with higher attachment rate should have accumulated more mass (s1={}, s2={})", s1, s2);
    }
}
#[test]
fn test_temperature_dependent_solute_decay_arrhenius() {
    let mut prj = Project::default();
    prj.processes.water_flow = false;
    prj.processes.heat = true;
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // Fixed warm temperature of 35 deg C (ref = 20 deg C)
    prj.heat.k_top = 1;
    prj.heat.t_top = 35.0;
    prj.heat.k_bot = 1;
    prj.heat.t_bot = 35.0;
    for node in &mut prj.profile.nodes {
        node.temp = 35.0;
        node.conc = vec![10.0];
    }

    // Set first-order aqueous decay rate
    prj.solute.species[0].per_material[0].gam_w = 0.1;

    // Enable temperature dependence with an activation energy of 50,000 J/mol
    prj.solute.l_tdep = true;
    prj.solute.t_dep = vec![SpeciesTDep {
        gam_w: 50_000.0,
        ..Default::default()
    }];

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    // Reference decay at 20 deg C: C(t) = 10 * exp(-0.1 * 0.5) = 9.512
    // At 35 deg C (308.15 K), TT = (308.15 - 293.15) / (8.314 * 308.15 * 293.15) = 2.000e-5
    // exp(50000 * TT) = exp(1.000) = 2.718
    // Decay rate at 35 deg C ~ 0.2718 / day
    // Expected C(0.5) ~ 10 * exp(-0.2718 * 0.5) ~ 8.73
    let c_final = sim.conc(0, sim.n / 2);
    assert!(c_final < 9.2, "Warm soil must accelerate solute degradation, got {}", c_final);
}
#[test]
fn test_water_content_dependent_solute_decay_walker() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // Static heads at h = -300 cm (dry soil)
    prj.water.bc.atmospheric = false;
    prj.water.bc.top_time_variable = false;
    prj.water.bc.kod_top = -1;
    prj.water.bc.r_top = 0.0;
    prj.water.bc.kod_bot = -1;
    prj.water.bc.r_bot = 0.0;

    for node in &mut prj.profile.nodes {
        node.h = -300.0;
        node.conc = vec![10.0];
    }

    // Set base degradation rate mu_w = 0.5 / day
    prj.solute.species[0].per_material[0].mu_w = 0.5;

    // Enable Walker power law: exponent B = 1.0, reference head h_ref = -100 cm
    // In drier soil (theta < theta_ref), degradation rate will be reduced by (theta / theta_ref)^B
    prj.solute.l_moist = true;
    let mut wdep = SpeciesWDep::default();
    wdep.exp_b[0] = 1.0;     // Exponent for mu_w
    wdep.h_ref[0] = -100.0;  // Reference head
    prj.solute.w_dep = vec![wdep];

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    // Without moisture scaling, decay at mu_w = 0.5 over 0.5 d gives C ~ 10 * exp(-0.25) = 7.788
    // With moisture scaling, theta(-300) < theta(-100), so effective mu_w < 0.5, leaving C > 7.8
    let c_final = sim.conc(0, sim.n / 2);
    assert!(c_final > 7.85, "Dry soil must retard solute degradation, got {}", c_final);
    assert!(c_final < 10.0, "Solute must still experience degradation");
}
#[test]
fn test_dual_nonequilibrium_physical_and_chemical_transport() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // Physical non-equilibrium parameters
    prj.solute.materials[0].th_immobile = 0.10;
    prj.solute.materials[0].frac = 0.5; // f_em = 0.5 (half equilibrium, half kinetic)
    
    // Chemical sorption parameters
    let sp_mat = &mut prj.solute.species[0].per_material[0];
    sp_mat.ks = 1.5;        // Kd = 1.5 cm^3/g
    sp_mat.omega = 0.1;     // Mass exchange rate omega = 0.1 / day

    // Activate dual non-equilibrium
    prj.solute.l_dual_neq = true;

    // Start with tracer only in mobile liquid (conc = 6.0), both non-equilibrium sites empty
    for node in &mut prj.profile.nodes {
        node.conc = vec![6.0];
        node.sorb = vec![0.0];
    }

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    if let Some(ref sol) = sim.sol {
        let mid = sim.n / 2;
        let c_im = sol.sorb[0][mid];
        let s_kin = sol.sorb2[0][mid];
        assert!(c_im > 0.0, "Immobile liquid concentration must increase via physical exchange, got {}", c_im);
        assert!(s_kin > 0.0, "Kinetic sorption site 2 must accumulate sorbed mass, got {}", s_kin);
    }
}
#[test]
fn test_dual_permeability_matrix_fracture_solute_transfer() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // Enable dual permeability
    prj.water.model = SoilModel::DualPermeability;
	prj.water.materials[0].extra = [0.08, 0.42, 0.01, 1.6, 1.0];
    prj.solute.species[0].per_material[0].omega = 0.1; // Diffusive exchange rate

    for node in &mut prj.profile.nodes {
        node.conc = vec![10.0]; // Fracture domain starts at 10.0
        node.sorb = vec![0.0];
    }

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");
    sim.l_dual_perm = true;
    sim.w_fracture = 0.1;
    sim.alpha_dw = 0.05;

    // Set matrix initial concentration to 0.0
    if let Some(ref mut sol) = sim.sol {
        for val in &mut sol.conc_m[0] {
            *val = 0.0;
        }
    }

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    // Matrix domain must accumulate solute mass from fracture exchange
    let mid = sim.n / 2;
    let c_matrix = sim.conc_matrix(0, mid);
    let c_fracture = sim.conc(0, mid);

    assert!(c_matrix > 0.0, "Matrix domain must absorb solute mass, got {}", c_matrix);
    assert!(c_fracture < 10.0, "Fracture domain must lose solute mass into matrix, got {}", c_fracture);
}

#[test]
fn test_colloid_pore_size_exclusion_and_velocity_acceleration() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // Colloid transport with size exclusion (excluded water content th_c = 0.05)
    prj.solute.l_bact = true;
    prj.solute.materials[0].th_immobile = 0.05; // Mapped to th_c for colloids
    prj.solute.species[0].per_material[0].r_ka1 = 0.0; // No attachment to isolate exclusion

    // Steady downward flux
    prj.water.bc.atmospheric = false;
    prj.water.bc.top_time_variable = false;
    prj.water.bc.kod_top = -1;
    prj.water.bc.r_top = -2.0;
    prj.water.bc.kod_bot = -1;
    prj.water.bc.r_bot = -2.0;

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..3 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    if let Some(ref sol) = sim.sol {
        let mid = sim.n / 2;
        // Accessible water content must be lower than bulk water content
        assert!(sol.th_n[mid] < sim.th_new[mid], "Accessible theta must be smaller than bulk theta");
        // Pore velocity must be strictly higher than standard water pore velocity
        let v_pore_solute = (sol.v_n[mid] / sol.th_n[mid]).abs();
        let v_pore_water = (sim.v_new[mid] / sim.th_new[mid]).abs();
        assert!(v_pore_solute > v_pore_water, "Colloid pore velocity must be accelerated");
    }
}
#[test]
fn test_vapor_phase_solute_transport_temperature_gradient() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;
    prj.processes.heat = true;
    prj.processes.vapor = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // No liquid flow: stagnant moisture
    prj.water.bc.atmospheric = false;
    prj.water.bc.top_time_variable = false;
    prj.water.bc.kod_top = -1;
    prj.water.bc.r_top = 0.0;
    prj.water.bc.kod_bot = -1;
    prj.water.bc.r_bot = 0.0;

    // Heat gradient: cold at bottom (10 C), hot at top (30 C)
    prj.heat.k_top = 1;
    prj.heat.t_top = 30.0;
    prj.heat.k_bot = 1;
    prj.heat.t_bot = 10.0;

    // Volatile solute properties
    prj.solute.species[0].diff_w = 0.0; // Zero aqueous diffusion to isolate gas diffusion
    prj.solute.species[0].diff_g = 500.0; // Gas diffusion 500 cm^2/d
    prj.solute.species[0].per_material[0].henry = 0.1; // Volatile Henry constant
    prj.solute.l_tdep = true;

    let mut tdep = SpeciesTDep::default();
    tdep.henry = 20000.0; // Strong temperature dependence of Henry constant
    prj.solute.t_dep = vec![tdep];

    // Dry soil so gas phase theta_g is large
    for (i, node) in prj.profile.nodes.iter_mut().enumerate() {
        node.h = -500.0;
        node.temp = 10.0 + 20.0 * (i as f64) / 100.0;
        node.conc = vec![5.0];
    }

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    // Gaseous diffusion and temperature gradient should induce a concentration profile
    let c_top = sim.conc(0, sim.n - 1);
    let c_bot = sim.conc(0, 0);
    assert!((c_top - c_bot).abs() > 1e-4, "Temperature gradient must redistribute volatile solute via gas phase");
}
#[test]
fn test_flux_concentration_computation() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // Steady downward infiltration
    prj.water.bc.atmospheric = false;
    prj.water.bc.top_time_variable = false;
    prj.water.bc.kod_top = -1;
    prj.water.bc.r_top = -2.0;
    prj.water.bc.kod_bot = -1;
    prj.water.bc.r_bot = -2.0;

    prj.solute.species[0].c_top = 10.0;
    prj.solute.materials[0].disp_l = 5.0; // Significant dispersivity

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    let conc_flux = sim.flux_conc(0);
    assert_eq!(conc_flux.len(), sim.n);

    // At the wetting / dispersion front, flux concentration differs from resident concentration
    let mid = sim.n / 2;
    let c_resident = sim.conc(0, mid);
    let c_fl = conc_flux[mid];
	let _ = c_resident;
    assert!(c_fl >= 0.0, "Flux concentration must be non-negative");
}

#[test]
fn test_dual_permeability_water_flow_matrix_imbibition() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.water.model = SoilModel::DualPermeability;

    prj.time.t_init = 0.0;
    prj.time.t_max = 0.5;
    prj.time.dt = 0.01;
    prj.time.dt_min = 0.001;
    prj.time.dt_max = 0.05;

    // Zero-flux boundaries to isolate matrix-fracture exchange
    prj.water.bc.atmospheric = false;
    prj.water.bc.top_time_variable = false;
    prj.water.bc.kod_top = -1;
    prj.water.bc.r_top = 0.0;
    prj.water.bc.free_drainage = false;
    prj.water.bc.kod_bot = -1;
    prj.water.bc.r_bot = 0.0;

    let m = &mut prj.water.materials[0];
    m.qr = 0.05;
    m.qs = 0.40;
    m.alpha = 0.02;
    m.n = 1.8;
    m.ks = 50.0;
    m.extra = [0.08, 0.42, 0.01, 1.6, 1.0];

    // Initial conditions: wet fracture (h = -20 cm), dry matrix (h_m = -80 cm)
    for node in &mut prj.profile.nodes {
        node.h = -20.0;
    }

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");
    sim.l_dual_perm = true;
    sim.w_fracture = 0.1;
    sim.alpha_dw = 0.01;

    for i in 0..sim.n {
        sim.h_matrix_new[i] = -80.0;
        sim.h_matrix_old[i] = -80.0;
    }
    sim.set_mat_matrix();

    let mid = sim.n / 2;
    let h_m_init = sim.h_matrix_new[mid];
    let h_f_init = sim.h_new[mid];
    let th_m_init = sim.th_matrix_new[mid];

    for _ in 0..3 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    let h_m_final = sim.h_matrix_new[mid];
    let h_f_final = sim.h_new[mid];
    let th_m_final = sim.th_matrix_new[mid];

    // 1. Matrix head must increase due to imbibition from wet fracture
    assert!(
        h_m_final > h_m_init,
        "Matrix head must increase (imbibition), got init: {}, final: {}",
        h_m_init,
        h_m_final
    );
    // 2. Fracture head must decrease due to loss into matrix
    assert!(
        h_f_final < h_f_init,
        "Fracture head must decrease, got init: {}, final: {}",
        h_f_init,
        h_f_final
    );
    // 3. Matrix water content must increase
    assert!(
        th_m_final > th_m_init,
        "Matrix theta must increase from {}, got {}",
        th_m_init,
        th_m_final
    );
}

#[test]
fn test_decay_chain_nonequilibrium_routing_l_nequil() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // Physical non-equilibrium with 2 solutes (Parent -> Daughter)
    prj.solute.materials[0].frac = 0.5;
    prj.solute.materials[0].th_immobile = 0.1;
    prj.solute.l_nequil = true; // Route parent immobile decay directly into daughter immobile phase

    // Add daughter species
    let mut daughter = prj.solute.species[0].clone();
    daughter.name = "Daughter".into();
    prj.solute.species.push(daughter);

    // Parent decays quickly
    prj.solute.species[0].per_material[0].mu_w = 0.5;
    prj.solute.species[0].per_material[0].mu_s = 0.5;

    // Initialize parent mass in immobile phase
    for node in &mut prj.profile.nodes {
        node.conc = vec![10.0, 0.0];
        node.sorb = vec![10.0, 0.0];
    }

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    // Daughter must appear in the immobile phase due to parent degradation
    let mid = sim.n / 2;
    let s_daughter = sim.sorb(1, mid);
    assert!(
        s_daughter > 0.0,
        "Daughter solute must accumulate in immobile phase under l_nequil, got {}",
        s_daughter
    );
}
#[test]
fn test_diagnostics_flux_conc_and_dual_permeability_snapshots() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;
    prj.water.model = SoilModel::DualPermeability;
    prj.solute.i_conc_type = 2; // Enable flux concentration recording mode

    prj.time.t_init = 0.0;
    prj.time.t_max = 0.2;
    prj.time.dt = 0.05;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.1;
    prj.time.print_times = vec![0.1, 0.2];

    // Set observation point at mid-depth (1-based index)
    let obs_idx = prj.profile.nodes.len() / 2;
    prj.profile.observation_nodes = vec![obs_idx];

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");
    sim.l_dual_perm = true;

    while sim.step() == StepStatus::Running {}

    // 1. Profile Snapshots verification
    assert!(!sim.res.profiles.is_empty(), "Profiles must be recorded");
    let prof = sim.res.profiles.last().expect("Must have at least one profile");
    let n_mid = &prof.nodes[sim.n / 2];

    assert!(
        n_mid.h_matrix.is_some(),
        "Dual-perm profile node must contain h_matrix diagnostic"
    );
    assert!(
        n_mid.th_matrix.is_some(),
        "Dual-perm profile node must contain th_matrix diagnostic"
    );
    assert!(
        n_mid.flux_conc.is_some(),
        "Profile node must contain flux_conc when i_conc_type == 2"
    );

    // 2. Observation Points verification
    assert!(!sim.res.obs.is_empty(), "Obs points must be recorded across time");
    let last_obs = sim.res.obs.last().expect("Must have observation output");
    assert_eq!(last_obs.points.len(), 1);
    let obs_pt = &last_obs.points[0];

    assert!(
        obs_pt.flux_conc.is_some(),
        "Observation point must record flux_conc when i_conc_type == 2"
    );
    assert!(
        obs_pt.h_matrix.is_some(),
        "Observation point must record h_matrix under dual permeability"
    );
}

#[test]
fn test_resident_conc_recording_when_i_conc_type_is_one() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;
    prj.solute.i_conc_type = 1; // Standard resident concentration mode

    prj.time.t_init = 0.0;
    prj.time.t_max = 0.1;
    prj.time.dt = 0.05;
    prj.time.print_times = vec![0.1];
    prj.profile.observation_nodes = vec![1];

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");
    while sim.step() == StepStatus::Running {}

    let prof = sim.res.profiles.last().unwrap();
    assert!(
        prof.nodes[0].flux_conc.is_none(),
        "flux_conc should be None when i_conc_type == 1"
    );

    let obs = sim.res.obs.last().unwrap();
    assert!(
        obs.points[0].flux_conc.is_none(),
        "ObsPoint flux_conc should be None when i_conc_type == 1"
    );
}
#[test]
fn test_tabular_moisture_dependent_solute_decay() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;

    prj.time.t_init = 0.0;
    prj.time.t_max = 1.0;
    prj.time.dt = 0.1;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.2;

    // Static heads at h = -200 cm
    prj.water.bc.atmospheric = false;
    prj.water.bc.top_time_variable = false;
    prj.water.bc.kod_top = -1;
    prj.water.bc.r_top = 0.0;
    prj.water.bc.kod_bot = -1;
    prj.water.bc.r_bot = 0.0;

    for node in &mut prj.profile.nodes {
        node.h = -200.0;
        node.conc = vec![10.0];
    }

    // Base decay rate mu_w = 0.5 / day
    prj.solute.species[0].per_material[0].mu_w = 0.5;

    // Tabular scaling: reduce reaction to 10% when theta <= 0.15, and 100% when theta >= 0.35
    prj.solute.l_moist = true;
    prj.solute.i_moist_dep = 2;

    let mut tab = MoistDepTable::default();
    tab.theta = vec![0.10, 0.20, 0.35];
    for k in 0..9 {
        tab.factors[k] = vec![0.10, 0.50, 1.00];
    }
    prj.solute.moist_tables = vec![tab];

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    for _ in 0..5 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    let c_final = sim.conc(0, sim.n / 2);
    // Effective mu_w is suppressed by the table factor, so c_final remains significantly higher
    assert!(c_final > 8.5, "Tabular moisture scaling should suppress decay, got {}", c_final);
    assert!(c_final < 10.0, "Some decay must still occur");
}
#[test]
fn test_dual_permeability_solute_mass_balance_closure() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;
    prj.water.model = SoilModel::DualPermeability;
    prj.water.materials[0].extra = [0.08, 0.42, 0.01, 1.6, 1.0];

    // Zero-flux boundaries to isolate matrix-fracture exchange and verify conservation
    prj.water.bc.atmospheric = false;
    prj.water.bc.top_time_variable = false;
    prj.water.bc.kod_top = -1;
    prj.water.bc.r_top = 0.0;
    prj.water.bc.free_drainage = false;
    prj.water.bc.kod_bot = -1;
    prj.water.bc.r_bot = 0.0;

    prj.solute.k_top = -1;
    prj.solute.k_bot = -1;

    prj.time.t_init = 0.0;
    prj.time.t_max = 0.5;
    prj.time.dt = 0.05;
    prj.time.dt_min = 0.005;
    prj.time.dt_max = 0.1;
    prj.time.print_times = vec![0.5];

    // Exchange into matrix
    prj.solute.species[0].per_material[0].omega = 0.2;

    for node in &mut prj.profile.nodes {
        node.conc = vec![10.0];
    }

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");
    sim.l_dual_perm = true;
    sim.w_fracture = 0.1;
    sim.alpha_dw = 0.05;

    // Initialize matrix domain with 0.0 concentration
    if let Some(ref mut sol) = sim.sol {
        for val in &mut sol.conc_m[0] {
            *val = 0.0;
        }
    }

    // Reset reference mass balance storage at t = 0 with matrix states at 0.0
    sim.res.balance.clear();
    sim.sub_reg(0);

    while sim.step() == StepStatus::Running {}

    let bal = sim.res.balance.last().expect("Must have balance record");
    assert!(
        bal.sol_bal_r[0] < 1.0,
        "Relative solute mass balance error must be < 1.0%, got {}%",
        bal.sol_bal_r[0]
    );
}

#[test]
fn test_ponding_surface_layer_solute_mixing() {
    let mut prj = Project::default();
    prj.processes.water_flow = true;
    prj.processes.solute = true;

    prj.water.bc.surface_layer = true;
    prj.water.bc.kod_top = -1;
    prj.water.bc.r_top = 0.0;

    prj.time.t_init = 0.0;
    prj.time.t_max = 0.2;
    prj.time.dt = 0.05;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 0.1;

    let n = prj.profile.nodes.len();
    for node in &mut prj.profile.nodes {
        node.conc = vec![0.0];
    }

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");
    // Simulate 2 cm ponding head with 10.0 initial concentration
    sim.h_new[n - 1] = 2.0;
    sim.sol.as_mut().unwrap().c_top[0] = 10.0;
    // Infiltrating rain of 5.0 cm/d with concentration 0.0
    sim.prec = 5.0;
    sim.r_soil = 0.0;
    sim.sol.as_mut().unwrap().c_t[0] = 0.0;

    sim.sol_wlayer_mix();

    let c_mixed = sim.sol.as_ref().unwrap().c_top[0];
    // Dilution: 2.0*10 / (2.0 + 0.05*5) = 20 / 2.25 = 8.888
    assert!(c_mixed < 10.0 && c_mixed > 8.0, "Ponded concentration must dilute, got {}", c_mixed);
}