use hydrus_core::*;

#[test]
fn test_dual_porosity_model_6_water_transfer() {
    let mut prj = Project::default();
    prj.water.model = SoilModel::DualPorosityW;
    prj.time.t_init = 0.0;
    prj.time.t_max = 600.0;
    prj.time.dt = 1.0;
    prj.time.dt_min = 0.01;
    prj.time.dt_max = 10.0;

    // Skaggs' Column Test 1 parameters
    let mut mat = SoilMaterial::default();
    mat.qr = 0.0;
    mat.qs = 0.20;
    mat.alpha = 0.041;
    mat.n = 1.964;
    mat.ks = 0.000722;
    mat.l = 0.5;
    // [thr_im, ths_im, omega, alpha_im, ks_im]
    mat.extra = [0.0, 0.15, 1e-4, 0.0, 0.0];
    prj.water.materials = vec![mat];

    // Constant ponding at top (Dirichlet head = 1.0), free drainage at bottom
    prj.water.bc.kod_top = 1;
    prj.water.bc.free_drainage = true;
    prj.water.bc.kod_bot = -5;

    // Start with relatively dry soil profile throughout, except the top surface
    for node in &mut prj.profile.nodes {
        node.h = -100.0;
    }
    // Set top boundary node head (index 0 in Profile is top surface):
    prj.profile.nodes[0].h = 1.0;

    let mut sim = Simulation::new(prj).expect("Simulation should initialize");

    // Initial immobile water content should start at ths_im
    assert_eq!(sim.i_dual_por, 1);
    assert_eq!(sim.th_old_im.len(), sim.n);

    // Advance 10 steps
    for _ in 0..10 {
        let st = sim.step();
        assert_eq!(st, StepStatus::Running);
    }

    // Verify that water transfer sink/source term has been evaluated
    assert!(sim.sink_im.iter().any(|&s| s.abs() > 0.0));
}