use hydrus_core::*;

fn run(p: Project) -> Simulation {
    let mut s = Simulation::new(p).expect("valid");
    let st = s.run(|_| true);
    assert_eq!(st, StepStatus::Finished, "{:?}", s.res.messages);
    s
}

#[test]
fn infiltration_example_balances() {
    let s = run(examples::infiltration_evaporation());
    let b = s.res.balance.last().unwrap();
    assert!(b.wat_bal_r < 1.0, "water balance error {}%", b.wat_bal_r);
    assert!(s.res.profiles.len() >= 7);
}

#[test]
fn solute_example_balances() {
    let s = run(examples::solute_pulse());
    let b = s.res.balance.last().unwrap();
    assert!(b.wat_bal_r < 0.5, "water balance error {}%", b.wat_bal_r);
    assert!(b.sol_bal_r[0] < 1.0, "solute balance error {}%", b.sol_bal_r[0]);
    // tracer must have arrived at 100 cm depth
    let last = s.res.tlevel.last().unwrap();
    assert!(last.solutes[0].cum_top > 0.0);
}

#[test]
fn hydraulic_functions_are_consistent() {
    use hydrus_core::material::*;
    let m = SoilMaterial::default();
    let p = par_of(&m, SoilModel::VanGenuchten, 100.0);
    // C = dtheta/dh
    let h = -50.0;
    let c = fc(SoilModel::VanGenuchten, h, &p);
    let num = (fq(SoilModel::VanGenuchten, h + 1e-4, &p) - fq(SoilModel::VanGenuchten, h - 1e-4, &p)) / 2e-4;
    assert!((c - num).abs() / c < 1e-5);
    // h(Se) inverse
    let se = fs(SoilModel::VanGenuchten, h, &p);
    assert!((fh(SoilModel::VanGenuchten, se, &p) - h).abs() < 1e-6);
}
