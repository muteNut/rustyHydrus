use hydrus_core::*;

#[test]
fn native_project_roundtrip_and_run() {
    let p = examples::solute_pulse();
    let dir = std::env::temp_dir().join("hydrus_rt_test");
    std::fs::create_dir_all(&dir).unwrap();
    let f = dir.join("p.h1dr");
    hydrus_io::save_project(&p, &f).unwrap();
    let q = hydrus_io::load_project(&f).unwrap();
    assert_eq!(q.profile.nodes.len(), p.profile.nodes.len());
    let mut s = Simulation::new(q).unwrap();
    assert_eq!(s.run(|_| true), StepStatus::Finished);
    hydrus_io::writer::write_all(&s.res, &dir.join("out")).unwrap();
    assert!(dir.join("out/T_LEVEL.csv").exists());
}
