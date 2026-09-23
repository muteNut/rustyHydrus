use hydrus_core::{Simulation, StepStatus};
use hydrus_io::read_legacy_project;
use std::fs;

#[test]
fn test_end_to_end_meteo_energy_balance() {
    let dir = std::env::temp_dir().join("hydrus_test_meteo_eb");
    let _ = fs::create_dir_all(&dir);
    let p = &dir;

    let selector = r#"Pcp_File_Version=4
*** BLOCK A: BASIC INFORMATION *****************************************
Heading
Welcome to HYDRUS-1D
LUnit  TUnit  MUnit  (indicated units are obligatory for all input data)
cm
days
mmol
lWat   lChem lTemp  lSink lRoot lShort lWDep lScreen lVariabBC lEquil lInverse
 t     f     t      f     f     f      t     t       t         t         f
lSnow  lHP1   lMeteo  lVapor  lDummy  lFluxes lDummy  lDummy  lDummy  lDummy
 f       f       t       t       f       f       f       f       f       f
NMat    NLay  CosAlpha
  1       1       1
*** BLOCK B: WATER FLOW INFORMATION ************************************
MaxIt   TolTh   TolH       (maximum number of iterations and tolerances)
  10    0.001      1
TopInf WLayer KodTop InitCond
 t     t      -1       t
BotInf qGWLF FreeD SeepF KodBot DrainF  hSeep
 f     f     t     f     -1      f      0
    hTab1   hTabN
    1e-006   10000
    Model   Hysteresis
      0          0
   thr     ths    Alfa      n         Ks       l
  0.011   0.445  0.0277    1.38       34.2     0.5 
*** BLOCK C: TIME INFORMATION ******************************************
        dt       dtMin       dtMax     DMul    DMul2  ItMin ItMax  MPL
      0.001      1e-005        0.01     1.3     0.7     3     7    24
      tInit        tMax
        328         340
  lPrintD  nPrintSteps tPrintInterval lEnter
     f           1             1       t
TPrint(1),TPrint(2),...,TPrint(MPL)
      328.5         329       329.5         330       330.5         331 
      331.5         332       332.5         333       333.5         334 
      334.5         335       335.5         336       336.5         337 
      337.5         338       338.5         339       339.5         340 
*** BLOCK E: HEAT TRANSPORT INFORMATION *********************************************************
    Qn      Qo    Disper.    B1          B2          B3          Cn          Co           Cw
  0.555       0       5 1.47054e+016 -1.5518e+017 3.16617e+017 1.43327e+014 1.8737e+014 3.12035e+014 
      tAmpl     tPeriod    Campbell   MeltConst  lDummy  lDummy  lDummy  lDummy  lDummy
          0           0          0       0.43       f       f       f       f       f
      kTopT       TTop      kBotT       TBot
         -1         20           0         20
*** END OF INPUT FILE 'SELECTOR.IN' ************************************
"#;

    let meteo = r#"Pcp_File_Version=4
* METEOROLOGICAL PARAMETERS AND INFORMATION |||||||||||||||||||||||||||||||
 MeteoRecords Radiation Penman-Hargreaves
            2        1       f
  lEnBal  lDaily  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy
       t       f       f       f       f       t       f       f       f       f
 Latitude  Altitude
    33.58        306
 ShortWaveRadA  ShortWaveRadB
          0.25            0.5
 LongWaveRadA   LongWaveRadB
           0.9            0.1
 LongWaveRadA1  LongWaveRadB1
          0.34         -0.139
 WindHeight     TempHeight
        200            150
 iCrop (=0: no crop, =1: constant, =2: table, =3: daily)  SunShine  RelativeHum
         0                                                2         0
    Albedo
      0.23
Daily values
       t        Rad        TMax        TMin     RHMean      Wind    SunHours CropHeight     Albedo   LAI(SCF)      rRoot
      [T]  [MJ/m2/d]       [C]         [C]       [%]     [km/d]     [hour]      [L]           [-]        [-]        [L]
   328.042          0       15.2       15.2         33     138.24      0.545 
   328.083          0       14.9       14.9         32     103.68      0.545 
end *** END OF INPUT FILE 'METEO.IN' **********************************
"#;

    let atmosph = r#"Pcp_File_Version=4
*** BLOCK I: ATMOSPHERIC INFORMATION  **********************************
   MaxAL                    (MaxAL = number of atmospheric data-records)
      2
 DailyVar  SinusVar  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy
       f       f       f       f       f       t       f       f       f       f
 hCritS                 (max. allowed pressure head at the soil surface)
      0
       tAtm        Prec       rSoil       rRoot      hCritA          rB          hB          ht        tTop        tBot        Ampl
        329           0           0           0      100000           0           0           0           0           0           0 
        340           0           0           0      100000           0           0           0           0           0           0 
end*** END OF INPUT FILE 'ATMOSPH.IN' **********************************
"#;

    let mut profile = String::from(
        r#"Pcp_File_Version=4
    2
    1  0.000000e+000  1.000000e+000  1.000000e+000
    2 -5.000000e+001  1.000000e+000  1.000000e+000
  101    1    0    1 x         h      Mat  Lay      Beta           Axz            Bxz            Dxz          Temp          Conc 
"#,
    );
    for i in 1..=101 {
        let x = -0.5 * (i - 1) as f64;
        profile.push_str(&format!(
            "{:5} {:14.6e}  1.300000e-001    1    1  0.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000  1.700000e+001\n",
            i, x
        ));
    }
    profile.push_str("    3\n    5   15   25\n");

    fs::write(p.join("SELECTOR.IN"), selector).unwrap();
    fs::write(p.join("METEO.IN"), meteo).unwrap();
    fs::write(p.join("ATMOSPH.IN"), atmosph).unwrap();
    fs::write(p.join("PROFILE.DAT"), profile).unwrap();

    let prj = read_legacy_project(p).expect("read legacy project");
    assert!(prj.processes.water_flow);
    assert!(prj.processes.heat);
    let meteo_cfg = prj.atmosphere.meteo.as_ref().expect("meteo parsed");
    assert!(meteo_cfg.l_en_bal);
    assert_eq!(meteo_cfg.records.len(), 2);

    let mut sim = Simulation::new(prj).expect("init simulation");
    assert!(sim.t >= 328.0);

    // Step forward past the first sub-daily hourly record
	while sim.t < 340.0 {
        let status = sim.step();
        assert_eq!(status, StepStatus::Running, "step failed at t = {}: {:?}", sim.t, status);
    }

    assert!((sim.t - 340.0).abs() < 1e-3);

    let _ = fs::remove_dir_all(&dir);
}