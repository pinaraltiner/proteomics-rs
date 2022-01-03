#[allow(dead_code)]

use crate::chemistry::constants::PROTON_MASS;

#[allow(non_camel_case_types)]
#[derive(Clone, Eq, PartialEq, Debug)]
enum MassTolUnit {
    Da,
    mmu,
    ppm
}

impl MassTolUnit {
    fn new(unit: &str) -> Option<MassTolUnit> {
        match unit {
            "Da"  => Some(MassTolUnit::Da),
            "mmu" => Some(MassTolUnit::mmu),
            "ppm" => Some(MassTolUnit::ppm),
            _     => None
        }
    }
}

impl std::fmt::Display for MassTolUnit {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

fn calc_mz_tol_in_daltons(mz: f64, mz_tol: f64, tol_unit: MassTolUnit) -> f64 {
    match tol_unit {
        MassTolUnit::Da => mz_tol,
        MassTolUnit::mmu => mz_tol / 1000.0,
        MassTolUnit::ppm => mz_tol * mz / 1000000.0
    }
}

pub fn mz_to_mass( mz: f64, charge: i32 ) -> f64 {
    let z = charge as f64;
    mz * z.abs() - z * PROTON_MASS
}
pub fn mass_to_mz( mass: f64, charge: i32 ) -> f64 {
    let z = charge as f64;
    (mass + z * PROTON_MASS) / z.abs()
}

