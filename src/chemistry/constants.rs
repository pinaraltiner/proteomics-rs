
// Source: http://pdg.lbl.gov/2012/reviews/rpp2012-rev-phys-constants.pdf
pub const AVERAGE_AA_MASS: f64 = 111.1254; // TODO: marco => why difference with 111.10523866044295 by computation
// TODO: refine this value and put a source reference here (publication ?)
pub const AVERAGE_PEPTIDE_ISOTOPE_MASS_DIFF: f64 = 1.0027;

pub const ELECTRON_MASS: f64 = 0.00054857990946; // Source: NIST 2010 CODATA
pub const PROTON_MASS: f64 = 1.007276466812; // Source: NIST 2010 CODATA

pub const CO_MONO_MASS: f64 = 27.99491461956;
pub const CO2_MONO_MASS: f64 = 0.0;  // FIXME
pub const H2O_MONO_MASS: f64 = 18.010565;
pub const NH3_MONO_MASS: f64 = 17.02654910101;

pub const WATER_MONO_MASS: f64 = 18.010565;
pub const WATER_AVERAGE_MASS: f64 = 18.01525697318;

#[allow(dead_code)]
pub mod aa {
    pub const A: char = 'A';
    pub const B: char = 'B';
    pub const C: char = 'C';
    pub const D: char = 'D';
    pub const E: char = 'E';
    pub const F: char = 'F';
    pub const G: char = 'G';
    pub const H: char = 'H';
    pub const J: char = 'J';
    pub const I: char = 'I';
    pub const K: char = 'K';
    pub const L: char = 'L';
    pub const M: char = 'M';
    pub const N: char = 'N';
    pub const O: char = 'O';
    pub const P: char = 'P';
    pub const Q: char = 'Q';
    pub const R: char = 'R';
    pub const S: char = 'S';
    pub const T: char = 'T';
    pub const U: char = 'U';
    pub const V: char = 'V';
    pub const W: char = 'W';
    pub const X: char = 'X';
    pub const Y: char = 'Y';
    pub const Z: char = 'Z';
}

#[allow(dead_code)]
pub mod atom {
    pub const C: &'static str = "C";
    pub const H: &'static str = "H";
    pub const O: &'static str = "O";
    pub const N: &'static str = "N";
    pub const P: &'static str = "P";
    pub const S: &'static str = "S";
}