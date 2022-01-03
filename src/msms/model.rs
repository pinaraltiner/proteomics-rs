//use std::fmt::Display;
// You need to bring the ToString trait into scope to use it
//use std::string::ToString;
//use strum_macros;

use crate::atom::H;
use crate::BIOMOLECULE_ATOM_TABLE;
use crate::chemistry::constants::*;
use crate::chemistry::model::Peptide;

#[derive(Clone, Eq, PartialEq, Debug)]
pub enum ActivationType {
    CID,
    ECD,
    ETD,
    HCD,
    PSD,
}

impl std::fmt::Display for ActivationType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub enum MsAnalyzer {
    FTMS,
    TRAP,
}

impl std::fmt::Display for MsAnalyzer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Eq, PartialEq, Debug)] //, strum_macros::Display
pub enum FragmentIonSeries {
    //M,
    //M_H2O,
    //M_NH3,
    a,
    //#[strum(serialize = "a-H2O")]
    a_H2O,
    //#[strum(serialize = "a-NH3")]
    a_NH3,
    b,
    //#[strum(serialize = "b-H2O")]
    b_H2O,
    //#[strum(serialize = "b-NH3")]
    b_NH3,
    c,
    c_dot,
    c_m1,
    c_p1,
    c_p2,
    c_H2O,
    c_NH3,
    d,
    v,
    w,
    x,
    x_H2O,
    x_NH3,
    y,
    //#[strum(serialize = "y-H2O")]
    y_H2O,
   //#[strum(serialize = "y-NH3")]
    y_NH3,
    ya,
    yb,
    z,
    z_H2O,
    z_NH3,
    z_dot,
    z_p1,
    z_p2,
    z_p3,
    immonium,
}


pub fn get_ion_mono_mass_shift(ion_type: FragmentIonSeries) -> f64 {

    // TODO: put in constants
    let h_mono_mass = BIOMOLECULE_ATOM_TABLE.atom_by_symbol.get("H").unwrap().mono_mass();
    let n_mono_mass = BIOMOLECULE_ATOM_TABLE.atom_by_symbol.get("N").unwrap().mono_mass();
    let o_mono_mass = BIOMOLECULE_ATOM_TABLE.atom_by_symbol.get("O").unwrap().mono_mass();

    use FragmentIonSeries::*;
    let mass_shift = match ion_type {
        //M => 0.0,
        //M_H2O => - WATER_MONO_MASS,
        //M_NH3 => - NH3_MONO_MASS,
        a =>  - (WATER_MONO_MASS + CO_MONO_MASS), //Composition(formula='H-2O-1' + 'C-1O-1'),
        a_H2O => - (2.0 * WATER_MONO_MASS + CO_MONO_MASS), //Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
        a_NH3 => - (CO_MONO_MASS + NH3_MONO_MASS + WATER_MONO_MASS), //Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
        b => - WATER_MONO_MASS, //Composition(formula='H-2O-1'),
        b_H2O => -2.0 * WATER_MONO_MASS, //Composition(formula='H-2O-1' + 'H-2O-1'),
        b_NH3 => - (WATER_MONO_MASS + NH3_MONO_MASS), //Composition(formula='H-2O-1' + 'N-1H-3'),
        c => - WATER_MONO_MASS + NH3_MONO_MASS,
        c_dot => - WATER_MONO_MASS + NH3_MONO_MASS + PROTON_MASS,
        c_m1 => - WATER_MONO_MASS + NH3_MONO_MASS - PROTON_MASS,
        c_p1 => - WATER_MONO_MASS + NH3_MONO_MASS + PROTON_MASS,
        c_p2 => - WATER_MONO_MASS + NH3_MONO_MASS + 2.0 * PROTON_MASS,
        c_H2O => - 2.0 * WATER_MONO_MASS + NH3_MONO_MASS,
        c_NH3 => - WATER_MONO_MASS,
        d => 0.0, // FIXME
        v => 0.0, // FIXME
        w => 0.0, // FIXME
        x => - WATER_MONO_MASS + CO2_MONO_MASS,
        x_H2O => - 2.0 * WATER_MONO_MASS + CO2_MONO_MASS,
        x_NH3 => - WATER_MONO_MASS - NH3_MONO_MASS + CO2_MONO_MASS,
        y => 0.0,
        y_H2O => - WATER_MONO_MASS,
        y_NH3 => - NH3_MONO_MASS,
        ya => 0.0, // FIXME
        yb => 0.0, // FIXME
        z => - WATER_MONO_MASS + o_mono_mass - n_mono_mass,
        z_H2O => - 2.0 * WATER_MONO_MASS + o_mono_mass - n_mono_mass - h_mono_mass,
        z_NH3 => - WATER_MONO_MASS - NH3_MONO_MASS + o_mono_mass - n_mono_mass - h_mono_mass,
        z_dot => - WATER_MONO_MASS + o_mono_mass - n_mono_mass, // TODO: check me
        z_p1 => - WATER_MONO_MASS + o_mono_mass - n_mono_mass + h_mono_mass,
        z_p2 => - WATER_MONO_MASS + o_mono_mass - n_mono_mass + 2.0 * h_mono_mass,
        z_p3 => - WATER_MONO_MASS + o_mono_mass - n_mono_mass + 3.0 * h_mono_mass,
        immonium => 0.0 // FIXME
    };


    mass_shift

    /*
        'z':        Composition(formula='H-2O-1' + 'ON-1H-1'),
    'z-dot':    Composition(formula='H-2O-1' + 'ON-1'),
    'z+1':      Composition(formula='H-2O-1' + 'ON-1H1'),
    'z+2':      Composition(formula='H-2O-1' + 'ON-1H2'),
    'z+3':      Composition(formula='H-2O-1' + 'ON-1H3'),
    'z-H2O':    Composition(formula='H-2O-1' + 'ON-1H-1' + 'H-2O-1'),
    'z-NH3':    Composition(formula='H-2O-1' + 'ON-1H-1' + 'N-1H-3'),
            'M':        Composition(formula=''),
    'M-H2O':    Composition(formula='H-2O-1'),
    'M-NH3':    Composition(formula='N-1H-3'),
    'a':        Composition(formula='H-2O-1' + 'C-1O-1'),
    'a-H2O':    Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
    'a-NH3':    Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
    'b':        Composition(formula='H-2O-1'),
    'b-H2O':    Composition(formula='H-2O-1' + 'H-2O-1'),
    'b-NH3':    Composition(formula='H-2O-1' + 'N-1H-3'),
    'c':        Composition(formula='H-2O-1' + 'NH3'),
    'c-1':      Composition(formula='H-2O-1' + 'NH3' + 'H-1'),
    'c-dot':    Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+1':      Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+2':      Composition(formula='H-2O-1' + 'NH3' + 'H2'),
    'c-H2O':    Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1'),
    'c-NH3':    Composition(formula='H-2O-1'),
    'x':        Composition(formula='H-2O-1' + 'CO2'),
    'x-H2O':    Composition(formula='H-2O-1' + 'CO2' + 'H-2O-1'),
    'x-NH3':    Composition(formula='H-2O-1' + 'CO2' + 'N-1H-3'),
    'y':        Composition(formula=''),
    'y-H2O':    Composition(formula='H-2O-1'),
    'y-NH3':    Composition(formula='N-1H-3'),

*/
}
// --- function determines ion_types (forward, reverse or none (to get rid of crashing) --- //
pub fn is_ion_forward(ion_type: FragmentIonSeries) -> Option<bool> {
    use FragmentIonSeriesDirection::*;

    match get_ion_series_direction(ion_type) {
        FORWARD => Some(true),
        REVERSE => Some(false),
        NONE => None
    }
}

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Eq, PartialEq, Debug)] //, strum_macros::Display
// --- create enumeration to define ion_types direction --- //
pub enum FragmentIonSeriesDirection {
    FORWARD,
    REVERSE,
    NONE,
}

//adv using match: if there is a change in fragment ion series, it warns you about it.
pub fn get_ion_series_direction(ion_type: FragmentIonSeries) -> FragmentIonSeriesDirection {

    use FragmentIonSeries::*;
    use FragmentIonSeriesDirection::*;

    let ion_series_direction = match ion_type {
        a => FORWARD,
        a_H2O => FORWARD,
        a_NH3 => FORWARD,
        b => FORWARD,
        b_H2O => FORWARD,
        b_NH3 => FORWARD,
        c => FORWARD,
        c_dot => FORWARD,
        c_m1 => FORWARD,
        c_p1 => FORWARD,
        c_p2 => FORWARD,
        c_H2O => FORWARD,
        c_NH3 => FORWARD,
        d => NONE, // FIXME
        v => NONE, // FIXME
        w => NONE, // FIXME
        x => REVERSE,
        x_H2O => REVERSE,
        x_NH3 => REVERSE,
        y => REVERSE,
        y_H2O => REVERSE,
        y_NH3 => REVERSE,
        ya => REVERSE,
        yb => REVERSE,
        z => REVERSE,
        z_H2O => REVERSE,
        z_NH3 => REVERSE,
        z_dot => REVERSE,
        z_p1 => REVERSE,
        z_p2 => REVERSE,
        z_p3 => REVERSE,
        immonium => NONE,
    };

    ion_series_direction
}

impl std::fmt::Display for FragmentIonSeries {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        use FragmentIonSeries::*;

        // TODO: update me
        match self {
            a => write!(f, "a"),
            a_NH3 => write!(f, "a-NH3"),
            a_H2O => write!(f, "a-H2O"),
            b_NH3 => write!(f, "b-NH3"),
            b_H2O => write!(f, "a-H2O"),
            c => write!(f, "c"),
            d => write!(f, "d"),
            v => write!(f, "v"),
            w => write!(f, "w"),
            x => write!(f, "x"),
            y => write!(f, "y"),
            y_NH3 => write!(f, "y-NH3"),
            y_H2O => write!(f, "y-H2O"),
            ya => write!(f, "ya"),
            yb => write!(f, "yb"),
            z => write!(f, "z"),
            z1 => write!(f, "z+1"),
            z2 => write!(f, "z+2"),
            immonium => write!(f, "immonium"),
        }


    }
}

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub enum FragmentType {
    IMMONIUM,
    INTERNAL,
    SATELLITE,
    SEQUENCE,
}

impl std::fmt::Display for FragmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        use FragmentType::*;

        match self {
            IMMONIUM => write!(f, "IM"),
            INTERNAL => write!(f, "IN"),
            SATELLITE => write!(f, "SAT"),
            SEQUENCE => write!(f, "SEQ"),
        }
    }
}

#[derive(Clone, PartialEq, Debug)]
pub struct FragmentIonType {
    pub ion_series: FragmentIonSeries,
    pub neutral_loss: Option<NeutralLoss>,
    pub is_forward_ion: bool
}

#[derive(Clone, Default, PartialEq, Debug)]
struct FragmentMatch {
    pub label: String,
    pub r#type: Option<String>, // = None,
    pub moz: f64,
    pub calculated_moz: f64,
    pub intensity: f32,
    pub neutral_loss_mass: Option<f64> // = None
}

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub enum NeutralLoss {
    CH4OS,
    H2O,
    H3PO4,
    HPO3,
    NH3,
}

impl std::fmt::Display for NeutralLoss {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Clone, PartialEq, Debug)]
struct TheoreticalFragmentIon {
    pub position: i32,
    pub moz: f64,
    pub nl_mass: f64, // = 0.0,
    pub charge: i8, //= 1,
    pub fragment_type: FragmentType, // = FragmentType.SEQUENCE,
    pub frag_series: Option<String>, // = None,
    pub is_reverse_series: bool,
    pub is_alternative_nl: bool // = false
}

#[derive(Clone, Default, PartialEq, Debug)]
struct FragmentIonTable {
    peptide: Peptide,
    charge_by_ion_type: std::collections::HashMap<char,i8>,
    //sequence: Option<Vec<char>>, // = None,
    ptm_neutral_losses: Option<std::collections::HashMap<i32,f64>> // NL mass by seq position
}
