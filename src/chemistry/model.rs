#![allow(dead_code, non_camel_case_types)]

use anyhow::*;

#[derive(Clone, Default, PartialEq, Debug)]
struct MolecularEntity {
    // require(StringUtils.isNotEmpty(symbol), "symbol is empty")
    pub symbol: String,
    // require(name != null, "name is null")
    pub name: String,
    // require(monoMass > 0, "mono_mass must be a strictly positive number")
    pub mono_mass: f64,
    // require(averageMass > 0, "average_mass must be a strictly positive number")
    pub average_mass: f64
}

impl MolecularEntity {
    pub fn new(symbol: &str, name: &str, mono_mass: f64, average_mass: f64) -> anyhow::Result<MolecularEntity> {
        if symbol.is_empty() { bail!("symbol is empty") }
        if name.is_empty() { bail!("name is empty") }
        if mono_mass <= 0.0 { bail!("mono_mass must be a strictly positive number") }
        if average_mass <= 0.0 { bail!("average_mass must be a strictly positive number") }

        Ok(MolecularEntity {
            symbol: symbol.to_string(),
            name: name.to_string(),
            mono_mass: mono_mass,
            average_mass: average_mass,
        })
    }
}
/*trait IMolecule extends IMolecularEntity {

val formula: String

// TODO: require it is matching the regex /^(\w{1}(\(-*\d+\))*\s*)+$/

def getAtomComposition(atomTable: AtomTableLike): AtomComposition = {
new AtomComposition(this.formula, atomTable)
}

}

trait IPolymer extends IMolecularEntity {
def sequence: String
}*/

#[derive(Clone, Default, PartialEq, Debug)]
pub struct AminoAcidResidue {
    pub code1: char,
    pub code3: String,
    pub name: String,
    pub formula: Option<String>,
    pub mono_mass: f64,
    pub average_mass: f64,
    pub occurrence: f32, // = 0f, // occurrence in human proteins
    pub pka1: f32, // = 0, // C-term pKa
    pub pka2: f32, // = 0, // N-term pKa
    pub pka3: f32, // = 0, // side chain pKa
    pub pi: f32, // = 0,   // isoelectric point
    pub codons: Vec<String>
}

const ACUG_STR: &str = "ACUG";

impl AminoAcidResidue {
    pub fn new(
        code1: char,
        code3: &str,
        name: &str,
        formula: &str,
        mono_mass: f64,
        average_mass: f64,
        codons: Vec<String>
    ) -> anyhow::Result<AminoAcidResidue> {

        if code3.len() < 3 { bail!("code3 must contain three characters") }
        if name.is_empty() { bail!("name is empty") }
        if formula.is_empty() { bail!("formula is empty") }

        for codon in &*codons {
            if codon.len() != 3 { bail!("a codon must contain three characters") }
            let only_acug_chars = codon.chars().all(|c| ACUG_STR.contains(c));
            if !only_acug_chars { bail!("a codon must only contain ACUG letters") }
        }

        if mono_mass <= 0.0 { bail!("mono_mass must be a strictly positive number") }
        if average_mass <= 0.0 { bail!("average_mass must be a strictly positive number") }

        Ok(AminoAcidResidue {
            code1: code1,
            code3: code3.to_string(),
            name: name.to_string(),
            formula: Some(formula.to_string()),
            mono_mass: mono_mass,
            average_mass: average_mass,
            occurrence: 0.0,
            pka1: 0.0,
            pka2: 0.0,
            pka3: 0.0,
            pi: 0.0,
            codons: codons //.iter().map(|s| *s.to_string()).collect()
        })
    }
}

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub enum PtmLocation {
    PROT_N_TERM, // = Value("Protein N-term")
    PROT_C_TERM, // = Value("Protein C-term")
    ANY_N_TERM, // = Value("Any N-term")
    ANY_C_TERM, // = Value("Any C-term")
    ANYWHERE, // = Value("Anywhere")
}

impl PtmLocation {
    fn new(location: &str) -> Option<PtmLocation> {
        match location {
            "Protein N-term" => Some(PtmLocation::PROT_N_TERM),
            "Protein C-term" => Some(PtmLocation::PROT_C_TERM),
            "Any N-term"     => Some(PtmLocation::ANY_N_TERM),
            "Any C-term"     => Some(PtmLocation::ANY_C_TERM),
            "Anywhere"       => Some(PtmLocation::ANYWHERE),
            _                => None
        }
    }
}

impl std::fmt::Display for PtmLocation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        use PtmLocation::*;

        match self {
            PROT_N_TERM => write!(f, "Protein N-term"),
            PROT_C_TERM => write!(f, "Protein C-term"),
            ANY_N_TERM  => write!(f, "Any N-term"),
            ANY_C_TERM  => write!(f, "Any C-term"),
            ANYWHERE    => write!(f, "Anywhere"),
        }
    }
}

#[derive(Clone, PartialEq, Debug)]
struct AminoAcidPtm {
    pub id: i64,
    pub name: String,
    pub formula: String,
    pub mono_mass: f64,
    pub average_mass: f64,
    pub position_constraint: PtmLocation, // any N-term, any C-term, protein N-term, protein C-term
    pub residue_constraint: char
}

impl AminoAcidPtm {
    pub fn new(id: i64, name: &str, formula: &str, mono_mass: f64, average_mass: f64, position_constraint: PtmLocation, residue_constraint: char) -> anyhow::Result<AminoAcidPtm> {
        if name.is_empty() { bail!("name is empty") }
        if formula.is_empty() { bail!("formula is empty") }

        if mono_mass <= 0.0 { bail!("mono_mass must be a strictly positive number") }
        if average_mass <= 0.0 { bail!("average_mass must be a strictly positive number") }

        Ok(AminoAcidPtm {
            id: id,
            name: name.to_string(),
            formula: formula.to_string(),
            mono_mass: mono_mass,
            average_mass: average_mass,
            position_constraint: position_constraint,
            residue_constraint: residue_constraint,
        })
    }
}

// The atomic_number uniquely identifies an element
#[derive(Clone, Default, PartialEq, Debug)]
pub struct Atom {
    pub atomic_number: u16,
    pub symbol: String,
    pub name: String,
    pub isotopes: Vec<Isotope>,
}

impl Atom {
    pub fn new(
        atomic_number: u16,
        symbol: &str,
        name: &str,
        isotopes: Vec<Isotope>,
    ) -> anyhow::Result<Atom> {
        if symbol.is_empty() { bail!("symbol is empty") }
        if name.is_empty() { bail!("name is empty") }
        if isotopes.is_empty() { bail!("isotopes is empty") }

        Ok(Atom {
            atomic_number: atomic_number,
            symbol: symbol.to_string(),
            name: name.to_string(),
            isotopes: isotopes,
        })
    }

    pub fn proton_number(&self) -> u16 { self.atomic_number }

    pub fn get_neutron_number(&self, isotope_idx: usize) -> u16 {
        self.isotopes[isotope_idx].get_neutron_number(self.proton_number())
    }

    pub fn mono_mass(&self) -> f64 { self.isotopes.first().unwrap().mass }

    pub fn calc_average_mass(&self) -> f64 {
        let mut weighted_mass_sum: f64 = 0.0;
        let mut weight_sum: f64 = 0.0;

        for iso in self.isotopes.iter() {
            let ab = iso.abundance as f64;
            weighted_mass_sum += iso.mass * ab;
            weight_sum += ab;
        }

        weighted_mass_sum / weight_sum
    }

    pub fn to_isotopic_variants(&self) -> Vec<AtomIsotopicVariant> {
        self.isotopes.iter().map(|isotope| AtomIsotopicVariant::new(self.clone(), isotope.clone())).collect()
    }

}

#[derive(Clone, Default, PartialEq, Debug)]
pub struct AtomIsotopicVariant { atom: Atom, isotope: Isotope }

impl AtomIsotopicVariant {
    pub fn new(atom: Atom, isotope: Isotope) -> AtomIsotopicVariant {
        AtomIsotopicVariant {
            atom: atom,
            isotope: isotope,
        }
    }
}

#[derive(Clone, Default, PartialEq, Debug)]
pub struct Isotope { mass_number: u16, mass: f64, abundance: f32 }

impl Isotope {
    pub fn new(mass_number: u16, mass: f64, abundance: f32) -> anyhow::Result<Isotope> {
        if mass_number <= 0 { bail!("mass_number must be a strictly positive number") }
        if mass <= 0.0 { bail!("mass must be a strictly positive number") }
        if abundance < 0.0 { bail!("abundance must be a positive number") }

        Ok(Isotope {
            mass_number: mass_number,
            mass: mass,
            abundance: abundance,
        })
    }

    pub fn nucleon_number(&self) -> u16 { self.mass_number }

    pub fn get_neutron_number(&self, proton_number: u16) -> u16 {
        self.mass_number - proton_number
    }

}

#[derive(Clone, Default, PartialEq, Debug)]
pub struct Peptide {
    sequence: String,
    mods: Vec<(i64, i32)>, // (ptm_id, seq_position)
    mono_mass: f64,
    average_mass: f64
}

impl Peptide {
    pub fn new(
        sequence: &str,
        mods: Vec<(i64, i32)>,
        mono_mass: f64,
        average_mass: f64
    ) -> anyhow::Result<Peptide> {

        if sequence.is_empty() { bail!("sequence is empty") }
        if mono_mass <= 0.0 { bail!("mono_mass must be a strictly positive number") }
        if average_mass <= 0.0 { bail!("average_mass must be a strictly positive number") }

        Ok(Peptide {
            sequence: sequence.to_string(),
            mods: mods,
            mono_mass: mono_mass,
            average_mass: average_mass,
        })
    }
}
