
use anyhow::*;
use std::collections::HashMap;

use crate::chemistry::composition::*;
use crate::chemistry::table::AminoAcidTable;
use crate::chemistry::constants::{WATER_MONO_MASS, WATER_AVERAGE_MASS};

pub fn calc_aa_seq_mass(aa_seq: &str, aa_table: &AminoAcidTable, mono_mass: bool) -> Result<f64> {

    let aa_comp = parse_aa_composition(aa_seq)?;

    let get_aa_mass = |aa_code1: char| -> Result<f64> {
        let aa = aa_table.aa_by_code1.get(&aa_code1).ok_or_else(
            || anyhow!("can't find amino acid '{}' in the provided table",aa_code1)
        )?;
        let m = if mono_mass {aa.mono_mass} else {aa.average_mass};
        Ok(m)
    };

    let seq_mass = _calc_mass(aa_comp, get_aa_mass)?;

    if mono_mass {
        Ok(seq_mass + WATER_MONO_MASS)
    } else {
        Ok(seq_mass + WATER_AVERAGE_MASS)
    }
}

fn _calc_mass<T,F>(abundance_map: HashMap<T, f32>, get_entity_mass: F) -> Result<f64> where F: Fn(T) -> Result<f64> {

    let mut mass: f64 = 0.0;
    for (entity, entity_ab) in abundance_map {
        let entity_mass = get_entity_mass(entity)?;
        mass += entity_ab as f64 * entity_mass;
    }

    Ok(mass)
}