use crate::msms::model::*;

use anyhow::*;
use std::collections::HashMap;
use crate::{AminoAcidTable, calc_aa_seq_mass, mass_to_mz, STANDARD_AMINO_ACID_TABLE};

#[derive(Clone, PartialEq, Debug)]
pub struct SimpleFragmentationRule {
    pub fragment_ion_type: FragmentIonType,
    pub fragment_ion_charge: i8, // = 1
}

#[derive(Clone, Debug)]
pub struct SimpleFragmentationConfig {
    pub frag_rules: Vec<SimpleFragmentationRule>,
    pub frag_rule_by_ion_series: HashMap<FragmentIonSeries,SimpleFragmentationRule>
}

impl SimpleFragmentationConfig {

    /*pub fn create_with_defaults() -> Result<SimpleFragmentationConfig> {
        create(1 as i8, ActivationType::HCD, MsAnalyzer::FTMS)
    }*/

    /*pub fn create(
        fragment_ion_charge: i8,
        activation_type: ActivationType,
        msn_analyzer: MsAnalyzer
    ) -> Result<SimpleFragmentationConfig> {

        use ActivationType::*;
        use FragmentIonSeries::*;
        use MsAnalyzer::*;

        let ion_series_list = match activation_type  {
            CID => vec!(b,b_NH3,b_H2O,y,y_NH3,y_H2O),
            ECD => vec!(c,y,z_p1,z_p2),
            ETD => if msn_analyzer == FTMS {vec!(c,y,z,z_p1,z_p2)} else {vec!(c,y,z_p1,z_p2)},
            HCD => vec!(a,a_NH3,a_H2O,b,b_NH3,b_H2O,y,y_NH3,y_H2O,ya,yb),
            PSD => vec!(a,a_NH3,a_H2O,b,b_NH3,b_H2O,y),
        };

        /*logger.debug(
            s"Fragmentation config for activationType=$activationType and msnAnalyzer=$msnAnalyzer will contain this list of fragment series: " + ionSeriesList.mkString(",")
        )*/

        // Map fragmentation rules by series
        let mut frag_rule_by_ion_series:  HashMap<FragmentIonSeries,SimpleFragmentationRule> = HashMap::new();

        let simple_frag_rules = ion_series_list.iter().map( |ion_series| {

            let frag_rule = SimpleFragmentationRule {
                fragment_ion_type: FragmentIonType { ion_series: (*ion_series).to_owned(), neutral_loss: None },
                fragment_ion_charge: fragment_ion_charge
            };

            let frag_ion_series = frag_rule.fragment_ion_type.ion_series;
            if frag_rule_by_ion_series.contains_key(&frag_ion_series) {
                bail!("Duplicated SimpleFragmentationRule found, the same ionSeries {} has been already defined",frag_ion_series)
            }

            frag_rule_by_ion_series.insert(frag_ion_series, &*frag_rule);

            frag_rule
        }).collect();

        Ok(SimpleFragmentationConfig {
            frag_rules: simple_frag_rules,
            frag_rule_by_ion_series: frag_rule_by_ion_series,
        })
    }*/

    /*pub fn  containsIonSeries(&self, ionSeries: FragmentIonSeries.Value) -> bool {
        frag_rule_by_ion_series.contains(ionSeries)
    }
    pub fn  getFragmentationRule(&self, ionSeries: FragmentIonSeries.Value) -> Option[SimpleFragmentationRule] {
    fragRuleByIonSeries.get(ionSeries)
    }
    pub fn  getRequiredFragmentIonCharge(&self, ionSeries: FragmentIonSeries.Value): Option[Int] = {
    fragRuleByIonSeries.get(ionSeries).map(_.fragmentIonCharge)
    }*/
}

// --- Compute m/z value of fragment ion series (it contains all ion_types into one vector)  --- //
pub fn compute_frag_series_mz_values(pep_seq: &str, ion_type: FragmentIonSeries, charge: i8, aa_table: &AminoAcidTable) -> Result<Vec<f64>> {

    let seq_len = pep_seq.len();
    let mut frag_series_mz_values = Vec::with_capacity(seq_len);
    //let t = &*STANDARD_AMINO_ACID_TABLE;


    // --- This function contains ion series directions itself --- //
    use FragmentIonSeriesDirection::*;
    match get_ion_series_direction(ion_type) {
        FORWARD => {
            // Forward ions loop
            for i in 1..seq_len { // starts at one because pep_seq slicing range will be between 0 and i-1
                let pep_subset_seq = &pep_seq[0..i];
                frag_series_mz_values.push(_calc_ion_mz(pep_subset_seq, aa_table, ion_type, charge));
                //println!("{} {:.4}", fion_aa, bmz);
            }

            Ok(frag_series_mz_values)
        },
        REVERSE => {
            // Reverse ions loop
            for i in (1..seq_len).rev() {
                // create slice (subset reference) of peptide sequence
                let pep_subset_seq = &pep_seq[i .. seq_len];
                frag_series_mz_values.push(_calc_ion_mz(pep_subset_seq, aa_table, ion_type, charge));
            }

            Ok(frag_series_mz_values)
        },
        NONE => {
            bail!("unsupported ion type")
        }
    }

}
// --- Calculate mz value of a given peptide sequence depending on ion type and charge state --- //
fn _calc_ion_mz(pep_subset_seq: &str, aa_table: &AminoAcidTable, ion_type: FragmentIonSeries, charge: i8) -> f64 {
    // calculate mass of peptide sequence subset
    let pep_subset_mass = calc_aa_seq_mass(pep_subset_seq, aa_table, true).unwrap();
    // calculate mass of fragment ion
    let ion_mass = pep_subset_mass + get_ion_mono_mass_shift(ion_type);
    // convert mass to m/z value
    let ion_mz = mass_to_mz(ion_mass, charge as i32);
    //println!("{} {:.4}", rion_aa, ymz);

    ion_mz
}

#[derive(Clone, PartialEq, Debug)]
pub struct TheoreticalFragmentIons {
    pub ion_type: FragmentIonSeries,
    pub charge: i8,
    pub mz_values: Vec<f64>
}

pub type FragmentationTable = Vec<TheoreticalFragmentIons>;
// --- by this function, we also add ion type, and charge state information additionally to this vector to distinguish between them. --- //
pub fn compute_fragmentation_table(pep_seq: &str, ion_types: &[FragmentIonSeries], frag_ion_charges: &Vec<i8>, aa_table: &AminoAcidTable) -> Result<FragmentationTable> {

    let mut frag_table: FragmentationTable = Vec::with_capacity(ion_types.len());
//for each charge state, ion_types should be recalculated so that "for loop of charge" into the "for loop of ion_type" --- //
// --- (e.g. b+1 for charge=1 and b+1 for charge=2) --- //
    for ion_type in ion_types {

        for charge in frag_ion_charges {

            let mz_values_res = compute_frag_series_mz_values(pep_seq, *ion_type, *charge, aa_table);
            let mz_values = mz_values_res?;

            // add to fragmentation table a new column containing different mz values for considered ion type and charge state
            frag_table.push(TheoreticalFragmentIons {
                ion_type: *ion_type,
                charge: *charge,
                mz_values: mz_values.clone()
            });
        }
    }
    //FragmentationTable
    Ok(frag_table)
}
