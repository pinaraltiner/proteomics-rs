
//use chemical_elements::{ChemicalComposition, ElementSpecification};

mod chemistry;
mod io;
mod ms;
mod msms;
mod phospho_test;


#[macro_use]
extern crate lazy_static;

use std::path::Path;
use anyhow::*;
use itertools::Itertools;
//use mzdata::{CentroidPeak, MGFReader, ParamDescribed};

use chemistry::constants::atom;
use chemistry::composition::*;
use chemistry::mass_calc::*;
use chemistry::table::*;
use msms::model::*;
//use test::bench::iter;
use crate::ms::utils::mass_to_mz;
use crate::msms::annotator::{annotate_spectrum, PeakSelectionStrategy};
use crate::msms::fragmentation::{compute_frag_series_mz_values, compute_fragmentation_table};
use crate::msms::io::read_mgf_file;
use crate::phospho_test::SpectrumPhosphoAnalysis;


fn main () {

    //phospho_test::test_phospho_analysis_on_example_data()?;
    phospho_test::test_phospho_analysis("./data/benchmarks/MGF/mgf_dataset_20220407.bin","./data/benchmarks/test_docker",["./data/benchmarks/phospho_evidence_matrices_list/results_correct_pos_from_eyers_pools/1comb_correct_pos_proline_no_isomer_pep_list_extracted_from_pool1-5.tsv"],true,1).expect("error");


}


