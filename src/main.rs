
//use chemical_elements::{ChemicalComposition, ElementSpecification};

mod chemistry;
mod io;
mod ms;
mod msms;
mod phospho_test;


#[macro_use]
extern crate lazy_static;

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
    phospho_test::test_phospho_analysis(true,1).expect("error");

}


