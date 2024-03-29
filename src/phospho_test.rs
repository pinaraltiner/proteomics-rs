use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;

use anyhow::*;
use itertools::Itertools;
use serde::{Serialize, Deserialize};

//use mzdata::{CentroidPeak, MGFReader, ParamDescribed};
use crate::io::mgf::*;

use crate::chemistry::constants::atom;
use crate::chemistry::composition::*;
use crate::chemistry::mass_calc::*;
use crate::chemistry::table::*;
use crate::msms::model::*;
//use test::bench::iter;
use crate::ms::utils::mass_to_mz;
use crate::msms::annotator::*;
use crate::msms::fragmentation::*;

#[derive(Clone, PartialEq, Debug, Deserialize)]
pub struct PsmRecord {
    pub peptide_id: u64,
    pub calculated_mass: f64,
    pub charge: u8,
    pub sequence: String,
    pub modifications: String,
    pub spectrum_title: String,
    pub peaklist_name: String,
    //pub phospho_pos: usize
}
// --- This function read the tsv file that contains pep_seq, modification (filtered only "phospho") and spectrum title. --- //

fn _read_psms_from_tsv_files(tsv_file_paths: &[&str]) -> Result<Vec<PsmRecord>> {
    let mut all_psms: Vec<PsmRecord> = Vec::with_capacity(1000);

    for tsv_file_path in tsv_file_paths {
        let mut tmp_psm_records = _read_psms_from_tsv(tsv_file_path)?;

        all_psms.append(&mut tmp_psm_records);
    }

    Ok(all_psms)
}

fn _read_psms_from_tsv(file_path: &str) -> Result<Vec<PsmRecord>> {

    //let file = File::open(file_path)?;
    //let mut rdr = csv::Reader::from_reader(file);

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .double_quote(false)
        .flexible(true)
        .escape(None)
        .from_path(file_path)?;

    let mut psms = Vec::with_capacity(100);

    for record_res in rdr.deserialize() {

        let psm: PsmRecord = record_res?;
        if psm.sequence.contains('X') == false { // TODO: make this configurable
            psms.push(psm);
        }

    }
    Ok(psms)
}

#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Spectrum {
    peaklist_name: String,
    title: String,
    precursor_mz: f64,
    charge: Option<i8>,
    //retention_time: f32, // TODO: add this
    peaks: Vec<[f64;2]>, // could also be Vec<(f64,f64)>
}

// --- This function reads the mgf file and we create a struct to put all (title, charge, pep_mass,etc) info into it. --- //
// --- Spectrum struct has MS/MS data (raw files means all spectra in one run (in one mgf file).) --- //
/*fn _load_spectra(file_path: &str) -> Result<Vec<Spectrum>> {

    // See: https://github.com/mobiusklein/mzdata/blob/4050dd35d53ac461f5c71f3f7839f88ef40eb506/src/io/mgf.rs#L728
    use std::fs;
    use std::path;
    use core::result::Result::Ok;

    let path = path::Path::new(file_path);
    let file = fs::File::open(path).expect("Test file doesn't exist");

    let reader = MGFReader::new(file); //MGFReaderType::<_, CentroidPeak, DeconvolutedPeak>::new(file);

    let mut spectra = Vec::with_capacity(1000);
    for scan in reader {

        let scan_description = scan.description;
        /*for p in scan_description.to_owned().params {
            println!("{}", p.name);
        }*/

        let scan_desc_ref = &scan_description;

        let title = (*scan_desc_ref).id.to_string();

        let prec_charge_opt = scan_desc_ref.get_param_by_name("charge").map(|charge_param| {
            let charge_str= charge_param.to_owned().value.chars().take_while(|char| char.is_ascii_digit()).collect::<String>();
            let charge = charge_str.parse::<i8>().unwrap();
            charge
        });

        let prec_mz_opt = (*scan_desc_ref).precursor.as_ref().map(|prec| {prec.ion.mz});

        let peaks = scan.peaks.unwrap().peaks;

        let peaks_as_arrays: Vec<[f64;2]> = peaks.iter().map(|peak| {
            //peaks_as_arrays.push([peak.mz,peak.intensity as f64])
            [peak.mz,peak.intensity as f64]
        }).collect();

        let spectrum = Spectrum {
            title: title,
            precursor_mz: prec_mz_opt.unwrap_or(0.0),
            charge: prec_charge_opt,
            peaks: peaks_as_arrays
        };

        spectra.push(spectrum);
    } // ends for scan in reader

    Ok(spectra)
}*/

fn _load_spectra(file_path: &str) -> Result<Vec<Spectrum>> {

    use std::fs;
    use std::path;

    let path = path::Path::new(file_path);
    let file = fs::File::open(path).expect("Test file doesn't exist");
    let peaklist_name = path.file_stem().unwrap().to_string_lossy();

    let mut spectra = Vec::with_capacity(1000);
    for_each_spectrum(&*path.to_string_lossy(), |spectrum| {

        let peaks_as_arrays = spectrum.data.mz_list.iter()
            .zip(spectrum.data.intensity_list)
            .map(|(&mz,int)| [mz, int as f64]).collect_vec();

        let spectrum = Spectrum {
            peaklist_name: peaklist_name.to_string(),
            title: spectrum.header.title,
            precursor_mz: spectrum.header.precursor_mz,//prec_mz_opt.unwrap_or(0.0),
            charge: spectrum.header.precursor_charge,
            peaks: peaks_as_arrays
        };

        spectra.push(spectrum);

        ()
    })?;

    Ok(spectra)
}
// --- This function is roughly done by each step in planned workflow --- //
    // From matching btw MS/MS output and their PSMs, calculating theo frag mz (with/out including phospho mass shift)
    // and comparing them to exp frag mz then, extracting intensity info of matched peaks. Finally writing into separate tsv file

pub fn test_phospho_analysis(mgf_dataset_path:&str ,final_path: &str,input_tsv_files: [&str;1], int_count: bool,phospho_num: usize) -> Result<()> {


    let aa_table = &*PROTEINOGENIC_AMINO_ACID_TABLE;

    //let output_path = Path::new("./data/benchmarks/example");
    let output_path = Path::new(&final_path);
    //phospho_evidence_matrices_list/results_correct_pos_from_eyers_pools/eliminate_duplicate_ions/two_int_count/
    if output_path.is_dir() == false {
        std::fs::create_dir_all(output_path)?;
    }

    let mz_error_tol = 0.02;
    let phospho_residues = ['S','T','Y'];
    // Only use for testing
    //let mgf_dataset_path = "./data/rim_dataset_NS_vs_45/mgf_dataset_20220215.bin"; //not exist!

    // --- Rim's input files --- //
        //Rim's Actual data

    //let mgf_dataset_path = "./data/rim_dataset_NS_vs_45/mgf_dataset_20220216.bin";
    /*
    let input_mgf_files = [
        "./data/rim_dataset_NS_vs_45/OXRIA211021_03_CV45_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211021_03_CV60_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211021_09_CV45_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211021_09_CV60_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211021_13_CV45_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211021_13_CV60_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211021_19_CV45_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211021_19_CV60_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211117_10_CV45_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211117_10_CV60_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211117_45_CV45_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211117_45_CV60_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211117_51_CV45_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211117_51_CV60_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211117_55_CV45_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211117_55_CV60_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211117_61_CV45_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA211117_61_CV60_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA241021_02_CV45_mzcal.mzDB.mgf",
        "./data/rim_dataset_NS_vs_45/OXRIA241021_02_CV60_mzcal.mzDB.mgf"
    ];
*/
/*
    let input_tsv_files = [
        "./data/tsv_processing/bio_interesting_psms_NS_vs_45.tsv"
        /*"./data/tsv_processing/OXRIA211021_03_CV45.tsv",
        "./data/tsv_processing/OXRIA211021_03_CV60.tsv",
        "./data/tsv_processing/OXRIA241021_02_CV45.tsv",
        "./data/tsv_processing/OXRIA241021_02_CV60.tsv"*/
    ];
*/
    // --- Eyer's input files --- //
    // Eyer's Actual data

    //let mgf_dataset_path = "./data/benchmarks/MGF/mgf_dataset_20220407.bin";

    let input_mgf_files = [
        "./data/benchmarks/MGF/SF_200217_pPeptideLibrary_pool1_HCDOT_rep1.mgf",
        "./data/benchmarks/MGF/SF_200217_pPeptideLibrary_pool1_HCDOT_rep2.mgf",
        "./data/benchmarks/MGF/SF_200217_pPeptideLibrary_pool2_HCDOT_rep1.mgf",
        "./data/benchmarks/MGF/SF_200217_pPeptideLibrary_pool2_HCDOT_rep2.mgf",
        "./data/benchmarks/MGF/SF_200217_pPeptideLibrary_pool3_HCDOT_rep1.mgf",
        "./data/benchmarks/MGF/SF_200217_pPeptideLibrary_pool3_HCDOT_rep2.mgf",
        "./data/benchmarks/MGF/SF_200217_pPeptideLibrary_pool4_HCDOT_rep1.mgf",
        "./data/benchmarks/MGF/SF_200217_pPeptideLibrary_pool4_HCDOT_rep2.mgf",
        "./data/benchmarks/MGF/SF_200217_pPeptideLibrary_pool5_HCDOT_rep1.mgf",
        "./data/benchmarks/MGF/SF_200217_pPeptideLibrary_pool5_HCDOT_rep2.mgf"
    ];

    //let input_tsv_files = ["./data/benchmarks/phospho_evidence_matrices_list/results_correct_pos_from_eyers_pools/1comb_correct_pos_proline_no_isomer_pep_list_extracted_from_pool1-5.tsv"];


    // ---Code used to cache the spectra --- //

    // In here, we created a vector of specific object (in this case spectra) and was written in a .bin file,
    // now we only use this bin file for playing the data

    /*let spectra = {
        let mut all_spectra: Vec<Spectrum> = Vec::with_capacity(10000);
        for mgf_file in input_mgf_files {
            let mut tmp_spectra = _load_spectra(mgf_file)?;
            println!("N spectra={} in {}", tmp_spectra.len(), mgf_file);

            all_spectra.append(&mut tmp_spectra);
        }

        let encoded_spectra: Vec<u8> = bincode::serialize(&all_spectra).unwrap();

        let mut file = File::create(mgf_dataset_path)?;

        // Write a slice of bytes to the file
        file.write_all(&encoded_spectra[..]);

        /*let mut spectrum_by_title = HashMap::new();
        for spectrum in all_spectra {
            spectrum_by_title.insert(spectrum.to_owned().title, spectrum.to_owned());
        }

        spectrum_by_title*/

        all_spectra
    };*/

    // --- Code used to read spectra from the binary cache --- //

    let spectra = {
        let mut file = File::open(mgf_dataset_path)?;

        // Read the same file back into a Vec of bytes
        let mut encoded_spectra = Vec::<u8>::new();
        file.read_to_end(&mut encoded_spectra)?;

        let decoded_spectra: Vec<Spectrum> = bincode::deserialize(&encoded_spectra[..]).unwrap();

        decoded_spectra
    };

    let mut spectrum_by_unique_title = HashMap::new();
    for spectrum in spectra.iter() {
        // TODO remove me
        //if spectrum.title == "2159" {
            let new_title = spectrum.title.to_owned() + " FROM " + &*spectrum.peaklist_name;
            //println!("new_title={}",new_title);
            spectrum_by_unique_title.insert(
                new_title,
                spectrum.to_owned()
            );
        //}
    }

    let psms = _read_psms_from_tsv_files(&input_tsv_files)?;

    // Group PSMs by sequence, calculated_mass and charge (because  pep_id would be changed for each phospho possibility.)
    // With the combination of all three features, each possibility could be specified uniquely.

    let psms_grouped_by_key = psms.iter()
        //.filter(|psm| { psm.sequence == "GGFFSSFMK"})  // TODO remove me
        .into_group_map_by(|psm| {
        //(psm.sequence.clone(),psm.calculated_mass.to_string(),psm.charge,psm.spectrum_title.clone())
        vec!(psm.sequence.clone(),psm.calculated_mass.to_string(),psm.charge.to_string()).join("|")
    });

    // See for group_by weird behavior (requires sorting):
    // - https://github.com/rust-itertools/itertools/issues/374
    // - https://stackoverflow.com/questions/69662321/how-do-i-perform-a-group-by-in-a-vector-of-structs

    //for ((pep_seq, pep_mass,charge, spectrum_title),peptide_psms_groups) in psms_seq_iter {
    for (group_by_key,peptide_psms_groups) in psms_grouped_by_key {

        let mut pep_ion_spectra = Vec::with_capacity(1000);
        let mut pep_ion_spectra_titles = Vec::with_capacity(1000);

        let mut pep_seq = String::new();
        let mut mods = String::new();
        let mut pep_mass = 0.0;
        let mut charge = 0;

        let mut seen_spectra = HashSet::new();

        // Iterate a group of PSMs that have the same peptide_id
        for (i,psm) in peptide_psms_groups.into_iter().enumerate() {
            /*if psm.sequence == "HTDDEMTGYVATR" && psm.charge == 2 {
                println!("psm={:?}",psm);
            }*/

            // if one psm has one spectrum, directly copied to PSMRecord
            if i == 0 {
                pep_seq = psm.sequence.clone();
                mods = psm.modifications.clone();
                pep_mass = psm.calculated_mass;
                charge = psm.charge;
            }

            // if one psm has either multiple spectra or zero spectrum (that's why option was generated)
            if seen_spectra.contains(&*psm.spectrum_title) == false {
                let unique_title = psm.spectrum_title.to_owned() + " FROM " +  &*psm.peaklist_name;

                let spectrum_opt = spectrum_by_unique_title.get(&*unique_title);
                if spectrum_opt.is_some() {
                    seen_spectra.insert(&*psm.spectrum_title);
                    pep_ion_spectra_titles.push(psm.spectrum_title.clone());

                    let spectrum = spectrum_opt.unwrap();
                    let peaks = spectrum.peaks.iter().map(|peak| *peak);

                    pep_ion_spectra.push(peaks.collect());
                } /*else {
                    println!("can't find spectrum for psm={:?}",psm);
                }*/
            }
        }

        // Splitting all modification using ";" and collect all modification that contains "Phospho"
        let phospho_mods: Vec<&str> = mods.split("; ").filter(|s| s.starts_with("Phospho")).collect();

        // The length of this object is number of phosphorylation for one pep seq
        let n_phosphates = phospho_mods.len();

        // By searching each peptides specifically, all possible phospho-positions (S,T,Y) were extracted.
        let mut phospho_positions = Vec::with_capacity(n_phosphates);
        for (char_idx, char) in pep_seq.chars().into_iter().enumerate() {
            if phospho_residues.contains(&char) {
                phospho_positions.push(char_idx + 1);
            }
        }

        // "phospho_num" was added as a new parameter to the main function (test_phospho_analysis)
        // To adjust which phosphorylation state of phospho-peptide we want to get
        if n_phosphates == phospho_num {

            let phospho_positions_combinations = phospho_positions.to_vec().into_iter().combinations(n_phosphates).collect_vec();
            let phospho_matched_peaks = _match_phospho_peaks(
                pep_seq.as_str(),
                &phospho_positions_combinations,
                n_phosphates,
                &pep_ion_spectra,
                mz_error_tol,
                aa_table
            )?;

            let phospho_evidence_matrix = _compute_phospho_evidence_matrix(
                &phospho_matched_peaks,
                &phospho_positions,
                pep_ion_spectra.len(),
                true
            );
            if pep_seq == "VGSLTPPSSPK" {
                println!("phospho_matched_peaks222={:?}",phospho_evidence_matrix);
            }

            // Verify that pep_ion_spectra and pep_ion_spectra_titles contain the same number of entries
            assert!(pep_ion_spectra.len() == pep_ion_spectra_titles.len());

            let matrix_name = format!("PEPTIDE={} MODS={} CHARGE={} MASS={}",pep_seq, mods, charge, pep_mass);
            let matrix_name_ref = &matrix_name;

            // "int_count" was added as a new parameter to the main function (test_phospho_analysis)
            //To adjust which matched intensity option (summation or count) we want to get

            let mut matrix_as_text_opt = _print_phospho_evidence_matrix(
                matrix_name_ref,
                &phospho_evidence_matrix,
                &phospho_positions_combinations,
                &pep_ion_spectra_titles,
                int_count
            );

            /*if pep_seq == "VGSLTPPSSPK" {
                println!("matrix_as_text_opt_222={:?}",matrix_as_text_opt);
            }*/

            for matrix_as_text in matrix_as_text_opt {
                let mut outfile = File::create(format!("{}/{}.tsv",output_path.to_string_lossy(),matrix_name_ref))?; //.expect("unable to create file");

                writeln!(&mut outfile, "{}", matrix_as_text.header)?; //.expect("unable to write");
                for line in matrix_as_text.lines.iter() {
                    writeln!(&mut outfile, "{}", line)?; //.expect("unable to write");

                }
            }
        }

    }

    Ok(())
}

fn _compute_phospho_evidence_matrix(
    phospho_matched_peaks: &Vec<PhosphoMatchedPeak>,
    phospho_positions: &Vec<usize>,
    n_spectra: usize,
    specific_peaks_only: bool
) -> Vec<SpectrumPhosphoAnalysis> {

    // Create Vector of Spectrum Phospho Analysis struct to collect all results
    let mut phospho_evidence_matrix = Vec::with_capacity(phospho_positions.len() * n_spectra);

    let filtered_phospho_matched_peaks = if !specific_peaks_only {
        phospho_matched_peaks.to_vec()
    } else {
        let pmps_by_ion_pos = phospho_matched_peaks.iter().into_group_map_by(|pmp| {
            vec!(pmp.matched_peak.ion_type.to_string(),pmp.matched_peak.aa_index.to_string()).join("|")
        });

        pmps_by_ion_pos.iter()
            // keep PhosphoMatchedPeak only if they are unique in the initial vector (phospho_matched_peaks)
            .filter( |(_key,pmps)| pmps.len() == 1)
            .map(|(_key,pmps)| (*pmps.first().unwrap()).to_owned()).collect_vec()
    };

    // First sort by ascending order and group all spectrum indices
    let spectra_pmps_iter = filtered_phospho_matched_peaks.into_iter()
        .sorted_by_key(|pmp| pmp.spectrum_index)
        .group_by(|pmp| pmp.spectrum_index);

    for (si, spectrum_pmps_groups) in &spectra_pmps_iter {
        // Iterate both spectrum and its index and change values to dereference each spectrum
        let spectrum_pmps: Vec<PhosphoMatchedPeak> = spectrum_pmps_groups.into_iter().map(|pmp_ref| (pmp_ref.clone())).collect();
        //let spectrum_pmps: Vec<PhosphoMatchedPeak> = spectrum_pmps_groups.into_iter().map(|pmp_ref| *pmp_ref.clone()).collect();

        // First sort by ascending order and group all spectrum_pmp groups depending on phospho_position
        /*let spectrum_pmps_groups = spectrum_pmps.into_iter()
            .sorted_by_key(|pmp| pmp.phospho_positions)
            .group_by(|pmp| pmp.phospho_positions);*/
        let spectrum_pmps_groups = spectrum_pmps.into_iter()
            .sorted_by_key(|pmp| pmp.phospho_positions.to_owned())
            .group_by(|pmp| pmp.phospho_positions.to_owned());


        for (phospho_positions, spectrum_pmp_group) in &spectrum_pmps_groups {

            let phospho_matched_peaks: Vec<PhosphoMatchedPeak> = spectrum_pmp_group.into_iter().map(|pmp| pmp).collect();
            let phospho_evidence_count = phospho_matched_peaks.iter().filter(|pmp| pmp.is_phospho_evidence).count() as u32;
            let non_phospho_evidence_count = phospho_matched_peaks.iter().filter(|pmp| !pmp.is_phospho_evidence).count() as u32;

            let phospho_evidence_intensity = phospho_matched_peaks.iter()
                .filter(|pmp| pmp.is_phospho_evidence)
                .fold(0.0, |sum, pmp| sum + pmp.matched_peak.peak_intensity);

            let non_phospho_evidence_intensity = phospho_matched_peaks.iter()
                .filter(|pmp| !pmp.is_phospho_evidence)
                .fold(0.0, |sum, pmp| sum + pmp.matched_peak.peak_intensity);

            let spa = SpectrumPhosphoAnalysis {
                phospho_matched_peaks: phospho_matched_peaks,
                spectrum_index: si,
                phospho_positions: phospho_positions,
                phospho_evidence_count: phospho_evidence_count,
                non_phospho_evidence_count: non_phospho_evidence_count,
                phospho_evidence_intensity: phospho_evidence_intensity,
                non_phospho_evidence_intensity: non_phospho_evidence_intensity,
            };

            phospho_evidence_matrix.push(spa);
        }
    }

    phospho_evidence_matrix
}
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct PhosphoEvidenceMatrixAsText {
    pub name: String,
    pub header: String,
    pub lines: Vec<String>
}
fn _print_phospho_evidence_matrix(
    name: &String,
    phospho_evidence_matrix: &Vec<SpectrumPhosphoAnalysis>,
    phospho_positions_combinations: &Vec<Vec<usize>>,
    spectra_titles: &Vec<String>,
    export_counts: bool
) -> Option<PhosphoEvidenceMatrixAsText> {
    let n_spectra = spectra_titles.len();
    if n_spectra == 0 { return None; }

    println!("\n{}",name);

    let mut header = Vec::with_capacity(2 * phospho_positions_combinations.len());
    header.push("SPECTRUM".to_string());
    for phospho_positions_combination in phospho_positions_combinations {
        //header.push(format!("P{}", phospho_positions_combination.join("&")));
        //header.push(format!("NP{}", phospho_positions_combination.join("&")));
        header.push(format!("P{:?}", phospho_positions_combination));
        header.push(format!("NP{:?}", phospho_positions_combination));
    }

    println!("{}",header.join("\t"));

    let mut lines = Vec::with_capacity(n_spectra);
    for spectrum_idx in 0..n_spectra {
        let spectrum_title = spectra_titles[spectrum_idx].clone();
        let spectrum_analyses: Vec<&SpectrumPhosphoAnalysis> = phospho_evidence_matrix.into_iter().filter(|spa| spa.spectrum_index == spectrum_idx).collect();
        let spectrum_analyses_ref = &spectrum_analyses;

        let mut values_as_strings: Vec<String> = Vec::new();
        values_as_strings.push(spectrum_title);
        //for phospho_positions_combination in &phospho_positions_combinations {
        for phospho_positions_combination in phospho_positions_combinations {
            let spectrum_analyses_for_pos: Vec<&&SpectrumPhosphoAnalysis> = spectrum_analyses_ref
                .into_iter().filter(|spa| spa.phospho_positions == phospho_positions_combination.clone()).collect();
            //.into_iter().filter(|spa| spa.phospho_positions == phospho_positions_combination).collect();
            if spectrum_analyses_for_pos.len() == 0 {
                values_as_strings.push("0".to_string());
                values_as_strings.push("0".to_string());
            } else {
                for spa in spectrum_analyses_for_pos {
                    if export_counts {
                        values_as_strings.push(spa.phospho_evidence_count.to_string());
                        values_as_strings.push(spa.non_phospho_evidence_count.to_string());
                    } else {
                        values_as_strings.push(spa.phospho_evidence_intensity.to_string());
                        values_as_strings.push(spa.non_phospho_evidence_intensity.to_string());
                    }

                }
            }
        }

        let line = values_as_strings.join("\t");
        println!("{}",line);

        lines.push(line);
    }

    let phospho_evidence_table = PhosphoEvidenceMatrixAsText {
        name: name.to_owned(),
        header: header.join("\t"),
        lines: lines
    };

    Some(phospho_evidence_table)
}

#[derive(Clone, PartialEq, Debug)]
pub struct SpectrumPhosphoAnalysis {
    pub phospho_matched_peaks: Vec<PhosphoMatchedPeak>,
    pub spectrum_index: usize,
    pub phospho_positions: Vec<usize>,
    pub phospho_evidence_count: u32,
    pub non_phospho_evidence_count: u32,
    pub phospho_evidence_intensity: f64,
    pub non_phospho_evidence_intensity: f64,
}

pub fn test_phospho_analysis_on_example_data() -> Result<()> {

    let aa_table = &*PROTEINOGENIC_AMINO_ACID_TABLE;

    let mz_error_tol = 0.02;
    let mut spectra = Vec::with_capacity(10);

    //type MatchedPeaksEachSpec = Vec<MatchedPeaks>
    //pub fn phospho_analysis(spectrum_peaks_array(list of all peaks),pep_seq,ion_type,charge,phospho_aa_pos) -> Result<Vec<MatchedPeaksEachSpec> {
/*    let spec_60736_peaks = [[110.071,68300.0],[112.087,8896.0],[115.087,9141.0],[120.081,41070.0],[121.028,6148.0],[126.055,57850.0],[127.039,9819.0],[127.086,7703.0],[129.066,11090.0],[129.102,75700.0],[130.086,18350.0],[130.107,5711.0],[133.043,10340.0],[133.061,7333.0],[136.063,9977.0],[136.076,63320.0],[138.055,108600.0],[139.059,7977.0],[141.066,7054.0],[141.103,11350.0],[143.118,12990.0],[144.065,5889.0],[145.098,6876.0],[147.113,34660.0],[155.081,8579.0],[157.061,7337.0],[157.097,5970.0],[157.134,6687.0],[159.077,11100.0],[159.092,10760.0],[159.113,9360.0],[166.061,6781.0],[167.081,7572.0],[168.066,46200.0],[169.098,14830.0],[169.134,8553.0],[171.077,9867.0],[171.114,12880.0],[173.092,12690.0],[173.13,7572.0],[175.119,60720.0],[183.114,7940.0],[183.15,19520.0],[185.091,15270.0],[185.165,11070.0],[186.076,19740.0],[186.235,6828.0],[187.073,13310.0],[187.108,11140.0],[187.144,10060.0],[195.077,7458.0],[197.127,7977.0],[198.089,8152.0],[199.072,11880.0],[199.108,12790.0],[200.14,6279.0],[201.087,12240.0],[201.123,12930.0],[203.066,9210.0],[204.087,72600.0],[209.09,6245.0],[211.145,15560.0],[212.103,9992.0],[215.139,32340.0],[216.097,7229.0],[217.082,13360.0],[218.15,8363.0],[226.119,13410.0],[226.155,7714.0],[227.067,10860.0],[227.103,15860.0],[228.135,6354.0],[229.121,8365.0],[230.079,5945.0],[230.114,9443.0],[231.061,7593.0],[234.144,6460.0],[235.119,40070.0],[240.103,7877.0],[242.114,7832.0],[243.135,21060.0],[244.094,9758.0],[244.136,8863.0],[244.167,16920.0],[245.077,10690.0],[251.154,7142.0],[256.128,7231.0],[258.146,7917.0],[259.103,5851.0],[272.173,28110.0],[274.093,21310.0],[286.677,6928.0],[296.123,7352.0],[314.099,7766.0],[315.205,17370.0],[327.129,8089.0],[341.145,11090.0],[348.203,13100.0],[356.154,7340.0],[358.163,9481.0],[366.14,17010.0],[394.202,16190.0],[394.7,7004.0],[457.222,6114.0],[501.314,23470.0],[502.32,7618.0],[505.732,5939.0],[538.222,6111.0],[545.277,7428.0],[554.269,7088.0],[556.244,7178.0],[572.355,31070.0],[574.249,7062.0],[592.27,5933.0],[600.777,8400.0],[601.276,9777.0],[646.292,8228.0],[658.288,12380.0],[679.305,16430.0],[687.377,20550.0],[689.325,9293.0],[689.838,8913.0],[690.334,22800.0],[697.826,15540.0],[715.801,7989.0],[720.33,7237.0],[721.337,12210.0],[722.334,7840.0],[741.349,11950.0],[772.34,8945.0],[787.398,46970.0],[788.399,20280.0],[802.403,26770.0],[803.406,11180.0],[808.356,25180.0],[822.377,8881.0],[855.395,10530.0],[856.41,8522.0],[868.488,15270.0],[871.405,12450.0],[875.387,11800.0],[877.423,14630.0],[877.921,14410.0],[878.424,8789.0],[880.444,9155.0],[886.413,11140.0],[887.395,13980.0],[892.088,11540.0],[892.417,21370.0],[892.748,13540.0],[893.409,11280.0],[899.408,9633.0],[903.418,21490.0],[904.435,8495.0],[914.403,8796.0],[921.427,27550.0],[925.44,53430.0],[925.938,12820.0],[926.41,18510.0],[926.908,12150.0],[931.443,26460.0],[932.438,15210.0],[941.459,9946.0],[943.435,9075.0],[976.489,13470.0],[1002.48,23210.0],[1003.507,13990.0],[1014.427,9292.0],[1046.529,25290.0],[1051.491,18220.0],[1052.472,12960.0],[1099.49,12900.0],[1100.502,12490.0],[1107.509,8872.0],[1122.0,11800.0],[1131.553,10180.0],[1132.551,13420.0],[1159.527,13010.0],[1170.554,11340.0],[1200.539,19820.0],[1201.529,11040.0],[1216.57,14030.0],[1217.552,10130.0],[1231.558,10440.0],[1430.612,18020.0],[1754.815,17790.0]];
    let spec_60840_peaks  = [[110.071,98170.0],[110.079,9392.0],[111.075,6742.0],[115.087,13700.0],[116.071,11290.0],[117.997,6826.0],[120.081,52280.0],[121.028,12780.0],[126.055,69300.0],[127.039,14950.0],[129.066,16010.0],[129.102,96450.0],[130.067,8469.0],[130.086,27080.0],[133.043,16350.0],[133.061,9314.0],[136.075,54410.0],[138.055,181300.0],[139.059,10730.0],[140.371,6315.0],[141.103,12650.0],[143.081,6460.0],[143.118,16040.0],[144.065,17470.0],[145.049,11030.0],[145.062,8809.0],[146.094,6761.0],[147.077,7427.0],[147.113,50240.0],[155.082,12310.0],[157.097,8346.0],[157.133,11810.0],[158.093,7546.0],[159.077,16360.0],[159.092,9463.0],[166.061,9633.0],[167.034,15000.0],[167.082,12290.0],[168.065,73600.0],[169.098,10490.0],[169.134,10790.0],[171.113,21960.0],[173.092,10070.0],[173.129,19370.0],[175.119,99060.0],[181.096,7545.0],[182.503,6934.0],[183.15,17850.0],[185.056,8414.0],[185.092,17020.0],[186.077,30850.0],[187.072,8806.0],[187.109,10500.0],[197.045,8464.0],[197.127,12350.0],[199.071,26820.0],[199.108,14500.0],[200.139,12610.0],[201.087,10980.0],[201.124,18430.0],[204.086,113700.0],[204.109,10670.0],[205.091,7114.0],[207.125,8975.0],[211.144,32480.0],[212.104,10490.0],[212.146,6854.0],[213.086,12800.0],[215.139,100500.0],[216.097,9158.0],[216.142,11400.0],[217.083,15730.0],[218.151,7473.0],[223.157,12280.0],[226.119,20730.0],[226.155,10590.0],[227.065,8963.0],[227.104,14990.0],[228.102,7147.0],[228.133,9221.0],[229.121,11880.0],[230.117,11580.0],[234.147,8683.0],[235.119,95490.0],[236.122,9867.0],[238.117,6727.0],[239.114,8271.0],[240.099,10150.0],[241.083,8163.0],[242.113,8714.0],[243.134,72270.0],[244.093,11450.0],[244.137,9908.0],[244.166,19750.0],[245.077,12040.0],[251.156,13940.0],[256.086,7782.0],[256.127,7507.0],[258.107,11200.0],[258.146,7517.0],[261.158,8279.0],[270.11,11900.0],[272.17,11360.0],[274.092,47170.0],[275.094,8401.0],[286.105,8949.0],[286.68,13970.0],[288.204,6687.0],[291.196,6696.0],[292.103,16660.0],[303.216,9979.0],[314.098,33680.0],[315.204,21440.0],[327.128,7299.0],[332.109,11570.0],[341.147,17230.0],[344.194,12440.0],[348.202,33640.0],[358.162,27940.0],[366.138,33530.0],[394.197,16090.0],[394.701,7454.0],[401.707,16110.0],[404.26,15460.0],[423.208,7572.0],[427.185,9823.0],[429.124,12720.0],[439.259,7594.0],[445.191,13730.0],[451.192,7090.0],[466.226,14170.0],[467.23,7618.0],[501.316,81300.0],[501.749,22070.0],[502.243,11080.0],[502.324,14160.0],[546.254,9038.0],[554.281,17150.0],[556.226,12550.0],[559.256,15050.0],[564.282,12320.0],[566.218,18290.0],[567.234,8376.0],[572.353,98100.0],[573.354,31350.0],[574.241,12100.0],[592.262,10210.0],[600.783,27960.0],[601.28,18410.0],[629.236,9042.0],[649.303,7998.0],[651.31,15440.0],[658.298,25630.0],[658.797,10840.0],[659.3,9638.0],[662.295,9542.0],[679.303,30650.0],[680.313,10680.0],[681.028,7601.0],[686.299,10800.0],[687.38,57070.0],[688.376,10520.0],[695.264,9243.0],[701.678,8658.0],[715.808,20870.0],[716.298,12590.0],[752.346,12440.0],[772.353,10250.0],[780.349,21280.0],[781.359,13570.0],[785.384,15120.0],[787.4,34570.0],[788.393,27590.0],[790.335,21670.0],[802.405,77320.0],[803.409,46480.0],[804.408,13690.0],[808.353,64020.0],[809.354,25460.0],[819.383,8704.0],[820.395,9393.0],[821.341,11380.0],[827.461,12800.0],[837.381,11520.0],[848.372,8581.0],[850.346,8923.0],[854.399,8922.0],[868.411,8627.0],[872.396,8007.0],[877.42,46200.0],[877.919,49260.0],[878.407,22490.0],[886.416,20160.0],[887.401,12390.0],[892.083,14560.0],[892.417,27410.0],[892.755,27820.0],[893.419,38680.0],[894.439,13350.0],[901.447,11080.0],[903.413,23140.0],[904.425,9545.0],[913.436,22650.0],[914.437,11170.0],[921.43,78570.0],[922.436,36430.0],[923.429,13060.0],[925.436,101100.0],[926.412,36820.0],[926.909,26240.0],[931.449,73370.0],[932.444,32890.0],[933.465,10770.0],[934.455,18820.0],[973.41,10730.0],[983.435,12610.0],[1000.476,15650.0],[1002.486,66760.0],[1003.5,42460.0],[1004.478,14670.0],[1007.452,10930.0],[1018.472,16430.0],[1030.478,11390.0],[1046.517,14810.0],[1051.012,21080.0],[1051.5,33060.0],[1052.016,18780.0],[1083.492,13770.0],[1100.024,11650.0],[1100.505,32560.0],[1103.544,13190.0],[1113.539,12180.0],[1117.521,43130.0],[1118.528,20070.0],[1131.562,36530.0],[1132.551,25240.0],[1135.501,11350.0],[1157.542,11350.0],[1200.556,20770.0],[1201.543,22280.0],[1217.058,16510.0],[1217.553,23380.0],[1244.637,19200.0],[1245.664,13500.0],[1281.579,18720.0],[1315.584,25450.0],[1360.684,13290.0],[1414.596,14850.0],[1430.606,40280.0],[1431.601,22370.0],[1753.815,36160.0],[1754.821,36220.0],[1755.811,19430.0]];
    let spec_61751_peaks = [[110.071,145100.0],[111.075,7668.0],[112.087,15430.0],[113.072,12220.0],[115.087,24770.0],[116.035,9839.0],[116.071,23190.0],[120.081,88920.0],[126.055,58280.0],[127.039,14100.0],[127.088,13940.0],[129.066,22260.0],[129.102,173900.0],[130.051,11460.0],[130.086,29120.0],[133.043,14070.0],[133.061,11120.0],[136.076,68840.0],[138.055,146600.0],[138.068,10260.0],[141.102,11380.0],[143.118,24030.0],[144.066,20980.0],[145.049,14860.0],[147.113,56200.0],[155.082,23250.0],[157.096,14910.0],[157.134,9146.0],[158.093,19240.0],[159.077,40040.0],[159.092,15110.0],[159.113,12580.0],[163.059,12970.0],[166.062,10130.0],[167.035,9209.0],[167.082,11180.0],[168.066,48860.0],[169.097,16280.0],[169.133,15630.0],[171.076,26650.0],[171.113,15130.0],[173.092,32800.0],[173.128,14580.0],[175.119,231600.0],[176.122,27750.0],[181.062,12230.0],[181.098,11010.0],[183.113,11510.0],[183.148,94710.0],[183.164,8898.0],[184.154,10830.0],[185.056,16340.0],[185.092,18280.0],[185.165,19040.0],[186.075,20760.0],[187.071,29220.0],[187.145,13240.0],[189.087,24200.0],[197.128,16390.0],[199.071,36180.0],[199.108,15080.0],[199.18,11440.0],[200.138,12800.0],[201.087,29080.0],[201.124,24960.0],[203.067,13420.0],[204.086,92690.0],[207.124,29460.0],[208.128,9927.0],[211.144,130000.0],[212.104,47560.0],[212.146,32880.0],[213.087,19920.0],[214.157,9264.0],[215.139,218900.0],[216.099,20200.0],[216.142,52900.0],[217.081,45590.0],[218.086,9223.0],[218.151,12290.0],[223.156,21120.0],[226.082,7755.0],[226.119,25690.0],[226.155,14840.0],[227.066,26020.0],[227.101,12070.0],[229.12,11740.0],[230.116,12980.0],[230.15,20390.0],[231.061,44990.0],[232.067,11270.0],[235.119,191500.0],[236.121,52090.0],[237.126,13930.0],[239.113,8753.0],[240.101,15710.0],[240.135,15930.0],[242.115,9830.0],[242.15,11910.0],[243.134,178100.0],[244.093,13990.0],[244.136,51070.0],[244.166,32100.0],[245.076,10240.0],[251.158,34200.0],[251.663,20660.0],[258.145,10430.0],[272.173,12940.0],[274.094,22650.0],[274.188,39890.0],[281.126,14280.0],[286.103,17480.0],[286.145,7506.0],[286.189,16150.0],[286.679,41370.0],[287.181,37310.0],[287.683,14430.0],[288.12,15950.0],[295.14,24420.0],[296.088,19740.0],[303.213,21440.0],[312.157,16980.0],[313.152,30740.0],[314.098,65070.0],[315.102,18680.0],[316.113,20290.0],[324.228,35030.0],[325.186,12350.0],[325.23,11140.0],[327.204,25720.0],[328.207,10250.0],[332.11,30370.0],[333.112,10810.0],[341.146,16690.0],[342.09,12220.0],[344.146,22310.0],[344.194,30430.0],[344.695,11380.0],[348.202,60390.0],[349.206,24240.0],[358.162,34600.0],[359.163,22850.0],[366.136,21970.0],[387.177,9479.0],[398.24,11550.0],[401.706,28670.0],[402.209,22000.0],[403.142,10280.0],[404.263,13850.0],[405.265,9335.0],[408.225,12260.0],[425.133,11040.0],[425.204,25030.0],[426.236,20750.0],[427.183,30020.0],[428.184,18480.0],[429.125,21600.0],[431.14,17730.0],[433.181,17060.0],[439.256,30790.0],[440.262,21310.0],[442.203,8958.0],[443.139,15440.0],[447.137,23050.0],[448.139,14620.0],[448.215,9492.0],[457.225,22260.0],[457.727,10220.0],[458.217,26450.0],[458.722,12470.0],[466.227,12510.0],[466.731,25760.0],[467.232,20290.0],[467.735,10720.0],[476.213,20390.0],[478.246,14740.0],[479.266,10600.0],[484.298,14460.0],[484.69,8368.0],[485.286,11420.0],[493.243,13500.0],[500.167,11750.0],[500.744,10740.0],[501.24,8708.0],[501.316,203600.0],[501.402,17300.0],[501.749,62220.0],[502.243,54210.0],[502.325,114800.0],[502.75,38310.0],[503.245,19110.0],[503.325,22980.0],[518.174,12170.0],[518.264,21130.0],[520.212,16130.0],[528.24,18760.0],[538.218,11000.0],[543.292,17530.0],[544.282,9537.0],[545.178,28040.0],[546.273,10340.0],[547.271,9482.0],[548.211,35010.0],[549.199,17950.0],[549.29,28650.0],[550.228,15490.0],[550.763,8806.0],[554.281,35000.0],[555.291,14660.0],[556.226,36550.0],[557.222,21780.0],[557.788,15340.0],[558.79,11130.0],[559.257,18880.0],[559.76,26180.0],[560.249,32000.0],[561.311,43290.0],[571.29,17710.0],[572.354,279900.0],[573.355,192200.0],[574.362,76400.0],[575.269,18190.0],[589.297,36480.0],[593.279,17060.0],[601.784,13570.0],[606.267,15300.0],[606.779,18710.0],[607.271,19840.0],[609.782,12830.0],[610.283,43990.0],[610.787,36360.0],[611.289,22720.0],[614.318,9001.0],[614.818,12490.0],[627.352,11920.0],[628.355,10060.0],[633.296,27690.0],[634.301,27020.0],[653.292,10580.0],[656.32,10910.0],[658.314,14930.0],[658.794,12290.0],[659.254,10690.0],[661.292,77660.0],[662.297,56670.0],[663.301,31490.0],[664.303,15880.0],[667.303,21380.0],[667.805,21450.0],[668.306,29180.0],[668.8,13750.0],[677.272,13810.0],[685.35,12620.0],[687.379,142100.0],[688.383,109200.0],[689.387,54640.0],[690.374,13000.0],[692.33,13090.0],[695.29,9639.0],[707.336,15460.0],[707.681,13340.0],[708.019,19550.0],[712.345,47130.0],[713.347,19300.0],[714.333,12880.0],[716.311,24610.0],[716.812,13660.0],[721.327,12070.0],[724.806,14680.0],[725.319,32380.0],[725.818,29010.0],[726.319,22600.0],[728.354,10350.0],[729.349,26650.0],[729.856,12650.0],[730.667,12080.0],[737.366,60530.0],[737.873,26590.0],[738.368,14980.0],[739.814,11400.0],[745.35,12500.0],[752.346,12290.0],[753.357,15800.0],[754.348,13190.0],[758.366,16680.0],[762.347,27090.0],[763.348,29960.0],[764.345,16010.0],[765.877,14160.0],[766.376,16520.0],[770.355,38340.0],[771.363,42630.0],[772.353,40550.0],[773.341,32440.0],[774.352,18280.0],[776.386,14200.0],[779.351,13570.0],[780.369,15470.0],[781.363,19870.0],[781.858,24360.0],[782.356,46200.0],[782.857,20000.0],[783.357,13800.0],[784.372,12690.0],[785.383,22150.0],[786.387,12860.0],[787.422,17350.0],[790.338,86570.0],[791.338,98630.0],[792.338,59740.0],[793.354,19510.0],[794.361,13380.0],[802.405,171900.0],[803.409,194500.0],[804.41,111400.0],[805.408,42930.0],[806.385,30850.0],[807.376,12780.0],[808.387,12430.0],[810.366,17830.0],[812.05,23680.0],[812.374,48280.0],[814.404,179000.0],[814.908,106400.0],[817.367,12420.0],[820.374,19810.0],[823.383,13260.0],[826.376,14200.0],[828.392,20990.0],[829.379,16710.0],[830.36,16740.0],[831.39,13120.0],[837.892,17650.0],[838.394,27840.0],[838.894,25470.0],[839.403,13960.0],[841.401,15770.0],[843.4,19300.0],[849.393,62490.0],[850.403,28560.0],[851.384,19690.0],[854.395,13470.0],[855.405,38870.0],[858.408,21710.0],[860.371,12860.0],[861.324,21250.0],[864.384,11390.0],[869.419,27950.0],[870.409,18850.0],[871.415,26900.0],[872.432,25120.0],[873.42,26460.0],[875.432,62900.0],[876.429,76530.0],[877.433,56160.0],[878.415,61930.0],[878.916,28590.0],[879.417,24940.0],[881.4,14650.0],[885.415,28670.0],[886.417,53230.0],[886.917,109900.0],[887.095,25160.0],[887.422,213000.0],[887.923,165000.0],[888.447,65160.0],[889.335,43530.0],[890.324,25430.0],[892.757,15850.0],[893.093,205200.0],[893.425,135700.0],[893.755,41070.0],[894.097,25220.0],[894.435,19790.0],[903.421,96950.0],[904.421,129900.0],[905.425,95170.0],[906.439,63680.0],[906.945,35880.0],[913.439,35730.0],[914.438,32910.0],[915.441,24190.0],[916.439,19890.0],[926.426,88570.0],[927.434,21590.0],[931.45,150200.0],[932.452,192100.0],[933.45,119000.0],[934.463,55540.0],[935.445,16400.0],[938.437,13310.0],[940.43,16880.0],[940.959,33160.0],[941.46,37860.0],[943.46,25680.0],[943.969,25560.0],[944.469,35740.0],[948.467,21170.0],[955.424,15940.0],[966.457,18920.0],[967.448,16870.0],[971.424,21950.0],[972.436,15820.0],[973.418,37030.0],[974.422,38890.0],[975.414,16720.0],[982.45,16480.0],[984.467,26820.0],[985.489,36870.0],[986.486,28510.0],[1000.457,19790.0],[1001.419,48150.0],[1002.479,153500.0],[1003.487,193100.0],[1004.491,116400.0],[1005.495,54780.0],[1012.482,27300.0],[1015.455,14670.0],[1020.442,13610.0],[1035.49,33170.0],[1038.011,18520.0],[1046.508,61900.0],[1046.988,43560.0],[1047.511,19260.0],[1060.523,16640.0],[1061.024,40220.0],[1061.518,73970.0],[1062.023,38690.0],[1085.49,16990.0],[1086.505,22090.0],[1087.537,19820.0],[1092.526,29160.0],[1095.5,39720.0],[1096.002,38850.0],[1096.533,18240.0],[1099.505,32760.0],[1100.529,41810.0],[1101.52,37900.0],[1102.482,21520.0],[1108.523,15370.0],[1113.545,41610.0],[1114.557,79320.0],[1115.565,49650.0],[1116.536,39550.0],[1117.521,91020.0],[1118.52,155200.0],[1119.521,119000.0],[1120.519,54700.0],[1121.527,17180.0],[1132.526,19560.0],[1133.477,33120.0],[1134.501,17950.0],[1149.545,69850.0],[1150.545,36190.0],[1160.553,21990.0],[1184.527,21610.0],[1201.521,20310.0],[1202.503,21070.0],[1209.058,22760.0],[1211.541,25090.0],[1212.54,39480.0],[1213.542,46550.0],[1214.534,23660.0],[1215.564,27870.0],[1216.55,27930.0],[1217.565,34030.0],[1218.068,71680.0],[1218.566,71740.0],[1219.547,55130.0],[1220.556,61410.0],[1221.562,32390.0],[1227.645,31060.0],[1228.634,30210.0],[1230.601,25940.0],[1243.568,17780.0],[1273.58,40470.0],[1282.581,38680.0],[1283.084,19270.0],[1315.582,23780.0],[1316.584,24600.0],[1318.594,19910.0],[1319.62,31860.0],[1326.613,30870.0],[1333.605,20070.0],[1334.583,46950.0],[1335.61,41280.0],[1336.61,27540.0],[1343.653,23460.0],[1344.632,19140.0],[1431.616,32130.0],[1432.633,37830.0],[1433.626,29910.0],[1448.611,33390.0],[1449.623,79310.0],[1450.618,98660.0],[1451.618,65350.0],[1452.626,20480.0],[1458.699,20090.0],[1531.708,17370.0],[1561.713,24130.0],[1562.707,32110.0],[1563.701,27370.0],[1627.798,33490.0],[1628.813,19120.0],[1675.748,21420.0],[1755.836,24940.0],[1756.845,21840.0],[1772.837,64820.0],[1773.842,90590.0],[1774.828,89010.0],[1775.858,25660.0]];
    let spec_61782_peaks = [[110.071,109400.0],[115.086,30520.0],[116.071,33450.0],[120.081,44250.0],[126.055,47510.0],[129.066,22390.0],[129.102,148300.0],[130.086,38660.0],[136.076,68660.0],[138.055,107000.0],[138.067,15170.0],[141.102,19710.0],[144.064,22740.0],[147.113,32430.0],[155.082,28520.0],[159.076,97380.0],[160.225,16510.0],[166.061,33230.0],[168.065,29600.0],[171.076,24250.0],[171.111,20120.0],[173.092,56440.0],[175.119,565000.0],[183.149,240000.0],[184.152,19650.0],[185.056,32810.0],[187.071,75640.0],[189.087,51690.0],[197.129,41850.0],[199.071,107300.0],[199.108,26780.0],[199.179,25330.0],[201.087,70590.0],[201.124,28610.0],[203.066,24670.0],[204.087,56410.0],[207.124,53230.0],[211.143,331900.0],[212.149,33210.0],[215.139,785500.0],[216.142,82370.0],[217.081,157400.0],[223.154,46330.0],[226.12,26220.0],[227.065,30930.0],[227.104,34840.0],[229.118,57100.0],[231.06,131100.0],[235.118,580900.0],[236.122,79560.0],[243.133,649100.0],[244.137,80690.0],[245.076,42240.0],[251.16,160400.0],[251.659,31870.0],[274.092,22600.0],[286.104,62370.0],[286.184,21880.0],[286.679,233200.0],[287.181,51610.0],[288.119,45370.0],[296.089,43180.0],[298.141,20360.0],[302.097,22640.0],[303.214,59930.0],[312.154,28650.0],[314.097,269800.0],[315.101,33770.0],[316.114,72320.0],[320.211,22250.0],[324.226,89170.0],[327.161,24160.0],[327.202,96780.0],[332.108,162200.0],[333.112,31240.0],[342.091,43920.0],[344.144,64600.0],[344.193,148800.0],[344.694,24270.0],[346.15,22980.0],[348.202,194100.0],[349.206,40480.0],[358.16,207400.0],[359.163,24080.0],[360.105,38710.0],[401.705,202500.0],[402.209,32110.0],[403.15,32710.0],[404.261,68830.0],[413.128,32390.0],[425.131,74170.0],[427.182,186500.0],[428.189,19850.0],[429.125,61820.0],[430.127,25030.0],[430.664,41850.0],[431.143,60810.0],[433.184,64170.0],[439.253,141200.0],[440.257,23910.0],[441.202,27970.0],[443.14,66390.0],[444.66,38510.0],[445.189,18700.0],[447.135,119400.0],[448.14,28440.0],[451.188,21210.0],[457.222,103700.0],[457.724,57180.0],[458.224,20610.0],[459.237,26750.0],[466.227,151700.0],[466.731,56370.0],[477.24,24490.0],[483.736,32960.0],[484.287,55580.0],[492.744,23580.0],[500.739,66320.0],[501.219,64790.0],[501.313,953000.0],[501.745,407000.0],[502.236,108400.0],[502.324,171500.0],[502.75,49250.0],[518.172,64350.0],[518.275,22400.0],[520.211,50190.0],[528.232,62210.0],[532.192,28180.0],[538.214,57210.0],[542.204,26910.0],[543.288,36950.0],[546.168,29140.0],[546.265,69820.0],[548.205,166800.0],[549.214,27560.0],[549.73,67670.0],[550.25,59360.0],[550.771,27490.0],[552.345,20410.0],[554.277,148800.0],[555.288,33860.0],[556.222,241600.0],[557.28,79470.0],[557.785,55060.0],[558.338,43670.0],[559.257,205700.0],[559.761,95150.0],[560.228,67610.0],[572.35,1388000.0],[573.354,299800.0],[574.36,37010.0],[592.271,52750.0],[592.773,25250.0],[600.777,37930.0],[601.275,24980.0],[606.266,114700.0],[606.774,44870.0],[607.276,30810.0],[608.296,28030.0],[609.78,206200.0],[610.281,111000.0],[610.784,36440.0],[613.824,39760.0],[614.324,34310.0],[624.288,22060.0],[626.202,24640.0],[627.353,36410.0],[629.199,25930.0],[631.248,21800.0],[633.297,126600.0],[634.299,29150.0],[644.26,24640.0],[647.216,54120.0],[648.829,30140.0],[649.283,35240.0],[649.794,26560.0],[653.282,38120.0],[654.202,26410.0],[655.328,34720.0],[658.305,72340.0],[658.804,21010.0],[659.244,24550.0],[661.293,323800.0],[662.299,76440.0],[663.307,40010.0],[667.294,158700.0],[667.795,91840.0],[668.303,28630.0],[669.366,35030.0],[670.359,45790.0],[671.346,74270.0],[671.841,19900.0],[677.254,51510.0],[687.376,746600.0],[688.38,208300.0],[689.389,28440.0],[706.827,31570.0],[707.018,64450.0],[707.347,95130.0],[713.823,27240.0],[715.814,64940.0],[720.332,41540.0],[720.822,39920.0],[722.333,23140.0],[722.835,28040.0],[724.812,178700.0],[725.306,98030.0],[725.814,39490.0],[728.848,83960.0],[729.351,64420.0],[729.853,23050.0],[744.338,34730.0],[752.347,53970.0],[757.218,25800.0],[758.353,33610.0],[759.269,28770.0],[760.302,30840.0],[762.352,112600.0],[763.344,51290.0],[767.371,35820.0],[770.356,185200.0],[771.356,80220.0],[772.332,131000.0],[773.033,39640.0],[773.354,81370.0],[773.703,33150.0],[777.839,58900.0],[778.334,34500.0],[779.371,69270.0],[779.877,59920.0],[780.84,27070.0],[781.352,163700.0],[781.857,75240.0],[782.36,59740.0],[784.395,58200.0],[785.389,35650.0],[786.452,32320.0],[787.394,59130.0],[790.335,478000.0],[791.338,165900.0],[792.364,44550.0],[802.404,1054000.0],[803.409,373200.0],[804.417,88740.0],[805.383,43120.0],[805.71,57480.0],[807.358,22560.0],[811.383,133800.0],[811.704,144900.0],[812.034,108400.0],[812.387,47640.0],[815.861,29850.0],[818.332,28780.0],[819.373,28920.0],[827.889,41160.0],[828.373,35930.0],[828.873,44730.0],[836.888,116300.0],[837.391,110000.0],[837.891,74170.0],[838.388,56800.0],[839.396,26250.0],[841.394,36000.0],[842.386,24660.0],[848.398,85230.0],[848.719,39810.0],[849.059,65530.0],[849.903,29790.0],[854.397,112300.0],[854.724,64330.0],[855.062,54840.0],[855.382,42950.0],[857.391,27190.0],[858.4,35490.0],[860.325,40220.0],[863.395,39220.0],[864.402,32960.0],[868.398,36260.0],[870.301,75910.0],[871.437,99050.0],[872.425,51200.0],[872.907,29420.0],[873.396,39400.0],[875.427,358500.0],[876.433,129200.0],[876.874,62950.0],[877.42,205200.0],[877.909,106300.0],[878.412,57400.0],[880.083,41360.0],[880.407,67950.0],[883.446,27280.0],[885.412,70430.0],[885.875,79000.0],[886.076,77880.0],[886.419,1123000.0],[886.916,640500.0],[887.07,52690.0],[887.417,305600.0],[887.922,84950.0],[888.302,257100.0],[889.308,89290.0],[892.088,542300.0],[892.419,563000.0],[892.751,391000.0],[893.081,176700.0],[901.471,73950.0],[903.42,743800.0],[904.424,288700.0],[905.427,66640.0],[913.433,157700.0],[914.432,74680.0],[915.428,32930.0],[921.42,34840.0],[924.449,80230.0],[931.446,1009000.0],[932.45,369800.0],[933.459,108100.0],[933.956,40730.0],[934.44,34960.0],[934.953,29660.0],[942.966,126800.0],[943.461,145300.0],[943.967,78920.0],[944.481,33480.0],[966.472,36150.0],[969.415,47440.0],[970.424,49650.0],[973.406,189600.0],[974.409,74980.0],[982.465,55560.0],[983.413,33040.0],[984.493,124400.0],[985.498,85520.0],[1000.473,77110.0],[1001.393,384400.0],[1002.485,961000.0],[1003.487,353300.0],[1004.499,85350.0],[1011.482,39720.0],[1011.987,41090.0],[1012.493,51840.0],[1043.445,29570.0],[1051.009,33110.0],[1051.511,44130.0],[1051.993,36290.0],[1052.457,33500.0],[1060.024,319200.0],[1060.515,297400.0],[1061.012,100400.0],[1061.502,49430.0],[1080.439,49650.0],[1081.48,38790.0],[1082.482,52930.0],[1083.004,28340.0],[1083.518,41340.0],[1085.518,65980.0],[1095.539,38630.0],[1096.544,36890.0],[1098.476,35250.0],[1099.529,181600.0],[1100.499,111500.0],[1113.554,362800.0],[1114.55,166700.0],[1115.544,92690.0],[1116.541,32010.0],[1117.511,833200.0],[1118.514,369600.0],[1119.524,95270.0],[1150.536,45730.0],[1183.531,36010.0],[1193.521,33080.0],[1197.535,42690.0],[1200.508,103400.0],[1201.513,39080.0],[1207.562,44870.0],[1208.052,69150.0],[1208.572,43690.0],[1210.573,31040.0],[1211.53,236400.0],[1212.542,122700.0],[1213.543,34180.0],[1214.557,99060.0],[1215.571,70620.0],[1216.562,217000.0],[1217.067,213800.0],[1217.562,104800.0],[1218.054,44690.0],[1218.558,296500.0],[1219.555,156600.0],[1220.559,50040.0],[1226.648,187600.0],[1227.629,105300.0],[1228.59,63090.0],[1229.583,35030.0],[1272.081,74690.0],[1272.558,106000.0],[1273.079,43860.0],[1273.583,30290.0],[1281.104,92660.0],[1281.582,136300.0],[1282.07,46900.0],[1282.596,36130.0],[1297.603,39520.0],[1312.547,35050.0],[1315.593,85550.0],[1316.594,57410.0],[1317.613,45480.0],[1324.615,127000.0],[1325.606,62490.0],[1326.615,45070.0],[1333.586,262700.0],[1334.584,122300.0],[1335.589,49890.0],[1341.661,96710.0],[1342.672,75120.0],[1343.645,43800.0],[1412.622,32910.0],[1413.609,38790.0],[1430.624,167700.0],[1431.624,126900.0],[1432.633,61490.0],[1439.633,44870.0],[1448.611,515900.0],[1449.613,274600.0],[1450.609,99130.0],[1456.697,105300.0],[1457.695,77320.0],[1559.688,34280.0],[1561.705,143200.0],[1562.689,93790.0],[1672.769,79910.0],[1673.757,66500.0],[1674.78,64740.0],[1753.813,78570.0],[1754.829,61680.0],[1771.83,448500.0],[1772.811,288400.0],[1773.779,119400.0],[1774.824,41300.0]];

    let pep_seq = "IEDSEPHIPLIDDTDAEDDAPTKR";
    spectra.push(spec_60736_peaks.to_vec());
    spectra.push(spec_60840_peaks.to_vec());
    spectra.push(spec_61751_peaks.to_vec());
    spectra.push(spec_61782_peaks.to_vec());

    let phospho_positions = [4, 14];
*/
//[[110.68624,757.5],[123.33151,900.7999877929688]
    let spec_peaks = [[110.68624,757.5],[123.33151,900.7999877929688],[129.10222,14770.0],[130.08601,1256.0],[143.1178,12020.0],[145.06033,3333.0],[147.11276,9363.0],[154.02582,3491.0],[155.1174,3975.0],[157.09706,33440.0],[171.1127,8172.0],[173.12842,15900.0],[173.6218,1622.0],[181.09642,3689.0],[183.11235,3803.0],[184.13239,1023.0],[185.09134,1368.0],[195.11304,1060.0],[197.12776,4835.0],[199.10709,3850.0],[201.12326,7243.0],[209.12831,1562.0],[212.1389,2029.0],[216.0416,3459.0],[216.13362,1728.0],[226.11836,16260.0],[226.15466,6481.0],[227.1385,1292.0],[228.1338,17490.0],[230.14963,1526.0],[240.13417,10150.0],[244.12907,30070.0],[244.16539,36530.0],[245.13211,3327.0],[245.16783,4331.0],[251.17473,2099.0],[258.14468,1412.0],[262.86487,1142.0],[266.14917,5422.0],[278.15009,3175.0],[282.1442,3805.0],[283.17599,7588.0],[284.16055,9656.0],[285.16351,1126.0],[292.12878,1805.0],[295.10553,1989.0],[306.67166,8577.0],[307.17169,1684.0],[308.19794,1116.0],[309.15527,5814.0],[311.20706,13440.0],[313.18625,1630.0],[323.17078,6798.0],[327.16464,7549.0],[329.21771,13180.0],[330.22125,1738.0],[331.19754,8413.0],[339.20291,6889.0],[340.22382,1203.0],[341.18149,30250.0],[342.18445,4789.0],[348.18933,32380.0],[348.69125,12530.0],[351.16537,3641.0],[357.21249,4380.0],[369.17609,7309.0],[377.21811,1238.0],[382.13617,7802.0],[394.24442,4755.0],[395.22791,4274.0],[396.18753,3746.0],[404.73123,20910.0],[405.23489,11090.0],[412.25415,4124.0],[418.22879,3563.0],[421.14661,1343.0],[422.23938,23190.0],[423.24277,7806.0],[439.24173,19080.0],[439.74179,10620.0],[440.24991,31210.0],[441.25323,9410.0],[448.24725,23500.0],[448.74832,13550.0],[453.71954,6797.0],[454.21924,1968.0],[458.74664,1302.0],[467.7532,8817.0],[468.25333,9659.0],[476.75775,44330.0],[477.25916,24920.0],[484.23868,4348.0],[488.22974,4590.0],[488.73282,4021.0],[496.67752,3846.0],[497.23672,9686.0],[497.27289,1563.0],[497.73529,5873.0],[508.28113,4846.0],[508.7822,4003.0],[509.23053,1534.0],[509.3078,1315.0],[515.28137,25320.0],[516.28375,7010.0],[516.73792,7028.0],[517.24188,3877.0],[517.28625,20300.0],[517.32996,1415.0],[517.78906,14030.0],[518.23578,1269.0],[520.21527,1475.0],[525.74683,5334.0],[526.24902,4902.0],[526.29144,59880.0],[526.79358,45800.0],[537.30182,4571.0],[565.77155,1795.0],[566.27441,10330.0],[566.77826,8106.0],[574.77765,12120.0],[576.31372,1838.0],[594.32416,11380.0],[595.32715,4880.0],[607.28436,3683.0],[610.31854,12100.0],[611.32312,15130.0],[612.33411,240000.0],[613.33728,77340.0],[634.31738,1273.0],[652.33044,1920.0],[677.36182,4607.0],[678.35889,1728.0],[679.34125,1268.0],[693.39246,13600.0],[694.39545,5577.0],[695.37122,65010.0],[696.37415,23900.0],[709.34973,4462.0],[710.35284,1398.0],[723.36993,6479.0],[724.36658,1330.0],[773.36328,1574.0],[775.33899,7572.0],[776.34351,1851.0],[790.40961,5293.0],[791.36584,20230.0],[792.37268,6816.0],[793.34729,28330.0],[794.35065,10310.0],[808.32831,4217.0],[808.46252,11500.0],[809.32837,6015.0],[809.42761,4686.0],[810.33368,1778.0],[881.45178,4629.0],[882.45996,4028.0],[895.48865,8745.0],[896.49042,4451.0],[906.43256,9387.0],[907.43427,3503.0],[934.5,1757.0],[935.49976,1806.0],[952.5097,23710.0],[953.50861,11960.0],[979.42664,9037.0],[980.43048,3411.0],[993.46515,12680.0],[994.46503,7237.0],[1033.47937,2059.0],[1050.48474,30580.0],[1051.48669,21580.0]];
    spectra.push(spec_peaks.to_vec());
    let pep_seq = String::from("VGSLTPPSSPK");
    let pep_seq1 = "VGSLTPPSSPK";
    let mut phospho_positions:Vec<usize> = Vec::new();
    phospho_positions.push(3);
    phospho_positions.push(5);
    phospho_positions.push(8);
    phospho_positions.push(9);

    let n_phosphates = phospho_positions.len();

    let spect_title = String::from("SF_200217_pPeptideLibrary_pool1_HCDOT_rep2");

    let mut spectra_titles: Vec<String> = Vec::with_capacity(10);
    spectra_titles.push(spect_title);

    let n_spectra = spectra.len();
    //output of phospho_matched_peaks= [PhosphoMatchedPeak {matched_peak: MatchedPeak {peak_mz,peak_intensity,theo_mz,mz_error,ion_type,charge,aa_index},spectrum_index, phospho_position, is_phospho_evidence}

    let phospho_positions_combinations = phospho_positions.to_vec().into_iter().combinations(1).collect_vec();
    //check it
    let phospho_matched_peaks = _match_phospho_peaks(
        pep_seq1,
        &phospho_positions_combinations,
        n_phosphates,
        &spectra,
        mz_error_tol,
        aa_table
    )?;
    //println!("{:?}",phospho_matched_peaks);

    let phospho_evidence_matrix = _compute_phospho_evidence_matrix(&phospho_matched_peaks, &phospho_positions, n_spectra, false);
    if pep_seq == "VGSLTPPSSPK" {
        println!("phospho_matched_peaks222={:?}", phospho_evidence_matrix);
        let x = 1;
    }


    /*
    let phospho_matched_peaks_60840 = annotate_spectrum(spec_60840_peaks.to_vec(), phospho_frag_table.clone(), mz_error_tol, PeakSelectionStrategy::HIGHEST_PEAK);

    println!("phospho_matched_peaks_60736={:?}", phospho_matched_peaks_60736);
    println!("phospho_matched_peaks_60840={:?}", phospho_matched_peaks_60840);

    println!("non_phospho_frag_table={:?}", non_phospho_frag_table[0]);

    //let non_phospho_matched_peaks_60840 = annotate_spectrum(spec_60840_peaks.to_vec(), non_phospho_frag_table.clone(), mz_error_tol, PeakSelectionStrategy::HIGHEST_PEAK);

    println!("non_phospho_matched_peaks_60736={:?}", non_phospho_matched_peaks_60736);
    println!("non_phospho_matched_peaks_60840={:?}", non_phospho_matched_peaks_60840);

    let count_matched_peaks_60736 = phospho_matched_peaks_60736.len();
    let count_unmatched_peaks_60736 = non_phospho_matched_peaks_60736.len();
    let sum_all_peaks_60736 = count_matched_peaks_60736 + count_unmatched_peaks_60736;
    let percentage_phospho = 100.0 * count_matched_peaks_60736 as f64 / sum_all_peaks_60736 as f64;
    let percentage_non_phospho = 100.0 * count_unmatched_peaks_60736 as f64 / sum_all_peaks_60736 as f64;
    println!("sum_all_peaks_60736={}",sum_all_peaks_60736);
    println!("percentage_phospho={}",percentage_phospho);
    println!("percentage_non_phospho={}",percentage_non_phospho);*/

    //pmp=phospho_matched_peak
    /*let mut final_result = Vec::with_capacity(phospho_positions.len() * n_spectra);
    let spectra_pmps_iter = phospho_matched_peaks.iter().into_group_map_by(|pmp| pmp.spectrum_index);
    for (si, spectrum_pmps) in spectra_pmps_iter {
        let spectrum_pmps_groups = spectrum_pmps.iter()
            .map(|pmp_ref| *pmp_ref).into_group_map_by(|pmp| pmp.phospho_positions.clone());
        //        let spectrum_pmps_groups = spectrum_pmps.iter()
        //             .map(|pmp_ref| *pmp_ref).into_group_map_by(|pmp| pmp.phospho_positions);
        //let spectrum_pmps_groups_vec: Vec<Vec<PhosphoMatchedPeak>> = spectrum_pmps_groups.into_iter().map(|(a,b)| b.collect()).collect();
        //println!("spectrum_pmps_groups_vec len={}",spectrum_pmps_groups_vec.len());

        for (phospho_positions, spectrum_pmp_group) in spectrum_pmps_groups {
            println!("phospho_position: {:?}", phospho_positions);
            //let phospho_matched_peaks2: Vec<PhosphoMatchedPeak> = spectrum_pmp_group.iter().map(|&pmp_ref| pmp_ref).collect();
            let phospho_matched_peaks2: Vec<PhosphoMatchedPeak> = spectrum_pmp_group.iter().map(|&pmp_ref| pmp_ref.clone()).collect();
            let phospho_evidence_count = phospho_matched_peaks2.iter().filter(|pmp| pmp.is_phospho_evidence).count() as u32;
            let non_phospho_evidence_count = phospho_matched_peaks2.iter().filter(|pmp| !pmp.is_phospho_evidence).count() as u32;

            let phospho_evidence_intensity = phospho_matched_peaks2.iter()
                .filter(|pmp| pmp.is_phospho_evidence)
                .fold(0.0, |sum, pmp| sum + pmp.matched_peak.peak_intensity);

            let non_phospho_evidence_intensity = phospho_matched_peaks2.iter()
                .filter(|pmp| !pmp.is_phospho_evidence)
                .fold(0.0, |sum, pmp| sum + pmp.matched_peak.peak_intensity);

            let spa = SpectrumPhosphoAnalysis {
                phospho_matched_peaks: phospho_matched_peaks2,
                spectrum_index: si,
                phospho_positions: phospho_positions,
                phospho_evidence_count: phospho_evidence_count,
                non_phospho_evidence_count: non_phospho_evidence_count,
                phospho_evidence_intensity: phospho_evidence_intensity,
                non_phospho_evidence_intensity: non_phospho_evidence_intensity,
            };

            final_result.push(spa);

        }

        println!("final_result_222={:?}",final_result);



    }*/

    //let header = ["P4","NP4","P14","NP14"];
    //let header = ["P3","NP3","P5","NP5","P8","NP8","P9","NP9"];
    //println!("{}",header.join("\t"));

/*
    let final_result_ref = &final_result;
    for spectrum_idx in (0..n_spectra) {
        let spectrum_analyses: Vec<&SpectrumPhosphoAnalysis> = final_result_ref.into_iter().filter(|spa| spa.spectrum_index == spectrum_idx).collect();
        let spectrum_analyses_ref = &spectrum_analyses;

        let mut values_as_strings: Vec<String> = Vec::new();
        for phospho_position in phospho_positions {
            let spectrum_analyses_for_pos: Vec<&&SpectrumPhosphoAnalysis> = spectrum_analyses_ref.into_iter().filter(|spa| spa.phospho_positions == phospho_positions).collect();
            if spectrum_analyses_for_pos.len() == 0 {
                values_as_strings.push("0".to_string());
                values_as_strings.push("0".to_string());
            } else {
                for spa in spectrum_analyses_for_pos {
                    values_as_strings.push(spa.phospho_evidence_intensity.to_string());
                    values_as_strings.push(spa.non_phospho_evidence_intensity.to_string());
                }
            }
        }
        println!("{}",values_as_strings.join("\t"));

    }
*/

    //println!("phospho count={}",phospho_matched_peaks.iter().filter(|p| p.is_phospho_evidence).count());
    //println!("non phospho count={}",phospho_matched_peaks.iter().filter(|p| p.is_phospho_evidence == false).count());
    let mut pem_in_text = _print_phospho_evidence_matrix(
        &pep_seq,
        &phospho_evidence_matrix,
        &phospho_positions_combinations,
        &spectra_titles,
        true
    );
    println!("PEM_in_text: {:?}",pem_in_text);

    Ok(())
}

#[derive(Clone, PartialEq, Debug)]
pub struct PhosphoMatchedPeak {
    matched_peak: MatchedPeak,
    spectrum_index: usize,
    phospho_positions: Vec<usize>,
    is_phospho_evidence: bool,
}

fn _match_phospho_peaks(
    pep_seq: &str,
    phospho_positions_combinations: &Vec<Vec<usize>>,
    n_phosphates: usize,
    spectra: &Vec<Vec<[f64;2]>>,
    mz_error_tol: f64,
    aa_table: &AminoAcidTable
) -> Result<Vec<PhosphoMatchedPeak>> {

    use FragmentIonSeries::*;

    // Put aa_pos as a variable in the compute_phospho_frag_table, and decide what else do you need as a variable to compute correctly your code
    let spectra_ref = &spectra;

    let expected_max_results_count = spectra.len() * 2 * phospho_positions_combinations.len() * pep_seq.len() * 2;
    let mut result = Vec::with_capacity(expected_max_results_count);

    // --- for each combination of phospho positions --- //
    for phospho_positions_combination in phospho_positions_combinations {
        println!("phospho_positions_combination={:?}",phospho_positions_combination);

        // --- 1. create phospho table --- //
        let phospho_frag_table = compute_phospho_frag_table(
            pep_seq,
            //&[a,a_H2O,a_NH3,b,b_H2O,b_NH3,y,y_H2O,y_NH3],
            &[b, y],
            &vec![1],
            &phospho_positions_combination,
            true,
            aa_table
        )?;

        // --- 2. create non-phospho table --- //
        let non_phospho_frag_table = compute_phospho_frag_table(
            pep_seq,
            //&[a,a_H2O,a_NH3,b,b_H2O,b_NH3,y,y_H2O,y_NH3],
            &[b, y],
            &vec![1],
            &phospho_positions_combination,
            false,
            aa_table
        )?;

        /*if pep_seq == "VGSLTPPSSPK" {
            println!("phospho_positions_combination={:?}",phospho_positions_combination);
            println!("phospho_frag_table={:?}",phospho_frag_table);
        }*/

        // --- 3. for loop for list of spectra -> annotate spectra  --- //
        for (spectrum_idx, spectrum_ref) in spectra_ref.into_iter().enumerate() {
            let phospho_matched_peaks = annotate_spectrum(spectrum_ref, &phospho_frag_table.clone(), mz_error_tol, PeakSelectionStrategy::HIGHEST_PEAK);
            let non_phospho_matched_peaks = annotate_spectrum(spectrum_ref, &non_phospho_frag_table.clone(), mz_error_tol, PeakSelectionStrategy::HIGHEST_PEAK);

            for matched_peak in phospho_matched_peaks {
                result.push(PhosphoMatchedPeak {
                    matched_peak: matched_peak,
                    spectrum_index: spectrum_idx,
                    phospho_positions: phospho_positions_combination.clone(),
                    is_phospho_evidence: true,
                });
            }

            for matched_peak in non_phospho_matched_peaks {
                result.push(PhosphoMatchedPeak {
                    matched_peak: matched_peak,
                    spectrum_index: spectrum_idx,
                    phospho_positions: phospho_positions_combination.clone(),
                    is_phospho_evidence: false,
                });
            }
        }
    }


   //println!("phospho_aa_frag_table={:?}", phospho_frag_table[0]);

    // TODO: check if we can replace multiple annotate_spectrum calls by a single one taking as input a Vector of many spectra

    Ok(result)
}

pub fn compute_phospho_frag_table(
    pep_seq: &str,
    ion_types: &[FragmentIonSeries],
    frag_ion_charges: &Vec<i8>,
    phospho_aa_positions: &Vec<usize>,
    apply_mass_inc: bool,
    aa_table: &AminoAcidTable
) -> Result<FragmentationTable> {

    let phospho_mass_inc = 79.966331; // H O(3) P => copy-pasted from https://www.unimod.org/modifications_view.php?editid1=21
    let pep_seq_len = pep_seq.len();
    let phospho_aa_forward_indices = phospho_aa_positions.iter().map(|&phospho_aa_pos| phospho_aa_pos - 1).collect_vec();
    let phospho_aa_reverse_indices = phospho_aa_positions.iter().map(|&phospho_aa_pos| pep_seq_len - phospho_aa_pos).collect_vec();
    let phospho_aa_forward_indices_ref = &phospho_aa_forward_indices;
    let phospho_aa_reverse_indices_ref = &phospho_aa_reverse_indices;

    use FragmentIonSeries::*;
    //the output type of compute_fragmentation_table is in Result<> because of that it did not use for determination of vector capacity
    let default_frag_table= compute_fragmentation_table(pep_seq, ion_types, frag_ion_charges, aa_table)?;
    //println!("frag_table_prec_3: {:?}",frag_table_prec_3)

    let mut updated_phospho_aa_frag_table = Vec::with_capacity(default_frag_table.len());
    // -- generation of multiple fragmentation table for each position's combination -- //
    // Create charge info by iterating fragment ions and create vectors to all frag_mz_values.
    for frag_series in default_frag_table.to_owned() {

        let frag_series_charge = frag_series.charge;
        let mut updated_frag_mz_values = Vec::with_capacity(frag_series.mz_values.len());
        let mut cloned_frag_series = frag_series.clone();

        // Iterate using both index and mz values of each fragment ions and determine phospho bound aminoacid position as index.
        // After determination of phospho bound aa indices, it searches and counts every frag_idx that has equal or lower than phospho index.
        // For more info look at your notebook.
        for (frag_idx, frag_mz) in frag_series.mz_values.into_iter().enumerate() {
            let phospho_aa_indices = if is_ion_forward(frag_series.ion_type).unwrap() { phospho_aa_forward_indices_ref } else { phospho_aa_reverse_indices_ref };
            let n_phosphates_at_cur_idx = phospho_aa_indices.iter()
                .filter(|&&phospho_aa_idx| phospho_aa_idx <= frag_idx).count();
                //.fold(0, |n_phosphates, &phospho_aa_idx| if phospho_aa_idx <= frag_idx { n_phosphates + 1 } else { n_phosphates });

            if n_phosphates_at_cur_idx == 0 {
                updated_frag_mz_values.push(0.0);
            } else {
                let updated_frag_mz: f64 = if apply_mass_inc {frag_mz + (n_phosphates_at_cur_idx as f64 * phospho_mass_inc) / frag_series_charge as f64} else {frag_mz};
                updated_frag_mz_values.push(updated_frag_mz);
            }
        }

        cloned_frag_series.mz_values = updated_frag_mz_values;
        updated_phospho_aa_frag_table.push(cloned_frag_series);
    }

    Ok(updated_phospho_aa_frag_table)
}

