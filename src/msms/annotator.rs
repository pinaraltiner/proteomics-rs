
use anyhow::*;
use itertools::Itertools;
use rstar::RTree;
//let mut tree = RTree::new();
//tree.insert([0.0]);

use rstar::primitives::Rectangle;
use rstar::PointDistance;
use rstar::primitives::GeomWithData;

use crate::FragmentIonSeries;
use crate::msms::fragmentation::FragmentationTable;

type CustomExpPeak = GeomWithData<[f64; 2], f64>;
type CustomTheoIon = GeomWithData<Rectangle<[f64; 2]>, (usize,usize)>;
// --- Creation of "struct" that contains the required info about the exp. and theo. data for their --- //
// comparison and specifications of each candidates. --- //
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct MatchedPeak {
    peak_mz: f64,
    peak_intensity: f64,
    theo_mz: f64,
    mz_error: f32,
    ion_type: FragmentIonSeries,
    charge: i8,
    aa_index: i32,
}

pub enum PeakSelectionStrategy {
    NEAREST_PEAK,
    HIGHEST_PEAK,
}

pub fn annotate_spectrum(spectrum_peaks: Vec<[f64;2]>, frag_table: FragmentationTable, mz_error_tol: f64, peak_sel_strategy: PeakSelectionStrategy) -> Vec<MatchedPeak> {

    // --- Compute R*Tree indexing of experimental data --- //
    let all_custom_points: Vec<CustomExpPeak> = spectrum_peaks.iter().map(|peak| {
        CustomExpPeak::new([peak[0], 0.0], peak[1])
    }).collect();

    let exp_tree = RTree::bulk_load(all_custom_points);

    // --- Compute R*Tree indexing of theoretical data --- //
    let mut all_theo_rects = Vec::with_capacity(frag_table.len() * frag_table.first().unwrap().mz_values.len());

    let mut frag_table_col_idx = 0;
    for frag_series_fragments in frag_table.clone() {
        let current_series = frag_series_fragments.mz_values;

        let mut frag_table_row_idx = 0;
        for mz_value in current_series {
            let rect = Rectangle::from_corners([mz_value - mz_error_tol, 0.0], [mz_value + mz_error_tol, 0.0]);
            all_theo_rects.push(CustomTheoIon::new(rect, (frag_table_col_idx, frag_table_row_idx) ) );

            frag_table_row_idx += 1
        }

        frag_table_col_idx += 1
    }
    let theo_tree = RTree::bulk_load(all_theo_rects);

    // --- Compute matched peaks using the intersection between the two R*Trees (experimental and theoretical) --- //
    let intersection = exp_tree.intersection_candidates_with_other_tree(&theo_tree);
    /*intersection.for_each(|(exp_point, theo_line)| {
        let peak_mz = exp_point.geom()[0];
        let peak_intensity = exp_point.data;
        let theo_mz = theo_line.geom().upper()[0] + 0.02;
        let (ion_type, charge, aa_index) = theo_line.data;
        println!("{} {} {} {} {}", peak_mz, theo_mz, ion_type, charge, aa_index)
    });*/

    let matched_peaks: Vec<MatchedPeak> = intersection.map(|(exp_point, theo_line)| {
        let peak_mz = exp_point.geom()[0];
        let peak_intensity = exp_point.data;
        let (frag_table_col_idx, frag_table_row_idx) = theo_line.data;
        let frag_table_col = frag_table.get(frag_table_col_idx).unwrap();
        let theo_mz = frag_table_col.mz_values.get(frag_table_row_idx).unwrap();
        //println!("{} {} {} {} {}", peak_mz, theo_mz, ion_type, charge, aa_index)

        let mz_error = peak_mz - theo_mz;

        MatchedPeak {
            peak_mz: peak_mz,
            peak_intensity: peak_intensity,
            theo_mz: *theo_mz,
            mz_error: mz_error as f32,
            ion_type: frag_table_col.ion_type,
            charge: frag_table_col.charge,
            aa_index: (frag_table_row_idx + 1) as i32,
        }
    }).collect();

    // --- Remove matched peaks redundancy by selection the best candidate for a given theoretical m/z value --- //
    let mut non_redundant_matched_peaks = Vec::with_capacity(matched_peaks.len());

    let groups_iter = &matched_peaks.into_iter().group_by(|matched_peak| (matched_peak.theo_mz,matched_peak.ion_type));
    for (_, grouped_items_iter) in groups_iter {
        let grouped_matched_peaks: Vec<MatchedPeak> = grouped_items_iter.collect();
        let n_items = grouped_matched_peaks.len();
        //println!("{} {}", theo_ion_mz, grouped_matched_peaks.first().unwrap().peak_mz);

        // if we have more than one items in grouped_matched_peaks, it means that we have ambiguous peak matching for considered theoretical value.
        // so we need to choose one of the matched peaks according to the given selection strategy.
        if n_items > 1 {
            //println!("found one duplicate");
            let best_matched_peak = match peak_sel_strategy {
                PeakSelectionStrategy::NEAREST_PEAK => {
                    let nearest_peak_opt = grouped_matched_peaks.into_iter().min_by(|x1_peak, x2_peak| {
                        // Use partial_cmp because cmp cannot deal with floating point numbers (either f32 or f64)
                        x1_peak.mz_error.abs().partial_cmp(&x2_peak.mz_error.abs()).unwrap_or(std::cmp::Ordering::Equal)
                    });

                    nearest_peak_opt.unwrap() // safe because n_items > 1
                }
                PeakSelectionStrategy::HIGHEST_PEAK => {
                    let highest_peak_opt = grouped_matched_peaks.into_iter().max_by(|x1_peak, x2_peak| {
                        // Use partial_cmp because cmp cannot deal with floating point numbers (either f32 or f64)
                        x1_peak.peak_intensity.partial_cmp(&x2_peak.peak_intensity).unwrap_or(std::cmp::Ordering::Equal)
                    });

                    highest_peak_opt.unwrap() // safe because n_items > 1
                }
            };
            non_redundant_matched_peaks.push(best_matched_peak);
        } else if n_items == 1 {
            let single_matched_peak = grouped_matched_peaks.first().unwrap(); // safe because n_items == 1
            non_redundant_matched_peaks.push(*single_matched_peak);
        }
    }

    /*non_redundant_matched_peaks.iter().for_each(|matched_peak| {
        let peak_mz = matched_peak.peak_mz;
        let peak_intensity = matched_peak.peak_intensity;
        let theo_mz = matched_peak.theo_mz;
        let (ion_type, charge, aa_index) = (matched_peak.ion_type, matched_peak.charge, matched_peak.aa_index);
        //println!("{} {} {} {} {}", peak_mz, theo_mz, ion_type, charge, aa_index);
        println!("{}", matched_peak.mz_error);
    });*/

    non_redundant_matched_peaks
}

/*
struct ObservableFragIonTypesConfig {
    pub ion_type_a: bool,
    pub ion_type_b: bool,
    pub ion_type_c: bool,
    pub ion_type_x: bool,
    pub ion_type_y: bool,
    pub ion_type_z: bool,
    pub ion_a_charge: i8,
    pub ion_type_b: i8,
    pub ion_type_c: i8,
    pub ion_type_x: i8,
    pub ion_type_y: i8,
    pub ion_type_z: i8,
}

// the user can say for each type if it is present or not
var ionTypeA: Boolean = false,
var ionTypeB: Boolean = false,
var ionTypeC: Boolean = false,
var ionTypeX: Boolean = false,
var ionTypeY: Boolean = false,
var ionTypeZ: Boolean = false,

// the default charge for each ion type
var chargeForIonsA: Int = 1,
var chargeForIonsB: Int = 1,
var chargeForIonsC: Int = 1,
var chargeForIonsX: Int = 1,
var chargeForIonsY: Int = 1,
var chargeForIonsZ: Int = 1) {

def setIonTypeAndCharge(ion: String, charge: Int) {
ion.toLowerCase() match {
case "a" => {
ionTypeA = true
chargeForIonsA = charge
}
case "b" => {
ionTypeB = true
chargeForIonsB = charge
}
case "c" => {
ionTypeC = true
chargeForIonsC = charge
}
case "x" => {
ionTypeX = true
chargeForIonsX = charge
}
case "y" => {
ionTypeY = true
chargeForIonsY = charge
}
case "z" => {
ionTypeZ = true
chargeForIonsZ = charge
}
case _  => {

}
}
}

def contains(ionType: String): Boolean = {
val ion = ionType.toLowerCase()
if (ion == "a" && this.ionTypeA) true
else if (ion == "b" && this.ionTypeB) true
else if (ion == "c" && this.ionTypeC) true
else if (ion == "x" && this.ionTypeX) true
else if (ion == "y" && this.ionTypeY) true
else if (ion == "z" && this.ionTypeZ) true
else false
}

def getCharge(ionType: String): Int = {
val ion = ionType.toLowerCase()
if (ion == "a" && this.ionTypeA) this.chargeForIonsA
else if (ion == "b" && this.ionTypeB) this.chargeForIonsB
else if (ion == "c" && this.ionTypeC) this.chargeForIonsC
else if (ion == "x" && this.ionTypeX) this.chargeForIonsX
else if (ion == "y" && this.ionTypeY) this.chargeForIonsY
else if (ion == "z" && this.ionTypeZ) this.chargeForIonsZ
else 0
}

}*/

/*pub fn compute_frag_table() {

}*/