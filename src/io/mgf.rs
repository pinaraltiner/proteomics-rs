use anyhow::*;
use serde::{Serialize, Deserialize};
use serde::de::StdError;

use crate::io::reader::TextReader;

#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct MgfSpectrum {
    pub header: MgfSpectrumHeader,
    pub data: SpectrumData,
}

#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct MgfSpectrumHeader {
    pub title: String,
    pub precursor_mz: f64,
    pub precursor_charge: Option<i8>,
    pub retention_time: Option<f32>, // TODO: add this
    //peaks: Vec<[f64;2]>, // could aslo be Vec<(f64,f64)>
}

#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct SpectrumData {
    pub mz_list: Vec<f64>,
    pub intensity_list: Vec<f32>,
}

pub fn for_each_spectrum<F>(path: &str, mut cb: F) -> Result<()> where F: FnMut(MgfSpectrum) -> () {

    let text_reader = TextReader::open(path, 10 * 1024 * 1024)?; // buffer capacity = 10MB

    let mut mz_list_buffer: Vec<f64> = Vec::with_capacity(10000);
    let mut intensity_list_buffer: Vec<f32> = Vec::with_capacity(10000);

    let mut title: String = "".to_string();
    let mut pep_mass: f64 = 0.0;
    let mut charge: Option<i8> = None;
    let mut rt: Option<f32> = None;

    let mut is_inside_spectrum_block = false;
    for line_res in text_reader {

        let line_rc = line_res?;

        //let line = line_rc.trim_end();
        if line_rc.len() <= 2 {
            continue;
        }

        let line = line_rc.as_str();
        let first_char = line.chars().next().unwrap();

        if first_char == 'B' && line.starts_with("BEGIN IONS") {
            // reset some variables (note: we don't reset the title because it should always be there and parsed)
            pep_mass = 0.0;
            charge = None;
            rt = None;

            is_inside_spectrum_block = true;
        } else if first_char == 'E' && line.starts_with("END IONS") {

            cb(MgfSpectrum {
                header: MgfSpectrumHeader {
                    title: title.trim_end().to_owned(),
                    precursor_mz: pep_mass,
                    precursor_charge: charge,
                    retention_time: rt
                },
                data: SpectrumData {
                    mz_list: mz_list_buffer.to_vec(),
                    intensity_list: intensity_list_buffer.to_vec(),
                }
            });

            mz_list_buffer.clear();
            intensity_list_buffer.clear();

            is_inside_spectrum_block = false;
        } else if is_inside_spectrum_block {

            // if peak line
            if first_char.is_numeric() {
                //line.split(' ');
                let mut parts = line.split_ascii_whitespace();
                let mz_str_opt = parts.next();
                let intensity_str_opt = parts.next();

                //let nparts = parts.len();
                if mz_str_opt.is_none() || intensity_str_opt.is_none() {
                    bail!("invalid number of columns for a peak line");
                }

                let mz: f64 = fast_float::parse(mz_str_opt.unwrap())?;
                let intensity: f32 = fast_float::parse(intensity_str_opt.unwrap())?;
                /*let mz: f64 = mz_str_opt.unwrap().parse()?;
                let intensity: f32 = intensity_str_opt.unwrap().parse()?;*/

                mz_list_buffer.push(mz);
                intensity_list_buffer.push(intensity);
            } else { // header line assumed

                match line.split_once('=') {
                    Some(("TITLE", value)) => {
                        title = value.to_string();
                    }
                    Some(("PEPMASS", value)) => {
                        let pep_mass_str_opt = value.split_ascii_whitespace().next();
                        pep_mass = pep_mass_str_opt.unwrap_or_else(|| "0.0").parse()?;
                    }
                    Some(("CHARGE", value)) => {
                        let charge_sign_idx = value.chars().position(|c| !c.is_numeric()).unwrap_or(value.len());
                        let charge_str = &value[0..charge_sign_idx];
                        // FIXME: parse the sign
                        charge = charge_str.parse::<i8>().ok();
                    }
                    Some(("RTINSECONDS", value)) => {
                        rt = value.parse::<f32>().ok();
                    }
                    _ => {}
                }
            }
        }
    }

    Ok(())
}