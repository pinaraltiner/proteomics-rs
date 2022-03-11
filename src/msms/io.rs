//use mzdata::MGFReader;

pub fn read_mgf_file() {
    // See: https://github.com/mobiusklein/mzdata/blob/4050dd35d53ac461f5c71f3f7839f88ef40eb506/src/io/mgf.rs#L728
    /*use std::fs;
    use std::path;
    use core::result::Result::Ok;

    let path = path::Path::new("./data/OXRIA211021_03_CV45_mzcal.mzDB.mgf");
    let file = fs::File::open(path).expect("Test file doesn't exist");

    let reader = MGFReader::new(file); //MGFReaderType::<_, CentroidPeak, DeconvolutedPeak>::new(file);
    for scan in reader {
        let data = scan.arrays.unwrap();
        let mzs = data.mzs();
        let intensities = data.intensities();

        let peaks = scan.peaks.unwrap().peaks;/*.map(|peak| {
            peak.
        });*/

        peaks.map(|peak|)

        panic!("stop here");
    }*/
}