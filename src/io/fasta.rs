
use anyhow::*;
use std::io;

use crate::io::reader::TextReader;

pub fn for_each_fasta_entry<F>(path: &str, mut cb: F) -> Result<()> where F: FnMut(&str,&str,u64,u32) -> () {

    //let file = File::open("./data/ups1.fasta")?;
    //let reader = BufReader::with_capacity(1024 * 1024, file);

    let mut header = String::with_capacity(512);
    let mut seq: String = String::with_capacity(100 * 1024);

    let mut begin_pos: u64 = 0;
    let mut end_pos: u64 = 0;

    let text_reader = TextReader::open(path, 10 * 1024 * 1024)?; // buffer capacity = 10MB

    for line in text_reader {
        let l = line?;
        //println!("Line len={}", l.len() );

        // When a new FASTA entry is met
        if l.as_str().starts_with(">") {
            if !seq.is_empty() {
                // Emit FASTA entry content as two strings (header and sequence data)
                cb(header.as_str(), seq.as_str(), begin_pos, (end_pos - begin_pos) as u32);

                begin_pos = end_pos;
            }

            let header_without_gt_char = l.char_indices()
                .nth(1)
                .map(|(i, _)| &l[i..])
                .unwrap();

            // See: https://github.com/hoodie/concatenation_benchmarks-rs
            header.clear();
            header.push_str(header_without_gt_char);

            seq.clear();
        } else {
            seq.push_str(l.as_str());
        }

        end_pos += l.len() as u64;
        //end_pos = reader.seek(SeekFrom::Current(0))?;
    }

    // Emit last entry
    if !seq.is_empty() {
        cb(header.as_str(),seq.as_str(), begin_pos, (end_pos - begin_pos) as u32 );
    }

    Ok(())
}