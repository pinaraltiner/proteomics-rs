
// Source: https://stackoverflow.com/questions/45882329/read-large-files-line-by-line-in-rust
use std::{
    fs::File,
    io::{self, prelude::*},
    rc::Rc,
};

pub struct TextReader {
    pub reader: io::BufReader<File>,
    buf: Rc<String>,
}

fn new_buf(capacity: usize) -> Rc<String> {
    Rc::new(String::with_capacity(capacity)) // Tweakable capacity
}

impl TextReader {
    pub fn open(path: impl AsRef<std::path::Path>, capacity: usize) -> io::Result<Self> {
        let file = File::open(path)?;
        let buf = new_buf(capacity);
        let reader = io::BufReader::with_capacity(buf.capacity(), file);

        Ok(Self { reader, buf })
    }
}

impl Iterator for TextReader {
    type Item = io::Result<Rc<String>>;

    fn next(&mut self) -> Option<Self::Item> {
        let buf = match Rc::get_mut(&mut self.buf) {
            Some(buf) => {
                buf.clear();
                buf
            }
            // TODO: try to understand why this could happen
            None => {
                self.buf = new_buf(1024 * 1024);
                Rc::make_mut(&mut self.buf)
            }
        };

        self.reader
            .read_line(buf)
            .map(|u| if u == 0 { None } else { Some(Rc::clone(&self.buf)) })
            .transpose()
    }
}
