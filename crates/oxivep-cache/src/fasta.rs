use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};

/// Indexed FASTA reader for reference sequence access.
///
/// Supports .fai index for random access to specific regions.
pub struct FastaReader {
    sequences: HashMap<String, Vec<u8>>,
}

impl FastaReader {
    /// Load a FASTA file entirely into memory.
    /// For large genomes, use `from_indexed` instead.
    pub fn from_reader<R: Read>(reader: R) -> Result<Self> {
        let buf = BufReader::new(reader);
        let mut sequences = HashMap::new();
        let mut current_name: Option<String> = None;
        let mut current_seq = Vec::new();

        for line in buf.lines() {
            let line = line?;
            let line = line.trim();
            if line.starts_with('>') {
                if let Some(name) = current_name.take() {
                    sequences.insert(name, std::mem::take(&mut current_seq));
                }
                let name = line[1..].split_whitespace().next().unwrap_or("").to_string();
                current_name = Some(name);
                current_seq.clear();
            } else if !line.is_empty() {
                // Store uppercase at load time to avoid repeated to_ascii_uppercase() on fetch
                for &b in line.as_bytes() {
                    current_seq.push(b.to_ascii_uppercase());
                }
            }
        }
        if let Some(name) = current_name {
            sequences.insert(name, current_seq);
        }

        Ok(Self { sequences })
    }

    /// Get the list of sequence names.
    pub fn sequence_names(&self) -> Vec<&str> {
        self.sequences.keys().map(|s| s.as_str()).collect()
    }

    /// Fetch a region as a borrowed slice (1-based, inclusive coordinates).
    /// Zero-allocation since data is already stored uppercase in memory.
    pub fn fetch_slice(&self, chrom: &str, start: u64, end: u64) -> Result<&[u8]> {
        let seq = self
            .sequences
            .get(chrom)
            .with_context(|| format!("Chromosome '{}' not found in FASTA", chrom))?;

        let start_idx = (start.saturating_sub(1)) as usize;
        let end_idx = (end as usize).min(seq.len());

        if start_idx >= seq.len() {
            anyhow::bail!(
                "Start position {} exceeds sequence length {} for {}",
                start,
                seq.len(),
                chrom
            );
        }

        Ok(&seq[start_idx..end_idx])
    }

    /// Fetch a region of a sequence (1-based, inclusive coordinates).
    /// Returns uppercase nucleotides as a new Vec. Prefer `fetch_slice` to avoid allocation.
    pub fn fetch(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        self.fetch_slice(chrom, start, end).map(|s| s.to_vec())
    }

    /// Get the length of a chromosome.
    pub fn sequence_length(&self, chrom: &str) -> Option<u64> {
        self.sequences.get(chrom).map(|s| s.len() as u64)
    }
}

/// Memory-mapped FASTA reader using .fai index.
/// Avoids loading the entire FASTA into RAM by memory-mapping the file.
pub struct MmapFastaReader {
    mmap: memmap2::Mmap,
    index: Vec<FaiEntry>,
}

impl MmapFastaReader {
    /// Create a memory-mapped FASTA reader from a file with a .fai index.
    pub fn open(fasta_path: &std::path::Path) -> Result<Self> {
        let fai_path = format!("{}.fai", fasta_path.display());
        let fai_contents = std::fs::read_to_string(&fai_path)
            .with_context(|| format!("Reading FASTA index: {}", fai_path))?;
        let index = parse_fai(&fai_contents)?;

        let file = std::fs::File::open(fasta_path)
            .with_context(|| format!("Opening FASTA: {}", fasta_path.display()))?;
        let mmap = unsafe { memmap2::Mmap::map(&file) }
            .with_context(|| "Memory-mapping FASTA file")?;

        Ok(Self { mmap, index })
    }

    /// Fetch a region as a new Vec (1-based, inclusive coordinates).
    pub fn fetch(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let entry = self.index.iter().find(|e| e.name == chrom)
            .with_context(|| format!("Chromosome '{}' not found in FASTA index", chrom))?;

        let start_0 = start.saturating_sub(1);
        let end_0 = end.min(entry.length).saturating_sub(1);

        if start_0 >= entry.length {
            anyhow::bail!("Start {} exceeds length {} for {}", start, entry.length, chrom);
        }

        let bases_needed = (end_0 - start_0 + 1) as usize;
        let mut result = Vec::with_capacity(bases_needed);

        // Calculate byte offset for start position
        let start_line = start_0 / entry.line_bases;
        let start_col = start_0 % entry.line_bases;
        let mut byte_offset = (entry.offset + start_line * entry.line_bytes + start_col) as usize;

        while result.len() < bases_needed && byte_offset < self.mmap.len() {
            let b = self.mmap[byte_offset];
            if b == b'\n' || b == b'\r' {
                byte_offset += 1;
                continue;
            }
            result.push(b.to_ascii_uppercase());
            byte_offset += 1;
        }

        Ok(result)
    }

    /// Get the length of a chromosome.
    pub fn sequence_length(&self, chrom: &str) -> Option<u64> {
        self.index.iter().find(|e| e.name == chrom).map(|e| e.length)
    }

    /// Get the list of sequence names.
    pub fn sequence_names(&self) -> Vec<&str> {
        self.index.iter().map(|e| e.name.as_str()).collect()
    }
}

/// Parsed FASTA index (.fai) entry.
#[derive(Debug, Clone)]
pub struct FaiEntry {
    pub name: String,
    pub length: u64,
    pub offset: u64,
    pub line_bases: u64,
    pub line_bytes: u64,
}

/// Parse a .fai index file.
pub fn parse_fai(contents: &str) -> Result<Vec<FaiEntry>> {
    let mut entries = Vec::new();
    for line in contents.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        entries.push(FaiEntry {
            name: fields[0].to_string(),
            length: fields[1].parse()?,
            offset: fields[2].parse()?,
            line_bases: fields[3].parse()?,
            line_bytes: fields[4].parse()?,
        });
    }
    Ok(entries)
}

/// Read a region from a FASTA file using a .fai index without loading the whole file.
/// This is the preferred method for large genomes.
pub fn fetch_with_index<RS: Read + Seek>(
    reader: &mut RS,
    fai_entries: &[FaiEntry],
    chrom: &str,
    start: u64,
    end: u64,
) -> Result<Vec<u8>> {
    let entry = fai_entries
        .iter()
        .find(|e| e.name == chrom)
        .with_context(|| format!("Chromosome '{}' not found in FASTA index", chrom))?;

    let start_0 = start.saturating_sub(1);
    let end_0 = (end.min(entry.length)).saturating_sub(1);

    if start_0 >= entry.length {
        anyhow::bail!(
            "Start position {} exceeds sequence length {} for {}",
            start,
            entry.length,
            chrom
        );
    }

    // Calculate byte offset for start position
    let start_line = start_0 / entry.line_bases;
    let start_col = start_0 % entry.line_bases;
    let start_byte = entry.offset + start_line * entry.line_bytes + start_col;

    let bases_needed = (end_0 - start_0 + 1) as usize;

    reader.seek(SeekFrom::Start(start_byte))?;
    let mut buf_reader = BufReader::new(reader);
    let mut result = Vec::with_capacity(bases_needed);

    while result.len() < bases_needed {
        let mut line = String::new();
        let bytes = buf_reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        for &b in line.trim_end().as_bytes() {
            if result.len() >= bases_needed {
                break;
            }
            result.push(b.to_ascii_uppercase());
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_reader() {
        let fasta = ">chr1\nACGTACGT\nAAAACCCC\n>chr2\nTTTTGGGG\n";
        let reader = FastaReader::from_reader(fasta.as_bytes()).unwrap();
        assert_eq!(reader.fetch("chr1", 1, 4).unwrap(), b"ACGT");
        assert_eq!(reader.fetch("chr1", 5, 8).unwrap(), b"ACGT");
        assert_eq!(reader.fetch("chr1", 9, 16).unwrap(), b"AAAACCCC");
        assert_eq!(reader.fetch("chr2", 1, 4).unwrap(), b"TTTT");
        assert_eq!(reader.sequence_length("chr1"), Some(16));
    }

    #[test]
    fn test_fasta_reader_missing_chrom() {
        let fasta = ">chr1\nACGT\n";
        let reader = FastaReader::from_reader(fasta.as_bytes()).unwrap();
        assert!(reader.fetch("chr99", 1, 4).is_err());
    }

    #[test]
    fn test_parse_fai() {
        let fai = "chr1\t248956422\t6\t70\t71\nchr2\t242193529\t252513167\t70\t71\n";
        let entries = parse_fai(fai).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].name, "chr1");
        assert_eq!(entries[0].length, 248956422);
        assert_eq!(entries[1].name, "chr2");
    }

    #[test]
    fn test_fetch_with_index() {
        use std::io::Cursor;
        // Simulate a FASTA file: header + sequence with 10 bases per line
        let fasta = ">chr1\nACGTACGTAC\nGGGGAAAACC\nTTTT\n";
        let fai_entries = vec![FaiEntry {
            name: "chr1".into(),
            length: 24,
            offset: 6, // byte offset after ">chr1\n"
            line_bases: 10,
            line_bytes: 11, // 10 bases + newline
        }];

        let mut cursor = Cursor::new(fasta.as_bytes().to_vec());
        let seq = fetch_with_index(&mut cursor, &fai_entries, "chr1", 1, 4).unwrap();
        assert_eq!(seq, b"ACGT");
    }
}
