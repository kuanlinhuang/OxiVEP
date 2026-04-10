//! Writer for .osa2 format (ZIP-based chunked annotation files).
//!
//! Organizes annotations into ~1MB genomic chunks with parallel u32 value
//! arrays, sorted Var32 keys, and delta encoding for efficient compression.

use crate::chunk::delta_encode;
use crate::fields::{Field, FieldType};
use crate::kmer16::{self, LongVariant};
use crate::var32;
use anyhow::Result;
use std::io::{Seek, Write};

/// Metadata for the .osa2 file.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Osa2Metadata {
    pub format_version: u32,
    pub name: String,
    pub version: String,
    pub assembly: String,
    pub json_key: String,
    pub match_by_allele: bool,
    pub is_array: bool,
    pub is_positional: bool,
    pub chunk_bits: u32,
    pub description: String,
}

/// A record to write into the .osa2 file.
pub struct Osa2Record {
    pub chrom: String,
    pub position: u32,
    pub ref_allele: Vec<u8>,
    pub alt_allele: Vec<u8>,
    /// Parallel field values (same order as the Field config).
    pub values: Vec<u32>,
    /// Optional JSON blob for JsonBlob fields.
    pub json_blob: Option<String>,
}

/// Builds an .osa2 file from sorted records.
pub struct Osa2Writer {
    metadata: Osa2Metadata,
    fields: Vec<Field>,
    /// Categorical string tables: field_idx -> list of unique strings.
    string_tables: Vec<Vec<String>>,
}

impl Osa2Writer {
    pub fn new(metadata: Osa2Metadata, fields: Vec<Field>) -> Self {
        let string_tables = fields.iter().map(|_| Vec::new()).collect();
        Self { metadata, fields, string_tables }
    }

    /// Set the string table for a categorical field.
    pub fn set_string_table(&mut self, field_idx: usize, strings: Vec<String>) {
        if field_idx < self.string_tables.len() {
            self.string_tables[field_idx] = strings;
        }
    }

    /// Write all records to a .osa2 ZIP file.
    ///
    /// Records MUST be sorted by (chrom, position).
    pub fn write_all<W: Write + Seek>(
        &self,
        writer: W,
        records: &[Osa2Record],
    ) -> Result<()> {
        let mut zip = zip::ZipWriter::new(writer);
        let options = zip::write::SimpleFileOptions::default()
            .compression_method(zip::CompressionMethod::Deflated);

        // Write metadata
        zip.start_file("fastsa/metadata.json", options)?;
        serde_json::to_writer_pretty(&mut zip, &self.metadata)?;

        // Write field config
        zip.start_file("fastsa/config.json", options)?;
        serde_json::to_writer_pretty(&mut zip, &self.fields)?;

        // Write string tables
        for (i, field) in self.fields.iter().enumerate() {
            if field.ftype == FieldType::Categorical && !self.string_tables[i].is_empty() {
                let path = format!("fastsa/strings/{}.txt", field.alias);
                zip.start_file(path, options)?;
                for s in &self.string_tables[i] {
                    writeln!(zip, "{}", s)?;
                }
            }
        }

        // Group records by (chrom, chunk_id)
        let chunk_bits = self.metadata.chunk_bits;
        let mut chunks: Vec<(String, u32, Vec<usize>)> = Vec::new();
        let mut current_key: Option<(String, u32)> = None;

        for (ri, record) in records.iter().enumerate() {
            let cid = record.position >> chunk_bits;
            let key = (record.chrom.clone(), cid);

            if current_key.as_ref() != Some(&key) {
                chunks.push((record.chrom.clone(), cid, Vec::new()));
                current_key = Some(key);
            }
            chunks.last_mut().unwrap().2.push(ri);
        }

        // Write each chunk
        for (chrom, chunk_id, record_indices) in &chunks {
            self.write_chunk(&mut zip, options, records, chrom, *chunk_id, chunk_bits, record_indices)?;
        }

        zip.finish()?;
        Ok(())
    }

    fn write_chunk<W: Write + Seek>(
        &self,
        zip: &mut zip::ZipWriter<W>,
        options: zip::write::SimpleFileOptions,
        records: &[Osa2Record],
        chrom: &str,
        chunk_id: u32,
        chunk_bits: u32,
        indices: &[usize],
    ) -> Result<()> {
        let chunk_mask = (1u32 << chunk_bits) - 1;

        // Build sorted Var32 keys and parallel value arrays
        let mut short_entries: Vec<(u32, usize)> = Vec::new(); // (var32_key, original_idx)
        let mut long_entries: Vec<(LongVariant, usize)> = Vec::new();

        for (sort_idx, &ri) in indices.iter().enumerate() {
            let record = &records[ri];
            let within_chunk_pos = record.position & chunk_mask;

            if var32::is_long(record.ref_allele.len(), record.alt_allele.len()) {
                long_entries.push((
                    LongVariant {
                        position: record.position,
                        idx: sort_idx as u32,
                        sequence: kmer16::encode_var(&record.ref_allele, &record.alt_allele),
                    },
                    ri,
                ));
            } else if let Some(key) = var32::encode(within_chunk_pos, &record.ref_allele, &record.alt_allele) {
                short_entries.push((key, ri));
            }
        }

        // Sort by Var32 key
        short_entries.sort_by_key(|(key, _)| *key);
        long_entries.sort_by(|(a, _), (b, _)| a.cmp(b));

        // Build sorted arrays
        let var32s: Vec<u32> = short_entries.iter().map(|(k, _)| *k).collect();
        let sorted_record_order: Vec<usize> = short_entries.iter().map(|(_, ri)| *ri).collect();

        // Delta-encode var32 keys
        let delta_var32s = delta_encode(&var32s);

        // Write var32 keys
        let prefix = format!("fastsa/{}/{}/", chrom, chunk_id);
        zip.start_file(format!("{}var32.bin", prefix), options)?;
        write_u32_array(zip, &delta_var32s)?;

        // Write long variants
        if !long_entries.is_empty() {
            let longs: Vec<&LongVariant> = long_entries.iter().map(|(lv, _)| lv).collect();
            zip.start_file(format!("{}too-long.enc", prefix), options)?;
            let data = bincode::serialize(&longs)?;
            zip.write_all(&data)?;
        }

        // Write parallel value arrays
        for (fi, field) in self.fields.iter().enumerate() {
            if field.ftype == FieldType::JsonBlob {
                continue;
            }

            let values: Vec<u32> = sorted_record_order
                .iter()
                .map(|&ri| {
                    if fi < records[ri].values.len() {
                        records[ri].values[fi]
                    } else {
                        field.missing_value
                    }
                })
                .collect();

            zip.start_file(format!("{}{}.bin", prefix, field.alias), options)?;
            write_u32_array(zip, &values)?;
        }

        // Write JSON blobs if any
        let has_blobs = self.fields.iter().any(|f| f.ftype == FieldType::JsonBlob);
        if has_blobs {
            let blobs: Vec<&str> = sorted_record_order
                .iter()
                .map(|&ri| {
                    records[ri].json_blob.as_deref().unwrap_or("")
                })
                .collect();
            if blobs.iter().any(|b| !b.is_empty()) {
                zip.start_file(format!("{}json_blobs.zst", prefix), options)?;
                let joined = blobs.join("\n");
                let compressed = zstd::encode_all(joined.as_bytes(), 3)?;
                zip.write_all(&compressed)?;
            }
        }

        Ok(())
    }
}

/// Write a u32 array as [4B count][4B * count values].
fn write_u32_array<W: Write>(writer: &mut W, values: &[u32]) -> Result<()> {
    writer.write_all(&(values.len() as u32).to_le_bytes())?;
    for &v in values {
        writer.write_all(&v.to_le_bytes())?;
    }
    Ok(())
}

/// Read a u32 array from [4B count][4B * count values].
pub fn read_u32_array(data: &[u8]) -> Result<Vec<u32>> {
    if data.len() < 4 {
        anyhow::bail!("u32 array too short");
    }
    let count = u32::from_le_bytes(data[0..4].try_into()?) as usize;
    let mut values = Vec::with_capacity(count);
    let mut offset = 4;
    for _ in 0..count {
        if offset + 4 > data.len() {
            anyhow::bail!("u32 array truncated");
        }
        values.push(u32::from_le_bytes(data[offset..offset + 4].try_into()?));
        offset += 4;
    }
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_u32_array_round_trip() {
        let values = vec![1, 2, 3, 100, 200];
        let mut buf = Vec::new();
        write_u32_array(&mut buf, &values).unwrap();
        let decoded = read_u32_array(&buf).unwrap();
        assert_eq!(decoded, values);
    }
}
