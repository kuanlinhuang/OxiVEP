//! Simple Bloom filter for fast negative lookups on genomic positions.
//!
//! Before loading and searching a chunk, check the Bloom filter to skip
//! decompression entirely when a position is definitely not present.
//! False positives are acceptable (they just cause an unnecessary chunk load);
//! false negatives never occur.

/// A compact Bloom filter backed by a bit array.
pub struct BloomFilter {
    bits: Vec<u64>,
    num_bits: usize,
    num_hashes: u32,
}

impl BloomFilter {
    /// Create a Bloom filter sized for `expected_items` with target false positive rate.
    pub fn new(expected_items: usize, fp_rate: f64) -> Self {
        let num_bits = optimal_num_bits(expected_items, fp_rate).max(64);
        let num_hashes = optimal_num_hashes(expected_items, num_bits);
        let num_words = (num_bits + 63) / 64;
        Self {
            bits: vec![0u64; num_words],
            num_bits,
            num_hashes,
        }
    }

    /// Insert a position into the filter.
    pub fn insert(&mut self, position: u32) {
        for i in 0..self.num_hashes {
            let bit = self.hash(position, i) % self.num_bits;
            self.bits[bit / 64] |= 1u64 << (bit % 64);
        }
    }

    /// Check if a position might be present (true = maybe, false = definitely not).
    #[inline]
    pub fn might_contain(&self, position: u32) -> bool {
        for i in 0..self.num_hashes {
            let bit = self.hash(position, i) % self.num_bits;
            if self.bits[bit / 64] & (1u64 << (bit % 64)) == 0 {
                return false;
            }
        }
        true
    }

    /// Simple double-hashing using mixing functions.
    #[inline]
    fn hash(&self, position: u32, i: u32) -> usize {
        let h1 = fmix32(position) as usize;
        let h2 = fmix32(position.wrapping_add(0x9e3779b9)) as usize;
        h1.wrapping_add(h2.wrapping_mul(i as usize))
    }
}

/// Murmur3 finalization mix for u32.
#[inline]
fn fmix32(mut h: u32) -> u32 {
    h ^= h >> 16;
    h = h.wrapping_mul(0x85ebca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2ae35);
    h ^= h >> 16;
    h
}

fn optimal_num_bits(n: usize, p: f64) -> usize {
    let ln2_sq = std::f64::consts::LN_2 * std::f64::consts::LN_2;
    (-(n as f64) * p.ln() / ln2_sq).ceil() as usize
}

fn optimal_num_hashes(n: usize, m: usize) -> u32 {
    let k = (m as f64 / n as f64) * std::f64::consts::LN_2;
    (k.ceil() as u32).max(1).min(16)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bloom_basic() {
        let mut bloom = BloomFilter::new(1000, 0.01);

        // Insert some positions
        for i in 0..100 {
            bloom.insert(i * 100);
        }

        // All inserted positions should be found
        for i in 0..100 {
            assert!(bloom.might_contain(i * 100), "Missing position {}", i * 100);
        }

        // Most non-inserted positions should NOT be found (allow ~1% FP)
        let mut false_positives = 0;
        for i in 0..1000 {
            let pos = i * 100 + 50; // Offsets that were never inserted
            if bloom.might_contain(pos) {
                false_positives += 1;
            }
        }
        // With 1% FP rate and 1000 queries, expect ~10 false positives
        assert!(false_positives < 50, "Too many false positives: {}", false_positives);
    }

    #[test]
    fn test_bloom_empty() {
        let bloom = BloomFilter::new(100, 0.01);
        assert!(!bloom.might_contain(12345));
    }
}
