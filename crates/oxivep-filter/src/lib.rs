//! Variant filtering for OxiVEP (filter_vep equivalent).
//!
//! Supports complex filter expressions on VEP annotation fields.
//!
//! # Syntax
//!
//! ```text
//! IMPACT is HIGH
//! Consequence is missense_variant
//! gnomAD_AF < 0.001
//! IMPACT is HIGH and gnomAD_AF < 0.001
//! IMPACT is HIGH or Consequence in missense_variant,stop_gained
//! not IMPACT is LOW
//! (IMPACT is HIGH) and (Consequence is missense_variant or Consequence is stop_gained)
//! ```
//!
//! # Operators
//! - `is` / `=` / `eq`: equality
//! - `!=` / `ne`: inequality
//! - `<`, `>`, `<=`, `>=`: numeric comparison
//! - `in`: value is one of a comma-separated list
//! - `match` / `~`: regex match
//! - `and`, `or`, `not`: logical operators

mod lexer;
mod parser;
mod evaluator;

pub use evaluator::FilterContext;
pub use parser::FilterExpr;

use anyhow::Result;

/// A compiled filter ready for evaluation.
pub struct Filter {
    expr: FilterExpr,
}

impl Filter {
    /// Parse a filter expression string.
    pub fn parse(input: &str) -> Result<Self> {
        let tokens = lexer::tokenize(input)?;
        let expr = parser::parse(tokens)?;
        Ok(Self { expr })
    }

    /// Evaluate this filter against a context of field values.
    /// Returns true if the variant passes the filter.
    pub fn matches(&self, ctx: &FilterContext) -> bool {
        evaluator::evaluate(&self.expr, ctx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ctx(fields: &[(&str, &str)]) -> FilterContext {
        let mut c = FilterContext::new();
        for (k, v) in fields {
            c.set(k, v);
        }
        c
    }

    #[test]
    fn test_simple_equality() {
        let f = Filter::parse("IMPACT is HIGH").unwrap();
        assert!(f.matches(&ctx(&[("IMPACT", "HIGH")])));
        assert!(!f.matches(&ctx(&[("IMPACT", "LOW")])));
    }

    #[test]
    fn test_numeric_comparison() {
        let f = Filter::parse("AF < 0.01").unwrap();
        assert!(f.matches(&ctx(&[("AF", "0.005")])));
        assert!(!f.matches(&ctx(&[("AF", "0.05")])));
    }

    #[test]
    fn test_and_or() {
        let f = Filter::parse("IMPACT is HIGH and AF < 0.01").unwrap();
        assert!(f.matches(&ctx(&[("IMPACT", "HIGH"), ("AF", "0.005")])));
        assert!(!f.matches(&ctx(&[("IMPACT", "HIGH"), ("AF", "0.05")])));
        assert!(!f.matches(&ctx(&[("IMPACT", "LOW"), ("AF", "0.005")])));

        let f = Filter::parse("IMPACT is HIGH or IMPACT is MODERATE").unwrap();
        assert!(f.matches(&ctx(&[("IMPACT", "HIGH")])));
        assert!(f.matches(&ctx(&[("IMPACT", "MODERATE")])));
        assert!(!f.matches(&ctx(&[("IMPACT", "LOW")])));
    }

    #[test]
    fn test_not() {
        let f = Filter::parse("not IMPACT is LOW").unwrap();
        assert!(f.matches(&ctx(&[("IMPACT", "HIGH")])));
        assert!(!f.matches(&ctx(&[("IMPACT", "LOW")])));
    }

    #[test]
    fn test_in_operator() {
        let f = Filter::parse("Consequence in missense_variant,stop_gained,frameshift_variant").unwrap();
        assert!(f.matches(&ctx(&[("Consequence", "missense_variant")])));
        assert!(f.matches(&ctx(&[("Consequence", "stop_gained")])));
        assert!(!f.matches(&ctx(&[("Consequence", "synonymous_variant")])));
    }

    #[test]
    fn test_parentheses() {
        let f = Filter::parse("(IMPACT is HIGH or IMPACT is MODERATE) and AF < 0.01").unwrap();
        assert!(f.matches(&ctx(&[("IMPACT", "HIGH"), ("AF", "0.005")])));
        assert!(f.matches(&ctx(&[("IMPACT", "MODERATE"), ("AF", "0.001")])));
        assert!(!f.matches(&ctx(&[("IMPACT", "LOW"), ("AF", "0.001")])));
        assert!(!f.matches(&ctx(&[("IMPACT", "HIGH"), ("AF", "0.05")])));
    }

    #[test]
    fn test_inequality() {
        let f = Filter::parse("IMPACT != LOW").unwrap();
        assert!(f.matches(&ctx(&[("IMPACT", "HIGH")])));
        assert!(!f.matches(&ctx(&[("IMPACT", "LOW")])));
    }

    #[test]
    fn test_missing_field() {
        // Missing fields are empty strings — comparisons should handle gracefully
        let f = Filter::parse("AF < 0.01").unwrap();
        assert!(!f.matches(&ctx(&[]))); // AF not present -> not < 0.01
    }
}
