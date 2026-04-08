//! Filter expression evaluator.
//!
//! Evaluates a FilterExpr AST against a context of field name -> value mappings.

use crate::parser::FilterExpr;
use std::collections::HashMap;

/// Context for filter evaluation: maps field names to their string values.
pub struct FilterContext {
    fields: HashMap<String, String>,
}

impl FilterContext {
    pub fn new() -> Self {
        Self {
            fields: HashMap::new(),
        }
    }

    pub fn set(&mut self, field: &str, value: &str) {
        self.fields.insert(field.to_string(), value.to_string());
    }

    pub fn get(&self, field: &str) -> &str {
        self.fields.get(field).map(|s| s.as_str()).unwrap_or("")
    }
}

impl Default for FilterContext {
    fn default() -> Self {
        Self::new()
    }
}

/// Evaluate a filter expression against a context.
pub fn evaluate(expr: &FilterExpr, ctx: &FilterContext) -> bool {
    match expr {
        FilterExpr::Eq { field, value } => {
            let actual = ctx.get(field);
            // Case-insensitive comparison for text fields
            actual.eq_ignore_ascii_case(value)
        }
        FilterExpr::Ne { field, value } => {
            let actual = ctx.get(field);
            !actual.eq_ignore_ascii_case(value)
        }
        FilterExpr::Lt { field, value } => {
            let actual = ctx.get(field);
            actual.parse::<f64>().map_or(false, |a| a < *value)
        }
        FilterExpr::Gt { field, value } => {
            let actual = ctx.get(field);
            actual.parse::<f64>().map_or(false, |a| a > *value)
        }
        FilterExpr::Le { field, value } => {
            let actual = ctx.get(field);
            actual.parse::<f64>().map_or(false, |a| a <= *value)
        }
        FilterExpr::Ge { field, value } => {
            let actual = ctx.get(field);
            actual.parse::<f64>().map_or(false, |a| a >= *value)
        }
        FilterExpr::In { field, values } => {
            let actual = ctx.get(field);
            // The actual value might be a single value or ampersand-separated list
            // (VEP uses & to separate multiple consequences)
            let actual_parts: Vec<&str> = actual.split('&').collect();
            actual_parts.iter().any(|a| {
                values.iter().any(|v| a.eq_ignore_ascii_case(v))
            })
        }
        FilterExpr::Match { field, pattern } => {
            let actual = ctx.get(field);
            // Simple glob-style matching (contains check for basic use)
            actual.contains(pattern)
        }
        FilterExpr::And(left, right) => {
            evaluate(left, ctx) && evaluate(right, ctx)
        }
        FilterExpr::Or(left, right) => {
            evaluate(left, ctx) || evaluate(right, ctx)
        }
        FilterExpr::Not(inner) => {
            !evaluate(inner, ctx)
        }
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
    fn test_eq() {
        let expr = FilterExpr::Eq {
            field: "IMPACT".into(),
            value: "HIGH".into(),
        };
        assert!(evaluate(&expr, &ctx(&[("IMPACT", "HIGH")])));
        assert!(evaluate(&expr, &ctx(&[("IMPACT", "high")]))); // case-insensitive
        assert!(!evaluate(&expr, &ctx(&[("IMPACT", "LOW")])));
    }

    #[test]
    fn test_numeric() {
        let expr = FilterExpr::Lt {
            field: "AF".into(),
            value: 0.01,
        };
        assert!(evaluate(&expr, &ctx(&[("AF", "0.005")])));
        assert!(!evaluate(&expr, &ctx(&[("AF", "0.05")])));
        assert!(!evaluate(&expr, &ctx(&[("AF", "not_a_number")]))); // graceful failure
    }

    #[test]
    fn test_in_with_ampersand_separated() {
        // VEP consequence fields use & for multiple values
        let expr = FilterExpr::In {
            field: "Consequence".into(),
            values: vec!["missense_variant".into(), "stop_gained".into()],
        };
        assert!(evaluate(&expr, &ctx(&[("Consequence", "missense_variant&splice_region_variant")])));
        assert!(!evaluate(&expr, &ctx(&[("Consequence", "synonymous_variant")])));
    }

    #[test]
    fn test_and_or_not() {
        let expr = FilterExpr::And(
            Box::new(FilterExpr::Eq { field: "IMPACT".into(), value: "HIGH".into() }),
            Box::new(FilterExpr::Not(
                Box::new(FilterExpr::Eq { field: "CANONICAL".into(), value: "NO".into() }),
            )),
        );
        assert!(evaluate(&expr, &ctx(&[("IMPACT", "HIGH"), ("CANONICAL", "YES")])));
        assert!(!evaluate(&expr, &ctx(&[("IMPACT", "HIGH"), ("CANONICAL", "NO")])));
    }
}
