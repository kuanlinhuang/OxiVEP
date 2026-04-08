//! Recursive descent parser for filter expressions.
//!
//! Produces a FilterExpr AST from a token stream.
//!
//! Grammar:
//!   expr     → or_expr
//!   or_expr  → and_expr ("or" and_expr)*
//!   and_expr → not_expr ("and" not_expr)*
//!   not_expr → "not" not_expr | primary
//!   primary  → "(" expr ")" | comparison
//!   comparison → FIELD operator value

use crate::lexer::Token;
use anyhow::{bail, Result};

/// Filter expression AST node.
#[derive(Debug, Clone)]
pub enum FilterExpr {
    /// field == value
    Eq { field: String, value: String },
    /// field != value
    Ne { field: String, value: String },
    /// field < number
    Lt { field: String, value: f64 },
    /// field > number
    Gt { field: String, value: f64 },
    /// field <= number
    Le { field: String, value: f64 },
    /// field >= number
    Ge { field: String, value: f64 },
    /// field is one of the comma-separated values
    In { field: String, values: Vec<String> },
    /// field matches regex pattern
    Match { field: String, pattern: String },
    /// expr AND expr
    And(Box<FilterExpr>, Box<FilterExpr>),
    /// expr OR expr
    Or(Box<FilterExpr>, Box<FilterExpr>),
    /// NOT expr
    Not(Box<FilterExpr>),
}

/// Parse a token stream into a FilterExpr.
pub fn parse(tokens: Vec<Token>) -> Result<FilterExpr> {
    let mut pos = 0;
    let expr = parse_or_expr(&tokens, &mut pos)?;
    if pos < tokens.len() {
        bail!("Unexpected tokens after expression at position {}", pos);
    }
    Ok(expr)
}

fn parse_or_expr(tokens: &[Token], pos: &mut usize) -> Result<FilterExpr> {
    let mut left = parse_and_expr(tokens, pos)?;
    while *pos < tokens.len() && tokens[*pos] == Token::Or {
        *pos += 1;
        let right = parse_and_expr(tokens, pos)?;
        left = FilterExpr::Or(Box::new(left), Box::new(right));
    }
    Ok(left)
}

fn parse_and_expr(tokens: &[Token], pos: &mut usize) -> Result<FilterExpr> {
    let mut left = parse_not_expr(tokens, pos)?;
    while *pos < tokens.len() && tokens[*pos] == Token::And {
        *pos += 1;
        let right = parse_not_expr(tokens, pos)?;
        left = FilterExpr::And(Box::new(left), Box::new(right));
    }
    Ok(left)
}

fn parse_not_expr(tokens: &[Token], pos: &mut usize) -> Result<FilterExpr> {
    if *pos < tokens.len() && tokens[*pos] == Token::Not {
        *pos += 1;
        let inner = parse_not_expr(tokens, pos)?;
        Ok(FilterExpr::Not(Box::new(inner)))
    } else {
        parse_primary(tokens, pos)
    }
}

fn parse_primary(tokens: &[Token], pos: &mut usize) -> Result<FilterExpr> {
    if *pos >= tokens.len() {
        bail!("Unexpected end of expression");
    }

    // Parenthesized expression
    if tokens[*pos] == Token::LParen {
        *pos += 1;
        let expr = parse_or_expr(tokens, pos)?;
        if *pos >= tokens.len() || tokens[*pos] != Token::RParen {
            bail!("Expected closing parenthesis");
        }
        *pos += 1;
        return Ok(expr);
    }

    // Comparison: FIELD op VALUE
    let field = match &tokens[*pos] {
        Token::Field(f) => f.clone(),
        other => bail!("Expected field name, got {:?}", other),
    };
    *pos += 1;

    if *pos >= tokens.len() {
        bail!("Expected operator after field '{}'", field);
    }

    let op = &tokens[*pos];
    *pos += 1;

    match op {
        Token::Is => {
            let value = consume_value(tokens, pos)?;
            Ok(FilterExpr::Eq { field, value })
        }
        Token::Ne => {
            let value = consume_value(tokens, pos)?;
            Ok(FilterExpr::Ne { field, value })
        }
        Token::Lt => {
            let num = consume_number(tokens, pos)?;
            Ok(FilterExpr::Lt { field, value: num })
        }
        Token::Gt => {
            let num = consume_number(tokens, pos)?;
            Ok(FilterExpr::Gt { field, value: num })
        }
        Token::Le => {
            let num = consume_number(tokens, pos)?;
            Ok(FilterExpr::Le { field, value: num })
        }
        Token::Ge => {
            let num = consume_number(tokens, pos)?;
            Ok(FilterExpr::Ge { field, value: num })
        }
        Token::In => {
            let list = consume_value(tokens, pos)?;
            let values: Vec<String> = list.split(',').map(|s| s.trim().to_string()).collect();
            Ok(FilterExpr::In { field, values })
        }
        Token::Match => {
            let pattern = consume_value(tokens, pos)?;
            Ok(FilterExpr::Match { field, pattern })
        }
        other => bail!("Expected operator, got {:?}", other),
    }
}

fn consume_value(tokens: &[Token], pos: &mut usize) -> Result<String> {
    if *pos >= tokens.len() {
        bail!("Expected value");
    }
    let val = match &tokens[*pos] {
        Token::Value(v) => v.clone(),
        Token::Number(n) => n.to_string(),
        Token::Field(f) => f.clone(), // Treat as value in this position
        other => bail!("Expected value, got {:?}", other),
    };
    *pos += 1;
    Ok(val)
}

fn consume_number(tokens: &[Token], pos: &mut usize) -> Result<f64> {
    if *pos >= tokens.len() {
        bail!("Expected number");
    }
    let num = match &tokens[*pos] {
        Token::Number(n) => *n,
        Token::Value(v) => v.parse::<f64>().map_err(|_| anyhow::anyhow!("Expected number, got '{}'", v))?,
        other => bail!("Expected number, got {:?}", other),
    };
    *pos += 1;
    Ok(num)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lexer::tokenize;

    #[test]
    fn test_parse_simple() {
        let tokens = tokenize("IMPACT is HIGH").unwrap();
        let expr = parse(tokens).unwrap();
        match expr {
            FilterExpr::Eq { field, value } => {
                assert_eq!(field, "IMPACT");
                assert_eq!(value, "HIGH");
            }
            _ => panic!("Expected Eq"),
        }
    }

    #[test]
    fn test_parse_and() {
        let tokens = tokenize("IMPACT is HIGH and AF < 0.01").unwrap();
        let expr = parse(tokens).unwrap();
        match expr {
            FilterExpr::And(_, _) => {}
            _ => panic!("Expected And"),
        }
    }

    #[test]
    fn test_parse_or() {
        let tokens = tokenize("IMPACT is HIGH or IMPACT is MODERATE").unwrap();
        let expr = parse(tokens).unwrap();
        match expr {
            FilterExpr::Or(_, _) => {}
            _ => panic!("Expected Or"),
        }
    }

    #[test]
    fn test_parse_not() {
        let tokens = tokenize("not IMPACT is LOW").unwrap();
        let expr = parse(tokens).unwrap();
        match expr {
            FilterExpr::Not(_) => {}
            _ => panic!("Expected Not"),
        }
    }

    #[test]
    fn test_parse_parens() {
        let tokens = tokenize("(IMPACT is HIGH or IMPACT is MODERATE) and AF < 0.01").unwrap();
        let expr = parse(tokens).unwrap();
        match expr {
            FilterExpr::And(left, _) => match *left {
                FilterExpr::Or(_, _) => {}
                _ => panic!("Expected Or inside And"),
            },
            _ => panic!("Expected And"),
        }
    }
}
