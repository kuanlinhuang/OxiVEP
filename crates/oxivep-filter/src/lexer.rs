//! Tokenizer for filter expressions.

use anyhow::{bail, Result};

#[derive(Debug, Clone, PartialEq)]
pub enum Token {
    Field(String),     // Field name (e.g., IMPACT, AF, Consequence)
    Value(String),     // String value (e.g., HIGH, missense_variant)
    Number(f64),       // Numeric value
    // Comparison operators
    Is,                // is, =, eq
    Ne,                // !=, ne
    Lt,                // <
    Gt,                // >
    Le,                // <=
    Ge,                // >=
    In,                // in
    Match,             // match, ~
    // Logical operators
    And,
    Or,
    Not,
    // Grouping
    LParen,
    RParen,
}

/// Tokenize a filter expression string.
pub fn tokenize(input: &str) -> Result<Vec<Token>> {
    let mut tokens = Vec::new();
    let mut chars = input.chars().peekable();

    while let Some(&ch) = chars.peek() {
        match ch {
            ' ' | '\t' | '\n' | '\r' => {
                chars.next();
            }
            '(' => {
                tokens.push(Token::LParen);
                chars.next();
            }
            ')' => {
                tokens.push(Token::RParen);
                chars.next();
            }
            '<' => {
                chars.next();
                if chars.peek() == Some(&'=') {
                    chars.next();
                    tokens.push(Token::Le);
                } else {
                    tokens.push(Token::Lt);
                }
            }
            '>' => {
                chars.next();
                if chars.peek() == Some(&'=') {
                    chars.next();
                    tokens.push(Token::Ge);
                } else {
                    tokens.push(Token::Gt);
                }
            }
            '!' => {
                chars.next();
                if chars.peek() == Some(&'=') {
                    chars.next();
                    tokens.push(Token::Ne);
                } else {
                    tokens.push(Token::Not);
                }
            }
            '=' => {
                chars.next();
                tokens.push(Token::Is);
            }
            '~' => {
                chars.next();
                tokens.push(Token::Match);
            }
            _ if ch.is_ascii_digit() || ch == '-' || ch == '.' => {
                // Could be a number or a negative number / field starting with digit
                let word = consume_word(&mut chars);
                if let Ok(n) = word.parse::<f64>() {
                    tokens.push(Token::Number(n));
                } else {
                    tokens.push(Token::Value(word));
                }
            }
            _ if is_word_char(ch) => {
                let word = consume_word(&mut chars);
                match word.to_lowercase().as_str() {
                    "and" => tokens.push(Token::And),
                    "or" => tokens.push(Token::Or),
                    "not" => tokens.push(Token::Not),
                    "is" | "eq" => tokens.push(Token::Is),
                    "ne" => tokens.push(Token::Ne),
                    "in" => tokens.push(Token::In),
                    "match" => tokens.push(Token::Match),
                    "lt" => tokens.push(Token::Lt),
                    "gt" => tokens.push(Token::Gt),
                    "le" => tokens.push(Token::Le),
                    "ge" => tokens.push(Token::Ge),
                    _ => {
                        // Determine if this is a field name or value based on context.
                        // If the previous token is a comparison operator, this is a value.
                        // Otherwise, it's a field name.
                        if is_value_position(&tokens) {
                            tokens.push(Token::Value(word));
                        } else {
                            tokens.push(Token::Field(word));
                        }
                    }
                }
            }
            '"' | '\'' => {
                // Quoted string value
                let quote = ch;
                chars.next();
                let mut s = String::new();
                while let Some(&c) = chars.peek() {
                    if c == quote {
                        chars.next();
                        break;
                    }
                    s.push(c);
                    chars.next();
                }
                tokens.push(Token::Value(s));
            }
            _ => {
                bail!("Unexpected character in filter expression: '{}'", ch);
            }
        }
    }

    Ok(tokens)
}

fn is_word_char(ch: char) -> bool {
    ch.is_alphanumeric() || ch == '_' || ch == '.' || ch == '-' || ch == '/' || ch == ':'
}

fn consume_word(chars: &mut std::iter::Peekable<std::str::Chars>) -> String {
    let mut word = String::new();
    while let Some(&ch) = chars.peek() {
        if is_word_char(ch) || ch == ',' {
            word.push(ch);
            chars.next();
        } else {
            break;
        }
    }
    word
}

/// Check if the current position in the token stream expects a value
/// (i.e., the last token was a comparison operator or 'in').
fn is_value_position(tokens: &[Token]) -> bool {
    matches!(
        tokens.last(),
        Some(Token::Is) | Some(Token::Ne) | Some(Token::Lt) | Some(Token::Gt)
            | Some(Token::Le) | Some(Token::Ge) | Some(Token::In) | Some(Token::Match)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_tokens() {
        let tokens = tokenize("IMPACT is HIGH").unwrap();
        assert_eq!(tokens, vec![
            Token::Field("IMPACT".into()),
            Token::Is,
            Token::Value("HIGH".into()),
        ]);
    }

    #[test]
    fn test_numeric_comparison() {
        let tokens = tokenize("AF < 0.01").unwrap();
        assert_eq!(tokens, vec![
            Token::Field("AF".into()),
            Token::Lt,
            Token::Number(0.01),
        ]);
    }

    #[test]
    fn test_complex_expression() {
        let tokens = tokenize("(IMPACT is HIGH or IMPACT is MODERATE) and AF < 0.01").unwrap();
        assert!(tokens.contains(&Token::LParen));
        assert!(tokens.contains(&Token::RParen));
        assert!(tokens.contains(&Token::Or));
        assert!(tokens.contains(&Token::And));
    }

    #[test]
    fn test_in_operator() {
        let tokens = tokenize("Consequence in missense_variant,stop_gained").unwrap();
        assert_eq!(tokens[0], Token::Field("Consequence".into()));
        assert_eq!(tokens[1], Token::In);
        assert_eq!(tokens[2], Token::Value("missense_variant,stop_gained".into()));
    }
}
