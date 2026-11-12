//! Error type for OrbitCore.
//!
//! One error enum for the whole workspace, so every crate reports failures
//! the same way. It is deliberately free of I/O: `orbitcore-core` contains
//! only math and physics, so there is no `std::io::Error` variant here —
//! callers that read files or the network map their own errors into
//! [`OrbitCoreError::InvalidValue`] with context.
//!
//! The type is `std`-based today; the variants carry no heap allocation
//! beyond one `String` detail, keeping a future `no_std` + `alloc` port
//! straightforward.

use std::fmt;

/// The single error type for all OrbitCore crates.
#[derive(Debug, Clone, PartialEq)]
pub enum OrbitCoreError {
    /// A numeric quantity fell outside its valid range
    /// (e.g. negative altitude, efficiency above 1).
    OutOfRange {
        /// What was out of range, e.g. `"altitude"`.
        quantity: &'static str,
        /// The offending value.
        value: f64,
        /// Minimum allowed value (inclusive, unless noted otherwise).
        min: f64,
        /// Maximum allowed value (inclusive, unless noted otherwise).
        max: f64,
    },
    /// A value is invalid in a way a range check cannot express
    /// (e.g. a TLE line with the wrong length).
    InvalidValue {
        /// What was invalid, e.g. `"TLE line 1"`.
        what: &'static str,
        /// Human-readable detail.
        detail: String,
    },
    /// Input text could not be parsed (TLE, OMM, mission file).
    Parse {
        /// What was being parsed, e.g. `"TLE"`.
        what: &'static str,
        /// Human-readable detail.
        detail: String,
    },
    /// A requested feature is not implemented yet.
    Unsupported {
        /// The missing feature, e.g. `"SGP4 propagation"`.
        feature: &'static str,
    },
}

impl OrbitCoreError {
    /// Convenience constructor for [`OrbitCoreError::OutOfRange`].
    pub fn out_of_range(quantity: &'static str, value: f64, min: f64, max: f64) -> Self {
        Self::OutOfRange {
            quantity,
            value,
            min,
            max,
        }
    }

    /// Convenience constructor for [`OrbitCoreError::InvalidValue`].
    pub fn invalid_value(what: &'static str, detail: impl Into<String>) -> Self {
        Self::InvalidValue {
            what,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`OrbitCoreError::Parse`].
    pub fn parse(what: &'static str, detail: impl Into<String>) -> Self {
        Self::Parse {
            what,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`OrbitCoreError::Unsupported`].
    pub fn unsupported(feature: &'static str) -> Self {
        Self::Unsupported { feature }
    }
}

impl fmt::Display for OrbitCoreError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OutOfRange {
                quantity,
                value,
                min,
                max,
            } => write!(f, "{quantity} = {value} is out of range [{min}, {max}]"),
            Self::InvalidValue { what, detail } => write!(f, "invalid {what}: {detail}"),
            Self::Parse { what, detail } => write!(f, "failed to parse {what}: {detail}"),
            Self::Unsupported { feature } => write!(f, "{feature} is not implemented yet"),
        }
    }
}

impl std::error::Error for OrbitCoreError {}

/// `Result` alias used across all OrbitCore crates.
pub type Result<T> = std::result::Result<T, OrbitCoreError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_out_of_range() {
        let err = OrbitCoreError::out_of_range("altitude", -100.0, 0.0, 2000.0);
        assert_eq!(err.to_string(), "altitude = -100 is out of range [0, 2000]");
    }

    #[test]
    fn display_invalid_value_and_parse() {
        assert_eq!(
            OrbitCoreError::invalid_value("TLE line 1", "line is 68 chars, expected 69")
                .to_string(),
            "invalid TLE line 1: line is 68 chars, expected 69"
        );
        assert_eq!(
            OrbitCoreError::parse("mission.toml", "unexpected token at line 3").to_string(),
            "failed to parse mission.toml: unexpected token at line 3"
        );
        assert_eq!(
            OrbitCoreError::unsupported("SGP4 propagation").to_string(),
            "SGP4 propagation is not implemented yet"
        );
    }

    #[test]
    fn error_is_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<OrbitCoreError>();
    }
}
