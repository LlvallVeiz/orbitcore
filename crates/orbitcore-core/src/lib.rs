//! `orbitcore-core`: shared core of OrbitCore.
//!
//! Pure math and physics only — no network, no files, no system clock — so
//! this crate compiles for the mission simulator on a laptop and for the
//! flight software on a microcontroller (`no_std`) alike.
//!
//! Day 2 adds the SI-units module ([`units`]) and the shared error type
//! ([`error`]); the mission-file types land in Fase 2 (Week 8) with the
//! `mission.toml` work.

pub mod error;
pub mod units;

pub use error::{OrbitCoreError, Result};

/// Crate version, kept in sync with the workspace version.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_set() {
        assert!(!VERSION.is_empty(), "crate version must not be empty");
        assert!(
            VERSION.starts_with("0.1.0"),
            "unexpected version: {VERSION}"
        );
    }
}
