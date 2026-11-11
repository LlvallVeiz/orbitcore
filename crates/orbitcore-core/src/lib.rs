//! `orbitcore-core`: shared core of OrbitCore.
//!
//! Pure math and physics only — no network, no files, no system clock — so
//! this crate compiles for the mission simulator on a laptop and for the
//! flight software on a microcontroller (`no_std`) alike.
//!
//! Day 1 only scaffolds the crate; the SI-units module and the error type
//! land on Day 2, followed by the mission-file types.

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
