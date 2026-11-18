//! `orbitcore-astro`: orbital mechanics for OrbitCore.
//!
//! This crate turns real tracking data into numbers the rest of OrbitCore
//! can use: it parses two-line element sets (TLEs) and orbit
//! mean-elements messages (OMMs, the JSON format CelesTrak serves) —
//! both into the same [`elements::OrbitalElements`] — and it converts
//! time and frames: Julian dates and Greenwich Mean Sidereal Time plus
//! the TEME → ECEF rotation in [`frames`], so a propagated state can be
//! placed over the Earth. Ground-station passes and eclipse modelling
//! land in the coming weeks of Fase 1.
//!
//! Like `orbitcore-core`, the parsing and the frame math here are pure
//! and deterministic — `tle::parse` and `omm::parse` take text and
//! return data, `frames::teme_to_ecef` takes a vector and an angle and
//! returns a vector, with no files and no network — so everything is
//! testable offline against snapshots.

pub mod elements;
pub mod frames;
pub mod omm;
pub mod tle;
