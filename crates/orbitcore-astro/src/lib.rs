//! `orbitcore-astro`: orbital mechanics for OrbitCore.
//!
//! This crate turns real tracking data into numbers the rest of OrbitCore
//! can use: it parses two-line element sets (TLEs) and orbit
//! mean-elements messages (OMMs, the JSON format CelesTrak serves) —
//! both into the same [`elements::OrbitalElements`] — and in the coming
//! days of Fase 1 it will propagate them with SGP4, transform frames
//! (TEME → ECEF → geodetic), predict ground-station passes, and model
//! eclipses and beta angles.
//!
//! Like `orbitcore-core`, the parsing here is pure and deterministic —
//! `tle::parse` and `omm::parse` take text and return data, with no
//! files and no network — so it is testable offline against snapshots.

pub mod elements;
pub mod omm;
pub mod tle;
