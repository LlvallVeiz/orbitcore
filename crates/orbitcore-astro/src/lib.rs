//! `orbitcore-astro`: orbital mechanics for OrbitCore.
//!
//! This crate turns real tracking data into numbers the rest of OrbitCore
//! can use: it parses two-line element sets (TLEs) from CelesTrak, and in
//! the coming days of Fase 1 it will propagate them with SGP4, transform
//! frames (TEME → ECEF → geodetic), predict ground-station passes, and
//! model eclipses and beta angles.
//!
//! Like `orbitcore-core`, the parsing here is pure and deterministic —
//! `tle::parse` takes text and returns a [`tle::TwoLineElement`], with no
//! files and no network — so it is testable offline against snapshot TLEs.

pub mod tle;
