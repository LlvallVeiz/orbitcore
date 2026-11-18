//! `orbitcore-sgp4`: SGP4 propagation for OrbitCore.
//!
//! This crate is the thin wrapper the project plan (§5.3, §7 Week 2 Day 1)
//! calls for: it turns the [`OrbitalElements`] that `orbitcore-astro`
//! parses from a TLE or OMM into propagated state vectors — position and
//! velocity in the TEME frame, kilometres and kilometres per second —
//! without re-implementing SGP4 itself. The propagation math lives in the
//! upstream `sgp4` crate (v2.4.0, neuromorphicsystems), validated against
//! CelesTrak to under 2×10⁻⁷ km in position; OrbitCore's rule is "don't
//! rewrite SGP4, the specification is its quirks", so this crate only
//! adapts types, maps errors into [`OrbitCoreError`], and adds the
//! conveniences the rest of OrbitCore needs (a range propagator, a
//! [`StateVector`] with radius/speed/altitude helpers).
//!
//! Two ways to build a [`Propagator`]:
//!
//! * [`Propagator::from_elements`] — from OrbitCore's own
//!   [`OrbitalElements`]. The fields a TLE/OMM carries that
//!   `OrbitalElements` does not (mean-motion derivatives, element-set
//!   number, classification) default to zero/unclassified; the drag term
//!   (BSTAR) is preserved. This is the path the rest of OrbitCore uses.
//! * [`Propagator::from_tle`] — parse a raw TLE with the upstream
//!   parser, so every field (including the mean-motion derivatives)
//!   survives exactly. Use it when the original TLE text is at hand and
//!   maximum fidelity matters.
//!
//! Both paths agree exactly over a day for the ISS (see
//! `tests/iss_24h.rs`): the only fields they differ in are the
//! mean-motion derivatives, and the upstream `sgp4` crate (v2.4.0)
//! parses those fields but never uses them in propagation — only BSTAR
//! feeds the model. The agreement is measured, not assumed.
//!
//! On top of the wrapper, the [`Ephemeris`] type (see
//! [`ephemeris`](crate::ephemeris)) turns a day of propagated states
//! into orbit facts: revolutions made, ascending-node crossings and
//! the nodal period, altitude and speed extremes. It is built from a
//! [`Propagator`] with no files and no network — the 24-hour checks in
//! `tests/iss_24h.rs` run offline against the CelesTrak snapshots in
//! `validation/data/`.
//!
//! A propagated state is in TEME; [`orbitcore_astro::frames`] rotates
//! it into the Earth-fixed frame (ECEF) — Week 2 Day 3 — with the
//! integration check in `tests/iss_frames.rs`.
//!
//! Like the rest of the workspace, propagation here is pure and
//! deterministic: `propagate` takes minutes since epoch and returns a
//! state vector, with no files and no network.
//!
//! [`OrbitalElements`]: orbitcore_astro::elements::OrbitalElements
//! [`OrbitCoreError`]: orbitcore_core::error::OrbitCoreError

pub mod ephemeris;
pub mod propagator;
pub mod state;

pub use ephemeris::Ephemeris;
pub use propagator::Propagator;
pub use state::StateVector;
