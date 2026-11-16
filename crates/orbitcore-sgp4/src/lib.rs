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
//! Both paths agree to well under a kilometre over a day for the ISS
//! (see `tests/iss_propagation.rs`), because the mean-motion derivatives
//! are tiny; the difference is documented, not hidden.
//!
//! Like the rest of the workspace, propagation here is pure and
//! deterministic: `propagate` takes minutes since epoch and returns a
//! state vector, with no files and no network. Day 2 builds the 24-hour
//! propagation and frame transforms on top of this wrapper.
//!
//! [`OrbitalElements`]: orbitcore_astro::elements::OrbitalElements
//! [`OrbitCoreError`]: orbitcore_core::error::OrbitCoreError

pub mod propagator;
pub mod state;

pub use propagator::Propagator;
pub use state::StateVector;
