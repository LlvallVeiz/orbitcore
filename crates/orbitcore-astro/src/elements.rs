//! Shared orbital elements: the common ground between TLE and OMM.
//!
//! A TLE (two-line element set) and an OMM (orbit mean-elements message,
//! the JSON format CelesTrak serves) describe the *same* orbit in two
//! different encodings — TLE with fixed-width fields and implied decimal
//! points, OMM with plain JSON numbers. Both parse into this struct, so
//! the rest of OrbitCore never has to care which format an element set
//! arrived in, and a TLE and an OMM for the same satellite at the same
//! epoch can be compared field by field (see the cross-format tests in
//! [`crate::omm`]).
//!
//! The derived quantities — orbital period, semi-major axis, altitude —
//! are computed here once, from the mean motion, with Kepler's third law,
//! and reused by both parsers.

use crate::tle::EpochDateTime;
use orbitcore_core::units::{EARTH_MEAN_RADIUS_KM, MU_EARTH_KM3_PER_S2};

/// The orbital elements common to a TLE and an OMM: one orbit at one epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct OrbitalElements {
    /// NORAD catalog number (satellite number).
    pub norad_cat_id: u32,
    /// Epoch as a calendar date/time (UTC).
    pub epoch: EpochDateTime,
    /// Epoch as a fractional day of year (1.0 = Jan 1, 00:00 UTC), the
    /// same convention a TLE epoch uses.
    pub epoch_day: f64,
    /// Inclination, degrees.
    pub inclination_deg: f64,
    /// Right ascension of the ascending node, degrees.
    pub raan_deg: f64,
    /// Eccentricity (dimensionless).
    pub eccentricity: f64,
    /// Argument of perigee, degrees.
    pub arg_perigee_deg: f64,
    /// Mean anomaly, degrees.
    pub mean_anomaly_deg: f64,
    /// Mean motion, revolutions per day.
    pub mean_motion_rev_per_day: f64,
    /// BSTAR drag term, 1/earth-radii. Zero when the source did not carry
    /// one (OMM may omit it).
    pub bstar: f64,
    /// Revolution number at epoch, when the source carries one.
    pub revolution_number: Option<u32>,
}

impl OrbitalElements {
    /// Orbital period in minutes, from mean motion: `T = 1440 / n`.
    pub fn orbital_period_minutes(&self) -> f64 {
        1440.0 / self.mean_motion_rev_per_day
    }

    /// Semi-major axis in km, from mean motion via Kepler's third law:
    /// `a = (μ / n²)^(1/3)` with `n` in rad/s.
    pub fn semi_major_axis_km(&self) -> f64 {
        let n_rad_per_s = self.mean_motion_rev_per_day * 2.0 * std::f64::consts::PI / 86_400.0;
        (MU_EARTH_KM3_PER_S2 / n_rad_per_s.powi(2)).cbrt()
    }

    /// Mean altitude above the Earth's surface in km (`a − R`).
    ///
    /// Uses the mean Earth radius from `orbitcore-core`; the value is
    /// ~7 km above the equatorial-radius figure, which is the convention
    /// difference, not an error.
    pub fn altitude_km(&self) -> f64 {
        self.semi_major_axis_km() - EARTH_MEAN_RADIUS_KM
    }
}
