//! Frame transforms: time (Julian date, GMST) and the TEME → ECEF rotation.
//!
//! SGP4 outputs state vectors in TEME — the True Equator, Mean Equinox
//! frame, an inertial frame tied to the epoch's equator and equinox
//! (see `orbitcore-sgp4`). Everything downstream that cares where the
//! satellite is *over the Earth* — ground-station passes (Week 3),
//! ground tracks, eclipse geometry — needs the Earth-fixed frame
//! (ECEF) instead, which rotates with the planet. This module is that
//! bridge, and it is deliberately small: two time conversions and one
//! rotation, each tested against a published anchor.
//!
//! The chain is
//!
//! ```text
//! EpochDateTime ──julian_date──▶ JD ──gmst_radians──▶ GMST
//! TEME position ──teme_to_ecef (rotate about z by GMST)──▶ ECEF position
//! ```
//!
//! Two approximations are baked in, and both are documented rather
//! than hidden:
//!
//! * **TEME → ECEF via a GMST rotation.** The rigorous IAU path
//!   (TEME → TIRS) adds precession, nutation, and the equation of the
//!   origins; the GMST rotation is the standard engineering
//!   approximation (the one satellite.js and most SGP4-based tools use),
//!   accurate to a few arcseconds — tens of metres at the surface.
//!   That is plenty for pass prediction and ground tracks, and not
//!   precision geodesy; the day the twin needs centimetres, this
//!   module is where the IAU terms land.
//! * **UT1 is taken as UTC.** The difference, |DUT1| ≤ 0.9 s by
//!   definition, shifts GMST by up to 0.0038° (~0.4 km at the
//!   equator). A DUT1 table (or an IERS bulletin) can tighten this
//!   later; for Week 2 it is a documented error, not a silent one.
//!
//! The geodetic step — ECEF → latitude/longitude/altitude on the
//! WGS-84 ellipsoid — is Week 2 Day 4, not this module: the rotation
//! is the frame transform, the ellipsoid is the datum.

use std::f64::consts::PI;

use crate::tle::EpochDateTime;

/// Julian date of the J2000.0 epoch (2000-01-01 12:00:00 UTC), the
/// reference the GMST polynomial is measured from.
pub const JULIAN_DATE_J2000: f64 = 2_451_545.0;

/// Days in a Julian century — the time unit of the GMST polynomial.
pub const DAYS_PER_JULIAN_CENTURY: f64 = 36_525.0;

/// GMST at the J2000.0 epoch, in degrees: 18h 41m 50.54841s. The
/// anchor every GMST implementation is checked against (Vallado,
/// *Fundamentals of Astrodynamics and Applications*).
pub const GMST_AT_J2000_DEG: f64 = 280.46061837;

/// Julian date (days since noon, 4713 BC) of a UTC calendar date/time.
///
/// The Gregorian-calendar algorithm from Meeus, *Astronomical
/// Algorithms* ch. 7: months January and February count as months 13
/// and 14 of the previous year, the century correction `B` handles the
/// Gregorian leap-year rule, and the day fraction is the time of day
/// over 86400 s. The result is a Julian *date* in UTC — the same
/// quantity the GMST polynomial takes as its argument, with UT1 ≈ UTC
/// (see the module docs).
///
/// # Examples
///
/// 2000-01-01 12:00:00 UTC is exactly JD 2451545.0 by definition;
/// 2000-01-01 00:00:00 UTC is half a day earlier, JD 2451544.5.
pub fn julian_date(epoch: &EpochDateTime) -> f64 {
    let (mut year, mut month) = (i64::from(epoch.year), i64::from(epoch.month));
    if month <= 2 {
        year -= 1;
        month += 12;
    }
    let a = year.div_euclid(100);
    // Gregorian century correction: 2 − A + ⌊A/4⌋.
    let b = 2 - a + a.div_euclid(4);
    let day_fraction =
        (f64::from(epoch.hour) * 3600.0 + f64::from(epoch.minute) * 60.0 + epoch.second) / 86_400.0;
    (365.25 * (year + 4716) as f64).floor()
        + (30.6001 * (month + 1) as f64).floor()
        + f64::from(epoch.day)
        + day_fraction
        + b as f64
        - 1524.5
}

/// Greenwich Mean Sidereal Time for a Julian date, in radians,
/// normalised to `[0, 2π)`.
///
/// Vallado's `gstime` polynomial (*Fundamentals of Astrodynamics and
/// Applications*): with `T` the Julian centuries since J2000.0,
/// `GMST = 67310.54841 + (876600·3600 + 8640184.812866)·T + 0.093104·T²
/// − 6.2×10⁻⁶·T³` seconds, converted to angle at 240 s per degree
/// (86400 s = 360°). The `876600·3600` term is the seconds in a Julian
/// century; the `8640184.812866` is the extra rotation the Earth makes
/// in that time because a sidereal day is shorter than a solar day
/// (360.98564736629° per solar day).
///
/// At J2000.0 this returns exactly [`GMST_AT_J2000_DEG`]; the unit
/// tests cross-check it against Meeus's independent polynomial at the
/// ISS epoch and in 2030, where the two agree to ~3×10⁻⁸ degrees.
pub fn gmst_radians(julian_date: f64) -> f64 {
    let t = (julian_date - JULIAN_DATE_J2000) / DAYS_PER_JULIAN_CENTURY;
    let gmst_seconds = 67310.54841 + (876_600.0 * 3600.0 + 8_640_184.812866) * t + 0.093104 * t * t
        - 6.2e-6 * t * t * t;
    // 240 seconds per degree: 86400 s = 360°.
    (gmst_seconds / 240.0).to_radians().rem_euclid(2.0 * PI)
}

/// Greenwich Mean Sidereal Time for a Julian date, in degrees,
/// normalised to `[0, 360)`.
pub fn gmst_degrees(julian_date: f64) -> f64 {
    gmst_radians(julian_date).to_degrees()
}

/// Rotate a TEME position into the Earth-fixed frame (ECEF): `r_ecef =
/// R₃(GMST) · r_teme`, a rotation about the z-axis by the Greenwich
/// Mean Sidereal Time.
///
/// `R₃(θ)` here is `[[cos θ, sin θ, 0], [−sin θ, cos θ, 0], [0, 0, 1]]`
/// — the matrix that *removes* the Earth's rotation: a point fixed in
/// inertial space at right ascension α ends up at ECEF longitude
/// `α − GMST` (the unit test checks exactly this, because getting the
/// sign wrong is the classic frame bug). The z-component is untouched,
/// so `|r|` and the equatorial-plane height are preserved exactly.
///
/// Positions are in kilometres; `gmst_rad` is in radians, as
/// [`gmst_radians`] returns.
pub fn teme_to_ecef(position_teme_km: [f64; 3], gmst_rad: f64) -> [f64; 3] {
    let (sin, cos) = gmst_rad.sin_cos();
    let [x, y, z] = position_teme_km;
    [cos * x + sin * y, -sin * x + cos * y, z]
}

/// The inverse of [`teme_to_ecef`]: rotate an ECEF position back into
/// TEME by `−GMST` (the transpose of the `R₃` matrix, so a round trip
/// recovers the input to floating-point precision).
pub fn ecef_to_teme(position_ecef_km: [f64; 3], gmst_rad: f64) -> [f64; 3] {
    let (sin, cos) = gmst_rad.sin_cos();
    let [x, y, z] = position_ecef_km;
    [cos * x - sin * y, sin * x + cos * y, z]
}

/// Convenience for the common case: rotate a TEME position (from
/// [`Propagator::propagate`](orbitcore_sgp4::Propagator::propagate) or
/// an [`Ephemeris`](orbitcore_sgp4::Ephemeris) state) into ECEF at the
/// instant it was propagated for — `epoch` plus `minutes_since_epoch`
/// minutes. The Julian date is `julian_date(epoch) + minutes / 1440`.
pub fn teme_to_ecef_at(
    position_teme_km: [f64; 3],
    epoch: &EpochDateTime,
    minutes_since_epoch: f64,
) -> [f64; 3] {
    let jd = julian_date(epoch) + minutes_since_epoch / 1440.0;
    teme_to_ecef(position_teme_km, gmst_radians(jd))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An `EpochDateTime` built by hand: the tests need exact calendar
    /// times, and the TLE parser's epoch conversion is tested
    /// elsewhere (`tle::tests`).
    fn epoch(year: u16, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> EpochDateTime {
        EpochDateTime {
            year,
            month,
            day,
            hour,
            minute,
            second,
        }
    }

    /// Meeus's GMST polynomial (*Astronomical Algorithms*, ch. 12) —
    /// an independent formula with different coefficients from
    /// Vallado's, used here only as a cross-check: if two published
    /// polynomials agree, the implementation is right, not just
    /// self-consistent.
    fn gmst_meeus_degrees(julian_date: f64) -> f64 {
        let t = (julian_date - JULIAN_DATE_J2000) / DAYS_PER_JULIAN_CENTURY;
        (GMST_AT_J2000_DEG
            + 360.98564736629 * (julian_date - JULIAN_DATE_J2000)
            + 0.000387933 * t * t
            - t * t * t / 38_710_000.0)
            .rem_euclid(360.0)
    }

    #[test]
    fn julian_date_of_j2000_is_2451545() {
        // By definition: 2000-01-01 12:00:00 UTC = JD 2451545.0.
        let jd = julian_date(&epoch(2000, 1, 1, 12, 0, 0.0));
        assert!((jd - 2_451_545.0).abs() < 1e-9, "JD = {jd}");
        // Midnight is half a day earlier.
        let jd_midnight = julian_date(&epoch(2000, 1, 1, 0, 0, 0.0));
        assert!(
            (jd_midnight - 2_451_544.5).abs() < 1e-9,
            "JD = {jd_midnight}"
        );
    }

    #[test]
    fn julian_date_of_iss_epoch() {
        // ISS snapshot epoch 2026-11-14 08:49:22.988352 UTC (TLE epoch
        // day 318.36762718). The JD was computed independently in
        // Python (datetime arithmetic from the J2000 anchor):
        // 2461358.86762718.
        let jd = julian_date(&epoch(2026, 11, 14, 8, 49, 22.988352));
        assert!(
            (jd - 2_461_358.86762718).abs() < 1e-6,
            "JD = {jd}, want 2461358.86762718"
        );
    }

    #[test]
    fn julian_date_handles_january_and_february() {
        // January/February count as months 13/14 of the previous year
        // in the algorithm — the branch that is easiest to get wrong.
        // 2000-02-29 (a leap day) 12:00 UTC is JD 2451604.0.
        let jd = julian_date(&epoch(2000, 2, 29, 12, 0, 0.0));
        assert!((jd - 2_451_604.0).abs() < 1e-9, "JD = {jd}");
        // 2026-01-01 12:00 UTC: 2451545.0 + 9497 days (26 years incl.
        // 7 leap days: 2000, 04, 08, 12, 16, 20, 24).
        let jd = julian_date(&epoch(2026, 1, 1, 12, 0, 0.0));
        assert!((jd - 2_461_042.0).abs() < 1e-9, "JD = {jd}");
    }

    #[test]
    fn gmst_at_j2000_is_the_published_value() {
        // 18h 41m 50.54841s = 280.46061837° at J2000.0.
        let gmst = gmst_degrees(JULIAN_DATE_J2000);
        assert!(
            (gmst - GMST_AT_J2000_DEG).abs() < 1e-6,
            "GMST = {gmst}°, want {GMST_AT_J2000_DEG}°"
        );
    }

    #[test]
    fn gmst_matches_meeus_at_iss_epoch_and_2030() {
        // Two published polynomials, different coefficients: they must
        // agree to ~3e-8° (checked in Python first); 1e-6° is the
        // assertion bound, ~0.1 m at the equator.
        for jd in [2_461_358.86762718, 2_462_668.0] {
            let vallado = gmst_degrees(jd);
            let meeus = gmst_meeus_degrees(jd);
            assert!(
                (vallado - meeus).abs() < 1e-6,
                "JD {jd}: Vallado {vallado}° vs Meeus {meeus}°"
            );
        }
    }

    #[test]
    fn gmst_at_iss_epoch_is_185_8192111_degrees() {
        // The absolute anchor for the ISS integration test: GMST at the
        // snapshot epoch, 2026-11-14 08:49:22.988352 UTC (both
        // polynomials agree to 3e-8° here, per the Python cross-check).
        let jd = julian_date(&epoch(2026, 11, 14, 8, 49, 22.988352));
        let gmst = gmst_degrees(jd);
        assert!(
            (gmst - 185.8192111).abs() < 1e-6,
            "GMST = {gmst}°, want 185.8192111°"
        );
    }

    #[test]
    fn gmst_is_normalised_to_a_full_circle() {
        for jd in [JULIAN_DATE_J2000, 2_461_358.86762718, 2_500_000.0] {
            let rad = gmst_radians(jd);
            assert!(
                (0.0..2.0 * PI).contains(&rad),
                "GMST = {rad} rad at JD {jd}"
            );
            assert!((gmst_degrees(jd) - rad.to_degrees()).abs() < 1e-12);
        }
    }

    #[test]
    fn teme_to_ecef_is_identity_at_zero_gmst() {
        let p = [1.0, 2.0, 3.0];
        let ecef = teme_to_ecef(p, 0.0);
        for (a, b) in ecef.iter().zip(p.iter()) {
            assert!((a - b).abs() < 1e-15);
        }
    }

    #[test]
    fn teme_to_ecef_rotates_90_degrees() {
        // GMST = +90°: the x-axis (RA 0°) moves to ECEF −y
        // (longitude −90°), and the y-axis moves to +x.
        let quarter = PI / 2.0;
        let ecef = teme_to_ecef([1.0, 0.0, 0.0], quarter);
        assert!((ecef[0]).abs() < 1e-15);
        assert!((ecef[1] + 1.0).abs() < 1e-15);
        assert!((ecef[2]).abs() < 1e-15);
        let ecef = teme_to_ecef([0.0, 1.0, 0.0], quarter);
        assert!((ecef[0] - 1.0).abs() < 1e-15);
        assert!((ecef[1]).abs() < 1e-15);
    }

    #[test]
    fn ecef_to_teme_is_the_inverse() {
        // A non-axis-aligned vector, at an arbitrary GMST: the round
        // trip must recover the input to floating-point precision, and
        // the rotation must preserve the radius exactly.
        let p = [5_000.0, -3_000.0, 2_000.0];
        let gmst = 3.245_678;
        let ecef = teme_to_ecef(p, gmst);
        let back = ecef_to_teme(ecef, gmst);
        for (a, b) in back.iter().zip(p.iter()) {
            assert!((a - b).abs() < 1e-9, "round trip: {back:?} vs {p:?}");
        }
        let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert!((norm(ecef) - norm(p)).abs() < 1e-9);
    }

    #[test]
    fn longitude_is_right_ascension_minus_gmst() {
        // The sign check that matters: a point on the equator at RA 0°
        // at J2000 (GMST 280.46061837°) must sit at ECEF longitude
        // 0 − 280.46061837 = −280.46061837 ≡ +79.53938163°. A flipped
        // rotation matrix would give −79.53938163° instead.
        let r = 6_796.0; // 425 km altitude, on the equator at RA 0°.
        let ecef = teme_to_ecef([r, 0.0, 0.0], gmst_radians(JULIAN_DATE_J2000));
        let longitude = ecef[1].atan2(ecef[0]).to_degrees();
        assert!(
            (longitude - 79.53938163).abs() < 1e-6,
            "longitude = {longitude}°, want 79.53938163°"
        );
    }

    #[test]
    fn teme_to_ecef_at_advances_with_earth_rotation() {
        // A fixed TEME vector's ECEF longitude must decrease as the
        // Earth turns under it: GMST advances 360.98564736629° per
        // solar day, i.e. 15.04106864° per hour.
        let iss_epoch = epoch(2026, 11, 14, 8, 49, 22.0);
        let p = [6_796.0, 0.0, 0.0];
        let lon_at = |minutes: f64| {
            let ecef = teme_to_ecef_at(p, &iss_epoch, minutes);
            ecef[1].atan2(ecef[0]).to_degrees()
        };
        let lon0 = lon_at(0.0);
        let lon60 = lon_at(60.0);
        // Unwrap the −180/180 jump before differencing.
        let mut delta = lon60 - lon0;
        if delta > 180.0 {
            delta -= 360.0;
        } else if delta < -180.0 {
            delta += 360.0;
        }
        assert!(
            (delta + 15.04106864).abs() < 1e-4,
            "longitude change over 1 h = {delta}°, want −15.04106864°"
        );
        // And `teme_to_ecef_at` at t = 0 agrees with the two-step call.
        let jd = julian_date(&iss_epoch);
        let direct = teme_to_ecef(p, gmst_radians(jd));
        let via_at = teme_to_ecef_at(p, &iss_epoch, 0.0);
        for (a, b) in direct.iter().zip(via_at.iter()) {
            assert!((a - b).abs() < 1e-9);
        }
    }
}
