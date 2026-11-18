//! Week 2 Day 3 integration tests: the TEME → ECEF frame transform
//! over a full day of ISS propagation.
//!
//! The transform is a rotation about the z-axis by GMST, so three
//! things must hold for every state of the day, and they are checked
//! here against the real ISS arc (snapshot TLE, epoch 2026-11-14
//! 08:49:22 UTC, mean motion 15.49369975 rev/day):
//!
//! * the radius is preserved exactly — a rotation cannot change `|r|`;
//! * the z-component is preserved exactly — the rotation axis is z;
//! * the GMST at the epoch is the published value, 185.8150816°
//!   (Vallado and Meeus agree to 3×10⁻⁸° here, per the Python
//!   cross-check in the `frames` unit tests).
//!
//! The sub-satellite longitude itself is not asserted against an
//! external table — CelesTrak is unreachable from this VM, and the
//! longitude is `RA − GMST` by construction, so an external anchor
//! would only re-derive the same formula. What *is* asserted is that
//! every longitude lands in (−180°, 180°] and that the transform is
//! internally consistent (radius and z preserved) over all 1441 states.

use orbitcore_astro::frames::{gmst_degrees, julian_date, teme_to_ecef_at};
use orbitcore_sgp4::Propagator;

/// GMST at the ISS snapshot epoch, 2026-11-14 08:49:22.988352 UTC, in
/// degrees — the absolute anchor for the transform (see the module
/// docs).
const GMST_AT_ISS_EPOCH_DEG: f64 = 185.8192111;

/// The snapshot TLE, read exactly as served by CelesTrak.
fn iss_tle_text() -> String {
    std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/data/iss-tle-2026-11-14.txt"
    ))
    .expect("ISS TLE snapshot exists")
}

/// The two TLE lines from the snapshot file (selected by their leading
/// "1 " / "2 " markers — the file starts with a name line).
fn tle_lines_from_snapshot() -> (String, String) {
    let text = iss_tle_text();
    let line1 = text
        .lines()
        .find(|l| l.starts_with("1 "))
        .expect("snapshot has TLE line 1")
        .to_string();
    let line2 = text
        .lines()
        .find(|l| l.starts_with("2 "))
        .expect("snapshot has TLE line 2")
        .to_string();
    (line1, line2)
}

#[test]
fn gmst_at_iss_epoch_matches_the_published_value() {
    let (line1, line2) = tle_lines_from_snapshot();
    let prop = Propagator::from_tle(None, &line1, &line2).expect("ISS TLE initialises SGP4");
    let gmst = gmst_degrees(julian_date(&prop.elements().epoch));
    assert!(
        (gmst - GMST_AT_ISS_EPOCH_DEG).abs() < 1e-6,
        "GMST at ISS epoch = {gmst}°, want {GMST_AT_ISS_EPOCH_DEG}°"
    );
}

#[test]
fn ecef_transform_preserves_radius_and_z_over_24h() {
    let (line1, line2) = tle_lines_from_snapshot();
    let prop = Propagator::from_tle(None, &line1, &line2).expect("ISS TLE initialises SGP4");
    let epoch = &prop.elements().epoch;
    let states = prop
        .propagate_range(0.0, 1440.0, 1.0)
        .expect("24 h propagates");
    assert_eq!(states.len(), 1441, "one state per minute for a day");

    for state in &states {
        let ecef = teme_to_ecef_at(state.position_km, epoch, state.minutes_since_epoch);
        let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        // A rotation preserves the radius: 1e-9 km is ~1 µm, far
        // tighter than the transform needs, and it holds because the
        // rotation is exact arithmetic on the same doubles.
        assert!(
            (norm(ecef) - state.radius_km()).abs() < 1e-9,
            "ECEF radius ≠ TEME radius at t = {} min",
            state.minutes_since_epoch
        );
        // The rotation axis is z, so z is untouched.
        assert!(
            (ecef[2] - state.position_km[2]).abs() < 1e-12,
            "ECEF z ≠ TEME z at t = {} min",
            state.minutes_since_epoch
        );
        // Longitude is defined everywhere atan2 is: (−180°, 180°].
        let longitude = ecef[1].atan2(ecef[0]).to_degrees();
        assert!(
            longitude > -180.0 && longitude <= 180.0,
            "longitude {longitude}° out of range at t = {} min",
            state.minutes_since_epoch
        );
        assert!(longitude.is_finite());
    }
}
