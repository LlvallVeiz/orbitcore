//! Day 2 integration tests: a full 24-hour ISS ephemeris, checked
//! against the orbit the CelesTrak element set itself describes.
//!
//! The snapshot is the live ISS element set fetched 2026-11-14, in two
//! formats served by CelesTrak for the *same* element set: the TLE
//! (`validation/data/iss-tle-2026-11-14.txt`, epoch 2026-11-14
//! 08:49:22 UTC, mean motion 15.49369975 rev/day, eccentricity
//! 0.0007385, NORAD 25544) and the OMM JSON
//! (`validation/data/iss-omm-2026-11-14.json`). Everything asserted
//! here follows from those numbers plus basic orbital mechanics — no
//! network, no CelesTrak call at test time. The "vs CelesTrak" check
//! is that the propagated arc reproduces what the published element
//! set implies: the revolution count matches the published mean
//! motion, the nodal period matches `1440 / n`, and the two CelesTrak
//! formats agree with each other over the whole day.

use orbitcore_astro::omm::Omm;
use orbitcore_astro::tle::TwoLineElement;
use orbitcore_core::units::EARTH_MEAN_RADIUS_KM;
use orbitcore_sgp4::{Ephemeris, Propagator};

/// ISS mean motion from the snapshot, rev/day — the anchor for the
/// revolution and period checks below.
const ISS_MEAN_MOTION_REV_PER_DAY: f64 = 15.49369975;

/// Nodal period the mean motion implies: `1440 / n` ≈ 92.94 min.
const ISS_PERIOD_MINUTES: f64 = 1440.0 / ISS_MEAN_MOTION_REV_PER_DAY;

/// The snapshot TLE, read exactly as served by CelesTrak.
fn iss_tle_text() -> String {
    std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/data/iss-tle-2026-11-14.txt"
    ))
    .expect("ISS TLE snapshot exists")
}

/// The snapshot OMM, read exactly as served by CelesTrak (a
/// one-element JSON array).
fn iss_omm_text() -> String {
    std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/data/iss-omm-2026-11-14.json"
    ))
    .expect("ISS OMM snapshot exists")
}

/// The two TLE lines from the snapshot file. The file starts with a
/// name line ("ISS (ZARYA)"), so the lines are selected by their
/// leading "1 " / "2 " markers, not by position.
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

fn propagator_from_tle() -> Propagator {
    let (line1, line2) = tle_lines_from_snapshot();
    Propagator::from_tle(None, &line1, &line2).expect("ISS TLE initialises SGP4")
}

fn propagator_from_tle_elements() -> Propagator {
    let parsed = TwoLineElement::parse(&iss_tle_text()).expect("ISS TLE parses");
    Propagator::from_elements(&parsed.elements()).expect("ISS elements initialise SGP4")
}

fn propagator_from_omm() -> Propagator {
    let omms = Omm::parse_many(&iss_omm_text()).expect("ISS OMM snapshot parses");
    assert_eq!(omms.len(), 1, "snapshot holds one element set");
    let elements = omms[0].elements().expect("OMM converts to elements");
    Propagator::from_elements(&elements).expect("OMM elements initialise SGP4")
}

/// The 24-hour ephemeris at a 1-minute step: 1441 states, 0..=1440 min.
fn iss_24h() -> Ephemeris {
    Ephemeris::propagate_24h(&propagator_from_tle(), 1.0).expect("24 h propagation works")
}

#[test]
fn ephemeris_24h_has_1441_states_over_1440_minutes() {
    let eph = iss_24h();
    assert_eq!(eph.len(), 1441, "1-minute samples over 24 h plus the end");
    assert!(!eph.is_empty());
    assert!((eph.duration_minutes() - 1440.0).abs() < 1e-9);
    assert!((eph.states().first().unwrap().minutes_since_epoch - 0.0).abs() < 1e-9);
    assert!((eph.states().last().unwrap().minutes_since_epoch - 1440.0).abs() < 1e-9);
}

#[test]
fn revolutions_in_24h_match_the_published_mean_motion() {
    // The "vs CelesTrak" check: the element set publishes a mean motion
    // of 15.4937 rev/day, so the propagated arc must make that many
    // revolutions in 24 h. The ±0.05 band absorbs the J2-driven
    // difference between the mean (Keplerian) rate and the actual
    // osculating arc — the count is measured from the propagated
    // states, not copied from the element set.
    let eph = iss_24h();
    let revolutions = eph.revolutions();
    assert!(
        (revolutions - ISS_MEAN_MOTION_REV_PER_DAY).abs() < 0.05,
        "ISS made {revolutions:.4} revolutions in 24 h, expected {ISS_MEAN_MOTION_REV_PER_DAY} ± 0.05"
    );
}

#[test]
fn nodal_period_matches_1440_over_mean_motion() {
    // Mean interval between ascending-node crossings ≈ 92.94 min.
    // The ±0.5 min band absorbs the nodal-vs-anomalistic period
    // difference (J2 node regression) and the interpolation error of
    // the 1-minute sampling.
    let eph = iss_24h();
    let crossings = eph.ascending_node_crossings();
    // About one crossing per orbit over 24 h, plus or minus the
    // partial orbit at the start and end of the arc.
    assert!(
        (15..=17).contains(&crossings.len()),
        "expected 15–17 node crossings in 24 h, got {}",
        crossings.len()
    );
    let period = eph
        .nodal_period_minutes()
        .expect("24 h arc has many crossings");
    assert!(
        (period - ISS_PERIOD_MINUTES).abs() < 0.5,
        "nodal period = {period:.3} min, expected {ISS_PERIOD_MINUTES:.2} ± 0.5 min"
    );
}

#[test]
fn altitude_stays_within_perigee_apogee_band() {
    // Perigee/apogee from the elements themselves: a = semi-major
    // axis from mean motion (Kepler), perigee = a(1 − e) − R, apogee
    // = a(1 + e) − R. The ISS eccentricity is 0.0007385, so the band
    // is only ~10 km wide (about 421–431 km) — but the *osculating*
    // radius SGP4 returns wobbles a few km beyond the mean orbit:
    // J2 short-periodic terms make the radius oscillate around the
    // mean ellipse, so the measured extremes (415.4–434.3 km over
    // this day) sit ~5 km outside the mean perigee/apogee. The ±8 km
    // margin absorbs that oscillation while still catching a unit or
    // frame slip.
    let elements = TwoLineElement::parse(&iss_tle_text())
        .expect("ISS TLE parses")
        .elements();
    let a = elements.semi_major_axis_km();
    let e = elements.eccentricity;
    let perigee = a * (1.0 - e) - EARTH_MEAN_RADIUS_KM;
    let apogee = a * (1.0 + e) - EARTH_MEAN_RADIUS_KM;

    let eph = iss_24h();
    let min = eph.altitude_min_km().expect("arc has states");
    let max = eph.altitude_max_km().expect("arc has states");
    assert!(
        min >= perigee - 8.0,
        "lowest altitude {min:.2} km fell below perigee − 8 km ({:.2} km)",
        perigee - 8.0
    );
    assert!(
        max <= apogee + 8.0,
        "highest altitude {max:.2} km rose above apogee + 8 km ({:.2} km)",
        apogee + 8.0
    );
}

#[test]
fn from_tle_and_from_elements_agree_over_the_whole_day() {
    // Same check as `iss_propagation.rs`, but at every minute of the
    // day instead of five samples. The measured difference is exactly
    // zero, not merely small: the two paths differ only in the
    // mean-motion derivatives, and the upstream `sgp4` crate (v2.4.0)
    // parses those fields but never uses them in propagation — only
    // BSTAR feeds the model. The < 5 km bound is kept as the
    // regression guard, so a future upstream version that does use
    // the derivatives still has to stay within it.
    let from_tle = Ephemeris::propagate_24h(&propagator_from_tle(), 1.0).unwrap();
    let from_elements = Ephemeris::propagate_24h(&propagator_from_tle_elements(), 1.0).unwrap();
    assert_eq!(from_tle.len(), from_elements.len());
    let max_distance = from_tle
        .states()
        .iter()
        .zip(from_elements.states())
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0, f64::max);
    assert!(
        max_distance < 5.0,
        "from_tle vs from_elements differ by up to {max_distance:.3} km over 24 h"
    );
}

#[test]
fn tle_and_omm_snapshots_agree_over_the_whole_day() {
    // The two CelesTrak formats carry the same element set (same
    // epoch, same catalog number), so propagators built from each must
    // agree over the whole day. The OMM path goes through
    // `Omm::parse_many` on the snapshot file, exactly as a CelesTrak
    // GP JSON response arrives.
    let from_tle = Ephemeris::propagate_24h(&propagator_from_tle(), 1.0).unwrap();
    let from_omm = Ephemeris::propagate_24h(&propagator_from_omm(), 1.0).unwrap();
    assert_eq!(from_tle.len(), from_omm.len());
    let max_distance = from_tle
        .states()
        .iter()
        .zip(from_omm.states())
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0, f64::max);
    assert!(
        max_distance < 5.0,
        "TLE vs OMM propagators differ by up to {max_distance:.3} km over 24 h"
    );
}

#[test]
fn speeds_stay_orbital_and_mean_altitude_matches_the_elements() {
    let eph = iss_24h();
    // Circular-orbit speed at ~420 km is ≈ 7.66 km/s; the band is the
    // same one the hourly Day 1 test uses, now checked at every minute.
    for state in eph.states() {
        let speed = state.speed_km_s();
        assert!(
            (7.55..=7.75).contains(&speed),
            "speed {speed:.3} km/s at t = {} min left the orbital range",
            state.minutes_since_epoch
        );
    }
    // Mean altitude over the day vs the altitude the elements imply
    // (a − R): the osculating radius wobbles around the mean orbit,
    // so ±10 km is a generous band that still catches a unit or
    // frame slip.
    let elements = TwoLineElement::parse(&iss_tle_text())
        .expect("ISS TLE parses")
        .elements();
    let mean_altitude = eph.mean_altitude_km().expect("arc has states");
    assert!(
        (mean_altitude - elements.altitude_km()).abs() < 10.0,
        "mean altitude {mean_altitude:.2} km vs elements altitude {:.2} km",
        elements.altitude_km()
    );
    // Speed extremes: the near-circular orbit keeps the speed in a
    // narrow band (measured 7.646–7.669 km/s over this day), so the
    // spread itself is the check — a wider spread would mean the
    // radius was swinging.
    let speed_spread = eph.speed_max_km_s().unwrap() - eph.speed_min_km_s().unwrap();
    assert!(
        speed_spread < 0.05,
        "speed spread over 24 h = {speed_spread:.4} km/s, expected < 0.05 km/s"
    );
}
