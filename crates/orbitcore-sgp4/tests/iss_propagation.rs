//! Integration tests for the SGP4 wrapper: propagate the real ISS element
//! set from the Day 5 CelesTrak snapshot (`validation/data/`) and check
//! the results against the physics the element set itself implies.
//!
//! The snapshot is the live ISS TLE fetched 2026-11-14 (epoch
//! 2026-11-14 08:49:22 UTC, mean motion 15.49369975 rev/day, inclination
//! 51.6323°, eccentricity 0.0007385). Everything asserted here follows
//! from those numbers plus basic orbital mechanics — no network, no
//! CelesTrak call at test time.

use orbitcore_astro::tle::TwoLineElement;
use orbitcore_sgp4::Propagator;

/// The snapshot TLE, read exactly as served by CelesTrak.
fn iss_tle_text() -> String {
    std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/data/iss-tle-2026-11-14.txt"
    ))
    .expect("ISS TLE snapshot exists")
}

/// ISS mean motion from the snapshot, rev/day — the anchor for the
/// period-based checks below.
const ISS_MEAN_MOTION_REV_PER_DAY: f64 = 15.49369975;

#[test]
fn iss_at_epoch_is_in_leo_at_orbital_speed() {
    let prop = propagator_from_snapshot();

    let state = prop.propagate(0.0).expect("propagation at epoch works");
    // ISS flies at ~420 km; allow a generous band — the point is that
    // the wrapper returns a real LEO state, not a unit or frame slip.
    let altitude = state.altitude_km();
    assert!(
        (380.0..=460.0).contains(&altitude),
        "ISS altitude at epoch = {altitude:.1} km, expected ~420 km"
    );
    // Circular-orbit speed at ~420 km: sqrt(mu / r) ≈ 7.66 km/s.
    let speed = state.speed_km_s();
    assert!(
        (7.60..=7.72).contains(&speed),
        "ISS speed at epoch = {speed:.3} km/s, expected ~7.66 km/s"
    );
    assert_eq!(state.minutes_since_epoch, 0.0);
}

#[test]
fn iss_returns_near_start_after_one_period() {
    let prop = propagator_from_snapshot();
    let period_min = 1440.0 / ISS_MEAN_MOTION_REV_PER_DAY; // ≈ 92.94 min
    let start = prop.propagate(0.0).unwrap();
    let after_one_orbit = prop.propagate(period_min).unwrap();

    // After one anomalistic period the satellite is back near the same
    // inertial point. J2 node regression (~5°/day) moves it ~0.3° in
    // 93 minutes, i.e. tens of km — 150 km is a generous bound that
    // still catches a broken propagator or a unit slip.
    let distance = start.distance_to(&after_one_orbit);
    assert!(
        distance < 150.0,
        "ISS moved {distance:.1} km after one period ({period_min:.2} min), expected < 150 km"
    );
    // The orbit itself is near-circular: the radius barely changes.
    assert!(
        (start.radius_km() - after_one_orbit.radius_km()).abs() < 20.0,
        "radius changed by {:.1} km over one orbit",
        (start.radius_km() - after_one_orbit.radius_km()).abs()
    );
}

#[test]
fn iss_24h_propagation_stays_in_leo() {
    let prop = propagator_from_snapshot();
    let states = prop
        .propagate_range(0.0, 1440.0, 60.0)
        .expect("24 h propagation works");
    assert_eq!(states.len(), 25, "hourly samples over 24 h plus the end");

    for state in &states {
        let radius = state.radius_km();
        assert!(
            (6700.0..=6900.0).contains(&radius),
            "radius {radius:.1} km at t = {} min left the LEO band",
            state.minutes_since_epoch
        );
        assert!(
            (7.55..=7.75).contains(&state.speed_km_s()),
            "speed {:.3} km/s at t = {} min left the orbital range",
            state.speed_km_s(),
            state.minutes_since_epoch
        );
    }
    // No secular drift over a day: the last radius is within a few km of
    // the first (drag over 24 h at ISS altitude is sub-kilometre).
    let drift = (states.first().unwrap().radius_km() - states.last().unwrap().radius_km()).abs();
    assert!(drift < 10.0, "radius drifted {drift:.2} km over 24 h");
}

#[test]
fn from_elements_agrees_with_from_tle_over_a_day() {
    // The OrbitCore-parser path drops the mean-motion derivatives
    // (OrbitCore's OrbitalElements does not carry them), while from_tle
    // keeps them. Both must still agree closely: the ISS mean-motion dot
    // is ~1e-4 rev/day^2, worth well under a kilometre over 12 h.
    let text = iss_tle_text();
    let (line1, line2) = tle_lines_from_snapshot();

    let from_tle = Propagator::from_tle(None, &line1, &line2).unwrap();
    let parsed = TwoLineElement::parse(&text).unwrap();
    let from_elements = Propagator::from_elements(&parsed.elements()).unwrap();

    for minutes in [0.0, 60.0, 360.0, 720.0, 1440.0] {
        let a = from_tle.propagate(minutes).unwrap();
        let b = from_elements.propagate(minutes).unwrap();
        let distance = a.distance_to(&b);
        assert!(
            distance < 5.0,
            "from_tle vs from_elements differ by {distance:.3} km at t = {minutes} min"
        );
    }
}

#[test]
fn from_tle_and_from_elements_agree_on_epoch_and_norad() {
    let text = iss_tle_text();
    let parsed = TwoLineElement::parse(&text).unwrap();
    let prop = Propagator::from_elements(&parsed.elements()).unwrap();
    assert_eq!(prop.norad_cat_id(), 25544);
    assert_eq!(prop.elements().epoch.year, 2026);
    assert_eq!(prop.elements().epoch.month, 11);
    assert_eq!(prop.elements().epoch.day, 14);
    assert!((prop.elements().epoch_day - 318.36762718).abs() < 1e-6);
}

fn propagator_from_snapshot() -> Propagator {
    let (line1, line2) = tle_lines_from_snapshot();
    Propagator::from_tle(None, &line1, &line2).expect("ISS TLE initialises SGP4")
}

/// The two TLE lines from the snapshot file. The file starts with a
/// name line ("ISS (ZARYA)"), so the lines are selected by their leading
/// "1 " / "2 " markers, not by position.
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
