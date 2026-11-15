//! Day 5 — validation of the TLE/OMM parsers against a fresh, live ISS
//! element set.
//!
//! The element set below was fetched live from CelesTrak on 2026-11-14
//! (`gp.php?CATNR=25544&FORMAT=tle` and `...&FORMAT=json`) and stored as
//! snapshots under `validation/data/` — the tests below read those exact
//! files, so the parser is checked against what CelesTrak actually
//! served that morning, not against a hand-typed fixture. Tests stay
//! offline: the network fetch happens once, at authoring time, and the
//! snapshot is what ships.
//!
//! What this validates that Days 3–4 did not:
//!
//! * The parser works on a *new* element set it was never tuned to — a
//!   different epoch, different drag terms, different revolution number.
//! * The two formats still agree on the fresh element set (TLE vs OMM).
//! * The numbers track the real satellite across days: the epoch advanced
//!   ~1.94 days and the revolution number ~30 revolutions since the Day 3
//!   snapshot (2026-11-12), exactly what ~15.49 rev/day predicts.
//!
//! Fresh element set (epoch 2026-11-14 08:49:22 UTC, element set 999):
//! mean motion 15.49369975 rev/day, inclination 51.6323°, eccentricity
//! 0.0007385, revolution 59503.

use orbitcore_astro::omm::Omm;
use orbitcore_astro::tle::{TwoLineElement, checksum_valid};

/// The live ISS TLE fetched from CelesTrak on 2026-11-14, exactly as
/// served (3-line form, name line included).
const LIVE_TLE: &str = include_str!("../../../validation/data/iss-tle-2026-11-14.txt");

/// The live ISS OMM fetched from CelesTrak on 2026-11-14, exactly as
/// served (a one-element JSON array, the GP endpoint's response shape).
const LIVE_OMM_JSON: &str = include_str!("../../../validation/data/iss-omm-2026-11-14.json");

/// The Day 3 snapshot (epoch 2026-11-12 10:22:57 UTC) — the previous
/// live element set — used for the cross-day continuity check.
const DAY3_TLE: &str = "1 25544U 98067A   26316.43260916  .00008191  00000+0  15481-3 0  9999
2 25544  51.6324 286.4897 0007327  17.4271 342.6968 15.49332729594739";

fn approx(what: &str, got: f64, want: f64, tol: f64) {
    assert!(
        (got - want).abs() <= tol,
        "{what} = {got}, expected {want} ± {tol}"
    );
}

#[test]
fn live_tle_checksums_are_valid() {
    let mut lines = LIVE_TLE.lines().filter(|l| !l.trim().is_empty());
    let name = lines.next().expect("name line");
    assert_eq!(name.trim(), "ISS (ZARYA)");
    let line1 = lines.next().expect("line 1");
    let line2 = lines.next().expect("line 2");
    assert!(checksum_valid(line1), "live line 1 checksum must be valid");
    assert!(checksum_valid(line2), "live line 2 checksum must be valid");
    assert_eq!(line1.chars().nth(68), Some('9'));
    assert_eq!(line2.chars().nth(68), Some('5'));
}

#[test]
fn live_tle_fields_match_what_celestrak_served() {
    let tle = TwoLineElement::parse(LIVE_TLE).expect("live ISS TLE must parse");
    assert_eq!(tle.name.as_deref(), Some("ISS (ZARYA)"));
    assert_eq!(tle.satellite_number, 25544);
    assert_eq!(tle.classification, 'U');
    assert_eq!(tle.international_designator, "98067A");
    assert_eq!(tle.epoch_year, 2026);
    approx("epoch day", tle.epoch_day, 318.36762718, 1e-8);
    approx(
        "first derivative of mean motion",
        tle.first_derivative_mean_motion,
        0.00009745,
        1e-10,
    );
    approx(
        "second derivative of mean motion",
        tle.second_derivative_mean_motion,
        0.0,
        1e-12,
    );
    // "18238-3" = 0.18238e-3.
    approx("BSTAR", tle.bstar, 0.00018238, 1e-10);
    assert_eq!(tle.ephemeris_type, 0);
    assert_eq!(tle.element_set_number, 999);
    approx("inclination", tle.inclination_deg, 51.6323, 1e-9);
    approx("RAAN", tle.raan_deg, 276.9113, 1e-9);
    approx("eccentricity", tle.eccentricity, 0.0007385, 1e-9);
    approx("argument of perigee", tle.arg_perigee_deg, 24.4215, 1e-9);
    approx("mean anomaly", tle.mean_anomaly_deg, 335.7122, 1e-9);
    approx(
        "mean motion",
        tle.mean_motion_rev_per_day,
        15.49369975,
        1e-8,
    );
    assert_eq!(tle.revolution_number, 59503);
}

#[test]
fn live_tle_epoch_is_2026_11_14() {
    // Day 318 of 2026 is Nov 14; 0.36762718 d = 08:49:22 UTC.
    let tle = TwoLineElement::parse(LIVE_TLE).expect("live ISS TLE must parse");
    let epoch = tle.epoch_datetime();
    assert_eq!((epoch.year, epoch.month, epoch.day), (2026, 11, 14));
    assert_eq!((epoch.hour, epoch.minute), (8, 49));
    approx("epoch second", epoch.second, 22.988, 0.01);
}

#[test]
fn live_tle_derived_orbit_is_still_the_iss() {
    // ISS, Nov 2026: ~92.9 min period, ~426 km altitude — the same
    // orbit as the Day 3 snapshot, two days later.
    let tle = TwoLineElement::parse(LIVE_TLE).expect("live ISS TLE must parse");
    approx(
        "ISS orbital period (min)",
        tle.orbital_period_minutes(),
        92.93,
        0.05,
    );
    approx(
        "ISS semi-major axis (km)",
        tle.semi_major_axis_km(),
        6796.9,
        2.0,
    );
    let alt = tle.altitude_km();
    assert!(
        (415.0..435.0).contains(&alt),
        "ISS altitude = {alt:.1} km, expected roughly 426 km"
    );
}

#[test]
fn live_omm_fields_match_what_celestrak_served() {
    let omms = Omm::parse_many(LIVE_OMM_JSON).expect("live OMM array must parse");
    assert_eq!(omms.len(), 1);
    let omm = &omms[0];
    assert_eq!(omm.object_name.as_deref(), Some("ISS (ZARYA)"));
    assert_eq!(omm.object_id.as_deref(), Some("1998-067A"));
    assert_eq!(omm.epoch, "2026-11-14T08:49:22.988352");
    approx("mean motion", omm.mean_motion, 15.49369975, 1e-9);
    approx("eccentricity", omm.eccentricity.unwrap(), 0.0007385, 1e-9);
    approx("inclination", omm.inclination.unwrap(), 51.6323, 1e-9);
    approx("RAAN", omm.ra_of_asc_node.unwrap(), 276.9113, 1e-9);
    approx(
        "argument of pericenter",
        omm.arg_of_pericenter.unwrap(),
        24.4215,
        1e-9,
    );
    approx("mean anomaly", omm.mean_anomaly.unwrap(), 335.7122, 1e-9);
    assert_eq!(omm.norad_cat_id, Some(25544));
    assert_eq!(omm.element_set_no, Some(999));
    assert_eq!(omm.rev_at_epoch, Some(59503));
    approx("BSTAR", omm.bstar.unwrap(), 0.00018238, 1e-10);
    approx(
        "mean motion dot",
        omm.mean_motion_dot.unwrap(),
        9.745e-5,
        1e-10,
    );
    approx(
        "mean motion ddot",
        omm.mean_motion_ddot.unwrap(),
        0.0,
        1e-15,
    );
}

#[test]
fn live_tle_and_live_omm_agree() {
    // The fresh element set in both encodings must describe the same
    // orbit — the cross-format check from Day 4, repeated on new data.
    let tle = TwoLineElement::parse(LIVE_TLE).expect("live ISS TLE must parse");
    let omms = Omm::parse_many(LIVE_OMM_JSON).expect("live OMM array must parse");
    let from_tle = tle.elements();
    let from_omm = omms[0].elements().expect("OMM elements must convert");
    assert_eq!(from_omm.norad_cat_id, from_tle.norad_cat_id);
    // Same epoch — compared field by field with tolerance on the
    // seconds, because the TLE epoch (day-of-year fraction × 86400)
    // carries float noise the ISO timestamp does not.
    assert_eq!(from_omm.epoch.year, from_tle.epoch.year);
    assert_eq!(from_omm.epoch.month, from_tle.epoch.month);
    assert_eq!(from_omm.epoch.day, from_tle.epoch.day);
    assert_eq!(from_omm.epoch.hour, from_tle.epoch.hour);
    assert_eq!(from_omm.epoch.minute, from_tle.epoch.minute);
    approx(
        "epoch second",
        from_omm.epoch.second,
        from_tle.epoch.second,
        1e-6,
    );
    approx("epoch day", from_omm.epoch_day, from_tle.epoch_day, 1e-6);
    approx(
        "inclination",
        from_omm.inclination_deg,
        from_tle.inclination_deg,
        1e-9,
    );
    approx("RAAN", from_omm.raan_deg, from_tle.raan_deg, 1e-9);
    approx(
        "eccentricity",
        from_omm.eccentricity,
        from_tle.eccentricity,
        1e-9,
    );
    approx(
        "argument of perigee",
        from_omm.arg_perigee_deg,
        from_tle.arg_perigee_deg,
        1e-9,
    );
    approx(
        "mean anomaly",
        from_omm.mean_anomaly_deg,
        from_tle.mean_anomaly_deg,
        1e-9,
    );
    approx(
        "mean motion",
        from_omm.mean_motion_rev_per_day,
        from_tle.mean_motion_rev_per_day,
        1e-8,
    );
    approx("BSTAR", from_omm.bstar, from_tle.bstar, 1e-10);
    assert_eq!(from_omm.revolution_number, from_tle.revolution_number);
    approx(
        "orbital period (min)",
        from_omm.orbital_period_minutes(),
        from_tle.orbital_period_minutes(),
        1e-6,
    );
}

#[test]
fn iss_tracked_across_days() {
    // Continuity against the Day 3 snapshot (epoch 2026-11-12
    // 10:22:57 UTC, rev 59473): the ISS really did fly ~30 more
    // revolutions in the ~1.94 days between the two element sets, and
    // its orbit barely changed — mean motion 15.4933 → 15.4937 rev/day,
    // inclination 51.6324° → 51.6323°. A parser that mis-reads the epoch
    // or the revolution number would fail this.
    let old = TwoLineElement::parse(DAY3_TLE).expect("Day 3 TLE must parse");
    let live = TwoLineElement::parse(LIVE_TLE).expect("live ISS TLE must parse");

    let epoch_gap_days = live.epoch_day - old.epoch_day;
    approx("epoch gap (days)", epoch_gap_days, 1.93501802, 1e-8);
    assert!(
        epoch_gap_days > 0.0,
        "live epoch must be after the Day 3 epoch"
    );

    let revs = live.revolution_number as i64 - old.revolution_number as i64;
    assert_eq!(
        revs, 30,
        "revolution number advanced by {revs}, expected 30"
    );
    // Revolutions flown ≈ gap × mean motion: 1.93501802 d × 15.4937
    // rev/day ≈ 29.98 — within one revolution of the observed 30.
    let expected_revs = epoch_gap_days * live.mean_motion_rev_per_day;
    assert!(
        (revs as f64 - expected_revs).abs() < 1.0,
        "revs flown = {revs}, mean motion predicts {expected_revs:.2}"
    );

    // The orbit itself is nearly unchanged over two days.
    approx(
        "mean motion drift (rev/day)",
        live.mean_motion_rev_per_day - old.mean_motion_rev_per_day,
        0.00037246,
        1e-8,
    );
    approx(
        "inclination drift (deg)",
        live.inclination_deg - old.inclination_deg,
        -0.0001,
        1e-9,
    );
    // The period shortened by ~0.002 min as drag lowered the orbit —
    // the sign and size a real decaying LEO orbit should show.
    let period_change = live.orbital_period_minutes() - old.orbital_period_minutes();
    assert!(
        (-0.01..0.0).contains(&period_change),
        "period change = {period_change:.4} min, expected a small decrease"
    );
}
