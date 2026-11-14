//! Orbit mean-elements message (OMM) parsing.
//!
//! OMM is the JSON format CelesTrak serves for the same element sets a
//! TLE carries (`gp.php?CATNR=25544&FORMAT=json` — note the explicit
//! `FORMAT=json`: since 9 May 2026 a GP query defaults to CSV). It is
//! the format that matters going forward, for two reasons:
//!
//! * A TLE's satellite number is five digits. NORAD catalog numbers ran
//!   past 100000 on 11 July 2026, so new objects simply cannot be
//!   represented as TLEs — OMM carries the full 1–9 digit
//!   `NORAD_CAT_ID`.
//! * OMM numbers are plain JSON numbers: eccentricity is `0.0007327`,
//!   not the TLE's implied-decimal `0007327`, and the epoch is an ISO
//!   8601 timestamp, not a two-digit year plus day-of-year.
//!
//! The parser is deliberately tolerant about *which* fields are present
//! — CelesTrak omits redundant fields (`CENTER_NAME`, `REF_FRAME`,
//! `TIME_SYSTEM`, `MEAN_ELEMENT_THEORY`) and any field it has no value
//! for, so every field except `EPOCH` and `MEAN_MOTION` is optional —
//! but strict about the values it does get: a missing `MEAN_MOTION` or
//! an unparseable epoch is an error, never a silent default, and
//! converting an OMM to [`OrbitalElements`] requires the full set of
//! Keplerian angles (an element set without them describes no orbit).
//!
//! Format reference: CCSDS 502.0-B-3, as documented by CelesTrak at
//! <https://celestrak.org/NORAD/documentation/gp-data-formats.php>.

use crate::elements::OrbitalElements;
use crate::tle::EpochDateTime;
use orbitcore_core::error::{OrbitCoreError, Result};
use serde::Deserialize;

/// One orbit mean-elements message: a single element set in JSON form.
///
/// Field names follow the CelesTrak JSON keys (`SCREAMING_SNAKE_CASE`).
/// Only `EPOCH` and `MEAN_MOTION` are required; everything else is
/// `Option` because CelesTrak omits fields it has no value for. Unknown
/// extra fields are ignored, so a future CelesTrak addition cannot break
/// parsing.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub struct Omm {
    /// Object name, e.g. `"ISS (ZARYA)"`.
    pub object_name: Option<String>,
    /// International designator, e.g. `"1998-067A"`.
    pub object_id: Option<String>,
    /// Epoch as an ISO 8601 timestamp, e.g. `"2026-11-12T10:22:57.431424"`.
    /// Always UTC; a trailing `Z` is accepted.
    pub epoch: String,
    /// Mean motion, revolutions per day.
    pub mean_motion: f64,
    /// Eccentricity (dimensionless, full decimal — unlike a TLE).
    pub eccentricity: Option<f64>,
    /// Inclination, degrees.
    pub inclination: Option<f64>,
    /// Right ascension of the ascending node, degrees.
    pub ra_of_asc_node: Option<f64>,
    /// Argument of pericenter, degrees.
    pub arg_of_pericenter: Option<f64>,
    /// Mean anomaly, degrees.
    pub mean_anomaly: Option<f64>,
    /// Ephemeris type (0 = SGP4).
    pub ephemeris_type: Option<i64>,
    /// Classification: `"U"` unclassified, `"C"` classified, `"S"` secret.
    pub classification_type: Option<String>,
    /// NORAD catalog number, 1–9 digits (this is the field a TLE cannot
    /// carry once the number exceeds 99999).
    pub norad_cat_id: Option<u32>,
    /// Element set number.
    pub element_set_no: Option<u32>,
    /// Revolution number at epoch.
    pub rev_at_epoch: Option<u32>,
    /// BSTAR drag term, 1/earth-radii.
    pub bstar: Option<f64>,
    /// First derivative of mean motion, rev/day² (may be scientific
    /// notation, e.g. `8.191e-5`).
    pub mean_motion_dot: Option<f64>,
    /// Second derivative of mean motion, rev/day³.
    pub mean_motion_ddot: Option<f64>,
}

impl Omm {
    /// Parse one OMM from a JSON object.
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError::Parse`] when the text is not valid JSON,
    /// is not a JSON object, or is missing `EPOCH` or `MEAN_MOTION`.
    pub fn parse(text: &str) -> Result<Self> {
        serde_json::from_str(text).map_err(|e| OrbitCoreError::parse("OMM", e.to_string()))
    }

    /// Parse the JSON array CelesTrak returns for a GP query
    /// (`[{...}, {...}]`), returning every element set in it.
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError::Parse`] when the text is not a JSON
    /// array of OMM objects, and [`OrbitCoreError::InvalidValue`] when
    /// the array is empty.
    pub fn parse_many(text: &str) -> Result<Vec<Self>> {
        let omms: Vec<Self> =
            serde_json::from_str(text).map_err(|e| OrbitCoreError::parse("OMM", e.to_string()))?;
        if omms.is_empty() {
            return Err(OrbitCoreError::invalid_value(
                "OMM",
                "JSON array contains no element sets",
            ));
        }
        Ok(omms)
    }

    /// Epoch as a calendar date/time (UTC).
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError::Parse`] when the `EPOCH` string is not
    /// an ISO 8601 timestamp.
    pub fn epoch_datetime(&self) -> Result<EpochDateTime> {
        parse_epoch_iso(&self.epoch).map(|(dt, _)| dt)
    }

    /// Epoch as a fractional day of year (1.0 = Jan 1, 00:00 UTC) — the
    /// same convention a TLE epoch uses, so OMM and TLE epochs compare
    /// directly.
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError::Parse`] when the `EPOCH` string is not
    /// an ISO 8601 timestamp.
    pub fn epoch_day(&self) -> Result<f64> {
        parse_epoch_iso(&self.epoch).map(|(_, day)| day)
    }

    /// Orbital period in minutes, from mean motion: `T = 1440 / n`.
    pub fn orbital_period_minutes(&self) -> f64 {
        1440.0 / self.mean_motion
    }

    /// The shared orbital elements this OMM describes, in the same form
    /// a TLE parses to — so a TLE and an OMM for the same satellite can
    /// be compared field by field.
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError::InvalidValue`] when any of the Keplerian
    /// angles (`INCLINATION`, `RA_OF_ASC_NODE`, `ECCENTRICITY`,
    /// `ARG_OF_PERICENTER`, `MEAN_ANOMALY`) is absent — an element set
    /// without them describes no orbit — or [`OrbitCoreError::Parse`]
    /// when the epoch is not a valid ISO 8601 timestamp.
    pub fn elements(&self) -> Result<OrbitalElements> {
        let (epoch, epoch_day) = parse_epoch_iso(&self.epoch)?;
        let angle = |field: &'static str, value: Option<f64>| {
            value.ok_or_else(|| {
                OrbitCoreError::invalid_value(
                    "OMM",
                    format!("cannot convert to orbital elements: {field} is missing"),
                )
            })
        };
        Ok(OrbitalElements {
            norad_cat_id: self.norad_cat_id.unwrap_or(0),
            epoch,
            epoch_day,
            inclination_deg: angle("INCLINATION", self.inclination)?,
            raan_deg: angle("RA_OF_ASC_NODE", self.ra_of_asc_node)?,
            eccentricity: angle("ECCENTRICITY", self.eccentricity)?,
            arg_perigee_deg: angle("ARG_OF_PERICENTER", self.arg_of_pericenter)?,
            mean_anomaly_deg: angle("MEAN_ANOMALY", self.mean_anomaly)?,
            mean_motion_rev_per_day: self.mean_motion,
            bstar: self.bstar.unwrap_or(0.0),
            revolution_number: self.rev_at_epoch,
        })
    }
}

/// Parse an ISO 8601 timestamp (`2026-11-12T10:22:57.431424`, optional
/// trailing `Z`, optional fractional seconds) into a calendar date/time
/// plus the fractional day of year.
fn parse_epoch_iso(text: &str) -> Result<(EpochDateTime, f64)> {
    let parse_err = |detail: String| OrbitCoreError::parse("OMM epoch", detail);
    let text = text.strip_suffix('Z').unwrap_or(text);
    let (date, time) = text
        .split_once(['T', ' '])
        .ok_or_else(|| parse_err(format!("no 'T' separator in epoch {text:?}")))?;

    let date_parts: Vec<&str> = date.split('-').collect();
    if date_parts.len() != 3 {
        return Err(parse_err(format!(
            "date part must be YYYY-MM-DD, got {date:?}"
        )));
    }
    let year: u16 = date_parts[0]
        .parse()
        .map_err(|_| parse_err(format!("bad epoch year in {date:?}")))?;
    let month: u8 = date_parts[1]
        .parse()
        .map_err(|_| parse_err(format!("bad epoch month in {date:?}")))?;
    let day: u8 = date_parts[2]
        .parse()
        .map_err(|_| parse_err(format!("bad epoch day in {date:?}")))?;
    if !(1..=12).contains(&month) {
        return Err(parse_err(format!("month {month} out of range in {date:?}")));
    }
    let days_in_month = days_in_month(year, month);
    if day == 0 || day > days_in_month {
        return Err(parse_err(format!(
            "day {day} out of range for {year}-{month:02} in {date:?}"
        )));
    }

    let time_parts: Vec<&str> = time.split(':').collect();
    if time_parts.len() != 3 {
        return Err(parse_err(format!(
            "time part must be HH:MM:SS[.fff], got {time:?}"
        )));
    }
    let hour: u8 = time_parts[0]
        .parse()
        .map_err(|_| parse_err(format!("bad epoch hour in {time:?}")))?;
    let minute: u8 = time_parts[1]
        .parse()
        .map_err(|_| parse_err(format!("bad epoch minute in {time:?}")))?;
    let second: f64 = time_parts[2]
        .parse()
        .map_err(|_| parse_err(format!("bad epoch second in {time:?}")))?;
    if hour > 23 {
        return Err(parse_err(format!("hour {hour} out of range in {time:?}")));
    }
    if minute > 59 {
        return Err(parse_err(format!(
            "minute {minute} out of range in {time:?}"
        )));
    }
    if !(0.0..61.0).contains(&second) {
        return Err(parse_err(format!(
            "second {second} out of range in {time:?}"
        )));
    }

    let day_of_year = day_of_year(year, month, day);
    let frac = (f64::from(hour) * 3600.0 + f64::from(minute) * 60.0 + second) / 86_400.0;
    Ok((
        EpochDateTime {
            year,
            month,
            day,
            hour,
            minute,
            second,
        },
        f64::from(day_of_year) + frac,
    ))
}

fn is_leap_year(year: u16) -> bool {
    year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400))
}

fn days_in_month(year: u16, month: u8) -> u8 {
    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 => {
            if is_leap_year(year) {
                29
            } else {
                28
            }
        }
        _ => 0,
    }
}

/// Day of year, 1-based (Jan 1 = 1).
fn day_of_year(year: u16, month: u8, day: u8) -> u16 {
    let mut doy: u16 = u16::from(day);
    for m in 1..month {
        doy += u16::from(days_in_month(year, m));
    }
    doy
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tle::TwoLineElement;

    /// The ISS OMM for the same element set Day 3 parsed as a TLE
    /// (epoch 2026-11-12 10:22:57 UTC), fetched live from CelesTrak on
    /// 2026-11-13 (`gp.php?CATNR=25544&FORMAT=json`). Tests run offline
    /// against this snapshot.
    pub const ISS_OMM_JSON: &str = r#"{"OBJECT_NAME":"ISS (ZARYA)","OBJECT_ID":"1998-067A","EPOCH":"2026-11-12T10:22:57.431424","MEAN_MOTION":15.49332729,"ECCENTRICITY":0.0007327,"INCLINATION":51.6324,"RA_OF_ASC_NODE":286.4897,"ARG_OF_PERICENTER":17.4271,"MEAN_ANOMALY":342.6968,"EPHEMERIS_TYPE":0,"CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":25544,"ELEMENT_SET_NO":999,"REV_AT_EPOCH":59473,"BSTAR":0.00015481,"MEAN_MOTION_DOT":8.191e-5,"MEAN_MOTION_DDOT":0}"#;

    /// The same element set as a TLE (Day 3's fixture) — the cross-format
    /// test below parses both and requires them to agree.
    const ISS_TLE: &str = "1 25544U 98067A   26316.43260916  .00008191  00000+0  15481-3 0  9999
2 25544  51.6324 286.4897 0007327  17.4271 342.6968 15.49332729594739";

    fn approx(what: &str, got: f64, want: f64, tol: f64) {
        assert!(
            (got - want).abs() <= tol,
            "{what} = {got}, expected {want} ± {tol}"
        );
    }

    #[test]
    fn parse_iss_omm_fields() {
        let omm = Omm::parse(ISS_OMM_JSON).expect("ISS OMM must parse");
        assert_eq!(omm.object_name.as_deref(), Some("ISS (ZARYA)"));
        assert_eq!(omm.object_id.as_deref(), Some("1998-067A"));
        assert_eq!(omm.epoch, "2026-11-12T10:22:57.431424");
        approx("mean motion", omm.mean_motion, 15.49332729, 1e-9);
        approx("eccentricity", omm.eccentricity.unwrap(), 0.0007327, 1e-9);
        approx("inclination", omm.inclination.unwrap(), 51.6324, 1e-9);
        approx("RAAN", omm.ra_of_asc_node.unwrap(), 286.4897, 1e-9);
        approx(
            "argument of pericenter",
            omm.arg_of_pericenter.unwrap(),
            17.4271,
            1e-9,
        );
        approx("mean anomaly", omm.mean_anomaly.unwrap(), 342.6968, 1e-9);
        assert_eq!(omm.ephemeris_type, Some(0));
        assert_eq!(omm.classification_type.as_deref(), Some("U"));
        assert_eq!(omm.norad_cat_id, Some(25544));
        assert_eq!(omm.element_set_no, Some(999));
        assert_eq!(omm.rev_at_epoch, Some(59473));
        approx("BSTAR", omm.bstar.unwrap(), 0.00015481, 1e-10);
        // Scientific notation for the mean motion derivatives.
        approx(
            "mean motion dot",
            omm.mean_motion_dot.unwrap(),
            8.191e-5,
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
    fn iss_epoch_converts_to_calendar() {
        // Same epoch as the Day 3 TLE: Nov 12, 10:22:57.431 UTC.
        let omm = Omm::parse(ISS_OMM_JSON).expect("ISS OMM must parse");
        let epoch = omm.epoch_datetime().expect("epoch must parse");
        assert_eq!((epoch.year, epoch.month, epoch.day), (2026, 11, 12));
        assert_eq!((epoch.hour, epoch.minute), (10, 22));
        approx("epoch second", epoch.second, 57.431424, 1e-6);
        // Day 316 of 2026 + 0.43260916 d — identical to the TLE epoch.
        approx("epoch day", omm.epoch_day().unwrap(), 316.43260916, 1e-7);
    }

    #[test]
    fn iss_derived_orbit_matches_tle() {
        let omm = Omm::parse(ISS_OMM_JSON).expect("ISS OMM must parse");
        approx(
            "ISS orbital period (min)",
            omm.orbital_period_minutes(),
            92.94,
            0.05,
        );
        let el = omm.elements().expect("elements must convert");
        approx(
            "ISS semi-major axis (km)",
            el.semi_major_axis_km(),
            6796.8,
            2.0,
        );
        let alt = el.altitude_km();
        assert!(
            (415.0..435.0).contains(&alt),
            "ISS altitude = {alt:.1} km, expected roughly 426 km"
        );
    }

    #[test]
    fn omm_and_tle_agree_on_same_element_set() {
        // The key cross-format check: the OMM and the TLE above are the
        // same element set (same epoch, same catalog number), so every
        // shared element must match.
        let omm = Omm::parse(ISS_OMM_JSON).expect("ISS OMM must parse");
        let tle = TwoLineElement::parse(ISS_TLE).expect("ISS TLE must parse");
        let from_omm = omm.elements().expect("OMM elements must convert");
        let from_tle = tle.elements();
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
    fn parse_celestrak_array_response() {
        // CelesTrak wraps the object in a one-element array.
        let array = format!("[{ISS_OMM_JSON}]");
        let omms = Omm::parse_many(&array).expect("array must parse");
        assert_eq!(omms.len(), 1);
        assert_eq!(omms[0].norad_cat_id, Some(25544));
        approx("mean motion", omms[0].mean_motion, 15.49332729, 1e-9);
    }

    #[test]
    fn parse_tolerates_missing_optional_fields() {
        // A bare OMM with only the two required fields still parses —
        // CelesTrak omits fields it has no value for.
        let minimal = r#"{"EPOCH":"2026-11-12T10:22:57.431424","MEAN_MOTION":15.49332729}"#;
        let omm = Omm::parse(minimal).expect("minimal OMM must parse");
        assert_eq!(omm.object_name, None);
        assert_eq!(omm.norad_cat_id, None);
        assert_eq!(omm.bstar, None);
        approx(
            "orbital period (min)",
            omm.orbital_period_minutes(),
            92.94,
            0.05,
        );
        // …but it describes no orbit without the Keplerian angles.
        let err = omm
            .elements()
            .expect_err("elements() without angles must fail");
        assert!(
            matches!(err, OrbitCoreError::InvalidValue { .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn parse_ignores_unknown_fields() {
        let extra = ISS_OMM_JSON.trim_end_matches('}');
        let with_extra = format!("{extra},\"FUTURE_FIELD\":42}}");
        let omm = Omm::parse(&with_extra).expect("unknown fields must be ignored");
        assert_eq!(omm.norad_cat_id, Some(25544));
    }

    #[test]
    fn parse_rejects_missing_mean_motion() {
        let no_mm = r#"{"EPOCH":"2026-11-12T10:22:57.431424","ECCENTRICITY":0.0007327}"#;
        let err = Omm::parse(no_mm).expect_err("missing MEAN_MOTION must fail");
        assert!(
            matches!(err, OrbitCoreError::Parse { .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn parse_rejects_invalid_epoch() {
        let bad_epoch = r#"{"EPOCH":"not-a-date","MEAN_MOTION":15.49332729}"#;
        let omm = Omm::parse(bad_epoch).expect("JSON itself is valid");
        let err = omm.epoch_datetime().expect_err("invalid epoch must fail");
        assert!(
            matches!(err, OrbitCoreError::Parse { .. }),
            "unexpected error: {err}"
        );
        assert!(omm.elements().is_err());
    }

    #[test]
    fn parse_rejects_non_object_json() {
        assert!(Omm::parse("[1, 2, 3]").is_err());
        assert!(Omm::parse("not json at all").is_err());
        assert!(Omm::parse("").is_err());
    }

    #[test]
    fn parse_many_rejects_empty_array() {
        let err = Omm::parse_many("[]").expect_err("empty array must fail");
        assert!(
            matches!(err, OrbitCoreError::InvalidValue { .. }),
            "unexpected error: {err}"
        );
        assert!(Omm::parse_many(ISS_OMM_JSON).is_err());
    }

    #[test]
    fn epoch_parsing_edge_cases() {
        // No fractional seconds, trailing Z.
        let (dt, day) = parse_epoch_iso("2026-01-01T00:00:00Z").unwrap();
        assert_eq!(
            dt,
            EpochDateTime {
                year: 2026,
                month: 1,
                day: 1,
                hour: 0,
                minute: 0,
                second: 0.0,
            }
        );
        assert!((day - 1.0).abs() < 1e-12, "day = {day}");
        // Dec 31 of a non-leap year, just before midnight.
        let (dt, day) = parse_epoch_iso("2025-12-31T23:59:59.999999").unwrap();
        assert_eq!((dt.year, dt.month, dt.day), (2025, 12, 31));
        assert_eq!((dt.hour, dt.minute), (23, 59));
        assert!((day - 366.0).abs() < 1e-5, "day = {day}");
        // Feb 29 exists in a leap year: day 60 of 2024.
        let (dt, day) = parse_epoch_iso("2024-02-29T12:00:00").unwrap();
        assert_eq!((dt.month, dt.day), (2, 29));
        assert!((day - 60.5).abs() < 1e-9, "day = {day}");
        // …and not in a non-leap year.
        assert!(parse_epoch_iso("2025-02-29T12:00:00").is_err());
        assert!(parse_epoch_iso("2026-13-01T00:00:00").is_err());
        assert!(parse_epoch_iso("2026-11-12T25:00:00").is_err());
        assert!(parse_epoch_iso("2026-11-12").is_err());
    }
}
