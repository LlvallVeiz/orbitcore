//! Two-line element set (TLE) parsing.
//!
//! A TLE is the standard NORAD format for sharing an orbit: two fixed-width
//! 69-character lines (optionally preceded by a name line, "line 0") that
//! fully describe a satellite's orbit at an epoch. OrbitCore reads TLEs
//! from CelesTrak and — from Week 2 — propagates them with SGP4.
//!
//! Field layout (1-based columns, per the NORAD / CelesTrak specification):
//!
//! ```text
//! Line 1: 1 NNNNNC U YYLLLPPP EEYYYDDD.DDDDDDDD +.DDDDDDDD +DDDDD-D +DDDDD-D 0 DDDD
//! Line 2: 2 NNNNNC III.IIII RRR.RRRR EEEEEEE AAA.AAAA BBB.BBBB MM.MMMMMMMMRRRRRR
//! ```
//!
//! Two quirks make the parser non-trivial and worth testing carefully:
//!
//! * **Implied decimal points.** The second derivative of mean motion,
//!   the BSTAR drag term, and the eccentricity are stored without a decimal
//!   point or (for the first two) with the exponent packed into the last two
//!   characters: `"-11606-4"` means `-0.11606e-4`, `"0006703"` means
//!   `0.0006703`.
//! * **Checksum.** Column 69 of each line is the modulo-10 sum of the
//!   digits in columns 1–68 (a `-` counts as 1, everything else as 0).
//!   A TLE whose checksum fails is corrupt data and is rejected, not
//!   silently trusted.
//!
//! The parser is strict about all of this — wrong line numbers, a wrong
//! length, mismatched satellite numbers, a bad checksum, or an
//! unparseable field all fail with [`OrbitCoreError::Parse`] or
//! [`OrbitCoreError::InvalidValue`] — because a silently mis-parsed orbit
//! is worse than no orbit at all.

use orbitcore_core::error::{OrbitCoreError, Result};
use orbitcore_core::units::{EARTH_MEAN_RADIUS_KM, MU_EARTH_KM3_PER_S2};

/// A parsed two-line element set.
#[derive(Debug, Clone, PartialEq)]
pub struct TwoLineElement {
    /// Optional name from line 0 (e.g. `"ISS (ZARYA)"`), if present.
    pub name: Option<String>,
    /// NORAD catalog number (satellite number).
    pub satellite_number: u32,
    /// Classification: `'U'` unclassified, `'C'` classified, `'S'` secret.
    pub classification: char,
    /// International designator, e.g. `"98067A"` (launch year + number + piece).
    pub international_designator: String,
    /// Epoch year (full four-digit year, expanded from the two-digit field).
    pub epoch_year: u16,
    /// Epoch day of year, with fractional part (1.0 = Jan 1, 00:00 UTC).
    pub epoch_day: f64,
    /// First derivative of mean motion, rev/day².
    pub first_derivative_mean_motion: f64,
    /// Second derivative of mean motion, rev/day³ (implied-decimal field).
    pub second_derivative_mean_motion: f64,
    /// BSTAR drag term, 1/earth-radii (implied-decimal field).
    pub bstar: f64,
    /// Ephemeris type (0 = SGP4).
    pub ephemeris_type: u8,
    /// Element set number.
    pub element_set_number: u32,
    /// Inclination, degrees.
    pub inclination_deg: f64,
    /// Right ascension of the ascending node, degrees.
    pub raan_deg: f64,
    /// Eccentricity (dimensionless, implied-decimal field).
    pub eccentricity: f64,
    /// Argument of perigee, degrees.
    pub arg_perigee_deg: f64,
    /// Mean anomaly, degrees.
    pub mean_anomaly_deg: f64,
    /// Mean motion, revolutions per day.
    pub mean_motion_rev_per_day: f64,
    /// Revolution number at epoch.
    pub revolution_number: u32,
}

/// A calendar date/time (UTC) derived from a TLE epoch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EpochDateTime {
    /// Full year, e.g. 2008.
    pub year: u16,
    /// Month, 1–12.
    pub month: u8,
    /// Day of month, 1–31.
    pub day: u8,
    /// Hour, 0–23.
    pub hour: u8,
    /// Minute, 0–59.
    pub minute: u8,
    /// Second, 0.0–60.0 (fractional).
    pub second: f64,
}

impl TwoLineElement {
    /// Parse a TLE from text: two lines, or three lines with a name line
    /// (line 0) first. Blank lines are ignored; trailing whitespace on a
    /// line is ignored. The checksum of both lines is verified.
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError::Parse`] or [`OrbitCoreError::InvalidValue`]
    /// when the text is not a well-formed TLE (wrong line count, wrong
    /// line numbers, wrong length, satellite-number mismatch, bad
    /// checksum, or an unparseable field).
    pub fn parse(text: &str) -> Result<Self> {
        let lines: Vec<&str> = text
            .lines()
            .map(str::trim_end)
            .filter(|line| !line.trim().is_empty())
            .collect();

        let (name, line1, line2) = match lines.as_slice() {
            [line1, line2] => (None, *line1, *line2),
            [name, line1, line2] => (Some((*name).trim().to_string()), *line1, *line2),
            _ => {
                return Err(OrbitCoreError::parse(
                    "TLE",
                    format!(
                        "expected 2 lines (or 3 with a name line), got {}",
                        lines.len()
                    ),
                ));
            }
        };

        let tle = Self::parse_lines(name, line1, line2)?;
        Ok(tle)
    }

    fn parse_lines(name: Option<String>, line1: &str, line2: &str) -> Result<Self> {
        for (n, line) in [(1, line1), (2, line2)] {
            if line.chars().count() < 69 {
                return Err(OrbitCoreError::invalid_value(
                    "TLE line",
                    format!("line {n} is {} chars, expected 69", line.chars().count()),
                ));
            }
        }
        if !line1.starts_with('1') {
            return Err(OrbitCoreError::invalid_value(
                "TLE line 1",
                "line 1 must start with '1'",
            ));
        }
        if !line2.starts_with('2') {
            return Err(OrbitCoreError::invalid_value(
                "TLE line 2",
                "line 2 must start with '2'",
            ));
        }
        for (n, line) in [(1, line1), (2, line2)] {
            if !checksum_valid(line) {
                return Err(OrbitCoreError::invalid_value(
                    "TLE line",
                    format!("line {n} checksum mismatch"),
                ));
            }
        }

        let sat1 = parse_u32(field(line1, 2, 7), "satellite number (line 1)")?;
        let sat2 = parse_u32(field(line2, 2, 7), "satellite number (line 2)")?;
        if sat1 != sat2 {
            return Err(OrbitCoreError::invalid_value(
                "TLE",
                format!("satellite number mismatch: line 1 has {sat1}, line 2 has {sat2}"),
            ));
        }

        let epoch_year_2digit = parse_u32(field(line1, 18, 20), "epoch year")?;
        let intl_year_2digit = parse_u32(field(line1, 9, 11), "international designator year")?;

        Ok(Self {
            name,
            satellite_number: sat1,
            classification: line1.chars().nth(7).unwrap_or('U'),
            international_designator: format!(
                "{:02}{:03}{}",
                intl_year_2digit,
                parse_u32(
                    field(line1, 11, 14),
                    "international designator launch number"
                )?,
                field(line1, 14, 17).trim(),
            ),
            epoch_year: full_year(epoch_year_2digit),
            epoch_day: parse_f64(field(line1, 20, 32), "epoch day")?,
            first_derivative_mean_motion: parse_f64(
                field(line1, 33, 43),
                "first derivative of mean motion",
            )?,
            second_derivative_mean_motion: parse_implied_decimal(
                field(line1, 44, 52),
                "second derivative of mean motion",
            )?,
            bstar: parse_implied_decimal(field(line1, 53, 61), "BSTAR drag term")?,
            ephemeris_type: parse_u32(field(line1, 62, 63), "ephemeris type")? as u8,
            element_set_number: parse_u32(field(line1, 64, 68), "element set number")?,
            inclination_deg: parse_f64(field(line2, 8, 16), "inclination")?,
            raan_deg: parse_f64(field(line2, 17, 25), "right ascension of ascending node")?,
            eccentricity: parse_implied_fraction(field(line2, 26, 33), "eccentricity")?,
            arg_perigee_deg: parse_f64(field(line2, 34, 42), "argument of perigee")?,
            mean_anomaly_deg: parse_f64(field(line2, 43, 51), "mean anomaly")?,
            mean_motion_rev_per_day: parse_f64(field(line2, 52, 63), "mean motion")?,
            revolution_number: parse_u32(field(line2, 63, 68), "revolution number")?,
        })
    }

    /// Epoch as a calendar date/time (UTC).
    pub fn epoch_datetime(&self) -> EpochDateTime {
        epoch_to_datetime(self.epoch_year, self.epoch_day)
    }

    /// Orbital period in minutes, from mean motion: `T = 1440 / n`.
    ///
    /// For the ISS example TLE this gives ≈ 91.6 min.
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

/// Compute the TLE checksum of a line: sum of the digits in the first 68
/// characters (a `-` counts as 1, anything else as 0), modulo 10.
pub fn checksum(line: &str) -> u8 {
    let mut sum: u32 = 0;
    for ch in line.chars().take(68) {
        match ch {
            '0'..='9' => sum += ch as u32 - '0' as u32,
            '-' => sum += 1,
            _ => {}
        }
    }
    (sum % 10) as u8
}

/// Check whether a TLE line's column-69 character equals its checksum.
pub fn checksum_valid(line: &str) -> bool {
    match line.chars().nth(68) {
        Some(digit) if digit.is_ascii_digit() => checksum(line) == digit as u8 - b'0',
        _ => false,
    }
}

/// Extract a fixed-width field by 0-based byte range `[start, end)`.
/// All TLE fields are ASCII, so byte indices are character indices.
fn field(line: &str, start: usize, end: usize) -> &str {
    line.get(start..end).unwrap_or("")
}

fn parse_f64(text: &str, what: &'static str) -> Result<f64> {
    text.trim().parse::<f64>().map_err(|_| {
        OrbitCoreError::parse("TLE", format!("field {what} is not a number: {text:?}"))
    })
}

fn parse_u32(text: &str, what: &'static str) -> Result<u32> {
    text.trim().parse::<u32>().map_err(|_| {
        OrbitCoreError::parse("TLE", format!("field {what} is not a number: {text:?}"))
    })
}

/// Parse an implied-decimal field such as `"-11606-4"` (BSTAR) or
/// `" 00000-0"` (second derivative): the leading digits are a mantissa
/// with an implied `0.` prefix, and the last two characters are a signed
/// single-digit exponent. `"-11606-4"` → `-0.11606e-4`.
fn parse_implied_decimal(text: &str, what: &'static str) -> Result<f64> {
    let s = text.trim();
    if s.len() < 3 {
        return Err(OrbitCoreError::parse(
            "TLE",
            format!("field {what} is too short for an implied-decimal value: {text:?}"),
        ));
    }
    let (mantissa, exp) = s.split_at(s.len() - 2);
    let exp: i32 = exp.parse().map_err(|_| {
        OrbitCoreError::parse("TLE", format!("field {what} has a bad exponent: {text:?}"))
    })?;
    let (sign, digits) = match mantissa.strip_prefix('-') {
        Some(rest) => (-1.0, rest),
        None => (1.0, mantissa.strip_prefix('+').unwrap_or(mantissa)),
    };
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Err(OrbitCoreError::parse(
            "TLE",
            format!("field {what} has a bad mantissa: {text:?}"),
        ));
    }
    let mantissa_value: f64 = digits.parse().map_err(|_| {
        OrbitCoreError::parse("TLE", format!("field {what} has a bad mantissa: {text:?}"))
    })?;
    let scale = 10f64.powi(-(digits.len() as i32));
    Ok(sign * mantissa_value * scale * 10f64.powi(exp))
}

/// Parse an implied-fraction field such as `"0006703"` (eccentricity):
/// the digits are the fractional part after an implied `0.`.
fn parse_implied_fraction(text: &str, what: &'static str) -> Result<f64> {
    let digits = text.trim();
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Err(OrbitCoreError::parse(
            "TLE",
            format!("field {what} is not a digit string: {text:?}"),
        ));
    }
    let value: f64 = digits.parse().map_err(|_| {
        OrbitCoreError::parse("TLE", format!("field {what} is not a number: {text:?}"))
    })?;
    Ok(value * 10f64.powi(-(digits.len() as i32)))
}

/// Expand a two-digit TLE year to a full year (NORAD convention:
/// 57–99 → 1957–1999, 00–56 → 2000–2056).
fn full_year(two_digit: u32) -> u16 {
    if two_digit < 57 {
        (2000 + two_digit) as u16
    } else {
        (1900 + two_digit) as u16
    }
}

/// Convert a TLE epoch (year + fractional day of year) to a calendar
/// date/time. Day 1.0 is Jan 1, 00:00 UTC.
fn epoch_to_datetime(year: u16, epoch_day: f64) -> EpochDateTime {
    let leap = year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400));
    let days_in_month = [
        31,
        if leap { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];

    let day_of_year = epoch_day.floor() as u32;
    let frac = epoch_day - epoch_day.floor();
    let total_seconds = frac * 86_400.0;
    let hour = (total_seconds / 3600.0).floor() as u8;
    let minute = ((total_seconds - f64::from(hour) * 3600.0) / 60.0).floor() as u8;
    let second = total_seconds - f64::from(hour) * 3600.0 - f64::from(minute) * 60.0;

    let mut month: u8 = 1;
    let mut day = day_of_year;
    for (i, days) in days_in_month.iter().enumerate() {
        if day <= *days {
            month = (i + 1) as u8;
            break;
        }
        day -= days;
    }

    EpochDateTime {
        year,
        month,
        day: day as u8,
        hour,
        minute,
        second,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The current ISS TLE, fetched live from CelesTrak on 2026-11-13
    /// (epoch 2026-11-12 10:22:57 UTC). Tests run offline against this
    /// snapshot; Day 5 re-validates the parser against a fresh live TLE.
    pub const ISS_TLE: &str =
        "1 25544U 98067A   26316.43260916  .00008191  00000+0  15481-3 0  9999
2 25544  51.6324 286.4897 0007327  17.4271 342.6968 15.49332729594739";

    /// The canonical ISS TLE from CelesTrak's documentation (epoch
    /// 2008-09-20) — a second, older element set proving the parser is
    /// not tuned to one TLE.
    pub const ISS_TLE_2008: &str =
        "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
2 25544  51.6416 247.4627 0006703  69.9862  25.2906 15.72193172510674";

    fn approx(what: &str, got: f64, want: f64, tol: f64) {
        assert!(
            (got - want).abs() <= tol,
            "{what} = {got}, expected {want} ± {tol}"
        );
    }

    #[test]
    fn checksums_of_iss_tle_are_valid() {
        let (line1, line2) = ISS_TLE.split_once('\n').unwrap();
        assert!(checksum_valid(line1), "line 1 checksum must be valid");
        assert!(checksum_valid(line2), "line 2 checksum must be valid");
        assert_eq!(checksum(line1), 9);
        assert_eq!(checksum(line2), 9);
        // The 2008 example is valid too — its line 2 checksum is 4, not 5.
        let (old1, old2) = ISS_TLE_2008.split_once('\n').unwrap();
        assert!(checksum_valid(old1), "2008 line 1 checksum must be valid");
        assert!(checksum_valid(old2), "2008 line 2 checksum must be valid");
        assert_eq!(checksum(old1), 7);
        assert_eq!(checksum(old2), 4);
    }

    #[test]
    fn parse_iss_line1_fields() {
        let tle = TwoLineElement::parse(ISS_TLE).expect("ISS TLE must parse");
        assert_eq!(tle.satellite_number, 25544);
        assert_eq!(tle.classification, 'U');
        assert_eq!(tle.international_designator, "98067A");
        assert_eq!(tle.epoch_year, 2026);
        approx("epoch day", tle.epoch_day, 316.43260916, 1e-8);
        approx(
            "first derivative mean motion",
            tle.first_derivative_mean_motion,
            0.00008191,
            1e-10,
        );
        approx(
            "second derivative mean motion",
            tle.second_derivative_mean_motion,
            0.0,
            1e-12,
        );
        approx("BSTAR", tle.bstar, 0.15481e-3, 1e-10);
        assert_eq!(tle.ephemeris_type, 0);
        assert_eq!(tle.element_set_number, 999);
        assert_eq!(tle.name, None);
    }

    #[test]
    fn parse_iss_line2_fields() {
        let tle = TwoLineElement::parse(ISS_TLE).expect("ISS TLE must parse");
        approx("inclination", tle.inclination_deg, 51.6324, 1e-9);
        approx("RAAN", tle.raan_deg, 286.4897, 1e-9);
        approx("eccentricity", tle.eccentricity, 0.0007327, 1e-9);
        approx("argument of perigee", tle.arg_perigee_deg, 17.4271, 1e-9);
        approx("mean anomaly", tle.mean_anomaly_deg, 342.6968, 1e-9);
        approx(
            "mean motion",
            tle.mean_motion_rev_per_day,
            15.49332729,
            1e-8,
        );
        assert_eq!(tle.revolution_number, 59473);
    }

    #[test]
    fn iss_epoch_converts_to_calendar() {
        // Day 316 of 2026 is Nov 12; 0.43260916 d = 10:22:57 UTC.
        let tle = TwoLineElement::parse(ISS_TLE).expect("ISS TLE must parse");
        let epoch = tle.epoch_datetime();
        assert_eq!((epoch.year, epoch.month, epoch.day), (2026, 11, 12));
        assert_eq!((epoch.hour, epoch.minute), (10, 22));
        approx("epoch second", epoch.second, 57.431, 0.01);
    }

    #[test]
    fn iss_derived_orbit_matches_reality() {
        // ISS, Nov 2026: ~92.9 min period, ~426 km altitude.
        let tle = TwoLineElement::parse(ISS_TLE).expect("ISS TLE must parse");
        approx(
            "ISS orbital period (min)",
            tle.orbital_period_minutes(),
            92.94,
            0.05,
        );
        approx(
            "ISS semi-major axis (km)",
            tle.semi_major_axis_km(),
            6796.8,
            2.0,
        );
        let alt = tle.altitude_km();
        assert!(
            (415.0..435.0).contains(&alt),
            "ISS altitude = {alt:.1} km, expected roughly 426 km"
        );
    }

    #[test]
    fn parse_2008_iss_tle_too() {
        // The older canonical example parses as well: epoch 2008-09-20,
        // mean motion 15.72193172 rev/day (~91.6 min period).
        let tle = TwoLineElement::parse(ISS_TLE_2008).expect("2008 ISS TLE must parse");
        assert_eq!(tle.epoch_year, 2008);
        assert_eq!(tle.international_designator, "98067A");
        approx("BSTAR 2008", tle.bstar, -0.11606e-4, 1e-10);
        approx(
            "mean motion 2008",
            tle.mean_motion_rev_per_day,
            15.72193172,
            1e-8,
        );
        assert_eq!(tle.revolution_number, 51067);
        let epoch = tle.epoch_datetime();
        assert_eq!((epoch.year, epoch.month, epoch.day), (2008, 9, 20));
        assert_eq!((epoch.hour, epoch.minute), (12, 25));
        approx("epoch second 2008", epoch.second, 40.104, 0.01);
        approx(
            "ISS period 2008 (min)",
            tle.orbital_period_minutes(),
            91.60,
            0.05,
        );
    }

    #[test]
    fn parse_with_name_line() {
        let with_name = format!("ISS (ZARYA)\n{ISS_TLE}");
        let tle = TwoLineElement::parse(&with_name).expect("3-line TLE must parse");
        assert_eq!(tle.name.as_deref(), Some("ISS (ZARYA)"));
        assert_eq!(tle.satellite_number, 25544);
    }

    #[test]
    fn parse_rejects_bad_checksum() {
        // Flip the last character of line 1 (the checksum) from 9 to 8.
        let (line1, line2) = ISS_TLE.split_once('\n').unwrap();
        let bad = format!("{}8\n{}", &line1[..line1.len() - 1], line2);
        let err = TwoLineElement::parse(&bad).expect_err("bad checksum must fail");
        assert!(
            matches!(err, OrbitCoreError::InvalidValue { .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn parse_rejects_wrong_line_numbers() {
        let (line1, line2) = ISS_TLE.split_once('\n').unwrap();
        let swapped = format!("2{}\n{}", &line1[1..], line2);
        assert!(TwoLineElement::parse(&swapped).is_err());
    }

    #[test]
    fn parse_rejects_short_lines() {
        assert!(TwoLineElement::parse("1 25544U\n2 25544").is_err());
    }

    #[test]
    fn parse_rejects_satellite_number_mismatch() {
        let (line1, line2) = ISS_TLE.split_once('\n').unwrap();
        // Change the satellite number on line 2, then repair its checksum
        // so the parser reaches the mismatch check.
        let mut line2 = line2.replacen("2 25544 ", "2 25545 ", 1);
        line2.truncate(68);
        line2.push(char::from_digit(u32::from(checksum(&line2)), 10).unwrap());
        let text = format!("{line1}\n{line2}");
        let err = TwoLineElement::parse(&text).expect_err("mismatch must fail");
        assert!(
            matches!(err, OrbitCoreError::InvalidValue { .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn parse_rejects_wrong_line_count() {
        assert!(TwoLineElement::parse(ISS_TLE.lines().next().unwrap()).is_err());
        assert!(TwoLineElement::parse("").is_err());
    }

    #[test]
    fn two_digit_year_mapping() {
        assert_eq!(full_year(98), 1998);
        assert_eq!(full_year(8), 2008);
        assert_eq!(full_year(56), 2056);
        assert_eq!(full_year(57), 1957);
        assert_eq!(full_year(99), 1999);
    }

    #[test]
    fn implied_decimal_parsing() {
        approx(
            "BSTAR-style negative",
            parse_implied_decimal("-11606-4", "test").unwrap(),
            -0.11606e-4,
            1e-12,
        );
        approx(
            "zero second derivative",
            parse_implied_decimal(" 00000-0", "test").unwrap(),
            0.0,
            1e-15,
        );
        approx(
            "positive exponent",
            parse_implied_decimal(" 12345+2", "test").unwrap(),
            12.345,
            1e-9,
        );
        assert!(parse_implied_decimal("abc", "test").is_err());
    }

    #[test]
    fn implied_fraction_parsing() {
        approx(
            "ISS eccentricity",
            parse_implied_fraction("0006703", "test").unwrap(),
            0.0006703,
            1e-12,
        );
        approx(
            "circular orbit",
            parse_implied_fraction("0000000", "test").unwrap(),
            0.0,
            1e-15,
        );
        assert!(parse_implied_fraction("12A4567", "test").is_err());
    }

    #[test]
    fn epoch_day_conversion_edge_cases() {
        // Jan 1, 00:00 exactly.
        let e = epoch_to_datetime(2026, 1.0);
        assert_eq!(
            e,
            EpochDateTime {
                year: 2026,
                month: 1,
                day: 1,
                hour: 0,
                minute: 0,
                second: 0.0,
            }
        );
        // Dec 31 of a non-leap year (day 365), 23:59:59.
        let e = epoch_to_datetime(2025, 365.99998843);
        assert_eq!(e.year, 2025);
        assert_eq!((e.month, e.day), (12, 31));
        assert_eq!((e.hour, e.minute), (23, 59));
        assert!((e.second - 59.0).abs() < 0.01, "second = {}", e.second);
        // Feb 29 exists in a leap year: day 60 of 2024.
        let e = epoch_to_datetime(2024, 60.5);
        assert_eq!((e.month, e.day), (2, 29));
        assert_eq!(e.hour, 12);
        // Day 60 of a non-leap year is Mar 1.
        let e = epoch_to_datetime(2025, 60.5);
        assert_eq!((e.month, e.day), (3, 1));
    }
}
