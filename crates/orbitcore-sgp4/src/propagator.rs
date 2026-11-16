//! The SGP4 propagator wrapper: build once from an element set, then
//! propagate to any time near the epoch.
//!
//! The wrapper is deliberately thin — the SGP4 math is the upstream
//! `sgp4` crate's, and this module only adapts OrbitCore's
//! [`OrbitalElements`] into it and the upstream [`Prediction`] out of it,
//! so a bug here can only be a conversion bug, never a propagator bug.
//! That is the point of the wrapper: one place where the two type systems
//! meet, tested against the ISS (see `tests/iss_propagation.rs`).
//!
//! [`OrbitalElements`]: orbitcore_astro::elements::OrbitalElements
//! [`Prediction`]: sgp4::Prediction

use chrono::{Datelike, NaiveDate, NaiveDateTime, Timelike};
use orbitcore_astro::elements::OrbitalElements;
use orbitcore_astro::tle::EpochDateTime;
use orbitcore_core::error::{OrbitCoreError, Result};

use crate::state::StateVector;

/// An SGP4 propagator built from one element set.
///
/// Build it with [`Propagator::from_elements`] (from OrbitCore's parsed
/// elements) or [`Propagator::from_tle`] (from raw TLE text, preserving
/// every field), then call [`Propagator::propagate`] with minutes since
/// the epoch. The propagator is cheap to clone in the sense that matters
/// here — it holds the precomputed SGP4 constants, so propagating many
/// times from one `Propagator` is the fast path; building a new one per
/// timestamp would redo the initialisation each time.
#[derive(Debug, Clone)]
pub struct Propagator {
    /// The element set this propagator was built from, in OrbitCore types.
    elements: OrbitalElements,
    /// Precomputed SGP4 constants (the expensive initialisation).
    constants: sgp4::Constants,
}

impl Propagator {
    /// Build a propagator from OrbitCore's parsed element set.
    ///
    /// Fields an `OrbitalElements` does not carry — the mean-motion
    /// derivatives, element-set number, classification, object name —
    /// default to zero / unclassified / absent, because neither the TLE
    /// parser nor the OMM parser keeps them (see
    /// `orbitcore-astro::elements`). The drag term (BSTAR) *is* carried
    /// and is preserved. If the original TLE text is available and those
    /// fields matter, use [`Propagator::from_tle`] instead.
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError::InvalidValue`] when the element set
    /// cannot initialise SGP4 (e.g. an eccentricity outside SGP4's valid
    /// range).
    pub fn from_elements(elements: &OrbitalElements) -> Result<Self> {
        let sgp_elements = sgp4_elements_from_orbital(elements)?;
        let constants = sgp4::Constants::from_elements(&sgp_elements)
            .map_err(|e| OrbitCoreError::invalid_value("SGP4 elements", e.to_string()))?;
        Ok(Self {
            elements: elements.clone(),
            constants,
        })
    }

    /// Build a propagator by parsing a raw TLE with the upstream parser,
    /// so every TLE field — including the mean-motion derivatives the
    /// OrbitCore parser does not keep — survives into the propagator.
    ///
    /// `name` is the optional third line (object name); `line1` and
    /// `line2` are the two 69-character TLE lines.
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError::Parse`] when the TLE text is malformed,
    /// and [`OrbitCoreError::InvalidValue`] when the parsed elements
    /// cannot initialise SGP4.
    pub fn from_tle(name: Option<String>, line1: &str, line2: &str) -> Result<Self> {
        let sgp_elements = sgp4::Elements::from_tle(name, line1.as_bytes(), line2.as_bytes())
            .map_err(|e| OrbitCoreError::parse("TLE", e.to_string()))?;
        let constants = sgp4::Constants::from_elements(&sgp_elements)
            .map_err(|e| OrbitCoreError::invalid_value("SGP4 elements", e.to_string()))?;
        Ok(Self {
            elements: orbital_elements_from_sgp4(&sgp_elements),
            constants,
        })
    }

    /// The element set this propagator was built from.
    pub fn elements(&self) -> &OrbitalElements {
        &self.elements
    }

    /// NORAD catalog number of the propagated object.
    pub fn norad_cat_id(&self) -> u32 {
        self.elements.norad_cat_id
    }

    /// Propagate to `minutes_since_epoch` minutes after the epoch
    /// (negative values propagate backwards) and return the state vector
    /// in the TEME frame.
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError::InvalidValue`] when the time is not
    /// finite, or when SGP4 itself reports divergence (the orbit decayed,
    /// or the time is so far from the epoch that the model breaks down —
    /// in practice, stay within a few days of the epoch for LEO).
    pub fn propagate(&self, minutes_since_epoch: f64) -> Result<StateVector> {
        if !minutes_since_epoch.is_finite() {
            return Err(OrbitCoreError::invalid_value(
                "minutes since epoch",
                format!("must be finite, got {minutes_since_epoch}"),
            ));
        }
        let prediction = self
            .constants
            .propagate(sgp4::MinutesSinceEpoch(minutes_since_epoch))
            .map_err(|e| {
                OrbitCoreError::invalid_value(
                    "SGP4 propagation",
                    format!("{e} (t = {minutes_since_epoch} min since epoch)"),
                )
            })?;
        Ok(StateVector {
            position_km: prediction.position,
            velocity_km_s: prediction.velocity,
            minutes_since_epoch,
        })
    }

    /// Propagate over a time range at a fixed step: from `start_minutes`
    /// to `end_minutes` inclusive, every `step_minutes`. Returns one
    /// [`StateVector`] per sampled time, in time order.
    ///
    /// The step is a sampling choice, not a model parameter — SGP4 is
    /// analytic, so a 1-minute step over 24 hours gives 1441 states with
    /// no accumulated integration error. Day 2 uses this for the 24-hour
    /// propagation test against CelesTrak.
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError::InvalidValue`] when the step is not
    /// positive and finite, when `end_minutes < start_minutes`, or when
    /// any single propagation in the range fails.
    pub fn propagate_range(
        &self,
        start_minutes: f64,
        end_minutes: f64,
        step_minutes: f64,
    ) -> Result<Vec<StateVector>> {
        if !step_minutes.is_finite() || step_minutes <= 0.0 {
            return Err(OrbitCoreError::invalid_value(
                "propagation step",
                format!("must be finite and > 0 minutes, got {step_minutes}"),
            ));
        }
        if !start_minutes.is_finite() || !end_minutes.is_finite() {
            return Err(OrbitCoreError::invalid_value(
                "propagation range",
                format!("start and end must be finite, got {start_minutes}..{end_minutes}"),
            ));
        }
        if end_minutes < start_minutes {
            return Err(OrbitCoreError::invalid_value(
                "propagation range",
                format!("end ({end_minutes} min) must not be before start ({start_minutes} min)"),
            ));
        }
        let count = ((end_minutes - start_minutes) / step_minutes).floor() as usize + 1;
        let mut states = Vec::with_capacity(count);
        for i in 0..count {
            let t = start_minutes + i as f64 * step_minutes;
            states.push(self.propagate(t)?);
        }
        // Include the exact end time when the step does not land on it.
        if states
            .last()
            .is_some_and(|last| (last.minutes_since_epoch - end_minutes).abs() > 1e-9)
        {
            states.push(self.propagate(end_minutes)?);
        }
        Ok(states)
    }
}

/// Convert OrbitCore's [`OrbitalElements`] into the upstream
/// `sgp4::Elements`. Fields `OrbitalElements` does not carry default to
/// zero / unclassified — see [`Propagator::from_elements`].
fn sgp4_elements_from_orbital(elements: &OrbitalElements) -> Result<sgp4::Elements> {
    Ok(sgp4::Elements {
        object_name: None,
        international_designator: None,
        norad_id: u64::from(elements.norad_cat_id),
        classification: sgp4::Classification::Unclassified,
        datetime: naive_datetime_from_epoch(&elements.epoch)?,
        mean_motion_dot: 0.0,
        mean_motion_ddot: 0.0,
        drag_term: elements.bstar,
        element_set_number: 0,
        inclination: elements.inclination_deg,
        right_ascension: elements.raan_deg,
        eccentricity: elements.eccentricity,
        argument_of_perigee: elements.arg_perigee_deg,
        mean_anomaly: elements.mean_anomaly_deg,
        mean_motion: elements.mean_motion_rev_per_day,
        revolution_number: u64::from(elements.revolution_number.unwrap_or(0)),
        ephemeris_type: 0,
    })
}

/// Convert the upstream `sgp4::Elements` back into OrbitCore's
/// [`OrbitalElements`], so a propagator built from raw TLE text still
/// exposes the element set in OrbitCore types. The epoch day is
/// recomputed from the calendar datetime (day-of-year + fraction of day,
/// the TLE epoch convention).
///
/// A free function rather than a `From` impl: the orphan rule forbids
/// implementing a foreign trait for a foreign type in this crate.
fn orbital_elements_from_sgp4(elements: &sgp4::Elements) -> OrbitalElements {
    let dt = elements.datetime;
    let epoch_day = dt.ordinal() as f64
        + (f64::from(dt.num_seconds_from_midnight())
            + f64::from(dt.nanosecond()) / 1_000_000_000.0)
            / 86_400.0;
    OrbitalElements {
        norad_cat_id: u32::try_from(elements.norad_id).unwrap_or(u32::MAX),
        epoch: EpochDateTime {
            year: u16::try_from(dt.year()).unwrap_or(0),
            month: u8::try_from(dt.month()).unwrap_or(0),
            day: u8::try_from(dt.day()).unwrap_or(0),
            hour: u8::try_from(dt.hour()).unwrap_or(0),
            minute: u8::try_from(dt.minute()).unwrap_or(0),
            second: f64::from(dt.second()) + f64::from(dt.nanosecond()) / 1_000_000_000.0,
        },
        epoch_day,
        inclination_deg: elements.inclination,
        raan_deg: elements.right_ascension,
        eccentricity: elements.eccentricity,
        arg_perigee_deg: elements.argument_of_perigee,
        mean_anomaly_deg: elements.mean_anomaly,
        mean_motion_rev_per_day: elements.mean_motion,
        bstar: elements.drag_term,
        revolution_number: Some(u32::try_from(elements.revolution_number).unwrap_or(0)),
    }
}

/// Convert an OrbitCore [`EpochDateTime`] (UTC calendar date/time) into a
/// `chrono::NaiveDateTime`, which is what `sgp4::Elements` stores as its
/// epoch.
fn naive_datetime_from_epoch(epoch: &EpochDateTime) -> Result<NaiveDateTime> {
    let date = NaiveDate::from_ymd_opt(
        i32::from(epoch.year),
        u32::from(epoch.month),
        u32::from(epoch.day),
    )
    .ok_or_else(|| {
        OrbitCoreError::invalid_value(
            "epoch date",
            format!(
                "{:04}-{:02}-{:02} is not a valid calendar date",
                epoch.year, epoch.month, epoch.day
            ),
        )
    })?;
    let whole_seconds = epoch.second.floor();
    let nanos = ((epoch.second - whole_seconds) * 1_000_000_000.0) as u32;
    date.and_hms_nano_opt(
        u32::from(epoch.hour),
        u32::from(epoch.minute),
        whole_seconds as u32,
        nanos,
    )
    .ok_or_else(|| {
        OrbitCoreError::invalid_value(
            "epoch time",
            format!(
                "{:02}:{:02}:{:06.3} is not a valid time of day",
                epoch.hour, epoch.minute, epoch.second
            ),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The ISS element set from the Day 5 snapshot, as OrbitCore parses
    /// it (epoch 2026-11-14 08:49:22 UTC).
    fn iss_elements() -> OrbitalElements {
        let text = include_str!("../../../validation/data/iss-tle-2026-11-14.txt");
        let tle = orbitcore_astro::tle::TwoLineElement::parse(text).expect("snapshot TLE parses");
        tle.elements()
    }

    #[test]
    fn epoch_converts_to_naive_datetime() {
        let elements = iss_elements();
        let dt = naive_datetime_from_epoch(&elements.epoch).expect("ISS epoch is valid");
        assert_eq!((dt.year(), dt.month(), dt.day()), (2026, 11, 14));
        assert_eq!((dt.hour(), dt.minute()), (8, 49));
        // Epoch second 22.4 (day fraction .36762718 × 86400 = 31762.4 s
        // into the day = 08:49:22.4).
        assert!((f64::from(dt.second()) - 22.0).abs() < 1.0);
    }

    #[test]
    fn invalid_epoch_date_is_rejected() {
        let bad = EpochDateTime {
            year: 2026,
            month: 13,
            day: 1,
            hour: 0,
            minute: 0,
            second: 0.0,
        };
        assert!(naive_datetime_from_epoch(&bad).is_err());
    }

    #[test]
    fn sgp4_elements_round_trip_through_orbital_elements() {
        // OrbitalElements -> sgp4::Elements -> OrbitalElements must
        // preserve every field the OrbitCore type carries.
        let original = iss_elements();
        let sgp = sgp4_elements_from_orbital(&original).expect("conversion works");
        let back = orbital_elements_from_sgp4(&sgp);
        assert_eq!(back.norad_cat_id, original.norad_cat_id);
        // Calendar fields round-trip exactly; the fractional second only
        // to nanosecond precision (chrono stores whole nanoseconds).
        assert_eq!(
            (
                back.epoch.year,
                back.epoch.month,
                back.epoch.day,
                back.epoch.hour,
                back.epoch.minute
            ),
            (
                original.epoch.year,
                original.epoch.month,
                original.epoch.day,
                original.epoch.hour,
                original.epoch.minute
            )
        );
        assert!(
            (back.epoch.second - original.epoch.second).abs() < 1e-6,
            "epoch second round-trip: {} vs {}",
            back.epoch.second,
            original.epoch.second
        );
        assert!((back.epoch_day - original.epoch_day).abs() < 1e-6);
        assert_eq!(back.inclination_deg, original.inclination_deg);
        assert_eq!(back.raan_deg, original.raan_deg);
        assert_eq!(back.eccentricity, original.eccentricity);
        assert_eq!(back.arg_perigee_deg, original.arg_perigee_deg);
        assert_eq!(back.mean_anomaly_deg, original.mean_anomaly_deg);
        assert_eq!(
            back.mean_motion_rev_per_day,
            original.mean_motion_rev_per_day
        );
        assert_eq!(back.bstar, original.bstar);
        assert_eq!(back.revolution_number, original.revolution_number);
    }

    #[test]
    fn propagate_rejects_non_finite_time() {
        let prop = Propagator::from_elements(&iss_elements()).expect("ISS initialises SGP4");
        assert!(prop.propagate(f64::NAN).is_err());
        assert!(prop.propagate(f64::INFINITY).is_err());
    }

    #[test]
    fn propagate_range_rejects_bad_inputs() {
        let prop = Propagator::from_elements(&iss_elements()).expect("ISS initialises SGP4");
        assert!(prop.propagate_range(0.0, 60.0, 0.0).is_err());
        assert!(prop.propagate_range(0.0, 60.0, -5.0).is_err());
        assert!(prop.propagate_range(60.0, 0.0, 1.0).is_err());
        assert!(prop.propagate_range(f64::NAN, 60.0, 1.0).is_err());
    }

    #[test]
    fn propagate_range_includes_end_time() {
        let prop = Propagator::from_elements(&iss_elements()).expect("ISS initialises SGP4");
        // Step 30 over 0..=100 lands on 0, 30, 60, 90, then 100 exactly.
        let states = prop
            .propagate_range(0.0, 100.0, 30.0)
            .expect("range propagates");
        assert_eq!(states.len(), 5);
        assert!((states.last().unwrap().minutes_since_epoch - 100.0).abs() < 1e-9);
    }
}
