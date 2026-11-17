//! A 24-hour (or arbitrary-span) ephemeris: a sampled arc of states
//! plus the analysis that turns it into orbit facts.
//!
//! [`Propagator::propagate`] answers "where is it at time *t*"; an
//! [`Ephemeris`] answers "what does the orbit look like over a day":
//! how many revolutions it makes, when it crosses the ascending node,
//! how the altitude and speed vary. The analysis is deliberately kept
//! separate from the propagation — SGP4 does the physics, this module
//! only measures the arc it returns — so each half can be tested on
//! its own: the math below is checked against synthetic states with
//! known answers (see the tests at the bottom), and the real-orbit
//! behaviour against the ISS snapshot in `tests/iss_24h.rs`.
//!
//! Everything here is offline and deterministic: an ephemeris is built
//! from a [`Propagator`], which needs no files and no network.
//!
//! [`Propagator::propagate`]: crate::propagator::Propagator::propagate

use std::f64::consts::PI;

use orbitcore_core::error::Result;

use crate::propagator::Propagator;
use crate::state::StateVector;

/// Minutes in a day — the span [`Ephemeris::propagate_24h`] covers.
pub const MINUTES_PER_DAY: f64 = 1440.0;

/// A sampled arc of propagated states, in time order, with the
/// analysis methods that turn the arc into orbit facts.
///
/// Build it with [`Ephemeris::propagate_24h`] (the Day 2 workhorse: a
/// full day from the epoch) or [`Ephemeris::propagate`] (any span), or
/// from states produced elsewhere with [`Ephemeris::from_states`].
/// All analysis methods are total over the stored states: an empty or
/// single-state ephemeris returns `0.0` / `None` instead of panicking,
/// because a one-sample arc simply contains no orbit information.
#[derive(Debug, Clone, PartialEq)]
pub struct Ephemeris {
    /// The sampled states, in ascending time order.
    states: Vec<StateVector>,
}

impl Ephemeris {
    /// Propagate a full day from the epoch: states at `0, step, 2·step,
    /// …, 1440` minutes, so a 1-minute step gives 1441 states.
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError`](orbitcore_core::error::OrbitCoreError)
    /// when the step is not positive and finite, or when any single
    /// propagation in the range fails (see
    /// [`Propagator::propagate_range`]).
    pub fn propagate_24h(propagator: &Propagator, step_minutes: f64) -> Result<Self> {
        Self::propagate(propagator, 0.0, MINUTES_PER_DAY, step_minutes)
    }

    /// Propagate over a time range at a fixed step: from `start_minutes`
    /// to `end_minutes` inclusive, every `step_minutes` (see
    /// [`Propagator::propagate_range`] for the sampling rules).
    ///
    /// # Errors
    ///
    /// Returns [`OrbitCoreError`](orbitcore_core::error::OrbitCoreError)
    /// when the step is not positive and finite, when `end_minutes <
    /// start_minutes`, or when any single propagation in the range
    /// fails.
    pub fn propagate(
        propagator: &Propagator,
        start_minutes: f64,
        end_minutes: f64,
        step_minutes: f64,
    ) -> Result<Self> {
        Ok(Self {
            states: propagator.propagate_range(start_minutes, end_minutes, step_minutes)?,
        })
    }

    /// Build an ephemeris from states produced elsewhere (another
    /// propagator, a file, a test fixture). The states are kept in the
    /// order given; callers are responsible for time order, since the
    /// analysis assumes ascending times.
    ///
    /// This is the constructor the unit tests use for synthetic arcs
    /// with known answers, and the one downstream code uses when it
    /// already holds states and only wants the analysis.
    pub fn from_states(states: Vec<StateVector>) -> Self {
        Self { states }
    }

    /// The sampled states, in time order.
    pub fn states(&self) -> &[StateVector] {
        &self.states
    }

    /// Number of sampled states.
    pub fn len(&self) -> usize {
        self.states.len()
    }

    /// Whether the ephemeris holds no states at all.
    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// Span covered by the arc in minutes: last time minus first time.
    /// `0.0` for an empty or single-state ephemeris, which covers no
    /// span.
    pub fn duration_minutes(&self) -> f64 {
        match (self.states.first(), self.states.last()) {
            (Some(first), Some(last)) => last.minutes_since_epoch - first.minutes_since_epoch,
            _ => 0.0,
        }
    }

    /// Complete revolutions made over the arc, as a fractional count.
    ///
    /// The angle measured is the *argument of latitude* — the angle of
    /// the position vector within the orbital plane, measured from the
    /// first state's position. The plane basis comes from the first
    /// state: `p = r₀/|r₀|` (in-plane, along the first position),
    /// `h = (r₀ × v₀)/|r₀ × v₀|` (the orbit normal), `q = h × p`
    /// (in-plane, 90° ahead of `p`). Each state's angle is
    /// `θ = atan2(r·q, r·p)`; the angles are unwrapped across samples
    /// (each step is assumed to advance by less than half a revolution
    /// — true for any sane sampling step, e.g. 1 minute ≈ 3.9° for the
    /// ISS) and the total unwrapped angle divided by `2π`.
    ///
    /// The result is the count of revolutions the satellite actually
    /// made, including the J2-driven drift of the real orbit — for the
    /// ISS over 24 h it lands within a few hundredths of the mean motion
    /// the element set publishes (see `tests/iss_24h.rs`). It is *not*
    /// simply `duration / period`: the period from the mean motion is
    /// the anomalistic period, while this counts plane crossings of the
    /// actual propagated arc.
    ///
    /// Returns `0.0` for an empty or single-state ephemeris, and for a
    /// degenerate first state (zero position or zero angular momentum),
    /// where no orbital plane exists to measure against.
    pub fn revolutions(&self) -> f64 {
        let Some(first) = self.states.first() else {
            return 0.0;
        };
        let Some((p, q)) = orbital_plane_basis(first) else {
            return 0.0;
        };
        let angle =
            |state: &StateVector| dot(&state.position_km, &q).atan2(dot(&state.position_km, &p));
        let mut total = 0.0;
        let mut prev = angle(first);
        for state in &self.states[1..] {
            let theta = angle(state);
            // Wrap the step into (−π, π]: the satellite cannot jump
            // more than half a revolution between samples.
            let delta = (theta - prev + PI).rem_euclid(2.0 * PI) - PI;
            total += delta;
            prev = theta;
        }
        total / (2.0 * PI)
    }

    /// Times (minutes since epoch) at which the satellite crosses the
    /// ascending node — the equatorial plane going from south to north,
    /// i.e. `z` going from negative to positive.
    ///
    /// Each crossing is located by linear interpolation between the
    /// last sample below the plane and the first at or above it, so the
    /// reported time is accurate to a fraction of the sampling step
    /// (seconds at a 1-minute step). A sample sitting exactly on the
    /// plane (`z = 0`) counts as the arrival side, so each crossing is
    /// reported exactly once.
    pub fn ascending_node_crossings(&self) -> Vec<f64> {
        let mut crossings = Vec::new();
        for pair in self.states.windows(2) {
            let (a, b) = (&pair[0], &pair[1]);
            let za = a.position_km[2];
            let zb = b.position_km[2];
            if za < 0.0 && zb >= 0.0 {
                // `zb − za > 0` here (za < 0 ≤ zb), so no division by zero.
                let frac = -za / (zb - za);
                crossings.push(
                    a.minutes_since_epoch + frac * (b.minutes_since_epoch - a.minutes_since_epoch),
                );
            }
        }
        crossings
    }

    /// Mean nodal period in minutes: the average interval between
    /// successive ascending-node crossings.
    ///
    /// The nodal period is the node-to-node period — the time between
    /// equator crossings — which for a J2-perturbed orbit differs
    /// slightly from the anomalistic period `1440 / n` the mean motion
    /// gives (for the ISS, by a fraction of a minute over a day). It is
    /// the period that matters for ground-track repetition and for
    /// eclipse timing.
    ///
    /// Returns `None` when fewer than two crossings were observed (an
    /// arc shorter than one orbit, or too few samples).
    pub fn nodal_period_minutes(&self) -> Option<f64> {
        let crossings = self.ascending_node_crossings();
        if crossings.len() < 2 {
            return None;
        }
        let total: f64 = crossings.windows(2).map(|w| w[1] - w[0]).sum();
        Some(total / (crossings.len() - 1) as f64)
    }

    /// Lowest altitude above the mean Earth radius over the arc, in km.
    /// `None` for an empty ephemeris.
    pub fn altitude_min_km(&self) -> Option<f64> {
        self.states
            .iter()
            .map(StateVector::altitude_km)
            .reduce(f64::min)
    }

    /// Highest altitude above the mean Earth radius over the arc, in
    /// km. `None` for an empty ephemeris.
    pub fn altitude_max_km(&self) -> Option<f64> {
        self.states
            .iter()
            .map(StateVector::altitude_km)
            .reduce(f64::max)
    }

    /// Mean altitude above the mean Earth radius over the arc, in km.
    /// `None` for an empty ephemeris.
    pub fn mean_altitude_km(&self) -> Option<f64> {
        if self.states.is_empty() {
            return None;
        }
        Some(
            self.states
                .iter()
                .map(StateVector::altitude_km)
                .sum::<f64>()
                / self.states.len() as f64,
        )
    }

    /// Lowest speed over the arc, in km/s. `None` for an empty
    /// ephemeris.
    pub fn speed_min_km_s(&self) -> Option<f64> {
        self.states
            .iter()
            .map(StateVector::speed_km_s)
            .reduce(f64::min)
    }

    /// Highest speed over the arc, in km/s. `None` for an empty
    /// ephemeris.
    pub fn speed_max_km_s(&self) -> Option<f64> {
        self.states
            .iter()
            .map(StateVector::speed_km_s)
            .reduce(f64::max)
    }
}

/// The in-plane basis (`p`, `q`) for the argument-of-latitude angle of
/// [`Ephemeris::revolutions`], from the first state: `p` along the
/// first position, `q` in-plane 90° ahead of `p` (via the orbit
/// normal `h = (r₀ × v₀)/|r₀ × v₀|`).
///
/// Returns `None` when the first state is degenerate — zero position,
/// or zero angular momentum (radial fall, no orbital plane) — so the
/// caller can report `0.0` revolutions instead of producing NaNs.
fn orbital_plane_basis(first: &StateVector) -> Option<([f64; 3], [f64; 3])> {
    let r = first.position_km;
    let v = first.velocity_km_s;
    let r_norm = norm(&r);
    if r_norm == 0.0 || !r_norm.is_finite() {
        return None;
    }
    let p = scale(&r, 1.0 / r_norm);
    let h = cross(&r, &v);
    let h_norm = norm(&h);
    if h_norm == 0.0 || !h_norm.is_finite() {
        return None;
    }
    let h = scale(&h, 1.0 / h_norm);
    let q = cross(&h, &p);
    Some((p, q))
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: &[f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn scale(a: &[f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A state on a circular equatorial orbit of radius `radius_km`,
    /// at angle `angle_rad` from the x-axis, moving counter-clockwise
    /// (`speed_km_s` tangential). The argument of latitude of this
    /// state is exactly `angle_rad`.
    fn circular_state(angle_rad: f64, minutes: f64) -> StateVector {
        let radius_km = 6796.0;
        let speed_km_s = 7.66;
        StateVector {
            position_km: [
                radius_km * angle_rad.cos(),
                radius_km * angle_rad.sin(),
                0.0,
            ],
            velocity_km_s: [
                -speed_km_s * angle_rad.sin(),
                speed_km_s * angle_rad.cos(),
                0.0,
            ],
            minutes_since_epoch: minutes,
        }
    }

    /// A state whose only meaningful component is its height above the
    /// equatorial plane, for node-crossing tests.
    fn state_at_z(z_km: f64, minutes: f64) -> StateVector {
        StateVector {
            position_km: [6796.0, 0.0, z_km],
            velocity_km_s: [0.0, 7.66, 0.0],
            minutes_since_epoch: minutes,
        }
    }

    #[test]
    fn revolutions_counts_two_full_circles() {
        // Eight samples per revolution, two revolutions: angles
        // 0, π/4, …, 4π. The last sample sits on top of the first
        // (angle 4π ≡ 0), so only the unwrapping recovers the count.
        let states: Vec<StateVector> = (0..=16)
            .map(|i| circular_state(i as f64 * PI / 4.0, i as f64))
            .collect();
        let eph = Ephemeris::from_states(states);
        assert!((eph.revolutions() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn revolutions_counts_a_partial_revolution() {
        // Angles 0, π/4, …, π: half a revolution.
        let states: Vec<StateVector> = (0..=4)
            .map(|i| circular_state(i as f64 * PI / 4.0, i as f64))
            .collect();
        let eph = Ephemeris::from_states(states);
        assert!((eph.revolutions() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn revolutions_counts_a_retrograde_arc_as_negative() {
        // Angles decreasing 0, −π/4, …, −π: half a revolution the
        // other way. The unwrapping is signed, so the count is −0.5.
        let states: Vec<StateVector> = (0..=4)
            .map(|i| circular_state(-(i as f64) * PI / 4.0, i as f64))
            .collect();
        let eph = Ephemeris::from_states(states);
        assert!((eph.revolutions() + 0.5).abs() < 1e-12);
    }

    #[test]
    fn revolutions_is_zero_without_an_arc() {
        assert_eq!(Ephemeris::from_states(vec![]).revolutions(), 0.0);
        assert_eq!(
            Ephemeris::from_states(vec![circular_state(0.0, 0.0)]).revolutions(),
            0.0
        );
    }

    #[test]
    fn revolutions_is_zero_for_a_degenerate_first_state() {
        // Zero velocity: no angular momentum, no orbital plane.
        let states = vec![
            StateVector {
                position_km: [6796.0, 0.0, 0.0],
                velocity_km_s: [0.0, 0.0, 0.0],
                minutes_since_epoch: 0.0,
            },
            circular_state(PI / 4.0, 1.0),
        ];
        let rev = Ephemeris::from_states(states).revolutions();
        assert_eq!(rev, 0.0);
        assert!(rev.is_finite());
    }

    #[test]
    fn ascending_node_crossing_is_interpolated() {
        // z: −2, −1, +1, +2 at t = 0..3 min. The crossing sits
        // between t = 1 (z = −1) and t = 2 (z = +1), halfway: t = 1.5.
        let states: Vec<StateVector> = [-2.0, -1.0, 1.0, 2.0]
            .iter()
            .enumerate()
            .map(|(i, &z)| state_at_z(z, i as f64))
            .collect();
        let crossings = Ephemeris::from_states(states).ascending_node_crossings();
        assert_eq!(crossings.len(), 1);
        assert!((crossings[0] - 1.5).abs() < 1e-12);
    }

    #[test]
    fn descending_crossings_are_not_counted() {
        // z goes +2, +1, −1, −2: a descending crossing only.
        let states: Vec<StateVector> = [2.0, 1.0, -1.0, -2.0]
            .iter()
            .enumerate()
            .map(|(i, &z)| state_at_z(z, i as f64))
            .collect();
        assert!(
            Ephemeris::from_states(states)
                .ascending_node_crossings()
                .is_empty()
        );
    }

    #[test]
    fn sample_exactly_on_the_plane_counts_once() {
        // z: −1, 0, +1 at t = 0, 1, 2. The t = 1 sample sits exactly
        // on the plane: the (0, 1) pair reports the crossing at t = 1,
        // and the (1, 2) pair must not report it again (its first
        // sample is not below the plane).
        let states: Vec<StateVector> = [-1.0, 0.0, 1.0]
            .iter()
            .enumerate()
            .map(|(i, &z)| state_at_z(z, i as f64))
            .collect();
        let crossings = Ephemeris::from_states(states).ascending_node_crossings();
        assert_eq!(crossings.len(), 1);
        assert!((crossings[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn nodal_period_is_the_mean_crossing_interval() {
        // Crossings at t = 0.5 and t = 2.5 (z −1→+1 twice, with a
        // descending leg between them), then a third at t = 4.5:
        // intervals 2.0 and 2.0, mean 2.0.
        let z = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0];
        let states: Vec<StateVector> = z
            .iter()
            .enumerate()
            .map(|(i, &z)| state_at_z(z, i as f64))
            .collect();
        let eph = Ephemeris::from_states(states);
        assert_eq!(eph.ascending_node_crossings().len(), 3);
        assert!((eph.nodal_period_minutes().unwrap() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn nodal_period_is_none_without_two_crossings() {
        assert!(
            Ephemeris::from_states(vec![])
                .nodal_period_minutes()
                .is_none()
        );
        // One crossing only.
        let states: Vec<StateVector> = [-1.0, 1.0]
            .iter()
            .enumerate()
            .map(|(i, &z)| state_at_z(z, i as f64))
            .collect();
        assert!(
            Ephemeris::from_states(states)
                .nodal_period_minutes()
                .is_none()
        );
    }

    #[test]
    fn altitude_and_speed_extremes() {
        // Radii 6796, 6800, 6790 km → altitudes 425, 429, 419 km;
        // speeds 7.66, 7.70, 7.60 km/s.
        let states = vec![
            StateVector {
                position_km: [6796.0, 0.0, 0.0],
                velocity_km_s: [0.0, 7.66, 0.0],
                minutes_since_epoch: 0.0,
            },
            StateVector {
                position_km: [6800.0, 0.0, 0.0],
                velocity_km_s: [0.0, 7.70, 0.0],
                minutes_since_epoch: 1.0,
            },
            StateVector {
                position_km: [6790.0, 0.0, 0.0],
                velocity_km_s: [0.0, 7.60, 0.0],
                minutes_since_epoch: 2.0,
            },
        ];
        let eph = Ephemeris::from_states(states);
        assert!((eph.altitude_min_km().unwrap() - 419.0).abs() < 1e-9);
        assert!((eph.altitude_max_km().unwrap() - 429.0).abs() < 1e-9);
        assert!((eph.mean_altitude_km().unwrap() - 424.333_333_333_333_3).abs() < 1e-9);
        assert!((eph.speed_min_km_s().unwrap() - 7.60).abs() < 1e-9);
        assert!((eph.speed_max_km_s().unwrap() - 7.70).abs() < 1e-9);
    }

    #[test]
    fn extremes_are_none_for_an_empty_ephemeris() {
        let eph = Ephemeris::from_states(vec![]);
        assert!(eph.altitude_min_km().is_none());
        assert!(eph.altitude_max_km().is_none());
        assert!(eph.mean_altitude_km().is_none());
        assert!(eph.speed_min_km_s().is_none());
        assert!(eph.speed_max_km_s().is_none());
        assert!(eph.is_empty());
        assert_eq!(eph.len(), 0);
        assert_eq!(eph.duration_minutes(), 0.0);
    }

    #[test]
    fn duration_is_last_minus_first_time() {
        let states = vec![
            circular_state(0.0, 10.0),
            circular_state(PI / 4.0, 100.0),
            circular_state(PI / 2.0, 1450.0),
        ];
        let eph = Ephemeris::from_states(states);
        assert_eq!(eph.len(), 3);
        assert!(!eph.is_empty());
        assert!((eph.duration_minutes() - 1440.0).abs() < 1e-12);
        assert_eq!(eph.states().len(), 3);
    }

    #[test]
    fn propagate_rejects_a_bad_step() {
        // The error path is the propagator's validation, surfaced
        // unchanged through the ephemeris constructor.
        let text = include_str!("../../../validation/data/iss-tle-2026-11-14.txt");
        let (line1, line2) = (
            text.lines().find(|l| l.starts_with("1 ")).unwrap(),
            text.lines().find(|l| l.starts_with("2 ")).unwrap(),
        );
        let prop = Propagator::from_tle(None, line1, line2).unwrap();
        assert!(Ephemeris::propagate_24h(&prop, 0.0).is_err());
        assert!(Ephemeris::propagate(&prop, 60.0, 0.0, 1.0).is_err());
    }
}
