//! SI units and physical constants for OrbitCore.
//!
//! Every physical quantity in `orbitcore-core` is a [`uom`] quantity, so a
//! meter/kilometer or watt/watt-hour mix-up fails at compile time instead of
//! silently corrupting a mission analysis — the classic unit bug in aerospace
//! software. Plain `f64` values are only used for dimensionless ratios
//! (efficiencies, duty cycles) and for constants whose unit is documented
//! in the name (e.g. `MU_EARTH_KM3_PER_S2`).
//!
//! Constants are grounded in the project plan (§2, §5.3):
//! * Earth's gravitational parameter μ = 398600.4418 km³/s²
//! * Earth's mean radius 6371 km
//! * Solar constant 1361 W/m²
//! * 500 km LEO: orbital period ≈ 94.6 min (94 min 37 s per the 1U EPS
//!   design paper), eclipse ≈ 35 min at β = 0°

use uom::si::{
    angle::degree,
    area::square_meter,
    energy::{joule, watt_hour},
    f64::{Angle, Area, Energy, Length, Power, Time},
    length::{kilometer, meter},
    power::watt,
    time::{minute, second},
};

use crate::error::{OrbitCoreError, Result};

/// Re-exported quantities, so downstream crates write `units::Length`
/// instead of the full `uom` path.
pub use uom::si::f64::{Mass, ThermodynamicTemperature, Velocity};

/// Earth's gravitational parameter, km³/s² (WGS-84 value used by SGP4).
pub const MU_EARTH_KM3_PER_S2: f64 = 398_600.441_8;

/// Earth's mean radius, in km.
pub const EARTH_MEAN_RADIUS_KM: f64 = 6_371.0;

/// Solar constant (solar irradiance at 1 AU), W/m².
pub const SOLAR_CONSTANT_W_PER_M2: f64 = 1_361.0;

/// Length in kilometers.
pub fn km(value: f64) -> Length {
    Length::new::<kilometer>(value)
}

/// Length in meters.
pub fn m(value: f64) -> Length {
    Length::new::<meter>(value)
}

/// Time in seconds.
pub fn seconds(value: f64) -> Time {
    Time::new::<second>(value)
}

/// Time in minutes.
pub fn minutes(value: f64) -> Time {
    Time::new::<minute>(value)
}

/// Power in watts.
pub fn watts(value: f64) -> Power {
    Power::new::<watt>(value)
}

/// Energy in watt-hours.
pub fn watt_hours(value: f64) -> Energy {
    Energy::new::<watt_hour>(value)
}

/// Energy in joules.
pub fn joules(value: f64) -> Energy {
    Energy::new::<joule>(value)
}

/// Angle in degrees.
pub fn degrees(value: f64) -> Angle {
    Angle::new::<degree>(value)
}

/// Area in square meters.
pub fn square_meters(value: f64) -> Area {
    Area::new::<square_meter>(value)
}

/// Semi-major axis of a circular orbit at the given altitude above the
/// Earth's surface: `a = R_earth + h`.
pub fn semi_major_axis(altitude: Length) -> Result<Length> {
    if altitude < Length::new::<kilometer>(0.0) {
        return Err(OrbitCoreError::out_of_range(
            "altitude",
            altitude.get::<kilometer>(),
            0.0,
            f64::MAX,
        ));
    }
    Ok(altitude + Length::new::<kilometer>(EARTH_MEAN_RADIUS_KM))
}

/// Orbital period of a circular orbit at the given altitude, from Kepler's
/// third law: `T = 2π√(a³/μ)`.
///
/// Validated against the literature: 500 km → ≈ 94.6 min (94 min 37 s,
/// 1U EPS design paper); 0 km (surface-skimming) → ≈ 84.4 min.
pub fn orbital_period(altitude: Length) -> Result<Time> {
    let a_km = semi_major_axis(altitude)?.get::<kilometer>();
    let period_s = 2.0 * std::f64::consts::PI * (a_km.powi(3) / MU_EARTH_KM3_PER_S2).sqrt();
    Ok(Time::new::<second>(period_s))
}

/// Power produced by a solar panel: `P = S × A × η × cos(θ)`, where `S` is
/// the solar constant (1361 W/m²), `A` the illuminated area, `η` the cell
/// efficiency, and `θ` the Sun incidence angle (`cos θ = 1` at normal
/// incidence, `0` edge-on).
///
/// `cosine` must be in [-1, 1] and `efficiency` in (0, 1]; anything else is
/// a caller bug and returns an error.
pub fn panel_power(area: Area, efficiency: f64, cosine: f64) -> Result<Power> {
    if !(0.0 < efficiency && efficiency <= 1.0) {
        return Err(OrbitCoreError::out_of_range(
            "panel efficiency",
            efficiency,
            0.0,
            1.0,
        ));
    }
    if !(-1.0..=1.0).contains(&cosine) {
        return Err(OrbitCoreError::out_of_range(
            "cosine of Sun incidence angle",
            cosine,
            -1.0,
            1.0,
        ));
    }
    let area_m2 = area.get::<square_meter>();
    Ok(Power::new::<watt>(
        SOLAR_CONSTANT_W_PER_M2 * area_m2 * efficiency * cosine.max(0.0),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::{energy::joule, length::meter, time::minute};

    #[test]
    fn unit_conversions_are_exact() {
        assert_eq!(km(1.0).get::<meter>(), 1000.0);
        assert_eq!(watt_hours(1.0).get::<joule>(), 3600.0);
        assert_eq!(minutes(1.0).get::<second>(), 60.0);
        assert_eq!(
            degrees(180.0).get::<uom::si::angle::radian>(),
            std::f64::consts::PI
        );
    }

    #[test]
    fn orbital_period_500_km_matches_literature() {
        // 1U EPS design paper: 94 min 37 s at 500 km.
        let period = orbital_period(km(500.0)).expect("500 km is a valid altitude");
        let minutes = period.get::<minute>();
        assert!(
            (minutes - 94.6).abs() < 0.5,
            "500 km period = {minutes:.2} min, expected ≈ 94.6 min"
        );
    }

    #[test]
    fn orbital_period_surface_skimming() {
        // Sanity anchor: a circular orbit at the surface takes ~84.4 min.
        let period = orbital_period(km(0.0)).expect("0 km is a valid altitude");
        let minutes = period.get::<minute>();
        assert!(
            (minutes - 84.4).abs() < 0.5,
            "0 km period = {minutes:.2} min, expected ≈ 84.4 min"
        );
    }

    #[test]
    fn orbital_period_grows_with_altitude() {
        let low = orbital_period(km(400.0)).unwrap();
        let high = orbital_period(km(800.0)).unwrap();
        assert!(high > low, "higher orbit must have a longer period");
    }

    #[test]
    fn orbital_period_rejects_negative_altitude() {
        let err = orbital_period(km(-100.0)).expect_err("negative altitude must fail");
        assert!(matches!(err, OrbitCoreError::OutOfRange { .. }));
    }

    #[test]
    fn panel_power_1u_reference() {
        // 1U face: 0.01 m² of 30% GaAs at normal incidence.
        // 1361 × 0.01 × 0.30 = 4.083 W.
        let p = panel_power(square_meters(0.01), 0.30, 1.0).expect("valid inputs");
        assert!(
            (p.get::<watt>() - 4.083).abs() < 1e-9,
            "panel power = {} W, expected 4.083 W",
            p.get::<watt>()
        );
    }

    #[test]
    fn panel_power_edge_on_is_zero() {
        let p = panel_power(square_meters(0.01), 0.30, 0.0).unwrap();
        assert_eq!(p.get::<watt>(), 0.0);
    }

    #[test]
    fn panel_power_rejects_bad_inputs() {
        assert!(panel_power(square_meters(0.01), 0.0, 1.0).is_err());
        assert!(panel_power(square_meters(0.01), 1.5, 1.0).is_err());
        assert!(panel_power(square_meters(0.01), 0.30, 2.0).is_err());
    }
}
