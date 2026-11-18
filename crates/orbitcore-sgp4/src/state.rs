//! A propagated state vector: where the satellite is and how fast it is
//! going, in the TEME (True Equator, Mean Equinox) frame that SGP4
//! outputs.
//!
//! TEME is the native frame of SGP4 — an Earth-centred inertial frame
//! whose equator and equinox are those of the epoch. The TEME → ECEF
//! rotation lives in `orbitcore-astro::frames` (Week 2 Day 3); the
//! geodetic step (ECEF → latitude/longitude/altitude) lands on Day 4.
//! Until then, radius, speed, and altitude above the mean Earth radius
//! are the quantities downstream code needs, and they are
//! frame-independent.

use orbitcore_core::units::EARTH_MEAN_RADIUS_KM;

/// Position and velocity of a satellite at one instant, in the TEME
/// frame: position in kilometres, velocity in kilometres per second.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StateVector {
    /// Position in the TEME frame, km: `[x, y, z]`.
    pub position_km: [f64; 3],
    /// Velocity in the TEME frame, km/s: `[vx, vy, vz]`.
    pub velocity_km_s: [f64; 3],
    /// Minutes since the element set's epoch that this state is for.
    /// 0 is the epoch itself; negative values propagate backwards.
    pub minutes_since_epoch: f64,
}

impl StateVector {
    /// Distance from Earth's centre in km: `|r|`.
    pub fn radius_km(&self) -> f64 {
        (self.position_km[0].powi(2) + self.position_km[1].powi(2) + self.position_km[2].powi(2))
            .sqrt()
    }

    /// Speed in km/s: `|v|`.
    pub fn speed_km_s(&self) -> f64 {
        (self.velocity_km_s[0].powi(2)
            + self.velocity_km_s[1].powi(2)
            + self.velocity_km_s[2].powi(2))
        .sqrt()
    }

    /// Altitude above the Earth's surface in km: `|r| − R`.
    ///
    /// Uses the same mean Earth radius (6371 km) as
    /// [`OrbitalElements::altitude_km`], so the epoch altitude from the
    /// elements and the propagated altitude agree at `t = 0`.
    ///
    /// [`OrbitalElements::altitude_km`]: orbitcore_astro::elements::OrbitalElements::altitude_km
    pub fn altitude_km(&self) -> f64 {
        self.radius_km() - EARTH_MEAN_RADIUS_KM
    }

    /// Straight-line distance to another state vector, in km.
    pub fn distance_to(&self, other: &StateVector) -> f64 {
        let dx = self.position_km[0] - other.position_km[0];
        let dy = self.position_km[1] - other.position_km[1];
        let dz = self.position_km[2] - other.position_km[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn radius_speed_altitude() {
        // A circular 500 km orbit: |r| = 6871 km, |v| = 7.6126 km/s
        // (sqrt(mu / r) with mu = 398600.4418 km^3/s^2).
        let sv = StateVector {
            position_km: [6871.0, 0.0, 0.0],
            velocity_km_s: [0.0, 7.6126, 0.0],
            minutes_since_epoch: 0.0,
        };
        assert!((sv.radius_km() - 6871.0).abs() < 1e-9);
        assert!((sv.speed_km_s() - 7.6126).abs() < 1e-9);
        assert!((sv.altitude_km() - 500.0).abs() < 1e-9);
    }

    #[test]
    fn distance_to_is_euclidean() {
        let a = StateVector {
            position_km: [0.0, 0.0, 0.0],
            velocity_km_s: [0.0, 0.0, 0.0],
            minutes_since_epoch: 0.0,
        };
        let b = StateVector {
            position_km: [3.0, 4.0, 0.0],
            velocity_km_s: [0.0, 0.0, 0.0],
            minutes_since_epoch: 10.0,
        };
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-12);
        assert!((b.distance_to(&a) - 5.0).abs() < 1e-12);
    }
}
