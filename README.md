# OrbitCore

> **Design, simulate, and fly a CubeSat from a single Rust codebase — orbit, power, link, and operations in one mission digital twin.**

OrbitCore is an open-source CubeSat mission platform. One `mission.toml` file
describes the whole mission — orbit, bus, solar panels, battery, radio,
payload, ground stations — and OrbitCore answers the question every small
university CubeSat team asks first: *"Can this satellite actually fly this
mission with the energy and link it has?"*

Today that answer is scattered across GMAT/STK, a power-budget spreadsheet, a
link-budget web calculator, and a thermal spreadsheet — none of which check
each other. OrbitCore puts orbit propagation, ground-station visibility,
eclipse modelling, power/energy budget, link budget, and an energy-aware
operations scheduler in one codebase, validated against real sky data
(CelesTrak TLEs and SatNOGS observations).

## Status

**Phase 1 — Orbit & Visibility (Weeks 1–4), in progress.** Day 1: workspace,
CI, and this README. Day 2: `orbitcore-core` gains type-safe SI units
(`uom`) — Kepler's third law (500 km → ≈ 94.6 min, matching the literature),
solar-panel power, and the shared `OrbitCoreError` type. Day 3: new crate
`orbitcore-astro` with a strict TLE parser (fixed-width fields,
implied-decimal BSTAR/eccentricity, column-69 checksum, epoch → calendar
conversion, period/semi-major axis from mean motion) — validated against
the live ISS TLE from CelesTrak (epoch 2026-11-12, 92.94 min period,
≈ 426 km) and the 2008 canonical example. Day 4: OMM parser
(`orbitcore-astro/src/omm.rs`, JSON via serde) for the same element
sets — tolerant of missing/extra fields, strict about epoch and mean
motion — plus a shared `OrbitalElements` type so a TLE and an OMM for
the same satellite are compared field by field and must agree (ISS,
2026-11-12 epoch: identical period 92.94 min, inclination, and
eccentricity in both formats). The goal of Phase 1
is a working `orbitcore pass` command driven by real TLE data, matching
CelesTrak/SatNOGS.

Roadmap (16 weeks, 4 phases):

| Phase | Focus | Deliverable |
|---|---|---|
| 1 | Orbit & visibility | SGP4 propagation, ground-station passes, eclipse & beta-angle model, `orbitcore pass` |
| 2 | Energy & link | Solar + Li-ion battery + EPS models, Friis link budget, data budget, `orbitcore budget` |
| 3 | Digital twin & scheduler | Discrete-event simulator, energy-aware scheduler (threshold policy), Monte Carlo, web dashboard v1 |
| 4 | Operations & flight | Axum + MQTT ops backend, 3D dashboard, Embassy flight-software skeleton on the same `orbitcore-core`, v0.1.0 |

## Quickstart

```bash
cargo build
cargo test
```

## Design principles

1. **Pure core, no I/O.** `orbitcore-core` is only math/physics — no network,
   no files, no system clock. It is deterministic, testable, `no_std`-capable,
   and shared by the simulator and the flight software (sim-to-flight).
2. **Validation is a feature.** Every release compares propagation against
   CelesTrak data and pass predictions against SatNOGS observations.
3. **One source of truth.** `mission.toml` is the single mission definition;
   the CLI, simulator, scheduler, and dashboard all read the same file.
4. **Physics first, features later.** No feature lands before its physics
   model has a validation test.

## License

Dual-licensed under MIT or Apache-2.0 — see `LICENSE-MIT` and
`LICENSE-APACHE`.
