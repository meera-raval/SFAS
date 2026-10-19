# S-FAS mixed-field localization — MATLAB implementation

Implements T. Wu et al., "Scalable Fluid Antenna Systems for Mixed-Field Source
Localization," IEEE ICC 2026 (DOI 10.1109/ICC59461.2026.11586773).
Runs in MATLAB (R2016b+, no toolboxes) and GNU Octave (tested on 8.4).

## Files
| File | Paper reference |
|---|---|
| `esg_steering.m`    | Eqs. (2), (4), (5) — exact spatial geometry steering vector |
| `ff_steering.m`     | far-field steering vector a_FF (Sec. II-C) |
| `coupling_matrix.m` | Eq. (6) — banded symmetric Toeplitz mutual coupling |
| `gen_snapshots.m`   | Eqs. (7), (8) — source/noise generation |
| `sfas_stage1.m`     | Algorithm 1, Eqs. (9)–(11) — central-subarray selection + MUSIC |
| `sfas_stage2.m`     | Algorithm 2, Eqs. (12)–(15) — 1-D range search + local 2-D MUSIC |
| `sfas_icc_demo.m`   | Sec. IV-A single-shot example (Figs. 2, 3) |
| `sfas_icc_mc.m`     | Monte-Carlo check of the same scenario (not in the paper) |

## Run
```matlab
sfas_icc_demo                                  % paper pipeline exactly as written
init_mode = 'coarse2d'; n_pass = 2; sfas_icc_demo   % robust variant (see below)
T = 100; init_mode = 'coarse2d'; n_pass = 2; sfas_icc_mc
```

## Assumptions (the paper does not specify these)
* Coupling: beta = 2 (stated), c0 = 0.33 so that |c1(alpha_c)| ≈ 0.30 (stated),
  phi0 = pi/4, coupling range p = 3 (central 26 of 32 elements kept).
  Coupling is applied in Stage 1 only; Stage 2 uses Eq. (8) with no coupling.
* Data: true channel always uses the exact ESG model (Stage 1 estimates with the
  far-field model, so the approximation error is included). Independent
  snapshots in each stage, unit-power circular Gaussian sources, SNR per source.
* Reference point: element 1 (p_1 = 0), as Eqs. (2) and (4) are written
  (the text says "array center"; the equations use element 1).
* Grids: Stage 1 0.1° over [-90°, 90°]; Stage 2 1-D range grid 3000 log-spaced
  points in [5λ, 10^4 λ]; local window Δθ = 5°, Δr = 0.2 r_max (stated), step
  0.01° in angle and 800 points uniform in 1/r.
* MUSIC denominators are normalised by ||a||² (`opt.normalize = true`), because
  the ESG amplitude term r/r_m makes unnormalised spectra favour small r.
  Set `opt.normalize = false` for Eqs. (13)/(15) exactly as written.

## Deviations available as options (NOT in the paper)
* `opt.init = 'coarse2d'`: replaces the fixed-angle 1-D range search of Eq. (14)
  with a coarse 2-D scan (0.25° steps) inside ±Δθ over the full range grid.
* `opt.n_pass = 2`: after the local 2-D search, redo Eq. (14) at the refined
  angle and run the local 2-D search again.

Why they exist: the Stage-1 aperture is only ~1.25λ, so Stage-1 DOA errors are
~0.5–1.5° at 10 dB. The fixed-angle range spectrum of Eq. (13) is very sensitive
to that error (0.2° shifts the 300λ peak to ~263λ; 1° shifts it to ~178λ, or
pushes it to r_max), and with Δr around a wrong r_init the true range is never
searched.

## Parameter consistency check: ICC paper vs. journal version (JSTSP 2026)
Both papers *state* alpha_c = 0.1, alpha_e = 1 (d_c = 0.05λ, d_e = 0.5λ), but the
journal's figures and text use d_c = 0.1λ, d_e = 1λ (alpha_c = 0.2, alpha_e = 2):
Sec. IV-A ("compressing the array to d_c = 0.1λ", "extended configuration (d_e = 1λ)"),
Fig. 4 and Fig. 5 captions, the "MUSIC/ESPRIT (d = 0.1λ)" legend entries in Fig. 7,
and CRB2 (d = 1λ) as the bound AAR "closely tracks". The ICC coupling claims also only
fit the second set: with beta = 2 and |c1(alpha_c)| = 0.30, |c1(alpha_e)| is 0.122 for
(0.1, 1) but 0.050 for (0.2, 2), which is the stated "< 0.05".

`journal_mc.m` runs both settings with the journal's Stage-2 grids
(Δθ = 5°, Δr = 60λ, δθ = 0.015°, δr = 0.1λ, G_r = 216), paper initialisation, no fixes:
* cfg 'B' (stated 0.05λ / 0.5λ): frequent Stage-2 failures for the 300λ–5000λ sources.
* cfg 'C' (0.1λ / 1λ): final DOA errors <= 0.05° for all four sources in nearly every
  trial; range error < 1% at 30λ and 300λ; 1000λ shows a constant ~5.7% bias (the
  216-point 1-D grid puts r_init > 60λ away, so the true range sits outside Δr);
  5000λ is mostly within 3.4%, with an occasional failure.
