# tract-reference-renderer

Small standalone Python helper package and localhost WebSocket service for rendering midsagittal tract SVGs and helper-owned 3D tract meshes from 19D tract parameter vectors.

## Scope

This repository contains GPL research/helper code, including copied or transliterated VocalTractLab-style geometry and rendering logic under `src/tract_reference_renderer/`. Treat it as a local research/helper service, not product-runtime code.

## What It Does

- exposes a Python API for rendering tract SVG from a 19-parameter vector
- exposes helper-owned Lab control drag solve requests for interactive 2D fitting workflows
- exposes `geometry_3d` WebSocket responses with surface meshes and lip paths for product UIs that must not use an internal geometry fallback
- runs a local WebSocket server on `127.0.0.1:8076`
- accepts JSON requests with:
  - `type?: "geometry_3d" | "health"`
  - `request_id: string`
  - `current_tract_params: array[19]`
  - `target_tract_params?: array[19]`
  - `render_target?: "current" | "target" | "both"`
  - `type?: "health" | "list_controls" | "drag_solve" | "multi_drag_solve"`
- returns JSON with:
  - `request_id`
  - `status: "ok" | "error"`
  - `current_svg?`
  - `target_svg?`
  - `surfaces?` and `paths?` for `type: "geometry_3d"`
  - `width_px`
  - `height_px`
  - `diagnostics?`
  - `solve?`
  - `controls?`
  - `error?`

## Install

```powershell
cd F:\tract-reference-renderer
py -3 -m pip install -e .
```

## Run The Server

```powershell
cd F:\tract-reference-renderer
py -3 -m tract_reference_renderer
```

Windows helper script:

```powershell
.\scripts\start_renderer.ps1 -HostName 127.0.0.1 -Port 8076
```

Optional host/port override:

```powershell
$env:TRACT_RENDERER_HOST = "127.0.0.1"
$env:TRACT_RENDERER_PORT = "8876"
py -3 -m tract_reference_renderer
```

## Example Request

```json
{
  "request_id": "demo-1",
  "current_tract_params": [1.0, -4.75, 0.0, -2.0, -0.07, 0.95, 0.0, -0.1, -0.4, -1.46, 3.5, -1.0, 2.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0],
  "render_target": "current"
}
```

## Smoke Test

```powershell
cd F:\tract-reference-renderer
py -3 scripts/smoke_render.py
```

## Health Request

Send this JSON over the local WebSocket:

```json
{"type": "health", "request_id": "health-1"}
```

The response includes helper version, protocol version, license status, IPC boundary, and truth-tier fields. The renderer is a reference visualization helper only; it does not release clinical or patient-truth claims.

## Drag Solve Request

Send this JSON over the local WebSocket:

```json
{
  "type": "drag_solve",
  "request_id": "drag-1",
  "current_tract_params": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0],
  "control_id": "tongue_body_handle",
  "x_cm": 0.52,
  "y_cm": 0.62,
  "z_cm": 0.0
}
```

The response returns `solve.parameter_values` keyed by the 19D parameter names (`HX`, `HY`, `JX`, ...), a `current_svg`, helper provenance, and `truth_tier: reference_interaction_not_patient_truth`. The solve is an interactive reference-helper projection, not patient anatomy, clinical truth, or a robot execution release.

## Release Package

```powershell
py -3 scripts/build_release_package.py --version 0.1.0
```

Artifacts are written under `dist/` and are intended for the helper repository GitHub release, not for bundling into the product runtime.

## Provenance And Licensing

See [LICENSE](LICENSE), [NOTICE.md](NOTICE.md), and [UPSTREAM_PROVENANCE.md](UPSTREAM_PROVENANCE.md).

Legal posture adopted for this repository:

- This repository is distributed under a conservative GPL helper-repository posture.
- The full text currently shipped in `LICENSE` is the GNU General Public License, Version 3.
- The repository's own evidence supports treating the code as GPL-derived helper code, but it does **not** prove that the exact upstream source-version mapping for every copied or transliterated segment has already been fully reconstructed.
- Until a more complete provenance review is documented, treat this repository as GPL-governed helper code for local research or tooling workflows, not as commercially cleared product-runtime code.
- The neutral runtime name does not change the underlying provenance or license obligations.
