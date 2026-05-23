# tract-reference-renderer

Small standalone Python helper package and localhost WebSocket service for rendering midsagittal tract SVGs from 19D tract parameter vectors.

## Scope

This repository contains GPL research/helper code, including copied or transliterated VocalTractLab-style geometry and rendering logic under `src/tract_reference_renderer/`. Treat it as a local research/helper service, not product-runtime code.

## What It Does

- exposes a Python API for rendering tract SVG from a 19-parameter vector
- runs a local WebSocket server on `127.0.0.1:8076`
- accepts JSON requests with:
  - `request_id: string`
  - `current_tract_params: array[19]`
  - `target_tract_params?: array[19]`
  - `render_target?: "current" | "target" | "both"`
- returns JSON with:
  - `request_id`
  - `status: "ok" | "error"`
  - `current_svg?`
  - `target_svg?`
  - `width_px`
  - `height_px`
  - `diagnostics?`
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

## Provenance And Licensing

See [LICENSE](LICENSE), [NOTICE.md](NOTICE.md), and [UPSTREAM_PROVENANCE.md](UPSTREAM_PROVENANCE.md).

Legal posture adopted for this repository:

- This repository is distributed under a conservative GPL helper-repository posture.
- The full text currently shipped in `LICENSE` is the GNU General Public License, Version 3.
- The repository's own evidence supports treating the code as GPL-derived helper code, but it does **not** prove that the exact upstream source-version mapping for every copied or transliterated segment has already been fully reconstructed.
- Until a more complete provenance review is documented, treat this repository as GPL-governed helper code for local research or tooling workflows, not as commercially cleared product-runtime code.
- The neutral runtime name does not change the underlying provenance or license obligations.
