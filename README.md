# tract-reference-renderer

Small standalone Python helper package and localhost WebSocket service for rendering midsagittal tract SVGs from 19D tract parameter vectors.

## Scope

This repository contains GPL research/helper code, including copied or transliterated VocalTractLab-style geometry and rendering logic under `src/vtl_renderer/`. Treat it as a local research/helper service, not product-runtime code.

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
py -3 -m vtl_renderer
```

Optional host/port override:

```powershell
$env:TRACT_RENDERER_HOST = "127.0.0.1"
$env:TRACT_RENDERER_PORT = "8876"
py -3 -m vtl_renderer
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

See [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md). The legal posture here is intentionally narrow and cautious: this repo should be handled as GPL-derived helper code for local research or tooling workflows. The neutral runtime name does not change the underlying provenance or license obligations.
