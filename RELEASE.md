# Release Process

This helper is released as a standalone GPL-governed local reference renderer. It must not be bundled into the commercial product runtime.

## Build

```powershell
cd F:\vocal-app-v3-authoritative\external_local\tract-reference-renderer
py -3 scripts\build_release_package.py --version 0.1.0
```

The build writes:

- `dist/tract-reference-renderer-0.1.0.zip`
- `dist/tract-reference-renderer-0.1.0.release.json`

## Verify

```powershell
$env:PYTHONPATH = "F:\vocal-app-v3-authoritative\external_local\tract-reference-renderer\src"
py -3 scripts\smoke_render.py
```

For WebSocket health:

```json
{"type":"health","request_id":"health-1"}
```

Expected fields include:

- `helper_version`
- `protocol_version`
- `license_status=external_gpl_helper`
- `truth_tier=reference_visualization_not_patient_truth`
- `clinical_truth_claim_allowed=false`

## Publish To GitHub

After `gh auth login` is valid:

```powershell
git tag v0.1.0
git push origin v0.1.0
gh release create v0.1.0 dist\tract-reference-renderer-0.1.0.zip dist\tract-reference-renderer-0.1.0.release.json --title "tract-reference-renderer v0.1.0" --notes-file RELEASE.md
```
