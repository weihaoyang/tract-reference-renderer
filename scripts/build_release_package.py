from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "dist"
PACKAGE_ROOT_NAME = "tract-reference-renderer"
INCLUDE_PATHS = [
    "LICENSE",
    "NOTICE.md",
    "README.md",
    "UPSTREAM_PROVENANCE.md",
    "pyproject.toml",
    "scripts/smoke_render.py",
    "scripts/smoke_ws_client.py",
    "scripts/start_renderer.ps1",
    "src",
]


def _git_value(args: list[str]) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()
    except Exception:
        return ""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_tree(staging: Path) -> None:
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    for relative in INCLUDE_PATHS:
        source = ROOT / relative
        target = staging / relative
        if source.is_dir():
            shutil.copytree(source, target, ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"))
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)


def _write_runtime_manifest(staging: Path, version: str) -> dict[str, object]:
    manifest = {
        "schema_version": "tract-reference-renderer.release.v1",
        "helper_name": "tract-reference-renderer",
        "helper_version": version,
        "protocol_version": "coach.reference_renderer.v1",
        "license_status": "external_gpl_helper",
        "ipc_boundary": "localhost_websocket",
        "default_ws_url": "ws://127.0.0.1:8076",
        "clinical_truth_claim_allowed": False,
        "truth_tier": "reference_visualization_not_patient_truth",
        "git_commit": _git_value(["rev-parse", "HEAD"]),
        "git_remote": _git_value(["remote", "get-url", "origin"]),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "entrypoints": {
            "python_module": "python -m tract_reference_renderer",
            "windows_script": "scripts/start_renderer.ps1",
        },
        "health_request": {"type": "health", "request_id": "health-1"},
    }
    manifest_path = staging / "release_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def build_package(version: str) -> Path:
    DIST.mkdir(exist_ok=True)
    staging = DIST / f"{PACKAGE_ROOT_NAME}-{version}"
    _copy_tree(staging)
    manifest = _write_runtime_manifest(staging, version)
    zip_path = DIST / f"{PACKAGE_ROOT_NAME}-{version}.zip"
    if zip_path.exists():
        zip_path.unlink()
    with ZipFile(zip_path, "w", ZIP_DEFLATED) as archive:
        for file_path in sorted(staging.rglob("*")):
            if file_path.is_file():
                archive.write(file_path, Path(staging.name) / file_path.relative_to(staging))
    manifest["artifact"] = {
        "file": zip_path.name,
        "sha256": _sha256(zip_path),
        "size_bytes": zip_path.stat().st_size,
    }
    (DIST / f"{PACKAGE_ROOT_NAME}-{version}.release.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return zip_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a standalone helper release package.")
    parser.add_argument("--version", default="0.1.0")
    args = parser.parse_args()
    zip_path = build_package(args.version)
    print(f"artifact={zip_path}")
    print(f"sha256={_sha256(zip_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
