from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Sequence

from .protocol import DEFAULT_IPC_BOUNDARY, HELPER_NAME, HELPER_VERSION, LICENSE_STATUS, PROTOCOL_VERSION
from .renderer import _apply_params, _coerce_param_vector
from .vocal_tract.geometry import VocalTract
from .vocal_tract import params as tract_params


@dataclass
class Mesh3DDiagnostics:
    parameter_space: str
    surface_count: int
    triangle_count: int
    path_count: int


def _surface_name_from_index(surface_index: int) -> str:
    try:
        return tract_params.SurfaceIndex(surface_index).name.lower()
    except ValueError:
        return f"surface_{surface_index}"


def _serialize_surface_mesh(surface_index: int, surface: Any) -> Dict[str, Any] | None:
    if surface is None:
        return None
    num_ribs = int(getattr(surface, "num_ribs", 0))
    num_points = int(getattr(surface, "num_points", 0))
    if num_ribs <= 0 or num_points <= 0:
        return None

    positions: List[float] = []
    for rib in range(num_ribs):
        for point_index in range(num_points):
            vertex = surface.get_vertex(rib, point_index)
            positions.extend([float(vertex.x), float(vertex.y), float(vertex.z)])

    indices: List[int] = []
    orientation_swapped = bool(getattr(surface, "_orientation_swapped", False))
    for rib in range(max(0, num_ribs - 1)):
        for point_index in range(max(0, num_points - 1)):
            p00 = rib * num_points + point_index
            p01 = rib * num_points + point_index + 1
            p10 = (rib + 1) * num_points + point_index
            p11 = (rib + 1) * num_points + point_index + 1
            if orientation_swapped:
                indices.extend([p00, p10, p01, p10, p11, p01])
            else:
                indices.extend([p00, p01, p10, p10, p01, p11])

    return {
        "surface_index": int(surface_index),
        "surface_name": _surface_name_from_index(surface_index),
        "num_ribs": num_ribs,
        "num_points": num_points,
        "num_vertices": num_ribs * num_points,
        "num_triangles": len(indices) // 3,
        "positions": positions,
        "indices": indices,
    }


def _serialize_polyline_3d(path_name: str, line_strip: Any) -> Dict[str, Any]:
    positions: List[float] = []
    for point in getattr(line_strip, "p", []) or []:
        positions.extend([float(point.x), float(point.y), float(getattr(point, "z", 0.0))])
    return {"name": path_name, "positions": positions}


def build_tract_mesh3d_payload(tract_params_vector: Sequence[float]) -> Dict[str, Any]:
    vector = _coerce_param_vector(tract_params_vector)
    vocal_tract = VocalTract()
    parameter_space = _apply_params(vocal_tract, vector)

    surfaces: List[Dict[str, Any]] = []
    for surface_index, surface in enumerate(getattr(vocal_tract, "surfaces_list", [])):
        payload = _serialize_surface_mesh(surface_index, surface)
        if payload is not None:
            surfaces.append(payload)

    paths = [
        _serialize_polyline_3d("wide_lip_corner_path", getattr(vocal_tract, "wide_lip_corner_path", None)),
        _serialize_polyline_3d("narrow_lip_corner_path", getattr(vocal_tract, "narrow_lip_corner_path", None)),
        _serialize_polyline_3d("lip_corner_path", getattr(vocal_tract, "lip_corner_path", None)),
    ]
    diagnostics = Mesh3DDiagnostics(
        parameter_space=parameter_space,
        surface_count=len(surfaces),
        triangle_count=sum(int(surface["num_triangles"]) for surface in surfaces),
        path_count=sum(1 for path in paths if len(path["positions"]) > 0),
    )
    return {
        "return_code": 0,
        "surfaces": surfaces,
        "paths": paths,
        "geometry_provenance": {
            "runtime": "tract_reference_renderer_ipc",
            "lineage": "standalone_vtl_reference_renderer_ipc_3d",
            "helper_name": HELPER_NAME,
            "helper_version": HELPER_VERSION,
            "protocol_version": PROTOCOL_VERSION,
            "parameter_schema": parameter_space,
            "ipc_boundary": DEFAULT_IPC_BOUNDARY,
            "external_reference_helper": True,
            "vtl_lineage": True,
            "fallback_allowed": False,
            "clinical_truth_claim_allowed": False,
        },
        "license_status": "external_reference_renderer_ipc",
        "helper_license_status": LICENSE_STATUS,
        "truth_tier": "reference_visualization_not_patient_truth",
        "clinical_truth_claim_allowed": False,
        "medical_anatomy_product_allowed": False,
        "medical_anatomy_gate": {
            "schema_version": "coach.reference_renderer.geometry_3d.v1",
            "status": "reference_only",
            "medical_anatomy_product_allowed": False,
            "clinical_truth_claim_allowed": False,
            "blockers": ["external_reference_visualization_only"],
        },
        "fem_case_handoff_manifest": {"schema_version": "coach.reference_renderer.geometry_3d.v1", "status": "not_applicable"},
        "fem_transfer_manifest": {"schema_version": "coach.reference_renderer.geometry_3d.v1", "status": "not_applicable"},
        "surrogate_manifest": {"schema_version": "coach.reference_renderer.geometry_3d.v1", "status": "not_applicable"},
        "reference_renderer_diagnostics": asdict(diagnostics),
    }


__all__ = ["build_tract_mesh3d_payload"]
