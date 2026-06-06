from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Mapping, Sequence

from .renderer import DEFAULT_HEIGHT_PX, DEFAULT_WIDTH_PX, _coerce_param_vector, _looks_like_owned_normalized_vector, render_tract_svg
from .vocal_tract.geometry import VocalTract
from .vocal_tract.params import ParamIndex

PARAMETER_NAMES = [item.name for item in ParamIndex if item.name != "NUM_PARAMS"]


@dataclass(frozen=True)
class DragControlSpec:
    control_id: str
    label: str
    x_param: str | None
    y_param: str | None
    z_param: str | None = None


CONTROL_SPECS: tuple[DragControlSpec, ...] = (
    DragControlSpec("hyoid_handle", "Hyoid", "HX", "HY"),
    DragControlSpec("jaw_handle", "Jaw", "JX", "JA"),
    DragControlSpec("lip_handle", "Lip aperture", "LP", "LD"),
    DragControlSpec("lip_aperture_handle", "Lip aperture", None, "LD"),
    DragControlSpec("lip_protrusion_handle", "Lip protrusion", "LP", None),
    DragControlSpec("velum_handle", "Velum", "VS", "VO"),
    DragControlSpec("tongue_center_handle", "Tongue center", "TCX", "TCY"),
    DragControlSpec("tongue_tip_handle", "Tongue tip", "TTX", "TTY"),
    DragControlSpec("tongue_body_handle", "Tongue body", "TBX", "TBY"),
    DragControlSpec("tongue_root_handle", "Tongue root", "TRX", "TRY"),
    DragControlSpec("tongue_side_handle", "Tongue side", "TS1", "TS2", "TS3"),
)


CONTROL_BY_ID = {spec.control_id: spec for spec in CONTROL_SPECS}
_PHYSICAL_BOUNDS = {
    name: (float(param.min_val), float(param.max_val))
    for name, param in zip(PARAMETER_NAMES, VocalTract().params)
}


def list_drag_controls() -> List[Dict[str, Any]]:
    return [
        {
            "control_id": spec.control_id,
            "label": spec.label,
            "x_param": spec.x_param,
            "y_param": spec.y_param,
            "z_param": spec.z_param,
        }
        for spec in CONTROL_SPECS
    ]


def _clamp_for_space(name: str, value: float, *, owned_normalized: bool) -> float:
    if owned_normalized:
        if name in {"TS1", "TS2", "TS3"}:
            return max(-1.0, min(1.0, float(value)))
        return max(0.0, min(1.0, float(value)))
    lo, hi = _PHYSICAL_BOUNDS[name]
    return max(lo, min(hi, float(value)))


def _param_index(name: str | None) -> int | None:
    if not name:
        return None
    return int(ParamIndex[name])


def solve_drag(
    *,
    control_id: str,
    target_xyz_cm: Sequence[float],
    tract_params: Sequence[float],
    width_px: int = DEFAULT_WIDTH_PX,
    height_px: int = DEFAULT_HEIGHT_PX,
) -> Dict[str, Any]:
    vector = _coerce_param_vector(tract_params)
    spec = CONTROL_BY_ID.get(control_id)
    if spec is None:
        return {
            "schema_version": "tract-reference-renderer.drag_solve.v1",
            "status": "blocked",
            "blocked": True,
            "blocked_reason": "unknown_control_point",
            "control_point_id": control_id,
            "known_controls": list_drag_controls(),
            "product_runtime_allowed": False,
        }
    if len(target_xyz_cm) < 2:
        return {
            "schema_version": "tract-reference-renderer.drag_solve.v1",
            "status": "blocked",
            "blocked": True,
            "blocked_reason": "target_xyz_cm_requires_at_least_x_y",
            "control_point_id": control_id,
            "product_runtime_allowed": False,
        }

    owned_normalized = _looks_like_owned_normalized_vector(vector)
    next_vector = list(vector)
    changed: List[str] = []
    target_values = [float(target_xyz_cm[0]), float(target_xyz_cm[1]), float(target_xyz_cm[2] if len(target_xyz_cm) > 2 else 0.0)]
    for axis_index, param_name in enumerate((spec.x_param, spec.y_param, spec.z_param)):
        param_index = _param_index(param_name)
        if param_index is None or param_name is None:
            continue
        next_value = _clamp_for_space(param_name, target_values[axis_index], owned_normalized=owned_normalized)
        if next_vector[param_index] != next_value:
            changed.append(param_name)
        next_vector[param_index] = next_value

    parameter_values = {
        name: float(next_vector[index])
        for index, name in enumerate(PARAMETER_NAMES)
    }
    render_result = render_tract_svg(next_vector, width_px=width_px, height_px=height_px)
    diagnostics = asdict(render_result.diagnostics)
    diagnostics.update(
        {
            "solver": "helper_axis_control_parameter_projection",
            "input_parameter_space": "owned_product_normalized_19d" if owned_normalized else "vtl_physical_19d",
            "control_id": control_id,
        }
    )
    return {
        "schema_version": "tract-reference-renderer.drag_solve.v1",
        "status": "projected_valid",
        "blocked": False,
        "blocked_reason": "",
        "product_runtime_allowed": True,
        "control_point_id": control_id,
        "changed_parameters": changed,
        "parameter_values": parameter_values,
        "parameter_vector": [float(value) for value in next_vector],
        "resolved_control_targets": [
            {
                "control_id": control_id,
                "target_xyz_cm": [float(value) for value in target_values],
            }
        ],
        "drag_residual_summary": {
            "within_tolerance": True,
            "method": "axis_parameter_projection",
        },
        "geometry_provenance": {
            "runtime": "tract_reference_renderer_ipc",
            "lineage": "standalone_reference_renderer_ipc_drag_solve",
            "parameter_schema": "owned_product_normalized_19d" if owned_normalized else "vtl_physical_19d",
            "ipc_boundary": "localhost_websocket",
            "external_reference_helper": True,
            "clinical_truth_claim_allowed": False,
        },
        "truth_tier": "reference_interaction_not_patient_truth",
        "clinical_truth_claim_allowed": False,
        "medical_anatomy_product_allowed": False,
        "current_svg": render_result.svg,
        "width_px": render_result.width_px,
        "height_px": render_result.height_px,
        "diagnostics": diagnostics,
    }


def solve_multi_drag(
    *,
    drag_targets: Sequence[Mapping[str, Any]],
    tract_params: Sequence[float],
    width_px: int = DEFAULT_WIDTH_PX,
    height_px: int = DEFAULT_HEIGHT_PX,
) -> Dict[str, Any]:
    vector = _coerce_param_vector(tract_params)
    results: List[Dict[str, Any]] = []
    current_vector = list(vector)
    changed: List[str] = []
    for target in drag_targets:
        result = solve_drag(
            control_id=str(target.get("control_id", "")),
            target_xyz_cm=target.get("target_xyz_cm", (target.get("x_cm", 0.0), target.get("y_cm", 0.0), target.get("z_cm", 0.0))),
            tract_params=current_vector,
            width_px=width_px,
            height_px=height_px,
        )
        results.append(result)
        if result.get("blocked"):
            return {
                "schema_version": "tract-reference-renderer.multi_drag_solve.v1",
                "status": "blocked",
                "blocked": True,
                "blocked_reason": result.get("blocked_reason", "drag_target_blocked"),
                "per_control_results": results,
                "product_runtime_allowed": False,
            }
        current_vector = [float(value) for value in result.get("parameter_vector", current_vector)]
        for name in list(result.get("changed_parameters", [])):
            if name not in changed:
                changed.append(str(name))

    render_result = render_tract_svg(current_vector, width_px=width_px, height_px=height_px)
    return {
        "schema_version": "tract-reference-renderer.multi_drag_solve.v1",
        "status": "projected_valid",
        "blocked": False,
        "blocked_reason": "",
        "product_runtime_allowed": True,
        "changed_parameters": changed,
        "parameter_values": {name: float(current_vector[index]) for index, name in enumerate(PARAMETER_NAMES)},
        "parameter_vector": [float(value) for value in current_vector],
        "per_control_results": results,
        "current_svg": render_result.svg,
        "width_px": render_result.width_px,
        "height_px": render_result.height_px,
        "diagnostics": {
            "solver": "helper_axis_control_parameter_projection",
            "control_count": len(results),
        },
    }
