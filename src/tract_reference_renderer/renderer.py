from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Sequence

from .vocal_tract.geometry import VocalTract
from .vocal_tract.params import ParamIndex

DEFAULT_WIDTH_PX = 640
DEFAULT_HEIGHT_PX = 640
_SVG_PADDING_PX = 34.0


@dataclass
class RenderDiagnostics:
    render_target: str
    param_count: int
    neutral_param_count: int
    parameter_space: str
    upper_outline_points: int
    lower_outline_points: int
    tongue_outline_points: int
    epiglottis_outline_points: int


@dataclass
class RenderResult:
    svg: str
    width_px: int
    height_px: int
    diagnostics: RenderDiagnostics


def _coerce_param_vector(values: Sequence[float]) -> List[float]:
    try:
        vector = [float(value) for value in values]
    except (TypeError, ValueError) as exc:
        raise ValueError("tract params must be numeric") from exc
    expected = int(ParamIndex.NUM_PARAMS)
    if len(vector) != expected:
        raise ValueError(f"tract params must have length {expected}")
    return vector


def _looks_like_owned_normalized_vector(vector: Sequence[float]) -> bool:
    if len(vector) != int(ParamIndex.NUM_PARAMS):
        return False
    for index, value in enumerate(vector):
        if index in (int(ParamIndex.TS1), int(ParamIndex.TS2), int(ParamIndex.TS3)):
            if not -1.0 <= float(value) <= 1.0:
                return False
        elif not 0.0 <= float(value) <= 1.0:
            return False
    return True


def _decode_owned_normalized_vector(vocal_tract: VocalTract, vector: Sequence[float]) -> List[float]:
    decoded: List[float] = []
    for index, value in enumerate(vector):
        param = vocal_tract.params[index]
        if index in (int(ParamIndex.TS1), int(ParamIndex.TS2), int(ParamIndex.TS3)):
            decoded.append(float(value))
            continue
        lo = float(param.min_val)
        neutral = float(param.neutral_val)
        hi = float(param.max_val)
        normalized = min(max(float(value), 0.0), 1.0)
        if normalized <= 0.5:
            decoded.append(lo + (normalized / 0.5) * (neutral - lo))
        else:
            decoded.append(neutral + ((normalized - 0.5) / 0.5) * (hi - neutral))
    return decoded


def _resolve_parameter_space(vocal_tract: VocalTract, vector: Sequence[float]) -> tuple[List[float], str]:
    if _looks_like_owned_normalized_vector(vector):
        return _decode_owned_normalized_vector(vocal_tract, vector), "owned_product_normalized_19d"
    return list(vector), "vtl_physical_19d"


def _iter_outline_points(vocal_tract: VocalTract) -> List[Any]:
    return [
        *vocal_tract.upper_outline.p,
        *vocal_tract.lower_outline.p,
        *vocal_tract.tongue_outline.p,
        *vocal_tract.epiglottis_outline.p,
    ]


def _build_view_transform(vocal_tract: VocalTract, width_px: int, height_px: int) -> tuple[float, float, float]:
    points = _iter_outline_points(vocal_tract)
    if not points:
        return 1.0, width_px * 0.5, height_px * 0.5
    min_x = min(float(point.x) for point in points)
    max_x = max(float(point.x) for point in points)
    min_y = min(float(point.y) for point in points)
    max_y = max(float(point.y) for point in points)
    span_x = max(max_x - min_x, 1.0e-6)
    span_y = max(max_y - min_y, 1.0e-6)
    drawable_width = max(float(width_px) - 2.0 * _SVG_PADDING_PX, 1.0)
    drawable_height = max(float(height_px) - 2.0 * _SVG_PADDING_PX, 1.0)
    scale = min(drawable_width / span_x, drawable_height / span_y)
    offset_x = (float(width_px) - span_x * scale) * 0.5 - min_x * scale
    offset_y = (float(height_px) - span_y * scale) * 0.5 + max_y * scale
    return scale, offset_x, offset_y


def _polyline_points(points: Sequence[Any], scale: float, offset_x: float, offset_y: float) -> str:
    return " ".join(
        f"{point.x * scale + offset_x:.2f},{offset_y - point.y * scale:.2f}"
        for point in points
    )


def _build_svg(vocal_tract: VocalTract, width_px: int, height_px: int) -> str:
    scale, offset_x, offset_y = _build_view_transform(vocal_tract, width_px, height_px)
    upper_points = _polyline_points(vocal_tract.upper_outline.p, scale, offset_x, offset_y)
    lower_points = _polyline_points(vocal_tract.lower_outline.p, scale, offset_x, offset_y)
    tongue_points = _polyline_points(vocal_tract.tongue_outline.p, scale, offset_x, offset_y)
    epiglottis_points = _polyline_points(vocal_tract.epiglottis_outline.p, scale, offset_x, offset_y)
    return (
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width_px}' height='{height_px}' "
        f"viewBox='0 0 {width_px} {height_px}'>"
        f"<rect x='0' y='0' width='{width_px}' height='{height_px}' fill='white'/>"
        f"<polyline fill='none' stroke='#1f77b4' stroke-width='2' points='{upper_points}'/>"
        f"<polyline fill='none' stroke='#d62728' stroke-width='2' points='{lower_points}'/>"
        f"<polyline fill='none' stroke='#2ca02c' stroke-width='1.5' points='{tongue_points}'/>"
        f"<polyline fill='none' stroke='#9467bd' stroke-width='1.5' points='{epiglottis_points}'/>"
        "</svg>"
    )


def _apply_params(vocal_tract: VocalTract, vector: Sequence[float]) -> str:
    resolved_vector, parameter_space = _resolve_parameter_space(vocal_tract, vector)
    for index, value in enumerate(vector):
        param = vocal_tract.params[index]
        bounded = min(max(float(resolved_vector[index]), float(param.min_val)), float(param.max_val))
        param.value = bounded
        param.limited_value = bounded
    vocal_tract.calculate_all()
    return parameter_space


def neutral_param_vector() -> List[float]:
    vocal_tract = VocalTract()
    return [float(param.neutral_val) for param in vocal_tract.params]


def render_tract_svg(
    tract_params: Sequence[float],
    *,
    width_px: int = DEFAULT_WIDTH_PX,
    height_px: int = DEFAULT_HEIGHT_PX,
    render_target: str = "current",
) -> RenderResult:
    vector = _coerce_param_vector(tract_params)
    vocal_tract = VocalTract()
    parameter_space = _apply_params(vocal_tract, vector)
    diagnostics = RenderDiagnostics(
        render_target=render_target,
        param_count=len(vector),
        neutral_param_count=sum(
            1 for index, value in enumerate(vector) if value == float(vocal_tract.params[index].neutral_val)
        ),
        parameter_space=parameter_space,
        upper_outline_points=len(vocal_tract.upper_outline.p),
        lower_outline_points=len(vocal_tract.lower_outline.p),
        tongue_outline_points=len(vocal_tract.tongue_outline.p),
        epiglottis_outline_points=len(vocal_tract.epiglottis_outline.p),
    )
    return RenderResult(
        svg=_build_svg(vocal_tract, width_px, height_px),
        width_px=width_px,
        height_px=height_px,
        diagnostics=diagnostics,
    )


def render_svg_pair(
    *,
    current_tract_params: Sequence[float],
    target_tract_params: Optional[Sequence[float]] = None,
    render_target: str = "current",
    width_px: int = DEFAULT_WIDTH_PX,
    height_px: int = DEFAULT_HEIGHT_PX,
) -> Dict[str, Any]:
    if render_target not in {"current", "target", "both"}:
        raise ValueError("render_target must be one of: current, target, both")
    current_result: Optional[RenderResult] = None
    target_result: Optional[RenderResult] = None
    if render_target in {"current", "both"}:
        current_result = render_tract_svg(
            current_tract_params,
            width_px=width_px,
            height_px=height_px,
            render_target="current",
        )
    if render_target in {"target", "both"}:
        if target_tract_params is None:
            raise ValueError("target_tract_params is required when render_target is 'target' or 'both'")
        target_result = render_tract_svg(
            target_tract_params,
            width_px=width_px,
            height_px=height_px,
            render_target="target",
        )
    diagnostics: Dict[str, Any] = {"render_target": render_target}
    if current_result is not None:
        diagnostics["current"] = asdict(current_result.diagnostics)
    if target_result is not None:
        diagnostics["target"] = asdict(target_result.diagnostics)
    return {
        "current_svg": current_result.svg if current_result is not None else None,
        "target_svg": target_result.svg if target_result is not None else None,
        "width_px": width_px,
        "height_px": height_px,
        "diagnostics": diagnostics,
    }
