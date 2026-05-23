from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Sequence

from .vocal_tract.geometry import VocalTract
from .vocal_tract.params import ParamIndex

DEFAULT_WIDTH_PX = 640
DEFAULT_HEIGHT_PX = 280
_DEFAULT_SCALE = 40.0
_DEFAULT_X_OFFSET = 40.0
_DEFAULT_Y_BASELINE = 200.0


@dataclass
class RenderDiagnostics:
    render_target: str
    param_count: int
    neutral_param_count: int
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


def _polyline_points(points: Sequence[Any]) -> str:
    return " ".join(
        f"{point.x * _DEFAULT_SCALE + _DEFAULT_X_OFFSET:.2f},{_DEFAULT_Y_BASELINE - point.y * _DEFAULT_SCALE:.2f}"
        for point in points
    )


def _build_svg(vocal_tract: VocalTract, width_px: int, height_px: int) -> str:
    upper_points = _polyline_points(vocal_tract.upper_outline.p)
    lower_points = _polyline_points(vocal_tract.lower_outline.p)
    tongue_points = _polyline_points(vocal_tract.tongue_outline.p)
    epiglottis_points = _polyline_points(vocal_tract.epiglottis_outline.p)
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


def _apply_params(vocal_tract: VocalTract, vector: Sequence[float]) -> None:
    for index, value in enumerate(vector):
        param = vocal_tract.params[index]
        bounded = min(max(float(value), float(param.min_val)), float(param.max_val))
        param.value = bounded
        param.limited_value = bounded
    vocal_tract.calculate_all()


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
    _apply_params(vocal_tract, vector)
    diagnostics = RenderDiagnostics(
        render_target=render_target,
        param_count=len(vector),
        neutral_param_count=sum(
            1 for index, value in enumerate(vector) if value == float(vocal_tract.params[index].neutral_val)
        ),
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
