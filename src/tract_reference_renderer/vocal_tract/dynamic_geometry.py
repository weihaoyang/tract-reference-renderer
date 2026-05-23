# -*- coding: utf-8 -*-
"""
动态声道几何计算

该模块实现 VocalTract 的动态更新主链路：
1. 由参考解剖表面组合得到当前动态表面
2. 由动态表面构建上下轮廓
3. 由轮廓求取中线与法线
4. 由中线切片计算横截面积与周长

当前实现属于项目自有运行时几何链，保持确定性与数值连续性。
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Iterable, List, Optional, Sequence, Tuple

from ..types import Articulator
from ..utils.geometry import Point2D, Point3D
from ..utils.splines import BezierCurve3D, LineStrip2D, LineStrip3D
from . import anatomy_geometry, params

if TYPE_CHECKING:
    from .geometry import Surface, VocalTract


def _clamp(value: float, min_value: float, max_value: float) -> float:
    return max(min_value, min(value, max_value))


def _get_circle_tangent(h: Point2D, c: Point2D, radius: float, clockwise: bool) -> float:
    p = h.x - c.x
    q = h.y - c.y
    radicand = p * p + q * q - radius * radius
    if radicand < 0.0:
        return 0.0
    n = math.sqrt(radicand)
    if clockwise:
        n = -n
    return math.atan2(p * n + q * radius, p * radius - q * n)


def _get_ellipse_tangent(h: Point2D, c: Point2D, a: float, b: float, clockwise: bool) -> float:
    r = c - h
    radicand = r.x * r.x * b * b - b * b * a * a + r.y * r.y * a * a
    if radicand < 0.0:
        radicand = 0.0
    root = math.sqrt(radicand)

    denominator = r.x * r.x * b * b + r.y * r.y * a * a
    if abs(denominator) < 1.0e-6:
        denominator = 1.0e-6

    alpha0 = math.atan2(
        -b * (r.y * a * a + r.x * root) / denominator,
        -a * (r.x * b * b - r.y * root) / denominator,
    )
    alpha1 = math.atan2(
        -b * (r.y * a * a - r.x * root) / denominator,
        -a * (r.x * b * b + r.y * root) / denominator,
    )

    cos0 = math.cos(alpha0)
    sin0 = math.sin(alpha0)
    v1 = c + Point2D(a * cos0, b * sin0) - h
    v2 = Point2D(-a * sin0, b * cos0)
    test = v1.x * v2.x + v1.y * v2.y

    if (test >= 0.0 and clockwise) or (test < 0.0 and not clockwise):
        return alpha1
    return alpha0


def _tongue_side_param_to_elevation_cm(param_value: float) -> float:
    if param_value >= 0.0:
        max_elevation_cm = 1.0
        elevation_cm = max_elevation_cm * (param_value / 0.3)
        return min(elevation_cm, max_elevation_cm)

    min_elevation_cm = -0.5
    elevation_cm = min_elevation_cm * (param_value / -0.15)
    return max(elevation_cm, min_elevation_cm)


def _limit_ellipse_pos(
    center: Point2D,
    rx: float,
    ry: float,
    border: LineStrip2D,
    anchor: Point2D,
) -> Point2D:
    epsilon = 1.0e-6
    num_points = border.get_num_points()
    if num_points <= 0:
        return Point2D(center.x, center.y)

    direction = center - anchor
    final_t = direction.magnitude()
    if final_t < epsilon:
        return Point2D(center.x, center.y)
    direction.normalize()

    rx2 = rx * rx
    ry2 = ry * ry
    for idx in range(num_points):
        p0 = border.get_control_point(idx)
        q = anchor - p0

        denominator = direction.x * direction.x * ry2 + direction.y * direction.y * rx2
        if abs(denominator) > epsilon:
            p_term = 2.0 * (q.x * direction.x * ry2 + q.y * direction.y * rx2) / denominator
            q_term = (q.x * q.x * ry2 + q.y * q.y * rx2 - rx2 * ry2) / denominator
            root = 0.25 * p_term * p_term - q_term
            if root >= 0.0:
                t = -0.5 * p_term - math.sqrt(root)
                if 0.0 <= t < final_t:
                    final_t = t

        if idx < num_points - 1:
            p1 = border.get_control_point(idx + 1)
            w = p1 - p0
            alpha = math.atan2(-ry * w.x, rx * w.y)
            u = Point2D(rx * math.cos(alpha), ry * math.sin(alpha))

            denominator = direction.x * w.y - direction.y * w.x
            if abs(denominator) > epsilon:
                q_plus = anchor - p0 + u
                s = (q_plus.y * direction.x - q_plus.x * direction.y) / denominator
                if 0.0 <= s <= 1.0:
                    t = (q_plus.y * w.x - q_plus.x * w.y) / denominator
                    if 0.0 <= t < final_t:
                        final_t = t

                q_minus = anchor - p0 - u
                s = (q_minus.y * direction.x - q_minus.x * direction.y) / denominator
                if 0.0 <= s <= 1.0:
                    t = (q_minus.y * w.x - q_minus.x * w.y) / denominator
                    if 0.0 <= t < final_t:
                        final_t = t

    return anchor + direction * final_t


def _get_param_value(vt: "VocalTract", index: int) -> float:
    param = vt.params[index]
    limited_value = float(param.limited_value)
    min_value = float(param.min_val)
    max_value = float(param.max_val)
    if limited_value < min_value or limited_value > max_value:
        limited_value = float(param.value)
    return _clamp(limited_value, min_value, max_value)


def _rebuild_cover_geometry_cpp(vt: "VocalTract") -> None:
    epsilon = 1.0e-6
    for param_def in vt.params:
        param_def.limited_value = _clamp(float(param_def.value), float(param_def.min_val), float(param_def.max_val))

    upper_cover = vt.surfaces_list[params.SurfaceIndex.UPPER_COVER]
    lower_cover = vt.surfaces_list[params.SurfaceIndex.LOWER_COVER]
    mandible = vt.surfaces_list[params.SurfaceIndex.MANDIBLE]
    palate = vt.surfaces_list[params.SurfaceIndex.PALATE]
    lower_teeth = vt.surfaces_list[params.SurfaceIndex.LOWER_TEETH]
    lower_teeth_original = vt.surfaces_list[params.SurfaceIndex.LOWER_TEETH_ORIGINAL]
    narrow_larynx_back = vt.surfaces_list[params.SurfaceIndex.NARROW_LARYNX_BACK]
    wide_larynx_back = vt.surfaces_list[params.SurfaceIndex.WIDE_LARYNX_BACK]
    narrow_larynx_front = vt.surfaces_list[params.SurfaceIndex.NARROW_LARYNX_FRONT]
    wide_larynx_front = vt.surfaces_list[params.SurfaceIndex.WIDE_LARYNX_FRONT]
    low_velum = vt.surfaces_list[params.SurfaceIndex.LOW_VELUM]
    mid_velum = vt.surfaces_list[params.SurfaceIndex.MID_VELUM]
    high_velum = vt.surfaces_list[params.SurfaceIndex.HIGH_VELUM]
    if (
        upper_cover is None
        or lower_cover is None
        or mandible is None
        or palate is None
        or lower_teeth is None
        or lower_teeth_original is None
        or narrow_larynx_back is None
        or wide_larynx_back is None
        or narrow_larynx_front is None
        or wide_larynx_front is None
        or low_velum is None
        or mid_velum is None
        or high_velum is None
    ):
        raise RuntimeError("[领域错误]: 声道主表面重建失败 - 关键参考表面缺失")

    jx = float(vt.params[int(params.ParamIndex.JX)].limited_value)
    ja = float(vt.params[int(params.ParamIndex.JA)].limited_value)
    hy = float(vt.params[int(params.ParamIndex.HY)].limited_value)
    hx = float(vt.params[int(params.ParamIndex.HX)].limited_value)
    angle_rad = ja * math.pi / 180.0
    cosine = math.cos(angle_rad)
    sine = math.sin(angle_rad)

    jaw_fulcrum = vt.anatomy.jaw_fulcrum
    jaw_rest_pos = vt.anatomy.jaw_rest_pos
    vertex = mandible.get_vertex(0, params.NUM_LOWER_COVER_POINTS - 1)
    dx = vertex.x + jaw_rest_pos.x + jx - jaw_fulcrum.x
    dy = vertex.y + jaw_rest_pos.y - jaw_fulcrum.y
    transformed_vertex = Point3D(
        cosine * dx - sine * dy + jaw_fulcrum.x,
        sine * dx + cosine * dy + jaw_fulcrum.y,
        vertex.z,
    )
    if hy > transformed_vertex.y:
        hy = transformed_vertex.y

    x_offset = anatomy_geometry.get_pharynx_back_x(vt, hy)
    p_last = narrow_larynx_front.get_vertex(params.NUM_LARYNX_RIBS - 1, params.NUM_LOWER_COVER_POINTS - 1)
    q_last = wide_larynx_front.get_vertex(params.NUM_LARYNX_RIBS - 1, params.NUM_LOWER_COVER_POINTS - 1)
    denominator = q_last.x - p_last.x
    if abs(denominator) < epsilon:
        denominator = epsilon
    max_hx = (transformed_vertex.x - x_offset - p_last.x) / denominator
    if hx > max_hx:
        hx = max_hx
    hx = _clamp(hx, float(vt.params[int(params.ParamIndex.HX)].min_val), float(vt.params[int(params.ParamIndex.HX)].max_val))

    vt.params[int(params.ParamIndex.HY)].limited_value = hy
    vt.params[int(params.ParamIndex.HX)].limited_value = hx

    rotation_rad = float(vt.anatomy.pharynx_rotation_angle_deg) * math.pi / 180.0
    sine_rotation = math.sin(rotation_rad)
    if abs(sine_rotation) < epsilon:
        sine_rotation = epsilon if sine_rotation >= 0.0 else -epsilon
    shear_coeff = math.cos(rotation_rad) / sine_rotation

    rib = 0
    for rib_idx in range(params.NUM_LARYNX_RIBS):
        for point_idx in range(params.NUM_UPPER_COVER_POINTS):
            narrow_vertex = narrow_larynx_back.get_vertex(rib_idx, point_idx)
            wide_vertex = wide_larynx_back.get_vertex(rib_idx, point_idx)
            mixed = (1.0 - hx) * narrow_vertex + hx * wide_vertex
            shear_x = mixed.y * shear_coeff
            mixed.x += x_offset + shear_x
            mixed.y += hy
            upper_cover.set_vertex(rib, point_idx, mixed)
        rib += 1

    pharynx_x_cm = [1.0, 0.69, 0.41, 0.19, 0.05, 0.0]
    pharynx_z_cm = [-1.0, -0.95, -0.81, -0.59, -0.31, 0.0]
    for rib_idx in range(params.NUM_PHARYNX_RIBS):
        ratio = float(rib_idx) / float(max(params.NUM_PHARYNX_RIBS - 1, 1))
        scale_z = 0.5 * (
            (1.0 - ratio) * float(vt.anatomy.pharynx_lower_depth_cm) + ratio * float(vt.anatomy.pharynx_upper_depth_cm)
        )
        scale_x = float(vt.anatomy.pharynx_back_width_cm)
        mix_ratio = float(rib_idx + 1) / float(params.NUM_PHARYNX_RIBS)
        y_value = (1.0 - mix_ratio) * hy + mix_ratio * float(vt.anatomy.pharynx_top_rib_y_cm)
        x_value = anatomy_geometry.get_pharynx_back_x(vt, y_value)
        for point_idx in range(params.NUM_UPPER_COVER_POINTS):
            upper_cover.set_vertex(
                rib,
                point_idx,
                Point3D(
                    x_value + scale_x * pharynx_x_cm[point_idx],
                    y_value,
                    scale_z * pharynx_z_cm[point_idx],
                ),
            )
        rib += 1

    vs = _get_param_value(vt, int(params.ParamIndex.VS))
    vo = _get_param_value(vt, int(params.ParamIndex.VO))
    vs_min = float(vt.params[int(params.ParamIndex.VS)].min_val)
    vs_max = float(vt.params[int(params.ParamIndex.VS)].max_val)
    vo_min = float(vt.params[int(params.ParamIndex.VO)].min_val)
    vo_max = float(vt.params[int(params.ParamIndex.VO)].max_val)
    s = 0.0 if abs(vs_max - vs_min) < epsilon else (vs - vs_min) / (vs_max - vs_min)
    t = 0.0 if abs(vo_max - vo_min) < epsilon else (vo - vo_min) / (vo_max - vo_min)
    t = _clamp(t, 0.0, 1.0)

    for rib_idx in range(params.NUM_VELUM_RIBS):
        for point_idx in range(params.NUM_UPPER_COVER_POINTS):
            high_vertex = high_velum.get_vertex(rib_idx, point_idx)
            mid_vertex = mid_velum.get_vertex(rib_idx, point_idx)
            low_vertex = low_velum.get_vertex(rib_idx, point_idx)
            mixed_closed = (1.0 - s) * high_vertex + s * mid_vertex
            mixed_open = (1.0 - t) * mixed_closed + t * low_vertex
            upper_cover.set_vertex(rib, point_idx, mixed_open)
        rib += 1

    last_velum_rib = rib - 1
    for rib_idx in range(params.NUM_JAW_RIBS):
        factor = float(vt.anatomy.palate_points[rib_idx].x) / 1.2
        for point_idx in range(params.NUM_UPPER_COVER_POINTS):
            palate_vertex = palate.get_vertex(rib_idx, point_idx)
            if factor < 1.0:
                velum_vertex = upper_cover.get_vertex(last_velum_rib, point_idx)
                merged_y = (1.0 - factor) * velum_vertex.y + factor * palate_vertex.y
                merged_z = velum_vertex.z if palate_vertex.x < 0.01 else palate_vertex.z
                upper_cover.set_vertex(rib, point_idx, Point3D(palate_vertex.x, merged_y, merged_z))
            else:
                upper_cover.set_vertex(rib, point_idx, Point3D(palate_vertex.x, palate_vertex.y, palate_vertex.z))
        rib += 1

    rib = 0
    for rib_idx in range(params.NUM_LARYNX_RIBS):
        for point_idx in range(params.NUM_LOWER_COVER_POINTS):
            narrow_vertex = narrow_larynx_front.get_vertex(rib_idx, point_idx)
            wide_vertex = wide_larynx_front.get_vertex(rib_idx, point_idx)
            mixed = (1.0 - hx) * narrow_vertex + hx * wide_vertex
            shear_x = mixed.y * shear_coeff
            mixed.x += x_offset + shear_x
            mixed.y += hy
            lower_cover.set_vertex(rib, point_idx, mixed)
        rib += 1

    rib += params.NUM_THROAT_RIBS
    for rib_idx in range(params.NUM_JAW_RIBS):
        for point_idx in range(params.NUM_LOWER_COVER_POINTS):
            mandible_vertex = mandible.get_vertex(rib_idx, point_idx)
            dx = mandible_vertex.x + jaw_rest_pos.x + jx - jaw_fulcrum.x
            dy = mandible_vertex.y + jaw_rest_pos.y - jaw_fulcrum.y
            lower_cover.set_vertex(
                rib,
                point_idx,
                Point3D(
                    cosine * dx - sine * dy + jaw_fulcrum.x,
                    sine * dx + cosine * dy + jaw_fulcrum.y,
                    mandible_vertex.z,
                ),
            )
        rib += 1

    for rib_idx in range(params.NUM_TEETH_RIBS):
        for point_idx in range(params.NUM_TEETH_POINTS):
            tooth_vertex = lower_teeth_original.get_vertex(rib_idx, point_idx)
            dx = tooth_vertex.x + jaw_rest_pos.x + jx - jaw_fulcrum.x
            dy = tooth_vertex.y + jaw_rest_pos.y - jaw_fulcrum.y
            lower_teeth.set_vertex(
                rib_idx,
                point_idx,
                Point3D(
                    cosine * dx - sine * dy + jaw_fulcrum.x,
                    sine * dx + cosine * dy + jaw_fulcrum.y,
                    tooth_vertex.z,
                ),
            )


def _copy_rib_by_ratio(dst: "Surface", dst_rib: int, src: "Surface", src_rib: int) -> None:
    if dst is None or src is None:
        return
    if dst.num_points <= 0 or src.num_points <= 0:
        return
    if src_rib < 0 or src_rib >= src.num_ribs:
        return
    if dst_rib < 0 or dst_rib >= dst.num_ribs:
        return

    if dst.num_points == 1:
        src_idx = 0
        src_point = src.get_vertex(src_rib, src_idx)
        dst.set_vertex(dst_rib, 0, Point3D(src_point.x, src_point.y, src_point.z))
        return

    for point_idx in range(dst.num_points):
        ratio = point_idx / (dst.num_points - 1)
        src_idx = int(round(ratio * (src.num_points - 1)))
        src_idx = _clamp(src_idx, 0, src.num_points - 1)
        source_point = src.get_vertex(src_rib, int(src_idx))
        dst.set_vertex(dst_rib, point_idx, Point3D(source_point.x, source_point.y, source_point.z))


def _fill_surface_from_sources(
    vt: "VocalTract",
    dst_surface_idx: int,
    source_specs: Sequence[Tuple[int, Sequence[int]]],
) -> None:
    dst = vt.surfaces_list[dst_surface_idx]
    if dst is None:
        return

    dst_rib = 0
    for source_surface_idx, ribs in source_specs:
        src = vt.surfaces_list[source_surface_idx]
        if src is None:
            continue
        for src_rib in ribs:
            if dst_rib >= dst.num_ribs:
                return
            safe_src_rib = int(_clamp(src_rib, 0, src.num_ribs - 1))
            _copy_rib_by_ratio(dst, dst_rib, src, safe_src_rib)
            dst_rib += 1

    if dst_rib == 0:
        return

    while dst_rib < dst.num_ribs:
        _copy_rib_by_ratio(dst, dst_rib, dst, dst_rib - 1)
        dst_rib += 1


def _mirror_point(point: Point3D) -> Point3D:
    return Point3D(point.x, point.y, -point.z)


def _set_twoside_surface_by_points(dst: "Surface", src: "Surface") -> None:
    if dst is None or src is None or dst.num_ribs != src.num_ribs or dst.num_points != src.num_points * 2 - 1:
        return
    for rib_idx in range(src.num_ribs):
        row = [src.get_vertex(rib_idx, point_idx) for point_idx in range(src.num_points)]
        mirrored = [_mirror_point(point) for point in row[-2::-1]]
        full_row = row + mirrored
        for point_idx, point in enumerate(full_row):
            dst.set_vertex(rib_idx, point_idx, point)


def _set_twoside_surface_by_ribs(dst: "Surface", src: "Surface") -> None:
    if dst is None or src is None or dst.num_points != src.num_points or dst.num_ribs != src.num_ribs * 2 - 1:
        return
    rows = [[src.get_vertex(rib_idx, point_idx) for point_idx in range(src.num_points)] for rib_idx in range(src.num_ribs)]
    mirrored_rows = [[_mirror_point(point) for point in row] for row in rows[-2::-1]]
    full_rows = rows + mirrored_rows
    for rib_idx, row in enumerate(full_rows):
        for point_idx, point in enumerate(row):
            dst.set_vertex(rib_idx, point_idx, point)


def _polyline_arc_lengths(points: Sequence[Point3D]) -> List[float]:
    if len(points) <= 0:
        return []
    arc = [0.0]
    total = 0.0
    for left, right in zip(points, points[1:]):
        total += (right - left).magnitude()
        arc.append(float(total))
    if total <= 1.0e-9:
        return [0.0 for _ in arc]
    return [float(value / total) for value in arc]


def _resample_point3d_polyline(points: Sequence[Point3D], count: int) -> List[Point3D]:
    if count <= 0:
        return []
    if len(points) <= 0:
        return [Point3D() for _ in range(count)]
    if len(points) == 1:
        point = points[0]
        return [Point3D(point.x, point.y, point.z) for _ in range(count)]
    params_arc = _polyline_arc_lengths(points)
    targets = [0.0] if count == 1 else [float(index / (count - 1)) for index in range(count)]
    sampled: List[Point3D] = []
    for target in targets:
        if target <= params_arc[0]:
            point = points[0]
            sampled.append(Point3D(point.x, point.y, point.z))
            continue
        if target >= params_arc[-1]:
            point = points[-1]
            sampled.append(Point3D(point.x, point.y, point.z))
            continue
        for index in range(len(params_arc) - 1):
            left = params_arc[index]
            right = params_arc[index + 1]
            if left <= target <= right:
                span = max(1.0e-9, right - left)
                ratio = float((target - left) / span)
                p0 = points[index]
                p1 = points[index + 1]
                sampled.append(
                    Point3D(
                        p0.x + (p1.x - p0.x) * ratio,
                        p0.y + (p1.y - p0.y) * ratio,
                        p0.z + (p1.z - p0.z) * ratio,
                    )
                )
                break
    return sampled


def _build_cover_side_surfaces(vt: "VocalTract") -> None:
    upper_cover = vt.surfaces_list[params.SurfaceIndex.UPPER_COVER]
    lower_cover = vt.surfaces_list[params.SurfaceIndex.LOWER_COVER]
    left_cover = vt.surfaces_list[params.SurfaceIndex.LEFT_COVER]
    right_cover = vt.surfaces_list[params.SurfaceIndex.RIGHT_COVER]
    if upper_cover is None or lower_cover is None or left_cover is None or right_cover is None:
        return

    upper_ribs = [5, 6, 7] + list(range(14, 23))
    lower_ribs = list(range(5, 17))
    upper_edge = [
        upper_cover.get_vertex(int(_clamp(rib_idx, 0, upper_cover.num_ribs - 1)), 0)
        for rib_idx in upper_ribs
    ]
    lower_edge = [
        lower_cover.get_vertex(int(_clamp(rib_idx, 0, lower_cover.num_ribs - 1)), 0)
        for rib_idx in lower_ribs
    ]
    sampled_upper = _resample_point3d_polyline(upper_edge, left_cover.num_ribs)
    sampled_lower = _resample_point3d_polyline(lower_edge, left_cover.num_ribs)
    for rib_idx in range(left_cover.num_ribs):
        lower_point = sampled_lower[rib_idx]
        upper_point = sampled_upper[rib_idx]
        bridge = [
            Point3D(lower_point.x, lower_point.y, lower_point.z),
            Point3D(
                0.67 * lower_point.x + 0.33 * upper_point.x,
                0.67 * lower_point.y + 0.33 * upper_point.y,
                0.67 * lower_point.z + 0.33 * upper_point.z,
            ),
            Point3D(
                0.33 * lower_point.x + 0.67 * upper_point.x,
                0.33 * lower_point.y + 0.67 * upper_point.y,
                0.33 * lower_point.z + 0.67 * upper_point.z,
            ),
            Point3D(upper_point.x, upper_point.y, upper_point.z),
        ]
        for point_idx, point in enumerate(bridge):
            left_cover.set_vertex(rib_idx, point_idx, point)
            right_cover.set_vertex(rib_idx, point_idx, _mirror_point(point))


def _copy_surface(dst: "Surface", src: "Surface") -> None:
    if dst is None or src is None:
        return
    for rib_idx in range(dst.num_ribs):
        mapped_rib = int(round(rib_idx * (src.num_ribs - 1) / max(dst.num_ribs - 1, 1)))
        _copy_rib_by_ratio(dst, rib_idx, src, mapped_rib)


def _apply_jaw_transform(vt: "VocalTract", surface_indices: Iterable[int]) -> None:
    jaw_angle_deg = _get_param_value(vt, params.ParamIndex.JA)
    jaw_angle_rad = jaw_angle_deg * math.pi / 180.0
    cosine = math.cos(jaw_angle_rad)
    sine = math.sin(jaw_angle_rad)
    jaw_shift_x = _get_param_value(vt, params.ParamIndex.JX)
    jaw_shift_y = float(vt.anatomy.jaw_rest_pos.y)
    fulcrum_x = float(vt.anatomy.jaw_fulcrum.x)
    fulcrum_y = float(vt.anatomy.jaw_fulcrum.y)

    for surface_idx in surface_indices:
        surface = vt.surfaces_list[surface_idx]
        if surface is None:
            continue
        for rib_idx in range(surface.num_ribs):
            for point_idx in range(surface.num_points):
                point = surface.get_vertex(rib_idx, point_idx)
                dx = point.x + jaw_shift_x - fulcrum_x
                dy = point.y + jaw_shift_y - fulcrum_y
                rotated_x = cosine * dx - sine * dy + fulcrum_x
                rotated_y = sine * dx + cosine * dy + fulcrum_y
                surface.set_vertex(rib_idx, point_idx, Point3D(rotated_x, rotated_y, point.z))


def _update_lower_gums_edges(vt: "VocalTract") -> None:
    jaw_angle_deg = _get_param_value(vt, params.ParamIndex.JA)
    jaw_angle_rad = jaw_angle_deg * math.pi / 180.0
    cosine = math.cos(jaw_angle_rad)
    sine = math.sin(jaw_angle_rad)
    jaw_shift_x = _get_param_value(vt, params.ParamIndex.JX)
    jaw_shift_y = float(vt.anatomy.jaw_rest_pos.y)
    fulcrum_x = float(vt.anatomy.jaw_fulcrum.x)
    fulcrum_y = float(vt.anatomy.jaw_fulcrum.y)

    num_ribs = min(
        len(vt.lower_gums_outer_edge_orig),
        len(vt.lower_gums_inner_edge_orig),
        len(vt.lower_gums_outer_edge),
        len(vt.lower_gums_inner_edge),
    )
    for rib_idx in range(num_ribs):
        outer_vertex = vt.lower_gums_outer_edge_orig[rib_idx]
        outer_dx = outer_vertex.x + jaw_shift_x - fulcrum_x
        outer_dy = outer_vertex.y + jaw_shift_y - fulcrum_y
        outer_x = cosine * outer_dx - sine * outer_dy + fulcrum_x
        outer_y = sine * outer_dx + cosine * outer_dy + fulcrum_y
        vt.lower_gums_outer_edge[rib_idx] = Point3D(outer_x, outer_y, outer_vertex.z)

        inner_vertex = vt.lower_gums_inner_edge_orig[rib_idx]
        inner_dx = inner_vertex.x + jaw_shift_x - fulcrum_x
        inner_dy = inner_vertex.y + jaw_shift_y - fulcrum_y
        inner_x = cosine * inner_dx - sine * inner_dy + fulcrum_x
        inner_y = sine * inner_dx + cosine * inner_dy + fulcrum_y
        vt.lower_gums_inner_edge[rib_idx] = Point3D(inner_x, inner_y, inner_vertex.z)


def _rebuild_throat_front(vt: "VocalTract") -> None:
    upper_cover = vt.surfaces_list[params.SurfaceIndex.UPPER_COVER]
    lower_cover = vt.surfaces_list[params.SurfaceIndex.LOWER_COVER]
    if upper_cover is None or lower_cover is None:
        raise RuntimeError("[领域错误]: 咽腔前壁重建失败 - 上下盖表面缺失")

    anchor = lower_cover.get_vertex(params.NUM_LARYNX_RIBS - 1, params.NUM_LOWER_COVER_POINTS - 1)
    control2 = Point3D(anchor.x, anchor.y, anchor.z)
    curve = BezierCurve3D()
    weights = [1.0, 0.71, 1.0]

    for throat_index in range(params.NUM_THROAT_RIBS):
        control0 = upper_cover.get_vertex(params.NUM_LARYNX_RIBS + throat_index, 0)
        control1 = Point3D(anchor.x, anchor.y, control0.z)
        curve.set_points([control0, control1, control2], weights)

        for point_idx in range(params.NUM_LOWER_COVER_POINTS):
            ratio = float(point_idx) / float(params.NUM_LOWER_COVER_POINTS - 1)
            lower_cover.set_vertex(params.NUM_LARYNX_RIBS + throat_index, point_idx, curve.get_point(ratio))


def _get_important_lip_points(vt: "VocalTract") -> Tuple[Point3D, Point3D, Point3D, Point3D, float]:
    upper_gums = vt.upper_gums_outer_edge
    lower_gums = vt.lower_gums_outer_edge
    narrow_path = vt.narrow_lip_corner_path
    wide_path = vt.wide_lip_corner_path

    if (
        len(upper_gums) < 9
        or len(lower_gums) < 9
        or narrow_path.num_points < 2
        or wide_path.num_points < 2
        or narrow_path.num_points != wide_path.num_points
    ):
        return Point3D(), Point3D(), Point3D(), Point3D(), 0.0

    angle_rad = 0.5 * _get_param_value(vt, params.ParamIndex.JA) * math.pi / 180.0
    cosine = math.cos(angle_rad)
    sine = math.sin(angle_rad)
    jaw_fulcrum = vt.anatomy.jaw_fulcrum
    jaw_rest_pos = vt.anatomy.jaw_rest_pos
    jx = _get_param_value(vt, params.ParamIndex.JX)

    lip_corner_path = vt.lip_corner_path
    lip_corner_path.reset(0)

    ld_param = vt.params[int(params.ParamIndex.LD)]
    ld_denom = max(1.0e-9, float(ld_param.max_val) - float(ld_param.min_val))
    ld_norm = (_get_param_value(vt, params.ParamIndex.LD) - float(ld_param.min_val)) / ld_denom
    ld_norm = _clamp(ld_norm, 0.0, 1.0)

    for index in range(narrow_path.num_points):
        q = (1.0 - ld_norm) * narrow_path.get_control_point(index) + ld_norm * wide_path.get_control_point(index)
        dx = q.x + float(jaw_rest_pos.x) + 0.5 * jx - float(jaw_fulcrum.x)
        dy = q.y + 0.5 * float(jaw_rest_pos.y) - float(jaw_fulcrum.y)
        r = Point3D(
            cosine * dx - sine * dy + float(jaw_fulcrum.x),
            sine * dx + cosine * dy + float(jaw_fulcrum.y),
            q.z,
        )
        lip_corner_path.add_point(r)

    lip_corner_separation_x = min(float(upper_gums[6].x), float(lower_gums[6].x))
    t0 = narrow_path.get_intersection(Point3D(lip_corner_separation_x, 0.0, 0.0), Point3D(1.0, 0.0, 0.0))
    t0 = _clamp(t0, 0.0, 1.0)

    lp_param = vt.params[int(params.ParamIndex.LP)]
    lp_denom = max(1.0e-9, float(lp_param.max_val) - float(lp_param.min_val))
    lp_norm = (_get_param_value(vt, params.ParamIndex.LP) - float(lp_param.min_val)) / lp_denom
    lp_norm = _clamp(lp_norm, 0.0, 1.0)

    corner = lip_corner_path.get_point(lp_norm)
    onset = lip_corner_path.get_point(t0) if lp_norm > t0 else Point3D(corner.x, corner.y, corner.z)

    t0_upper = upper_gums[8]
    t1_lower = lower_gums[8]
    x_offset = 0.1
    y_offset = 0.35 - lp_norm * 0.3
    l_max = 1.0
    l_min = 0.3
    l_value = (1.0 - lp_norm) * l_max + lp_norm * l_min

    f0 = Point3D()
    f1 = Point3D()
    f0.x = corner.x + l_value + 0.3
    f1.x = corner.x + l_value
    if f0.x < t0_upper.x + x_offset:
        f0.x = t0_upper.x + x_offset
    if f1.x < t1_lower.x + x_offset:
        f1.x = t1_lower.x + x_offset

    min_t = -0.05
    y_close = 0.5 * (t0_upper.y + t1_lower.y) + (1.0 - lp_norm) * y_offset
    lip_delta = 0.5 * _get_param_value(vt, params.ParamIndex.LD)
    if lip_delta < min_t:
        lip_delta = min_t

    f0.y = y_close + lip_delta - _get_param_value(vt, params.ParamIndex.JX)
    f1.y = y_close - lip_delta - _get_param_value(vt, params.ParamIndex.JX)
    if f0.y > t0_upper.y:
        f0.y = t0_upper.y
    if f1.y > t0_upper.y:
        f1.y = t0_upper.y
    if f0.y < t1_lower.y:
        f0.y = t1_lower.y
    if f1.y < t1_lower.y:
        f1.y = t1_lower.y

    lip_tangent = upper_gums[5] - upper_gums[0]
    tangent_x = lip_tangent.x if abs(lip_tangent.x) > 1.0e-9 else 1.0e-9
    f0.z = corner.z + (lip_tangent.z / tangent_x) * (f0.x - corner.x)
    f1.z = corner.z + (lip_tangent.z / tangent_x) * (f1.x - corner.x)
    return onset, corner, f0, f1, y_close


def _calc_radiation(vt: "VocalTract", lip_corner: Point3D) -> None:
    upper_lip = vt.surfaces_list[params.SurfaceIndex.UPPER_LIP]
    lower_lip = vt.surfaces_list[params.SurfaceIndex.LOWER_LIP]
    radiation = vt.surfaces_list[params.SurfaceIndex.RADIATION]
    if upper_lip is None or lower_lip is None or radiation is None:
        return

    lip_point = params.NUM_INNER_LIP_POINTS
    x_values = [0.0] * params.NUM_RADIATION_POINTS
    y_values = [0.0] * params.NUM_RADIATION_POINTS
    for point_index in range(params.NUM_RADIATION_POINTS):
        angle_rad = -0.5 * math.pi + math.pi * float(point_index) / float(params.NUM_RADIATION_POINTS - 1)
        x_values[point_index] = math.cos(angle_rad)
        y_values[point_index] = 0.5 * math.sin(angle_rad) + 0.5

    l = Point3D()
    direction = Point3D(0.0, 1.0, 0.0)
    normal = Point3D(0.0, 0.0, -1.0)
    min_z = 0.0
    h = 0.0

    for rib_index in range(params.NUM_RADIATION_RIBS):
        if rib_index < params.NUM_LIP_RIBS:
            normal = Point3D(0.0, 0.0, -1.0)
            u = upper_lip.get_vertex(rib_index, lip_point)
            l = lower_lip.get_vertex(rib_index, lip_point)
            mid = 0.5 * (u + l)
            direction = u - l
            h = u.y - l.y
            if h <= 0.0:
                h = 0.0
            if u.x <= lip_corner.x:
                h = 0.0
            else:
                new_z = mid.z - h
                if new_z > min_z:
                    new_z = min_z
                else:
                    min_z = new_z
                h = mid.z - min_z
        else:
            angle_rad = 0.5 * math.pi - 0.5 * math.pi * float(rib_index - params.NUM_LIP_RIBS + 1) / float(
                params.NUM_RADIATION_RIBS - params.NUM_LIP_RIBS
            )
            normal = Point3D(math.cos(angle_rad), 0.0, -math.sin(angle_rad))

        for point_index in range(params.NUM_RADIATION_POINTS):
            point = l + y_values[point_index] * direction + x_values[point_index] * h * normal
            radiation.set_vertex(rib_index, point_index, point)


def _build_lips_geometry(vt: "VocalTract") -> None:
    upper_lip = vt.surfaces_list[params.SurfaceIndex.UPPER_LIP]
    lower_lip = vt.surfaces_list[params.SurfaceIndex.LOWER_LIP]
    if upper_lip is None or lower_lip is None:
        return
    onset, corner, f0, f1, _y_close = _get_important_lip_points(vt)

    if len(vt.upper_gums_outer_edge) < params.NUM_JAW_RIBS or len(vt.lower_gums_outer_edge) < params.NUM_JAW_RIBS:
        return

    epsilon = 1.0e-6
    lip_radius = max(0.4, float(vt.anatomy.lips_width_cm) - 0.4)
    curve = BezierCurve3D()
    upper_path = LineStrip3D()
    lower_path = LineStrip3D()

    upper_path.reset(0)
    lower_path.reset(0)
    n_path = vt.lip_corner_path.num_points
    for index in range(max(0, n_path - 1)):
        p0 = vt.lip_corner_path.get_control_point(index)
        p1 = vt.lip_corner_path.get_control_point(index + 1)
        length = p1.x - p0.x
        if length < epsilon:
            length = epsilon

        if p0.x <= onset.x <= p1.x:
            upper_path.add_point(Point3D(onset.x, onset.y, onset.z))
            lower_path.add_point(Point3D(onset.x, onset.y, onset.z))
        if onset.x < p0.x < corner.x:
            upper_path.add_point(Point3D(p0.x, p0.y, p0.z))
            lower_path.add_point(Point3D(p0.x, p0.y, p0.z))
        if p0.x <= corner.x <= p1.x:
            upper_path.add_point(Point3D(corner.x, corner.y, corner.z))
            lower_path.add_point(Point3D(corner.x, corner.y, corner.z))

    curve.set_points([corner, f0, Point3D(f0.x, f0.y, 0.0)], [1.0, 1.5, 1.0])
    for index in range(1, 8):
        upper_path.add_point(curve.get_point(float(index) / 7.0))

    curve.set_points([corner, f1, Point3D(f1.x, f1.y, 0.0)], [1.0, 1.5, 1.0])
    for index in range(1, 8):
        lower_path.add_point(curve.get_point(float(index) / 7.0))

    def _build_lip_ribs(gums_outer: List[Point3D], gums_inner: List[Point3D], path: LineStrip3D, is_upper: bool) -> None:
        outer_origin = [Point3D() for _ in range(params.NUM_LIP_RIBS)]
        inner_origin = [Point3D() for _ in range(params.NUM_LIP_RIBS)]
        temp_outer = [Point3D() for _ in range(params.NUM_LIP_RIBS)]
        temp_inner = [Point3D() for _ in range(params.NUM_LIP_RIBS)]
        num_divisions = [1 for _ in range(params.NUM_LIP_RIBS)]

        num_ribs = 0
        for jaw_index in range(1, params.NUM_JAW_RIBS):
            p0 = gums_outer[jaw_index - 1]
            p1 = gums_outer[jaw_index]
            q0 = gums_inner[jaw_index - 1]
            q1 = gums_inner[jaw_index]
            length = p1.x - p0.x
            if length < epsilon:
                length = epsilon

            if p0.x <= onset.x <= p1.x:
                t = (onset.x - p0.x) / length
                outer_origin[num_ribs] = p0 + t * (p1 - p0)
                inner_origin[num_ribs] = q0 + t * (q1 - q0)
                num_ribs += 1
            if p0.x <= corner.x <= p1.x:
                t = (corner.x - p0.x) / length
                outer_origin[num_ribs] = p0 + t * (p1 - p0)
                inner_origin[num_ribs] = q0 + t * (q1 - q0)
                num_ribs += 1
            if p1.x > onset.x:
                outer_origin[num_ribs] = p1
                inner_origin[num_ribs] = q1
                num_ribs += 1

        if num_ribs < 2:
            return
        num_fixed = num_ribs
        while num_ribs < params.NUM_LIP_RIBS:
            best_segment = -1
            max_length = -1_000_000.0
            for index in range(num_fixed - 1):
                seg = (outer_origin[index + 1].z - outer_origin[index].z) / float(num_divisions[index])
                if seg > max_length:
                    max_length = seg
                    best_segment = index
            if best_segment < 0:
                break
            num_divisions[best_segment] += 1
            num_ribs += 1

        num_ribs = 0
        for fixed_index in range(num_fixed):
            if fixed_index < num_fixed - 1:
                p0 = outer_origin[fixed_index]
                p1 = outer_origin[fixed_index + 1]
                q0 = inner_origin[fixed_index]
                q1 = inner_origin[fixed_index + 1]
                for split in range(num_divisions[fixed_index]):
                    t = float(split) / float(num_divisions[fixed_index])
                    temp_outer[num_ribs] = p0 + (p1 - p0) * t
                    temp_inner[num_ribs] = q0 + (q1 - q0) * t
                    num_ribs += 1
            else:
                temp_outer[num_ribs] = outer_origin[fixed_index]
                temp_inner[num_ribs] = inner_origin[fixed_index]
                num_ribs += 1

        for index in range(num_ribs):
            outer_origin[index] = temp_outer[index]
            inner_origin[index] = temp_inner[index]

        max_x = path.get_point(1.0).x
        denom_outer = max(epsilon, outer_origin[params.NUM_LIP_RIBS - 1].x - outer_origin[0].x)

        target_surface = upper_lip if is_upper else lower_lip
        for rib_index in range(num_ribs):
            p0 = outer_origin[rib_index]
            t = float(rib_index) / float(max(1, num_ribs - 1))
            p1 = path.get_point(t)

            c0 = p0
            c2 = p1
            q = Point3D(p1.x, p0.y, p1.z)
            q = 0.3 * q + 0.7 * p1
            r = Point3D(p0.x, p1.y, p0.z)
            c1_mix = (p0.x - outer_origin[0].x) / denom_outer
            c1_mix = c1_mix * c1_mix * c1_mix
            c1 = (1.0 - c1_mix) * q + c1_mix * r
            curve.set_points([c0, c1, c2], [1.0, 2.0, 1.0])

            target_surface.set_vertex(rib_index, 0, inner_origin[rib_index])
            for point_index in range(params.NUM_INNER_LIP_POINTS - 1):
                t_inner = float(point_index) / float(max(1, params.NUM_INNER_LIP_POINTS - 2))
                target_surface.set_vertex(rib_index, point_index + 1, curve.get_point(t_inner))

            radius = lip_radius * (p1.x - corner.x) / max(epsilon, max_x - corner.x)
            if radius < 0.0:
                radius = 0.0
            q_outer = Point3D(p1.x, p1.y + radius if is_upper else p1.y - radius, p1.z)
            for point_index in range(params.NUM_OUTER_LIP_POINTS):
                if is_upper:
                    angle = -0.5 * math.pi + 0.5 * math.pi * float(point_index + 1) / float(params.NUM_OUTER_LIP_POINTS)
                else:
                    angle = 0.5 * math.pi - 0.5 * math.pi * float(point_index + 1) / float(params.NUM_OUTER_LIP_POINTS)
                point = Point3D(
                    q_outer.x + radius * math.cos(angle),
                    q_outer.y + radius * math.sin(angle),
                    q_outer.z,
                )
                target_surface.set_vertex(rib_index, params.NUM_INNER_LIP_POINTS + point_index, point)

    _build_lip_ribs(vt.upper_gums_outer_edge, vt.upper_gums_inner_edge, upper_path, is_upper=True)
    _build_lip_ribs(vt.lower_gums_outer_edge, vt.lower_gums_inner_edge, lower_path, is_upper=False)
    _calc_radiation(vt, corner)


def _build_tongue_geometry(vt: "VocalTract") -> None:
    tongue_surface = vt.surfaces_list[params.SurfaceIndex.TONGUE]
    lower_cover = vt.surfaces_list[params.SurfaceIndex.LOWER_COVER]
    if tongue_surface is None or lower_cover is None:
        return

    tcx = _get_param_value(vt, int(params.ParamIndex.TCX))
    tcy = _get_param_value(vt, int(params.ParamIndex.TCY))
    ttx = _get_param_value(vt, int(params.ParamIndex.TTX))
    tty = _get_param_value(vt, int(params.ParamIndex.TTY))
    tbx = _get_param_value(vt, int(params.ParamIndex.TBX))
    tby = _get_param_value(vt, int(params.ParamIndex.TBY))
    r0x = float(vt.anatomy.tongue_center_radius_x_cm)
    r0y = float(vt.anatomy.tongue_center_radius_y_cm)
    r1 = float(vt.anatomy.tongue_tip_radius_cm)

    upper_cover = vt.surfaces_list[params.SurfaceIndex.UPPER_COVER]
    upper_teeth = vt.surfaces_list[params.SurfaceIndex.UPPER_TEETH]
    lower_teeth = vt.surfaces_list[params.SurfaceIndex.LOWER_TEETH]
    if upper_cover is not None and upper_teeth is not None and lower_teeth is not None:
        delta = 0.3
        q0_anchor = upper_cover.get_vertex(
            params.NUM_UPPER_COVER_RIBS - params.NUM_PALATE_RIBS + 1,
            params.NUM_UPPER_COVER_POINTS - 1,
        ).to_point2d()
        q1_anchor = upper_cover.get_vertex(
            params.NUM_UPPER_COVER_RIBS - 1,
            params.NUM_UPPER_COVER_POINTS - 1,
        ).to_point2d()
        upper_anchor = Point2D(q0_anchor.x, q1_anchor.y - (q1_anchor.x - q0_anchor.x))

        upper_border = LineStrip2D()
        upper_border.reset(0)
        for rib_idx in range(params.NUM_LARYNX_RIBS, params.NUM_UPPER_COVER_RIBS):
            upper_border.add_point(upper_cover.get_vertex(rib_idx, params.NUM_UPPER_COVER_POINTS - 1).to_point2d())
        for point_idx in range(2):
            upper_border.add_point(upper_teeth.get_vertex(params.NUM_TEETH_RIBS - 1, point_idx).to_point2d())
        upper_tail = upper_border.get_control_point(upper_border.get_num_points() - 1)
        upper_border.add_point(Point2D(upper_tail.x, upper_tail.y - 10.0))

        for idx in range(upper_border.get_num_points()):
            border_point = upper_border.get_control_point(idx)
            radial = border_point - upper_anchor
            if radial.magnitude() > 1.0e-12:
                radial.normalize()
                border_point = border_point + delta * radial
            if border_point.x > upper_tail.x:
                border_point.x = upper_tail.x
            upper_border.set_point(idx, border_point)

        lower_border_tt = LineStrip2D()
        lower_border_tt.reset(1)
        for rib_idx in range(params.NUM_LARYNX_RIBS - 1, params.NUM_LOWER_COVER_RIBS):
            lower_border_tt.add_point(lower_cover.get_vertex(rib_idx, params.NUM_LOWER_COVER_POINTS - 1).to_point2d())
        for point_idx in range(2):
            lower_border_tt.add_point(lower_teeth.get_vertex(params.NUM_TEETH_RIBS - 1, point_idx).to_point2d())
        lower_tt_last = lower_border_tt.get_control_point(lower_border_tt.get_num_points() - 1)
        upper_gum = vt.upper_gums_outer_edge[8].to_point2d()
        if upper_gum.x < lower_tt_last.x + 0.1:
            upper_gum.x = lower_tt_last.x + 0.1
        lower_border_tt.add_point(upper_gum)
        lower_tt_shift = lower_border_tt.get_control_point(1)
        lower_tt_shift.x -= 10.0
        lower_border_tt.set_point(0, lower_tt_shift)

        lower_border_tb = LineStrip2D()
        lower_border_tb.reset(1)
        p0 = lower_cover.get_vertex(params.NUM_LARYNX_RIBS - 1, params.NUM_LOWER_COVER_POINTS - 1).to_point2d()
        p2 = lower_cover.get_vertex(params.NUM_LOWER_COVER_RIBS - 1, params.NUM_LOWER_COVER_POINTS - 1).to_point2d()
        p1 = Point2D(p2.x - 1.0, p0.y)
        if p1.x < p0.x:
            p1.x = p0.x
        lower_border_tb.add_point(p0)
        lower_border_tb.add_point(p1)
        lower_border_tb.add_point(p2)
        for point_idx in range(3):
            lower_border_tb.add_point(lower_teeth.get_vertex(params.NUM_TEETH_RIBS - 1, point_idx).to_point2d())
        lower_tb_shift = lower_border_tb.get_control_point(1)
        lower_tb_shift.x -= 10.0
        lower_border_tb.set_point(0, lower_tb_shift)

        if ttx < tcx + r0x + r1:
            ttx = tcx + r0x + r1

        tip_limited = _limit_ellipse_pos(Point2D(ttx, tty), r1, r1, upper_border, upper_anchor)
        ttx, tty = tip_limited.x, tip_limited.y

        if tcx > ttx - r0x - r1:
            tcx = ttx - r0x - r1

        body_limited = _limit_ellipse_pos(Point2D(tcx, tcy), r0x, r0y, upper_border, upper_anchor)
        tcx, tcy = body_limited.x, body_limited.y

        if ttx < tcx + r0x + r1:
            ttx = tcx + r0x + r1

        lower_anchor = Point2D(upper_anchor.x, 1000.0)
        tip_limited = _limit_ellipse_pos(Point2D(ttx, tty), r1, r1, lower_border_tt, lower_anchor)
        ttx, tty = tip_limited.x, tip_limited.y
        body_limited = _limit_ellipse_pos(Point2D(tcx, tcy), r0x, r0y, lower_border_tb, lower_anchor)
        tcx, tcy = body_limited.x, body_limited.y

    if bool(vt.anatomy.automatic_tongue_root_calc):
        hyoid = lower_cover.get_vertex(params.NUM_LARYNX_RIBS - 1, params.NUM_LOWER_COVER_POINTS - 1).to_point2d()
        c0_for_root = Point2D(tcx, tcy)
        distance = (c0_for_root - hyoid).magnitude()
        trx = vt.anatomy.tongue_root_trx_slope * distance + vt.anatomy.tongue_root_trx_intercept
        try_ = vt.anatomy.tongue_root_try_slope * tcx + vt.anatomy.tongue_root_try_intercept
    else:
        trx = _get_param_value(vt, int(params.ParamIndex.TRX))
        try_ = _get_param_value(vt, int(params.ParamIndex.TRY))

    # C++ restrictTongueParams() 中后半段约束（tongue root / blade）。
    tongue_radius_max = max(r0x, r0y)
    alpha_rad = math.atan2(r0y, r0x) + math.pi
    c_ref = Point2D(
        tcx + math.cos(alpha_rad) * (r0x + 1.0e-6),
        tcy + math.sin(alpha_rad) * (r0y + 1.0e-6),
    )
    hyoid = lower_cover.get_vertex(params.NUM_LARYNX_RIBS - 1, params.NUM_LOWER_COVER_POINTS - 1).to_point2d()

    if try_ < hyoid.y:
        try_ = hyoid.y
    if try_ > tcy:
        try_ = tcy

    if trx - c_ref.x > c_ref.y - try_:
        trx = c_ref.x + c_ref.y - try_

    min_x = tcx - r0x - (tcy - hyoid.y)
    if trx < min_x:
        trx = min_x

    max_x = tcx - tongue_radius_max
    if hyoid.x + 0.5 > max_x:
        max_x = hyoid.x + 0.5
    if trx > max_x:
        trx = max_x

    if tbx < tcx:
        tbx = tcx
    if tbx > ttx:
        tbx = ttx

    radius_ref = tongue_radius_max + 1.0e-6
    c_tb = Point2D(tcx, tcy + 1.415 * radius_ref)
    if tby - c_tb.y > tbx - c_tb.x:
        tby = c_tb.y + tbx - c_tb.x
    if tby - c_tb.y < -(tbx - c_tb.x):
        tby = c_tb.y - (tbx - c_tb.x)

    tip_radius_ref = r1 + 1.0e-6
    c_tt = Point2D(ttx, tty + 1.415 * tip_radius_ref)
    if tby - c_tt.y < tbx - c_tt.x:
        tby = c_tt.y + tbx - c_tt.x

    ts1 = _get_param_value(vt, int(params.ParamIndex.TS1))
    ts2 = _get_param_value(vt, int(params.ParamIndex.TS2))
    ts3 = _get_param_value(vt, int(params.ParamIndex.TS3))

    tongue_indices_and_values = (
        (int(params.ParamIndex.TCX), float(tcx)),
        (int(params.ParamIndex.TCY), float(tcy)),
        (int(params.ParamIndex.TTX), float(ttx)),
        (int(params.ParamIndex.TTY), float(tty)),
        (int(params.ParamIndex.TBX), float(tbx)),
        (int(params.ParamIndex.TBY), float(tby)),
        (int(params.ParamIndex.TRX), float(trx)),
        (int(params.ParamIndex.TRY), float(try_)),
        (int(params.ParamIndex.TS1), float(ts1)),
        (int(params.ParamIndex.TS2), float(ts2)),
        (int(params.ParamIndex.TS3), float(ts3)),
    )
    for param_index, param_value in tongue_indices_and_values:
        param_def = vt.params[param_index]
        param_def.limited_value = _clamp(param_value, float(param_def.min_val), float(param_def.max_val))

    c0 = Point2D(tcx, tcy)
    c1 = Point2D(ttx, tty)
    epsilon = 1.0e-6

    hyoid_point = lower_cover.get_vertex(params.NUM_LARYNX_RIBS - 1, params.NUM_LOWER_COVER_POINTS - 1)

    root_curve = BezierCurve3D()
    blade_curve = BezierCurve3D()
    weights = [1.0, 2.0, 1.0]

    q0 = Point3D(hyoid_point.x, hyoid_point.y, 0.0)
    q1 = Point3D(trx, try_, 0.0)
    alpha0 = _get_ellipse_tangent(q1.to_point2d(), c0, r0x, r0y, True)
    q2 = Point3D(c0.x + r0x * math.cos(alpha0), c0.y + r0y * math.sin(alpha0), 0.0)
    root_curve.set_points([q0, q1, q2], weights)

    blade_mid = Point3D(tbx, tby, 0.0)
    alpha1 = _get_ellipse_tangent(blade_mid.to_point2d(), c0, r0x, r0y, False)
    if (alpha1 > alpha0) and (alpha1 - alpha0 < 0.25 * math.pi):
        alpha1 = alpha0
    blade_start = Point3D(c0.x + r0x * math.cos(alpha1), c0.y + r0y * math.sin(alpha1), 0.0)

    alpha2 = _get_circle_tangent(blade_mid.to_point2d(), c1, r1, True)
    blade_end = Point3D(c1.x + r1 * math.cos(alpha2), c1.y + r1 * math.sin(alpha2), 0.0)
    blade_curve.set_points([blade_start, blade_mid, blade_end], weights)

    target_curve = LineStrip2D()
    target_curve.reset(0)
    tongue_side_index = [0] * 5
    tongue_side_index[0] = 0

    for i in range(32):
        target_curve.add_point(root_curve.get_point(float(i) / 32.0).to_point2d())
    tongue_side_index[1] = target_curve.get_num_points() - 8

    alpha0_unwrapped = alpha0 + (2.0 * math.pi if alpha0 < 0.0 else 0.0)
    delta = alpha1 - alpha0_unwrapped
    if delta > -epsilon:
        delta = -epsilon
    for i in range(32):
        angle = alpha0_unwrapped + delta * float(i) / 32.0
        target_curve.add_point(Point2D(c0.x + r0x * math.cos(angle), c0.y + r0y * math.sin(angle)))
    tongue_side_index[2] = target_curve.get_num_points() - 8

    for i in range(32):
        target_curve.add_point(blade_curve.get_point(float(i) / 32.0).to_point2d())
    tongue_side_index[3] = target_curve.get_num_points() - 1

    alpha2_unwrapped = alpha2 + (2.0 * math.pi if alpha2 < 0.0 else 0.0)
    delta = 0.0 - alpha2_unwrapped
    if delta > -epsilon:
        delta = -epsilon
    for i in range(8):
        angle = alpha2_unwrapped + delta * float(i) / 7.0
        target_curve.add_point(Point2D(c1.x + r1 * math.cos(angle), c1.y + r1 * math.sin(angle)))
    tongue_side_index[4] = target_curve.get_num_points() - 1

    tongue_side_elevation_cm = [
        1.0,
        _tongue_side_param_to_elevation_cm(ts1),
        _tongue_side_param_to_elevation_cm(ts2),
        _tongue_side_param_to_elevation_cm(ts3),
        _tongue_side_param_to_elevation_cm(ts3),
    ]
    if tongue_side_elevation_cm[4] > 0.0:
        tongue_side_elevation_cm[4] = 0.0

    left_side_height = LineStrip2D()
    left_side_height.reset(0)
    right_side_height = LineStrip2D()
    right_side_height.reset(0)
    for i in range(5):
        t = target_curve.get_curve_param(tongue_side_index[i])
        sample = Point2D(t, tongue_side_elevation_cm[i])
        left_side_height.add_point(sample)
        right_side_height.add_point(sample)

    dynamic_ribs = params.NUM_DYNAMIC_TONGUE_RIBS
    rib_points: List[Point2D] = [Point2D() for _ in range(dynamic_ribs)]
    rib_normals: List[Point2D] = [Point2D(0.0, 1.0) for _ in range(dynamic_ribs)]
    rib_left_height = [0.0 for _ in range(dynamic_ribs)]
    rib_right_height = [0.0 for _ in range(dynamic_ribs)]

    for i in range(dynamic_ribs):
        t = float(i) / float(max(dynamic_ribs - 1, 1))
        rib_points[i] = target_curve.get_point(t)
        rib_left_height[i] = left_side_height.get_function_value(t)
        rib_right_height[i] = right_side_height.get_function_value(t)

    filter_coeff = 0.5
    for i in range(1, dynamic_ribs):
        rib_left_height[i] += filter_coeff * (rib_left_height[i - 1] - rib_left_height[i])
        rib_right_height[i] += filter_coeff * (rib_right_height[i - 1] - rib_right_height[i])

    sector = 0
    for i in range(dynamic_ribs):
        point = rib_points[i]
        if sector == 0 and point.y > c0.y:
            sector = 1
        if sector == 1 and point.x > c0.x:
            sector = 2

        if sector == 0:
            rib_normals[i] = Point2D(-1.0, 0.0)
        elif sector == 1:
            normal = point - c0
            if normal.magnitude() < 1.0e-12:
                normal = Point2D(-1.0, 0.0)
            else:
                normal.normalize()
            rib_normals[i] = normal
        else:
            rib_normals[i] = Point2D(0.0, 1.0)

    num_tongue_ribs = params.NUM_TONGUE_RIBS
    tongue_ribs: List[params.TongueRib] = [params.TongueRib() for _ in range(num_tongue_ribs)]
    for rib_idx in range(dynamic_ribs):
        tongue_ribs[rib_idx].point = Point2D(rib_points[rib_idx].x, rib_points[rib_idx].y)
        tongue_ribs[rib_idx].left_side_height = float(rib_left_height[rib_idx])
        tongue_ribs[rib_idx].right_side_height = float(rib_right_height[rib_idx])
        tongue_ribs[rib_idx].normal = Point2D(rib_normals[rib_idx].x, rib_normals[rib_idx].y)

    lower_boundary = LineStrip2D()
    lower_boundary.reset(0)
    boundary_start = lower_cover.get_vertex(
        params.NUM_LARYNX_RIBS + params.NUM_THROAT_RIBS,
        params.NUM_LOWER_COVER_POINTS - 1,
    ).to_point2d()
    boundary_start.y += 0.5
    lower_boundary.add_point(boundary_start)
    found, _, boundary_hit = vt.lower_outline.get_closest_intersection(boundary_start, Point2D(1.0, 0.0))
    if found:
        lower_boundary.add_point(boundary_hit)
    else:
        lower_boundary.add_point(boundary_start)

    min_boundary_y = lower_boundary.get_control_point(1).y
    for rib_idx in range(params.NUM_LARYNX_RIBS + params.NUM_THROAT_RIBS, params.NUM_LOWER_COVER_RIBS):
        point = lower_cover.get_vertex(rib_idx, params.NUM_LOWER_COVER_POINTS - 1).to_point2d()
        if point.y > min_boundary_y:
            lower_boundary.add_point(point)
    for point_idx in range(3):
        lower_boundary.add_point(lower_teeth.get_vertex(params.NUM_TEETH_RIBS - 1, point_idx).to_point2d())
    boundary_tail = lower_boundary.get_control_point(lower_boundary.get_num_points() - 1)
    boundary_tail.x += 10.0
    lower_boundary.add_point(boundary_tail)

    final_dynamic_idx = dynamic_ribs - 1
    min_side = tongue_ribs[final_dynamic_idx].left_side_height
    if tongue_ribs[final_dynamic_idx].right_side_height < min_side:
        min_side = tongue_ribs[final_dynamic_idx].right_side_height
    final_point = tongue_ribs[final_dynamic_idx].point + tongue_ribs[final_dynamic_idx].normal * min_side
    if tty - r1 < final_point.y:
        final_point.y = tty - r1

    first_static_idx = dynamic_ribs
    tongue_ribs[first_static_idx].point = Point2D(final_point.x, final_point.y)
    tongue_ribs[first_static_idx].normal = Point2D(0.0, 1.0)
    tongue_ribs[first_static_idx].left_side_height = 0.0
    tongue_ribs[first_static_idx].right_side_height = 0.0

    upward = Point2D(0.0, 1.0)
    for rib_idx in range(num_tongue_ribs // 2, dynamic_ribs + 1):
        found, hit_t, hit_point = lower_boundary.get_closest_intersection(tongue_ribs[rib_idx].point, upward)
        if found and hit_t > 0.0:
            tongue_ribs[rib_idx].point = Point2D(hit_point.x, hit_point.y)

        found, hit_t, hit_point = lower_boundary.get_closest_intersection(
            tongue_ribs[rib_idx].point,
            tongue_ribs[rib_idx].normal,
        )
        if found:
            if tongue_ribs[rib_idx].left_side_height < hit_t:
                tongue_ribs[rib_idx].left_side_height = hit_t
            if tongue_ribs[rib_idx].right_side_height < hit_t:
                tongue_ribs[rib_idx].right_side_height = hit_t

    pivot = Point2D(
        tongue_ribs[num_tongue_ribs - 4].point.x,
        tongue_ribs[num_tongue_ribs - 4].point.y,
    )
    anterior = Point2D(pivot.x - 2.0 * r1, pivot.y - 2.0 * r1)
    found, hit_t, hit_point = lower_boundary.get_closest_intersection(pivot, anterior - pivot)
    if found and hit_t < 1.0:
        anterior = Point2D(hit_point.x, hit_point.y)

    tongue_ribs[num_tongue_ribs - 3].point = Point2D(anterior.x, anterior.y)
    tongue_ribs[num_tongue_ribs - 3].normal = Point2D(0.0, 1.0)
    tongue_ribs[num_tongue_ribs - 3].left_side_height = 0.0
    tongue_ribs[num_tongue_ribs - 3].right_side_height = 0.0

    tail_source = tongue_ribs[num_tongue_ribs - 3].point
    found, _, hit_point = lower_boundary.get_closest_intersection(tail_source, Point2D(0.0, -1.0))
    if found:
        tongue_ribs[num_tongue_ribs - 2].point = Point2D(hit_point.x, hit_point.y)
    else:
        tongue_ribs[num_tongue_ribs - 2].point = Point2D(tail_source.x, tail_source.y)

    if tongue_ribs[num_tongue_ribs - 2].point.x < lower_boundary.get_control_point(1).x:
        tongue_ribs[num_tongue_ribs - 1].point = Point2D(
            lower_boundary.get_control_point(1).x,
            lower_boundary.get_control_point(1).y,
        )
    else:
        tongue_ribs[num_tongue_ribs - 1].point = Point2D(
            tongue_ribs[num_tongue_ribs - 2].point.x,
            tongue_ribs[num_tongue_ribs - 2].point.y,
        )

    for rib_idx in (num_tongue_ribs - 2, num_tongue_ribs - 1):
        tongue_ribs[rib_idx].normal = Point2D(0.0, 1.0)
        tongue_ribs[rib_idx].left_side_height = 0.0
        tongue_ribs[rib_idx].right_side_height = 0.0

    if tongue_surface.num_points != 11 or lower_cover.num_points != 5:
        raise RuntimeError("[领域错误]: calcTongueRibs - 网格点数量与运行时预期不一致")

    base_rib = params.NUM_LARYNX_RIBS - 1
    point_mid = lower_cover.get_vertex(base_rib, 2)
    tongue_surface.set_vertex(0, 0, Point3D(point_mid.x, point_mid.y, point_mid.z))
    tongue_surface.set_vertex(0, 1, Point3D(point_mid.x, point_mid.y, point_mid.z))
    tongue_surface.set_vertex(0, 9, Point3D(point_mid.x, point_mid.y, -point_mid.z))
    tongue_surface.set_vertex(0, 10, Point3D(point_mid.x, point_mid.y, -point_mid.z))
    point_side = lower_cover.get_vertex(base_rib, 3)
    tongue_surface.set_vertex(0, 2, Point3D(point_side.x, point_side.y, point_side.z))
    tongue_surface.set_vertex(0, 3, Point3D(point_side.x, point_side.y, point_side.z))
    tongue_surface.set_vertex(0, 7, Point3D(point_side.x, point_side.y, -point_side.z))
    tongue_surface.set_vertex(0, 8, Point3D(point_side.x, point_side.y, -point_side.z))
    point_tip = lower_cover.get_vertex(base_rib, 4)
    tongue_surface.set_vertex(0, 4, Point3D(point_tip.x, point_tip.y, point_tip.z))
    tongue_surface.set_vertex(0, 5, Point3D(point_tip.x, point_tip.y, point_tip.z))
    tongue_surface.set_vertex(0, 6, Point3D(point_tip.x, point_tip.y, point_tip.z))
    tongue_ribs[0].min_x = tongue_surface.get_vertex(0, 0).z
    tongue_ribs[0].max_x = tongue_surface.get_vertex(0, tongue_surface.num_points - 1).z

    invalid = params.INVALID_PROFILE_SAMPLE
    standard_half_width = 1.75
    min_half_width = 0.75
    min_half_width_samples = int(min_half_width / params.PROFILE_SAMPLE_LENGTH)
    half_profile = params.NUM_PROFILE_SAMPLES // 2
    a2_left = [0.0] * num_tongue_ribs
    a3_left = [0.0] * num_tongue_ribs
    a2_right = [0.0] * num_tongue_ribs
    a3_right = [0.0] * num_tongue_ribs
    upper_profiles = [[invalid for _ in range(params.NUM_PROFILE_SAMPLES)] for _ in range(num_tongue_ribs)]
    lower_profiles = [[invalid for _ in range(params.NUM_PROFILE_SAMPLES)] for _ in range(num_tongue_ribs)]

    for rib_idx in range(1, num_tongue_ribs):
        tongue_ribs[rib_idx].left = Point2D(-standard_half_width, tongue_ribs[rib_idx].left_side_height)
        tongue_ribs[rib_idx].right = Point2D(standard_half_width, tongue_ribs[rib_idx].right_side_height)
        upper_profile, lower_profile, _ = get_cross_profiles(
            vt,
            tongue_ribs[rib_idx].point,
            tongue_ribs[rib_idx].normal,
            consider_tongue=False,
        )
        upper_profiles[rib_idx] = upper_profile
        lower_profiles[rib_idx] = lower_profile

        left_edge = tongue_ribs[rib_idx].left
        a2_left[rib_idx] = 3.0 * left_edge.y / max(left_edge.x * left_edge.x, 1.0e-9)
        a3_left[rib_idx] = -2.0 * left_edge.y / max(left_edge.x * left_edge.x * left_edge.x, 1.0e-9)
        tongue_ribs[rib_idx].min_x = -min_half_width
        last_inside = False
        for sample_idx in range(0, half_profile - min_half_width_samples):
            x_value = -0.5 * params.PROFILE_LENGTH + float(sample_idx) * params.PROFILE_SAMPLE_LENGTH
            is_inside = False
            if (
                upper_profile[sample_idx] != invalid
                and lower_profile[sample_idx] != invalid
                and x_value >= left_edge.x
            ):
                y_value = a2_left[rib_idx] * x_value * x_value + a3_left[rib_idx] * x_value * x_value * x_value
                if lower_profile[sample_idx] <= y_value <= upper_profile[sample_idx]:
                    is_inside = True
            if (not last_inside) and is_inside:
                tongue_ribs[rib_idx].min_x = x_value
            last_inside = is_inside

        right_edge = tongue_ribs[rib_idx].right
        a2_right[rib_idx] = 3.0 * right_edge.y / max(right_edge.x * right_edge.x, 1.0e-9)
        a3_right[rib_idx] = -2.0 * right_edge.y / max(right_edge.x * right_edge.x * right_edge.x, 1.0e-9)
        tongue_ribs[rib_idx].max_x = min_half_width
        last_inside = False
        for sample_idx in range(params.NUM_PROFILE_SAMPLES - 1, half_profile + min_half_width_samples, -1):
            x_value = float(sample_idx) * params.PROFILE_SAMPLE_LENGTH - 0.5 * params.PROFILE_LENGTH
            is_inside = False
            if (
                upper_profile[sample_idx] != invalid
                and lower_profile[sample_idx] != invalid
                and x_value <= right_edge.x
            ):
                y_value = a2_right[rib_idx] * x_value * x_value + a3_right[rib_idx] * x_value * x_value * x_value
                if lower_profile[sample_idx] <= y_value <= upper_profile[sample_idx]:
                    is_inside = True
            if (not last_inside) and is_inside:
                tongue_ribs[rib_idx].max_x = x_value
            last_inside = is_inside

        if tongue_ribs[rib_idx].max_x < min_half_width:
            tongue_ribs[rib_idx].max_x = min_half_width
        if tongue_ribs[rib_idx].min_x > -min_half_width:
            tongue_ribs[rib_idx].min_x = -min_half_width

    lowpass_coeff = 0.5
    for rib_idx in range(1, num_tongue_ribs):
        z_value = tongue_ribs[rib_idx - 1].max_x + lowpass_coeff * (
            tongue_ribs[rib_idx].max_x - tongue_ribs[rib_idx - 1].max_x
        )
        if z_value <= tongue_ribs[rib_idx].max_x:
            tongue_ribs[rib_idx].max_x = z_value

        z_value = tongue_ribs[rib_idx - 1].min_x + lowpass_coeff * (
            tongue_ribs[rib_idx].min_x - tongue_ribs[rib_idx - 1].min_x
        )
        if z_value >= tongue_ribs[rib_idx].min_x:
            tongue_ribs[rib_idx].min_x = z_value

    for rib_idx in range(num_tongue_ribs - 2, -1, -1):
        z_value = tongue_ribs[rib_idx + 1].max_x + lowpass_coeff * (
            tongue_ribs[rib_idx].max_x - tongue_ribs[rib_idx + 1].max_x
        )
        if z_value <= tongue_ribs[rib_idx].max_x:
            tongue_ribs[rib_idx].max_x = z_value

        z_value = tongue_ribs[rib_idx + 1].min_x + lowpass_coeff * (
            tongue_ribs[rib_idx].min_x - tongue_ribs[rib_idx + 1].min_x
        )
        if z_value >= tongue_ribs[rib_idx].min_x:
            tongue_ribs[rib_idx].min_x = z_value

    safety_margin_cm = 0.02
    for rib_idx in range(1, num_tongue_ribs):
        for point_idx in range(params.NUM_TONGUE_POINTS):
            z_value = tongue_ribs[rib_idx].min_x + (
                (tongue_ribs[rib_idx].max_x - tongue_ribs[rib_idx].min_x)
                * float(point_idx)
                / float(params.NUM_TONGUE_POINTS - 1)
            )
            if point_idx < params.NUM_TONGUE_POINTS // 2:
                offset_t = a2_left[rib_idx] * z_value * z_value + a3_left[rib_idx] * z_value * z_value * z_value
            else:
                offset_t = a2_right[rib_idx] * z_value * z_value + a3_right[rib_idx] * z_value * z_value * z_value

            sample_index = int((z_value + 0.5 * params.PROFILE_LENGTH) / params.PROFILE_SAMPLE_LENGTH)
            ratio = (
                z_value
                + 0.5 * params.PROFILE_LENGTH
                - float(sample_index) * params.PROFILE_SAMPLE_LENGTH
            ) / params.PROFILE_SAMPLE_LENGTH
            if sample_index < 0:
                sample_index = 0
                ratio = 0.0
            if sample_index > params.NUM_PROFILE_SAMPLES - 2:
                sample_index = params.NUM_PROFILE_SAMPLES - 2
                ratio = 1.0

            if (
                upper_profiles[rib_idx][sample_index] != invalid
                and upper_profiles[rib_idx][sample_index + 1] != invalid
            ):
                y_upper = upper_profiles[rib_idx][sample_index] + ratio * (
                    upper_profiles[rib_idx][sample_index + 1] - upper_profiles[rib_idx][sample_index]
                ) - safety_margin_cm
                if offset_t > y_upper:
                    offset_t = y_upper

            if rib_idx < num_tongue_ribs - 1:
                if (
                    lower_profiles[rib_idx][sample_index] != invalid
                    and lower_profiles[rib_idx][sample_index + 1] != invalid
                ):
                    y_lower = lower_profiles[rib_idx][sample_index] + ratio * (
                        lower_profiles[rib_idx][sample_index + 1] - lower_profiles[rib_idx][sample_index]
                    ) + safety_margin_cm
                    if offset_t < y_lower:
                        offset_t = y_lower

            tongue_surface.set_vertex(
                rib_idx,
                point_idx,
                Point3D(
                    tongue_ribs[rib_idx].point.x + tongue_ribs[rib_idx].normal.x * offset_t,
                    tongue_ribs[rib_idx].point.y + tongue_ribs[rib_idx].normal.y * offset_t,
                    z_value,
                ),
            )

    center_point_index = tongue_surface.num_points // 2
    if vt.anatomy.tongue_tonsil_length_cm < epsilon:
        vt.anatomy.tongue_tonsil_length_cm = epsilon
    rib_pos_cm = [0.0 for _ in range(num_tongue_ribs)]
    last_relevant_rib = 0
    for rib_idx in range(1, num_tongue_ribs):
        point_prev = tongue_surface.get_vertex(rib_idx - 1, center_point_index)
        point_curr = tongue_surface.get_vertex(rib_idx, center_point_index)
        delta_x = point_curr.x - point_prev.x
        delta_y = point_curr.y - point_prev.y
        rib_pos_cm[rib_idx] = rib_pos_cm[rib_idx - 1] + math.sqrt(delta_x * delta_x + delta_y * delta_y)
        if rib_pos_cm[rib_idx] > vt.anatomy.tongue_tonsil_length_cm:
            last_relevant_rib = rib_idx - 1
            break

    for rib_idx in range(0, last_relevant_rib + 1):
        normal = tongue_ribs[rib_idx].normal
        pos_cm = rib_pos_cm[rib_idx]
        for point_idx in range(tongue_surface.num_points):
            hump = (
                vt.anatomy.tongue_tonsil_height_cm
                * (0.5 - 0.5 * math.cos(2.0 * math.pi * pos_cm / vt.anatomy.tongue_tonsil_length_cm))
                * (0.5 - 0.5 * math.cos(2.0 * math.pi * float(point_idx) / float(tongue_surface.num_points - 1)))
            )
            point = tongue_surface.get_vertex(rib_idx, point_idx)
            point.x += normal.x * hump
            point.y += normal.y * hump
            tongue_surface.set_vertex(rib_idx, point_idx, point)


def _merge_outline_by_x(points: List[Point2D], upper: bool) -> List[Point2D]:
    if not points:
        return []
    points_sorted = sorted(points, key=lambda point: point.x)
    merged: List[Point2D] = []
    bucket: List[Point2D] = [points_sorted[0]]
    tol = 1e-5

    for point in points_sorted[1:]:
        if abs(point.x - bucket[-1].x) <= tol:
            bucket.append(point)
            continue

        y_value = max(p.y for p in bucket) if upper else min(p.y for p in bucket)
        x_value = sum(p.x for p in bucket) / len(bucket)
        merged.append(Point2D(x_value, y_value))
        bucket = [point]

    y_value = max(p.y for p in bucket) if upper else min(p.y for p in bucket)
    x_value = sum(p.x for p in bucket) / len(bucket)
    merged.append(Point2D(x_value, y_value))
    return merged


def _build_outline_points(surface: "Surface", point_idx: int) -> List[Point2D]:
    if surface is None:
        return []
    result: List[Point2D] = []
    safe_point_idx = int(_clamp(point_idx, 0, surface.num_points - 1))
    for rib_idx in range(surface.num_ribs):
        vertex = surface.get_vertex(rib_idx, safe_point_idx)
        result.append(Point2D(vertex.x, vertex.y))
    return result


def _sample_outline_y(outline: "LineStrip2D", x: float, upper: bool) -> Optional[float]:
    if outline.num_points == 0:
        return None
    if outline.num_points == 1:
        return outline.p[0].y

    candidates: List[float] = []
    for idx in range(1, outline.num_points):
        point0 = outline.p[idx - 1]
        point1 = outline.p[idx]
        x0 = point0.x
        x1 = point1.x
        dx = x1 - x0

        if abs(dx) < 1e-9:
            if abs(x - x0) <= 1e-6:
                candidates.extend([point0.y, point1.y])
            continue

        if (x0 <= x <= x1) or (x1 <= x <= x0):
            ratio = (x - x0) / dx
            y = point0.y + ratio * (point1.y - point0.y)
            candidates.append(y)

    if candidates:
        return max(candidates) if upper else min(candidates)

    nearest = min(outline.p, key=lambda point: abs(point.x - x))
    return nearest.y


def _line_strip_intersections(origin: Point2D, direction: Point2D, outline: "LineStrip2D") -> List[float]:
    if outline.num_points < 2:
        return []

    t_values: List[float] = []
    for idx in range(1, outline.num_points):
        point0 = outline.p[idx - 1]
        point1 = outline.p[idx]
        segment = point1 - point0

        denom = direction.x * segment.y - direction.y * segment.x
        if abs(denom) < 1e-9:
            continue

        diff = point0 - origin
        t_value = (diff.x * segment.y - diff.y * segment.x) / denom
        u_value = (diff.x * direction.y - diff.y * direction.x) / denom

        if -1e-6 <= u_value <= 1.0 + 1e-6:
            t_values.append(t_value)

    return t_values


def calc_surfaces(vt: "VocalTract") -> None:
    _rebuild_cover_geometry_cpp(vt)
    _copy_surface(vt.surfaces_list[params.SurfaceIndex.UVULA], vt.surfaces_list[params.SurfaceIndex.UVULA_ORIGINAL])
    _copy_surface(vt.surfaces_list[params.SurfaceIndex.EPIGLOTTIS], vt.surfaces_list[params.SurfaceIndex.EPIGLOTTIS_ORIGINAL])

    _update_lower_gums_edges(vt)
    _rebuild_throat_front(vt)

    _build_lips_geometry(vt)
    upper_cover = vt.surfaces_list[params.SurfaceIndex.UPPER_COVER]
    lower_cover = vt.surfaces_list[params.SurfaceIndex.LOWER_COVER]
    upper_teeth = vt.surfaces_list[params.SurfaceIndex.UPPER_TEETH]
    lower_teeth = vt.surfaces_list[params.SurfaceIndex.LOWER_TEETH]
    upper_lip = vt.surfaces_list[params.SurfaceIndex.UPPER_LIP]
    lower_lip = vt.surfaces_list[params.SurfaceIndex.LOWER_LIP]
    if (
        upper_cover is None
        or lower_cover is None
        or upper_teeth is None
        or lower_teeth is None
        or upper_lip is None
        or lower_lip is None
    ):
        raise RuntimeError("[领域错误]: 声道轮廓构建失败 - 关键表面缺失")

    upper_cover_point = params.NUM_UPPER_COVER_POINTS - 1
    lower_cover_point = params.NUM_LOWER_COVER_POINTS - 1
    lip_rib = params.NUM_LIP_RIBS - 1
    teeth_rib = params.NUM_TEETH_RIBS - 1

    upper_outline_points: List[Point2D] = []
    for rib_idx in range(params.NUM_UPPER_COVER_RIBS):
        vertex = upper_cover.get_vertex(rib_idx, upper_cover_point)
        upper_outline_points.append(Point2D(vertex.x, vertex.y))
    for point_idx in range(params.NUM_TEETH_POINTS - 1):
        vertex = upper_teeth.get_vertex(teeth_rib, point_idx)
        upper_outline_points.append(Point2D(vertex.x, vertex.y))
    for point_idx in range(params.NUM_INNER_LIP_POINTS + 2):
        vertex = upper_lip.get_vertex(lip_rib, point_idx)
        upper_outline_points.append(Point2D(vertex.x, vertex.y))
    while len(upper_outline_points) >= 2 and upper_outline_points[-1].x <= upper_outline_points[-2].x:
        upper_outline_points.pop()

    lower_outline_points: List[Point2D] = []
    for rib_idx in range(params.NUM_LOWER_COVER_RIBS):
        vertex = lower_cover.get_vertex(rib_idx, lower_cover_point)
        lower_outline_points.append(Point2D(vertex.x, vertex.y))
    for point_idx in range(params.NUM_TEETH_POINTS - 1):
        vertex = lower_teeth.get_vertex(teeth_rib, point_idx)
        lower_outline_points.append(Point2D(vertex.x, vertex.y))
    for point_idx in range(params.NUM_INNER_LIP_POINTS + 2):
        vertex = lower_lip.get_vertex(lip_rib, point_idx)
        lower_outline_points.append(Point2D(vertex.x, vertex.y))
    while len(lower_outline_points) >= 2 and lower_outline_points[-1].x <= lower_outline_points[-2].x:
        lower_outline_points.pop()

    if len(upper_outline_points) < 2 or len(lower_outline_points) < 2:
        raise RuntimeError("[领域错误]: 声道轮廓构建失败 - 有效轮廓点不足")

    vt.upper_outline.set_points(upper_outline_points)
    vt.lower_outline.set_points(lower_outline_points)

    _build_tongue_geometry(vt)

    tongue_surface = vt.surfaces_list[params.SurfaceIndex.TONGUE]
    epiglottis = vt.surfaces_list[params.SurfaceIndex.EPIGLOTTIS]
    epiglottis_original = vt.surfaces_list[params.SurfaceIndex.EPIGLOTTIS_ORIGINAL]
    uvula = vt.surfaces_list[params.SurfaceIndex.UVULA]
    uvula_original = vt.surfaces_list[params.SurfaceIndex.UVULA_ORIGINAL]
    if tongue_surface is None or epiglottis is None or epiglottis_original is None or uvula is None or uvula_original is None:
        raise RuntimeError("[领域错误]: 舌体/会厌/悬雍垂计算失败 - 表面缺失")

    anchor_p = lower_cover.get_vertex(params.NUM_LARYNX_RIBS - 2, params.NUM_LOWER_COVER_POINTS - 1)
    anchor_q = lower_cover.get_vertex(params.NUM_LARYNX_RIBS - 1, params.NUM_LOWER_COVER_POINTS - 1)
    direction = anchor_q - anchor_p
    direction.normalize()
    epiglottis_anchor = anchor_p + float(vt.anatomy.epiglottis_width_cm) * direction

    for point_idx in range(params.NUM_EPIGLOTTIS_POINTS):
        original = epiglottis_original.get_vertex(0, point_idx)
        rotated_x = direction.x * original.x - direction.y * original.y
        rotated_y = direction.y * original.x + direction.x * original.y
        epiglottis.set_vertex(
            0,
            point_idx,
            Point3D(
                rotated_x + epiglottis_anchor.x,
                rotated_y + epiglottis_anchor.y,
                original.z + epiglottis_anchor.z,
            ),
        )

    min_angle_rad = float(vt.anatomy.epiglottis_angle_deg) * math.pi / 180.0
    tongue_mid_point = int(params.NUM_TONGUE_POINTS / 2)
    max_tongue_rib = min(10, params.NUM_TONGUE_RIBS)
    for rib_idx in range(1, max_tongue_rib):
        tongue_offset = tongue_surface.get_vertex(rib_idx, tongue_mid_point) - epiglottis_anchor
        if tongue_offset.y < float(vt.anatomy.epiglottis_height_cm):
            if tongue_offset.y < 0.01:
                tongue_offset.y = 0.01
            angle_rad = math.atan2(tongue_offset.y, tongue_offset.x)
            if angle_rad > min_angle_rad:
                min_angle_rad = angle_rad

    min_angle_rad -= 0.5 * math.pi
    if min_angle_rad > 0.25 * math.pi:
        min_angle_rad = 0.25 * math.pi
    sine = math.sin(min_angle_rad)
    cosine = math.cos(min_angle_rad)
    for rib_idx in range(1, params.NUM_EPIGLOTTIS_RIBS):
        for point_idx in range(params.NUM_EPIGLOTTIS_POINTS):
            original = epiglottis_original.get_vertex(rib_idx, point_idx)
            rotated_x = cosine * original.x - sine * original.y
            rotated_y = sine * original.x + cosine * original.y
            epiglottis.set_vertex(
                rib_idx,
                point_idx,
                Point3D(
                    rotated_x + epiglottis_anchor.x,
                    rotated_y + epiglottis_anchor.y,
                    original.z + epiglottis_anchor.z,
                ),
            )

    uvula_anchor = upper_cover.get_vertex(
        params.NUM_LARYNX_RIBS + params.NUM_PHARYNX_RIBS + 1,
        params.NUM_UPPER_COVER_POINTS - 1,
    )
    for rib_idx in range(params.NUM_UVULA_RIBS):
        for point_idx in range(params.NUM_UVULA_POINTS):
            original = uvula_original.get_vertex(rib_idx, point_idx)
            uvula.set_vertex(rib_idx, point_idx, original + uvula_anchor)

    tongue_outline_points: List[Point2D] = []
    for rib_idx in range(params.NUM_TONGUE_RIBS):
        vertex = tongue_surface.get_vertex(rib_idx, tongue_mid_point)
        tongue_outline_points.append(Point2D(vertex.x, vertex.y))
    vt.tongue_outline.set_points(tongue_outline_points)

    epiglottis_outline_points: List[Point2D] = []
    for rib_idx in range(params.NUM_EPIGLOTTIS_RIBS):
        top_vertex = epiglottis.get_vertex(rib_idx, params.NUM_EPIGLOTTIS_POINTS - 1)
        epiglottis_outline_points.append(Point2D(top_vertex.x, top_vertex.y))
    for rib_idx in range(params.NUM_EPIGLOTTIS_RIBS - 2, -1, -1):
        bottom_vertex = epiglottis.get_vertex(rib_idx, 0)
        epiglottis_outline_points.append(Point2D(bottom_vertex.x, bottom_vertex.y))
    vt.epiglottis_outline.set_points(epiglottis_outline_points)

    _set_twoside_surface_by_points(
        vt.surfaces_list[params.SurfaceIndex.UPPER_COVER_TWOSIDE],
        vt.surfaces_list[params.SurfaceIndex.UPPER_COVER],
    )
    _set_twoside_surface_by_points(
        vt.surfaces_list[params.SurfaceIndex.LOWER_COVER_TWOSIDE],
        vt.surfaces_list[params.SurfaceIndex.LOWER_COVER],
    )
    _set_twoside_surface_by_ribs(
        vt.surfaces_list[params.SurfaceIndex.UPPER_TEETH_TWOSIDE],
        vt.surfaces_list[params.SurfaceIndex.UPPER_TEETH],
    )
    _set_twoside_surface_by_ribs(
        vt.surfaces_list[params.SurfaceIndex.LOWER_TEETH_TWOSIDE],
        vt.surfaces_list[params.SurfaceIndex.LOWER_TEETH],
    )
    _set_twoside_surface_by_ribs(
        vt.surfaces_list[params.SurfaceIndex.UPPER_LIP_TWOSIDE],
        vt.surfaces_list[params.SurfaceIndex.UPPER_LIP],
    )
    _set_twoside_surface_by_ribs(
        vt.surfaces_list[params.SurfaceIndex.LOWER_LIP_TWOSIDE],
        vt.surfaces_list[params.SurfaceIndex.LOWER_LIP],
    )
    _set_twoside_surface_by_points(
        vt.surfaces_list[params.SurfaceIndex.UVULA_TWOSIDE],
        vt.surfaces_list[params.SurfaceIndex.UVULA],
    )
    _set_twoside_surface_by_points(
        vt.surfaces_list[params.SurfaceIndex.EPIGLOTTIS_TWOSIDE],
        vt.surfaces_list[params.SurfaceIndex.EPIGLOTTIS],
    )
    _build_cover_side_surfaces(vt)

    # 与 C++ `intersectionsPrepared[] = false` 保持一致：
    # 每次几何更新后都需要重新准备切割加速结构，不能复用旧几何缓存。
    for surface in vt.surfaces_list:
        if surface is not None:
            setattr(surface, "_dvts_intersections_prepared", False)



def _calc_center_line_simple(vt: "VocalTract") -> None:
    if vt.upper_outline.num_points < 2 or vt.lower_outline.num_points < 2:
        raise RuntimeError("[领域错误]: 中线计算失败 - 上下轮廓未就绪")

    upper_min_x = min(point.x for point in vt.upper_outline.p)
    upper_max_x = max(point.x for point in vt.upper_outline.p)
    lower_min_x = min(point.x for point in vt.lower_outline.p)
    lower_max_x = max(point.x for point in vt.lower_outline.p)

    x_min = max(upper_min_x, lower_min_x)
    x_max = min(upper_max_x, lower_max_x)

    if x_max <= x_min + 1e-6:
        x_min = min(upper_min_x, lower_min_x)
        x_max = max(upper_max_x, lower_max_x)

    if x_max <= x_min + 1e-6:
        raise RuntimeError("[领域错误]: 中线计算失败 - 声道长度为零")

    center_points: List[Point2D] = []
    center_normals: List[Point2D] = []
    center_line: List[params.CenterLinePoint] = []

    for idx in range(params.NUM_CENTERLINE_POINTS):
        ratio = idx / (params.NUM_CENTERLINE_POINTS - 1)
        x_coord = x_min + ratio * (x_max - x_min)
        upper_y = _sample_outline_y(vt.upper_outline, x_coord, upper=True)
        lower_y = _sample_outline_y(vt.lower_outline, x_coord, upper=False)

        if upper_y is None or lower_y is None:
            raise RuntimeError("[领域错误]: 中线采样失败 - 无法匹配上下轮廓")
        if upper_y < lower_y:
            upper_y, lower_y = lower_y, upper_y

        center_y = 0.5 * (upper_y + lower_y)
        center_point = Point2D(x_coord, center_y)
        center_points.append(center_point)
        center_normals.append(Point2D(0.0, 1.0))
        center_line.append(params.CenterLinePoint(point=center_point, normal=Point2D(0.0, 1.0), min_v=lower_y, max_v=upper_y))

    for idx in range(params.NUM_CENTERLINE_POINTS):
        if idx == 0:
            tangent = center_points[1] - center_points[0]
        elif idx == params.NUM_CENTERLINE_POINTS - 1:
            tangent = center_points[idx] - center_points[idx - 1]
        else:
            tangent = center_points[idx + 1] - center_points[idx - 1]

        if tangent.magnitude() < 1e-9:
            normal = Point2D(0.0, 1.0)
        else:
            normal = Point2D(-tangent.y, tangent.x).normalize()
        center_normals[idx] = normal
        center_line[idx].normal = Point2D(normal.x, normal.y)

    cumulative_pos = 0.0
    center_line[0].pos = 0.0
    for idx in range(1, params.NUM_CENTERLINE_POINTS):
        cumulative_pos += (center_points[idx] - center_points[idx - 1]).magnitude()
        center_line[idx].pos = cumulative_pos

    vt.center_line_points = center_points
    vt.center_line_normals = center_normals
    vt.center_line = center_line
    vt.center_line_length = cumulative_pos


def _calc_center_line_cpp(vt: "VocalTract") -> None:
    if vt.upper_outline.num_points < 3 or vt.lower_outline.num_points < 3:
        raise RuntimeError("[领域错误]: 中线计算失败 - 上下轮廓点数不足")

    skipped_lip_points = 1
    num_points = params.NUM_CENTERLINE_POINTS
    exp = params.NUM_CENTERLINE_POINTS_EXPONENT

    tongue_center = Point2D(
        _get_param_value(vt, int(params.ParamIndex.TCX)),
        _get_param_value(vt, int(params.ParamIndex.TCY)),
    )
    tongue_radius_x = float(vt.anatomy.tongue_center_radius_x_cm)
    tongue_radius_y = float(vt.anatomy.tongue_center_radius_y_cm)

    min_y = vt.upper_outline.p[0].y
    lower0 = vt.lower_outline.p[0]
    if lower0.y > min_y:
        min_y = lower0.y

    max_upper_index = vt.upper_outline.num_points - 1 - skipped_lip_points
    max_lower_index = vt.lower_outline.num_points - 1 - skipped_lip_points
    if max_upper_index < 0 or max_lower_index < 0:
        raise RuntimeError("[领域错误]: 中线计算失败 - 轮廓边界索引非法")
    max_x = vt.upper_outline.p[max_upper_index].x
    lower_edge = vt.lower_outline.p[max_lower_index]
    if lower_edge.x < max_x:
        max_x = lower_edge.x

    mu_section_length_0 = tongue_center.y - min_y
    mu_section_length_1 = 0.5 * math.pi * 0.5 * (tongue_radius_x + tongue_radius_y)
    mu_section_length_2 = max_x - tongue_center.x
    mu_length = mu_section_length_0 + mu_section_length_1 + mu_section_length_2
    if mu_length == 0.0:
        raise RuntimeError("[领域错误]: 中线计算失败 - 初始中线长度为零")

    nu_line: List[Point2D] = [Point2D() for _ in range(num_points + 2)]
    rough_points: List[Point2D] = [Point2D() for _ in range(num_points)]
    rough_normals: List[Point2D] = [Point2D() for _ in range(num_points)]

    for idx in range(num_points):
        pos = mu_length * float(idx) / float(num_points - 1)

        if pos <= mu_section_length_0:
            point = Point2D(tongue_center.x - tongue_radius_x, tongue_center.y - mu_section_length_0 + pos)
            normal = Point2D(-1.0, 0.0)
        elif pos <= mu_section_length_0 + mu_section_length_1:
            angle = math.pi - 0.5 * math.pi * (pos - mu_section_length_0) / mu_section_length_1
            cosine = math.cos(angle)
            sine = math.sin(angle)
            point = Point2D(tongue_center.x + tongue_radius_x * cosine, tongue_center.y + tongue_radius_y * sine)
            normal = Point2D(tongue_radius_x * cosine, tongue_radius_y * sine).normalize()
        else:
            point = Point2D(
                tongue_center.x + pos - mu_section_length_0 - mu_section_length_1,
                tongue_center.y + tongue_radius_y,
            )
            normal = Point2D(0.0, 1.0)

        ok_top, _t_top, intersection_top = vt.upper_outline.get_special_intersection(point, normal)
        ok_bottom, t_bottom, intersection_bottom = vt.lower_outline.get_intersection_with_greatest_t(point, normal)
        ok_tongue, t_tongue, intersection_tongue = vt.tongue_outline.get_intersection_with_greatest_t(point, normal)
        ok_epi, t_epi, intersection_epi = vt.epiglottis_outline.get_intersection_with_greatest_t(point, normal)

        intersections = [intersection_top, intersection_bottom, intersection_tongue, intersection_epi]
        top_index = 0
        bottom_index = 1
        if ok_top:
            t_max = -1_000_000.0
            if ok_bottom and t_bottom > t_max:
                t_max = t_bottom
                bottom_index = 1
            if ok_tongue and t_tongue > t_max:
                t_max = t_tongue
                bottom_index = 2
            if ok_epi and t_epi > t_max:
                t_max = t_epi
                bottom_index = 3
        center_point = 0.5 * (intersections[top_index] + intersections[bottom_index])

        rough_points[idx] = center_point
        rough_normals[idx] = Point2D(normal.x, normal.y)
        nu_line[idx + 1] = Point2D(center_point.x, center_point.y)

    nu_pos = [0.0 for _ in range(num_points + 2)]
    nu_line_length = 0.0
    for idx in range(1, num_points):
        nu_pos[idx] = nu_line_length
        nu_line_length += (nu_line[idx + 1] - nu_line[idx]).magnitude()
    nu_pos[num_points] = nu_line_length

    range_length = 2.0
    half_range = 0.5 * range_length
    nu_line[0] = nu_line[1] + Point2D(0.0, -half_range)
    nu_line[num_points + 1] = nu_line[num_points] + Point2D(half_range, 0.0)
    nu_pos[0] = -half_range
    nu_pos[num_points + 1] = nu_pos[num_points] + half_range

    center_line: List[params.CenterLinePoint] = [params.CenterLinePoint() for _ in range(num_points)]
    center_points: List[Point2D] = [Point2D() for _ in range(num_points)]
    center_normals: List[Point2D] = [Point2D() for _ in range(num_points)]

    range_index = [0, 0]
    for idx in range(num_points):
        range_pos = [
            nu_line_length * float(idx) / float(num_points - 1) - half_range,
            0.0,
        ]
        range_pos[1] = range_pos[0] + range_length
        range_point = [Point2D(), Point2D()]

        for side in (0, 1):
            while (range_index[side] < num_points) and (range_pos[side] > nu_pos[range_index[side] + 1]):
                range_index[side] += 1
            denominator = nu_pos[range_index[side] + 1] - nu_pos[range_index[side]]
            if abs(denominator) < 1.0e-12:
                d = 0.0
            else:
                d = (range_pos[side] - nu_pos[range_index[side]]) / denominator
            range_point[side] = nu_line[range_index[side]] + d * (nu_line[range_index[side] + 1] - nu_line[range_index[side]])

        center_point = Point2D(0.0, 0.0)
        if range_index[1] > range_index[0]:
            seg_len = (nu_line[range_index[0] + 1] - range_point[0]).magnitude()
            center_point += seg_len * 0.5 * (range_point[0] + nu_line[range_index[0] + 1])
            for seg_idx in range(range_index[0] + 1, range_index[1]):
                seg_len = (nu_line[seg_idx + 1] - nu_line[seg_idx]).magnitude()
                center_point += seg_len * 0.5 * (nu_line[seg_idx] + nu_line[seg_idx + 1])
            seg_len = (range_point[1] - nu_line[range_index[1]]).magnitude()
            center_point += seg_len * 0.5 * (range_point[1] + nu_line[range_index[1]])
        else:
            seg_len = (range_point[1] - range_point[0]).magnitude()
            center_point += seg_len * 0.5 * (range_point[0] + range_point[1])

        center_point = center_point / range_length
        center_points[idx] = center_point
        center_line[idx].point = Point2D(center_point.x, center_point.y)

    q0 = vt.upper_outline.p[0]
    q1 = vt.lower_outline.p[0]
    x0 = center_line[0].point.x
    denominator = q1.x - q0.x
    if abs(denominator) < 1.0e-12:
        denominator = 1.0e-12
    y0 = q0.y + (q1.y - q0.y) * (x0 - q0.x) / denominator
    center_line[0].point.y = y0
    for idx in range(min(4, num_points)):
        if center_line[idx].point.y < y0:
            center_line[idx].point.y = y0
    normal0 = (q0 - q1).normalize()
    center_line[0].normal = Point2D(normal0.x, normal0.y)

    center_line[num_points - 1].normal = Point2D(0.0, 1.0)
    xl = center_line[num_points - 1].point.x
    for idx in range(max(0, num_points - 4), num_points):
        if center_line[idx].point.x > xl:
            center_line[idx].point.x = xl

    center_line_length = 0.0
    for idx in range(num_points - 1):
        center_line[idx].pos = center_line_length
        center_line_length += (center_line[idx + 1].point - center_line[idx].point).magnitude()
    center_line[num_points - 1].pos = center_line_length

    for idx in range(1, num_points - 1):
        q_left = center_line[idx].point - center_line[idx - 1].point
        q_right = center_line[idx + 1].point - center_line[idx].point
        normal = Point2D(-q_left.y - q_right.y, q_left.x + q_right.x).normalize()
        center_line[idx].normal = normal

    reserved = [0.0 for _ in range(num_points)]

    def _recalc_minmax(point_index: int) -> None:
        point = center_line[point_index].point
        normal = center_line[point_index].normal
        ok_top, t_top, _s_top = vt.upper_outline.get_special_intersection(point, normal)
        ok_lower, t_lower, _s_lower = vt.lower_outline.get_intersection_with_greatest_t(point, normal)
        ok_tongue, t_tongue, _s_tongue = vt.tongue_outline.get_intersection_with_greatest_t(point, normal)

        if ok_top:
            center_line[point_index].max_v = t_top
        else:
            center_line[point_index].max_v = 3.0
            reserved[point_index] = 1.0

        if ok_lower and ok_tongue:
            center_line[point_index].min_v = t_lower if t_lower > t_tongue else t_tongue
        elif ok_tongue:
            center_line[point_index].min_v = t_tongue
        elif ok_lower:
            center_line[point_index].min_v = t_lower
        else:
            center_line[point_index].min_v = -3.0
            reserved[point_index] = 1.0

    for idx in range(num_points):
        reserved[idx] = 0.0
        _recalc_minmax(idx)

    def _verify_center_line_normal(left_idx: int, mid_idx: int, right_idx: int) -> None:
        point = center_line[mid_idx].point
        normal = center_line[mid_idx].normal

        left_base = center_line[left_idx].point + center_line[left_idx].min_v * center_line[left_idx].normal
        left_vec = (
            center_line[left_idx].point + center_line[left_idx].max_v * center_line[left_idx].normal - left_base
        )
        delta = point - left_base
        denominator = normal.x * left_vec.y - normal.y * left_vec.x
        if denominator != 0.0:
            s = (normal.x * delta.y - delta.x * normal.y) / denominator
            if 0.0 <= s <= 1.0:
                t = (left_vec.x * delta.y - delta.x * left_vec.y) / denominator
                if center_line[mid_idx].min_v <= t <= 0.0:
                    new_normal = (center_line[mid_idx].point - left_base).normalize()
                    center_line[mid_idx].normal = new_normal
                    reserved[mid_idx] = 1.0
                elif 0.0 <= t <= center_line[mid_idx].max_v:
                    new_normal = (left_base + left_vec - center_line[mid_idx].point).normalize()
                    center_line[mid_idx].normal = new_normal
                    reserved[mid_idx] = 1.0

        point = center_line[mid_idx].point
        normal = center_line[mid_idx].normal
        right_base = center_line[right_idx].point + center_line[right_idx].min_v * center_line[right_idx].normal
        right_vec = (
            center_line[right_idx].point + center_line[right_idx].max_v * center_line[right_idx].normal - right_base
        )
        delta = point - right_base
        denominator = normal.x * right_vec.y - normal.y * right_vec.x
        if denominator != 0.0:
            s = (normal.x * delta.y - delta.x * normal.y) / denominator
            if 0.0 <= s <= 1.0:
                t = (right_vec.x * delta.y - delta.x * right_vec.y) / denominator
                if center_line[mid_idx].min_v <= t <= 0.0:
                    new_normal = (center_line[mid_idx].point - right_base).normalize()
                    center_line[mid_idx].normal = new_normal
                    reserved[mid_idx] = 1.0
                elif 0.0 <= t <= center_line[mid_idx].max_v:
                    new_normal = (right_base + right_vec - center_line[mid_idx].point).normalize()
                    center_line[mid_idx].normal = new_normal
                    reserved[mid_idx] = 1.0

    for level in range(exp):
        num_to_check = 1 << level
        distance = 1 << (exp - level)
        index = distance // 2
        for _ in range(num_to_check):
            _verify_center_line_normal(index - distance // 2, index, index + distance // 2)
            if reserved[index] != 0.0:
                reserved[index] = 0.0
                _recalc_minmax(index)
            index += distance

    vt.center_line = center_line
    vt.center_line_points = [Point2D(item.point.x, item.point.y) for item in center_line]
    vt.center_line_normals = [Point2D(item.normal.x, item.normal.y) for item in center_line]
    vt.center_line_length = center_line_length
    vt.rough_center_line = [
        params.CenterLinePoint(point=rough_points[idx], normal=rough_normals[idx]) for idx in range(num_points)
    ]


def calc_center_line(vt: "VocalTract") -> None:
    _calc_center_line_cpp(vt)


def _insert_upper_profile_line(
    point0: Point2D,
    point1: Point2D,
    surface_index: int,
    upper_profile: List[float],
    upper_profile_surface: List[int],
) -> None:
    if point0.x == point1.x:
        return

    p0 = Point2D(point0.x + 0.5 * params.PROFILE_LENGTH, point0.y)
    p1 = Point2D(point1.x + 0.5 * params.PROFILE_LENGTH, point1.y)

    if p0.x > p1.x:
        p0, p1 = p1, p0

    vector = p1 - p0
    vector.normalize()
    p0 = p0 - vector * 0.01
    p1 = p1 + vector * 0.01

    first_sample = int(p0.x / params.PROFILE_SAMPLE_LENGTH)
    last_sample = int(p1.x / params.PROFILE_SAMPLE_LENGTH)
    if first_sample == last_sample:
        return

    dx = params.PROFILE_SAMPLE_LENGTH
    dy = (p1.y - p0.y) * dx / (p1.x - p0.x)
    y_val = p0.y + (((first_sample + 1.0) * dx - p0.x) * dy / dx)

    for sample_idx in range(first_sample + 1, last_sample + 1):
        if 0 <= sample_idx < params.NUM_PROFILE_SAMPLES:
            if (params.MIN_PROFILE_VALUE <= y_val <= params.MAX_PROFILE_VALUE) and (y_val <= upper_profile[sample_idx]):
                upper_profile[sample_idx] = y_val
                upper_profile_surface[sample_idx] = surface_index
        y_val += dy



def _insert_lower_profile_line(
    point0: Point2D,
    point1: Point2D,
    surface_index: int,
    lower_profile: List[float],
    lower_profile_surface: List[int],
) -> None:
    if point0.x == point1.x:
        return

    p0 = Point2D(point0.x + 0.5 * params.PROFILE_LENGTH, point0.y)
    p1 = Point2D(point1.x + 0.5 * params.PROFILE_LENGTH, point1.y)

    if p0.x > p1.x:
        p0, p1 = p1, p0

    vector = p1 - p0
    vector.normalize()
    p0 = p0 - vector * 0.01
    p1 = p1 + vector * 0.01

    first_sample = int(p0.x / params.PROFILE_SAMPLE_LENGTH)
    last_sample = int(p1.x / params.PROFILE_SAMPLE_LENGTH)
    if first_sample == last_sample:
        return

    dx = params.PROFILE_SAMPLE_LENGTH
    dy = (p1.y - p0.y) * dx / (p1.x - p0.x)
    y_val = p0.y + (((first_sample + 1.0) * dx - p0.x) * dy / dx)

    for sample_idx in range(first_sample + 1, last_sample + 1):
        if 0 <= sample_idx < params.NUM_PROFILE_SAMPLES:
            if (params.MIN_PROFILE_VALUE <= y_val <= params.MAX_PROFILE_VALUE) and (y_val >= lower_profile[sample_idx]):
                lower_profile[sample_idx] = y_val
                lower_profile_surface[sample_idx] = surface_index
        y_val += dy



def _insert_lower_cover_profile_line(
    point0: Point2D,
    point1: Point2D,
    surface_index: int,
    upper_profile: List[float],
    upper_profile_surface: List[int],
    lower_profile: List[float],
    lower_profile_surface: List[int],
) -> None:
    if point0.x == point1.x:
        return

    p0 = Point2D(point0.x + 0.5 * params.PROFILE_LENGTH, point0.y)
    p1 = Point2D(point1.x + 0.5 * params.PROFILE_LENGTH, point1.y)

    if p0.x > p1.x:
        p0, p1 = p1, p0

    vector = p1 - p0
    vector.normalize()
    p0 = p0 - vector * 0.01
    p1 = p1 + vector * 0.01

    first_sample = int(p0.x / params.PROFILE_SAMPLE_LENGTH)
    last_sample = int(p1.x / params.PROFILE_SAMPLE_LENGTH)
    if first_sample == last_sample:
        return

    dx = params.PROFILE_SAMPLE_LENGTH
    dy = (p1.y - p0.y) * dx / (p1.x - p0.x)
    y_val = p0.y + (((first_sample + 1.0) * dx - p0.x) * dy / dx)

    for sample_idx in range(first_sample + 1, last_sample + 1):
        if 0 <= sample_idx < params.NUM_PROFILE_SAMPLES and params.MIN_PROFILE_VALUE <= y_val <= params.MAX_PROFILE_VALUE:
            if y_val >= lower_profile[sample_idx]:
                lower_profile[sample_idx] = y_val
                lower_profile_surface[sample_idx] = surface_index
            if upper_profile[sample_idx] == params.EXTREME_PROFILE_VALUE and y_val <= upper_profile[sample_idx]:
                upper_profile[sample_idx] = y_val
                upper_profile_surface[sample_idx] = surface_index
        y_val += dy



def _interpolate_invalid_samples(profile: List[float], invalid_value: float) -> None:
    left = 0
    right = 0
    sample_count = len(profile)

    for sample_idx in range(sample_count):
        if profile[sample_idx] == invalid_value:
            if sample_idx <= left or sample_idx >= right:
                left = sample_idx
                right = sample_idx
                while left > 0 and profile[left] == invalid_value:
                    left -= 1
                while right < sample_count - 1 and profile[right] == invalid_value:
                    right += 1

            if right > left:
                if profile[left] != invalid_value and profile[right] != invalid_value:
                    profile[sample_idx] = profile[left] + ((sample_idx - left) * (profile[right] - profile[left]) / (right - left))
                elif profile[left] != invalid_value:
                    profile[sample_idx] = profile[left]
                elif profile[right] != invalid_value:
                    profile[sample_idx] = profile[right]



def get_cross_profiles(
    vt: "VocalTract",
    point: Point2D,
    normal: Point2D,
    consider_tongue: bool = True,
) -> Tuple[List[float], List[float], Articulator]:
    min_squared_normal_length = 1e-7
    invalid = params.INVALID_PROFILE_SAMPLE
    sample_count = params.NUM_PROFILE_SAMPLES
    half_sample_count = params.NUM_PROFILE_SAMPLES // 2

    profile_surface_indices = [
        params.SurfaceIndex.UPPER_COVER,
        params.SurfaceIndex.UPPER_TEETH,
        params.SurfaceIndex.UPPER_LIP,
        params.SurfaceIndex.UVULA,
        params.SurfaceIndex.LOWER_COVER,
        params.SurfaceIndex.LOWER_TEETH,
        params.SurfaceIndex.LOWER_LIP,
        params.SurfaceIndex.EPIGLOTTIS,
        params.SurfaceIndex.LEFT_COVER,
        params.SurfaceIndex.RADIATION,
    ]

    upper_profile = [params.EXTREME_PROFILE_VALUE] * sample_count
    lower_profile = [-params.EXTREME_PROFILE_VALUE] * sample_count
    upper_profile_surface = [-1] * sample_count
    lower_profile_surface = [-1] * sample_count

    cuts: List[Tuple[Point2D, Point2D, Point2D, int, int]] = []
    max_list_entries = 2048

    for local_idx, global_idx in enumerate(profile_surface_indices):
        if (not consider_tongue) and global_idx in (params.SurfaceIndex.UVULA, params.SurfaceIndex.EPIGLOTTIS):
            continue

        surface = vt.surfaces_list[int(global_idx)]
        if surface is None:
            continue

        if not getattr(surface, "_dvts_intersections_prepared", False):
            surface.prepare_intersections()
            setattr(surface, "_dvts_intersections_prepared", True)
        surface.prepare_intersection(point, normal)
        triangle_indices = surface.get_triangle_list(max_entries=max_list_entries)

        for triangle_idx in triangle_indices:
            has_intersection, p0, p1, n = surface.get_triangle_intersection(triangle_idx, point, normal)
            if not has_intersection:
                continue
            if not (params.MIN_PROFILE_VALUE < p0.y < params.MAX_PROFILE_VALUE):
                continue
            if not (params.MIN_PROFILE_VALUE < p1.y < params.MAX_PROFILE_VALUE):
                continue
            cuts.append((p0, p1, n, int(global_idx), local_idx))

    upper_teeth_min_y = params.EXTREME_PROFILE_VALUE
    lower_teeth_max_y = -params.EXTREME_PROFILE_VALUE
    upper_teeth_max_y = -params.EXTREME_PROFILE_VALUE
    lower_teeth_min_y = params.EXTREME_PROFILE_VALUE
    upper_lip_min_y = params.EXTREME_PROFILE_VALUE
    lower_lip_max_y = -params.EXTREME_PROFILE_VALUE

    for p0, p1, _n, global_idx, _local_idx in cuts:
        if global_idx == int(params.SurfaceIndex.UPPER_TEETH):
            upper_teeth_min_y = min(upper_teeth_min_y, p0.y, p1.y)
            upper_teeth_max_y = max(upper_teeth_max_y, p0.y, p1.y)
        elif global_idx == int(params.SurfaceIndex.LOWER_TEETH):
            lower_teeth_max_y = max(lower_teeth_max_y, p0.y, p1.y)
            lower_teeth_min_y = min(lower_teeth_min_y, p0.y, p1.y)
        elif global_idx == int(params.SurfaceIndex.UPPER_LIP):
            upper_lip_min_y = min(upper_lip_min_y, p0.y, p1.y)
        elif global_idx == int(params.SurfaceIndex.LOWER_LIP):
            lower_lip_max_y = max(lower_lip_max_y, p0.y, p1.y)

    both_teeth_cut = upper_teeth_max_y > lower_teeth_min_y
    both_lips_cut = (upper_lip_min_y != params.EXTREME_PROFILE_VALUE) and (lower_lip_max_y != -params.EXTREME_PROFILE_VALUE)

    for p0, p1, n, global_idx, _local_idx in cuts:
        if global_idx in (
            int(params.SurfaceIndex.UPPER_COVER),
            int(params.SurfaceIndex.UPPER_TEETH),
            int(params.SurfaceIndex.UPPER_LIP),
            int(params.SurfaceIndex.UVULA),
        ):
            if n.y < 0.0:
                _insert_upper_profile_line(p0, p1, global_idx, upper_profile, upper_profile_surface)
            continue

        if global_idx == int(params.SurfaceIndex.LOWER_COVER):
            if n.y > 0.0:
                _insert_lower_cover_profile_line(
                    p0,
                    p1,
                    global_idx,
                    upper_profile,
                    upper_profile_surface,
                    lower_profile,
                    lower_profile_surface,
                )
            continue

        if global_idx in (
            int(params.SurfaceIndex.LOWER_TEETH),
            int(params.SurfaceIndex.LOWER_LIP),
            int(params.SurfaceIndex.EPIGLOTTIS),
        ):
            if n.y > 0.0:
                _insert_lower_profile_line(p0, p1, global_idx, lower_profile, lower_profile_surface)
            continue

        right_orientation = not (normal.x < 0.0 and normal.y < 0.0)
        if (n.x * n.x + n.y * n.y > min_squared_normal_length) and right_orientation:
            if n.y <= 0.0:
                _insert_upper_profile_line(p0, p1, global_idx, upper_profile, upper_profile_surface)
            else:
                _insert_lower_profile_line(p0, p1, global_idx, lower_profile, lower_profile_surface)

    for sample_idx in range(half_sample_count):
        if upper_profile[sample_idx] == params.EXTREME_PROFILE_VALUE:
            upper_profile[sample_idx] = invalid
        if lower_profile[sample_idx] == -params.EXTREME_PROFILE_VALUE:
            lower_profile[sample_idx] = invalid

    upper_limit = 1_000_000.0
    lower_limit = -1_000_000.0
    for sample_idx in range(half_sample_count - 1, -1, -1):
        if upper_profile[sample_idx] != invalid:
            if upper_profile_surface[sample_idx] == int(params.SurfaceIndex.UPPER_TEETH):
                upper_limit = min(upper_limit, upper_profile[sample_idx])
            if upper_profile[sample_idx] > upper_limit:
                upper_profile[sample_idx] = upper_limit

        if lower_profile[sample_idx] != invalid:
            if lower_profile_surface[sample_idx] == int(params.SurfaceIndex.LOWER_TEETH):
                lower_limit = max(lower_limit, lower_profile[sample_idx])
            if lower_profile[sample_idx] < lower_limit:
                lower_profile[sample_idx] = lower_limit

    last_value = invalid
    for sample_idx in range(half_sample_count):
        if upper_profile[sample_idx] == invalid:
            upper_profile[sample_idx] = last_value
        else:
            last_value = upper_profile[sample_idx]

    first_upper_valid = -1
    first_lower_valid = -1
    for sample_idx in range(half_sample_count):
        if first_lower_valid == -1 and lower_profile[sample_idx] != invalid:
            first_lower_valid = sample_idx
        if first_upper_valid == -1 and upper_profile[sample_idx] != invalid:
            first_upper_valid = sample_idx

    if first_lower_valid == -1:
        if first_upper_valid == -1:
            first_upper_valid = 0
        for sample_idx in range(first_upper_valid, half_sample_count):
            lower_profile[sample_idx] = params.MIN_PROFILE_VALUE
    else:
        sample_idx = half_sample_count - 1
        while sample_idx > 0 and lower_profile[sample_idx] == invalid:
            lower_profile[sample_idx] = params.MIN_PROFILE_VALUE
            sample_idx -= 1

    left = 0
    right = 0
    for sample_idx in range(half_sample_count):
        if lower_profile[sample_idx] == invalid:
            if sample_idx <= left or sample_idx >= right:
                left = sample_idx
                right = sample_idx
                while left > 0 and lower_profile[left] == invalid:
                    left -= 1
                while right < half_sample_count - 1 and lower_profile[right] == invalid:
                    right += 1

            if left < sample_idx < right:
                if lower_profile[left] != invalid and lower_profile[right] != invalid:
                    lower_profile[sample_idx] = lower_profile[left] + (
                        (sample_idx - left) * (lower_profile[right] - lower_profile[left]) / (right - left)
                    )
                elif lower_profile[left] != invalid:
                    lower_profile[sample_idx] = lower_profile[left]

    for sample_idx in range(half_sample_count):
        if lower_profile[sample_idx] != invalid and upper_profile[sample_idx] != invalid and lower_profile[sample_idx] > upper_profile[sample_idx]:
            lower_profile[sample_idx] = upper_profile[sample_idx]

    upper_left = 0
    while upper_left < half_sample_count - 1 and upper_profile[upper_left] == invalid:
        upper_left += 1

    lower_left = 0
    while lower_left < half_sample_count - 1 and lower_profile[lower_left] == invalid:
        lower_left += 1

    leftmost = half_sample_count - 1
    if upper_profile[upper_left] != invalid and upper_left < leftmost:
        leftmost = upper_left
    if lower_profile[lower_left] != invalid and lower_left < leftmost:
        leftmost = lower_left

    if upper_left < lower_left:
        for sample_idx in range(upper_left, lower_left):
            lower_profile[sample_idx] = upper_profile[sample_idx]
    else:
        for sample_idx in range(lower_left, upper_left):
            upper_profile[sample_idx] = lower_profile[sample_idx]

    if both_lips_cut and (not both_teeth_cut) and upper_profile[half_sample_count - 1] == lower_profile[half_sample_count - 1]:
        for sample_idx in range(half_sample_count):
            if upper_profile[sample_idx] != invalid and lower_profile[sample_idx] != invalid:
                merged_value = 0.5 * (upper_profile[sample_idx] + lower_profile[sample_idx])
                upper_profile[sample_idx] = merged_value
                lower_profile[sample_idx] = merged_value

    for sample_idx in range(half_sample_count + 1):
        mirror_idx = sample_count - 1 - sample_idx
        upper_profile[mirror_idx] = upper_profile[sample_idx]
        lower_profile[mirror_idx] = lower_profile[sample_idx]

    rightmost = sample_count - 1 - leftmost

    if consider_tongue:
        tongue_profile = [-params.EXTREME_PROFILE_VALUE] * sample_count
        tongue_profile_surface = [-1] * sample_count
        tongue_surface = vt.surfaces_list[int(params.SurfaceIndex.TONGUE)]

        if tongue_surface is not None:
            if not getattr(tongue_surface, "_dvts_intersections_prepared", False):
                tongue_surface.prepare_intersections()
                setattr(tongue_surface, "_dvts_intersections_prepared", True)
            tongue_surface.prepare_intersection(point, normal)
            triangle_indices = tongue_surface.get_triangle_list(max_entries=max_list_entries)
            for triangle_idx in triangle_indices:
                has_intersection, p0, p1, n = tongue_surface.get_triangle_intersection(triangle_idx, point, normal)
                if has_intersection and n.y >= 0.0:
                    _insert_lower_profile_line(
                        p0,
                        p1,
                        int(params.SurfaceIndex.TONGUE),
                        tongue_profile,
                        tongue_profile_surface,
                    )

        for sample_idx in range(sample_count):
            if tongue_profile[sample_idx] == -params.EXTREME_PROFILE_VALUE:
                tongue_profile[sample_idx] = invalid

        max_range_cm = 1.5
        highest_sample = max([value for value in tongue_profile if value != invalid], default=-1_000_000.0)
        for sample_idx in range(sample_count):
            if tongue_profile[sample_idx] != invalid and tongue_profile[sample_idx] < highest_sample - max_range_cm:
                tongue_profile[sample_idx] = invalid

        _interpolate_invalid_samples(tongue_profile, invalid)

        for sample_idx in range(leftmost, rightmost + 1):
            if (
                tongue_profile[sample_idx] != invalid
                and lower_profile[sample_idx] != invalid
                and tongue_profile[sample_idx] > lower_profile[sample_idx]
            ):
                lower_profile[sample_idx] = tongue_profile[sample_idx]
                lower_profile_surface[sample_idx] = tongue_profile_surface[sample_idx]

        while leftmost < half_sample_count and lower_profile[leftmost] > upper_profile[leftmost] - 0.1:
            leftmost += 1
        while rightmost > half_sample_count and lower_profile[rightmost] > upper_profile[rightmost] - 0.1:
            rightmost -= 1

        for sample_idx in range(sample_count):
            if leftmost <= sample_idx <= rightmost:
                if lower_profile[sample_idx] > upper_profile[sample_idx]:
                    lower_profile[sample_idx] = upper_profile[sample_idx]
            else:
                lower_profile[sample_idx] = invalid
                upper_profile[sample_idx] = invalid

    open_threshold = 0.1
    cutoff_threshold = 0.01

    sample_idx = half_sample_count - 1
    while sample_idx > 0 and upper_profile[sample_idx] - lower_profile[sample_idx] < open_threshold:
        sample_idx -= 1
    while sample_idx > 0 and upper_profile[sample_idx] - lower_profile[sample_idx] >= cutoff_threshold:
        sample_idx -= 1
    while sample_idx >= 0:
        upper_profile[sample_idx] = invalid
        lower_profile[sample_idx] = invalid
        sample_idx -= 1

    sample_idx = half_sample_count - 1
    while sample_idx < sample_count - 1 and upper_profile[sample_idx] - lower_profile[sample_idx] < open_threshold:
        sample_idx += 1
    while sample_idx < sample_count - 1 and upper_profile[sample_idx] - lower_profile[sample_idx] >= cutoff_threshold:
        sample_idx += 1
    while sample_idx <= sample_count - 1:
        upper_profile[sample_idx] = invalid
        lower_profile[sample_idx] = invalid
        sample_idx += 1

    leftmost = 0
    while leftmost < half_sample_count - 1 and (
        upper_profile[leftmost] == invalid or lower_profile[leftmost] == invalid
    ):
        leftmost += 1
    if upper_profile[leftmost] != invalid and lower_profile[leftmost] != invalid:
        merged_value = 0.5 * (upper_profile[leftmost] + lower_profile[leftmost])
        upper_profile[leftmost] = merged_value
        lower_profile[leftmost] = merged_value

    rightmost = sample_count - 1
    while rightmost > half_sample_count and (
        upper_profile[rightmost] == invalid or lower_profile[rightmost] == invalid
    ):
        rightmost -= 1
    if upper_profile[rightmost] != invalid and lower_profile[rightmost] != invalid:
        merged_value = 0.5 * (upper_profile[rightmost] + lower_profile[rightmost])
        upper_profile[rightmost] = merged_value
        lower_profile[rightmost] = merged_value

    articulator = Articulator.OTHER_ARTICULATOR
    has_tongue = False
    has_lower_lip = False
    has_lower_teeth = False

    check_range_cm = 0.5
    num_check_samples = int(check_range_cm / params.PROFILE_SAMPLE_LENGTH)
    for check_idx in range(1, max(2, num_check_samples)):
        sample_idx = half_sample_count - check_idx
        if sample_idx < 0:
            break
        profile_surface = lower_profile_surface[sample_idx]

        if profile_surface == int(params.SurfaceIndex.TONGUE):
            has_tongue = True
        elif profile_surface in (int(params.SurfaceIndex.LOWER_LIP), int(params.SurfaceIndex.RADIATION)):
            has_lower_lip = True
        elif profile_surface == int(params.SurfaceIndex.LOWER_TEETH):
            has_lower_teeth = True

    if has_lower_teeth:
        articulator = Articulator.LOWER_INCISORS
    elif has_lower_lip:
        articulator = Articulator.LOWER_LIP
    elif has_tongue:
        articulator = Articulator.TONGUE

    return upper_profile, lower_profile, articulator


def _get_center_line_pos(vt: "VocalTract", query_point: Point2D) -> float:
    center_line = getattr(vt, "center_line", None)
    if center_line is None or len(center_line) < 2:
        return 0.0

    epsilon = params.CENTERLINE_EPSILON
    min_dist = 1_000_000.0
    best_index = -1
    best_t = 0.0

    for index in range(len(center_line) - 1):
        p0 = Point2D(center_line[index].point.x, center_line[index].point.y)
        n0 = Point2D(center_line[index].normal.x, center_line[index].normal.y)
        p1 = Point2D(center_line[index + 1].point.x, center_line[index + 1].point.y)
        n1 = Point2D(center_line[index + 1].normal.x, center_line[index + 1].normal.y)

        separation = (p1 - p0) * epsilon
        p0 = p0 - separation
        p1 = p1 + separation

        condition_left = (query_point.x - p0.x) * n0.y - n0.x * (query_point.y - p0.y) >= 0.0
        condition_right = (query_point.x - p1.x) * n1.y - n1.x * (query_point.y - p1.y) <= 0.0
        if not (condition_left and condition_right):
            continue

        segment = p1 - p0
        relative = query_point - p0
        normal_delta = n1 - n0

        denominator = segment.y * normal_delta.x - segment.x * normal_delta.y
        if denominator == 0.0:
            denominator = epsilon

        p_term = (
            relative.x * normal_delta.y
            - relative.y * normal_delta.x
            + segment.y * n0.x
            - segment.x * n0.y
        ) / denominator
        q_term = (relative.x * n0.y - relative.y * n0.x) / denominator
        discriminant = 0.25 * p_term * p_term - q_term
        if discriminant < 0.0:
            discriminant = 0.0
        root = math.sqrt(discriminant)
        t0 = -0.5 * p_term + root
        t1 = -0.5 * p_term - root
        t = t0 if (-epsilon < t0 < 1.0 + epsilon) else t1

        projected = p0 + segment * t
        dist = (query_point - projected).magnitude()
        if dist < min_dist:
            min_dist = dist
            best_index = index
            best_t = t

    if best_index == -1:
        best_index = 0
        for index, center_point in enumerate(center_line):
            dist = (center_point.point - query_point).magnitude()
            if dist < min_dist:
                min_dist = dist
                best_index = index
                best_t = 0.0

    if best_index >= len(center_line) - 1:
        best_index = len(center_line) - 2
        best_t = 1.0

    left_pos = center_line[best_index].pos
    right_pos = center_line[best_index + 1].pos
    return left_pos + best_t * (right_pos - left_pos)


def _postprocess_cross_sections(vt: "VocalTract") -> None:
    cross_sections = getattr(vt, "cross_sections", None)
    center_line = getattr(vt, "center_line", None)
    if cross_sections is None or center_line is None or len(cross_sections) == 0 or len(center_line) == 0:
        vt.nasal_port_area_cm2 = 0.0
        vt.nasal_port_pos_cm = 0.0
        vt.incisor_pos_cm = 0.0
        vt.tongue_tip_side_elevation = 0.0
        return

    vt.nasal_port_area_cm2 = max(0.0, _get_param_value(vt, int(params.ParamIndex.VO)))
    vt.tongue_tip_side_elevation = _get_param_value(vt, int(params.ParamIndex.TS3))

    upper_cover = vt.surfaces_list[int(params.SurfaceIndex.UPPER_COVER)]
    port_rib = params.NUM_LARYNX_RIBS + params.NUM_PHARYNX_RIBS
    if upper_cover is not None and 0 <= port_rib < upper_cover.num_ribs and upper_cover.num_points > 0:
        port_point = upper_cover.get_vertex(port_rib, upper_cover.num_points // 2).to_point2d()
        vt.nasal_port_pos_cm = _get_center_line_pos(vt, port_point)
    else:
        vt.nasal_port_pos_cm = 0.0

    incisor_pos_cm = float(getattr(vt, "incisor_pos_cm", 0.0))
    upper_teeth = vt.surfaces_list[int(params.SurfaceIndex.UPPER_TEETH)]
    if upper_teeth is not None and upper_teeth.num_ribs >= params.NUM_TEETH_RIBS and upper_teeth.num_points > 2:
        teeth_x = upper_teeth.get_vertex(params.NUM_TEETH_RIBS - 1, 2).x
        for index in range(len(center_line) - 1):
            point_left = center_line[index].point
            point_right = center_line[index + 1].point
            if point_left.x < teeth_x <= point_right.x:
                delta_x = point_right.x - point_left.x
                if delta_x < params.CENTERLINE_EPSILON:
                    delta_x = params.CENTERLINE_EPSILON
                incisor_pos_cm = center_line[index].pos + ((teeth_x - point_left.x) / delta_x) * (
                    point_right - point_left
                ).magnitude()
    vt.incisor_pos_cm = incisor_pos_cm

    tongue_tip_radius = float(vt.anatomy.tongue_tip_radius_cm)
    tongue_tip_center = Point2D(
        _get_param_value(vt, int(params.ParamIndex.TTX)),
        _get_param_value(vt, int(params.ParamIndex.TTY)),
    )
    tongue_center_x = _get_param_value(vt, int(params.ParamIndex.TCX))

    for index, cross_section in enumerate(cross_sections):
        point = center_line[index].point
        normal = center_line[index].normal
        relative = tongue_tip_center - point
        is_right_of_tongue_tip = (relative.y * normal.x - relative.x * normal.y > tongue_tip_radius) and (
            point.x > tongue_center_x
        )
        if is_right_of_tongue_tip and cross_section.articulator == int(Articulator.TONGUE):
            cross_section.articulator = int(Articulator.OTHER_ARTICULATOR)

    min_incisor_circ = 2.0 * math.sqrt(params.MIN_INCISOR_AREA_CM2 * math.pi)
    for cross_section in cross_sections:
        if (
            cross_section.pos >= vt.incisor_pos_cm - params.LEFT_INCISOR_MARGIN_CM
            and cross_section.pos <= vt.incisor_pos_cm + params.RIGHT_INCISOR_MARGIN_CM
            and cross_section.area < params.MIN_INCISOR_AREA_CM2
        ):
            cross_section.area = params.MIN_INCISOR_AREA_CM2
            cross_section.circ = min_incisor_circ

    ts2 = _get_param_value(vt, int(params.ParamIndex.TS2))
    ts3 = _get_param_value(vt, int(params.ParamIndex.TS3))
    min_area_back = params.tongue_side_param_to_min_area_cm2(ts2)
    if ts2 < 0.0:
        min_area_back = 0.0
    min_area_tip = params.tongue_side_param_to_min_area_cm2(ts3)
    min_circ_back = 2.0 * math.sqrt(min_area_back * math.pi)
    min_circ_tip = 2.0 * math.sqrt(min_area_tip * math.pi)

    tongue_tip_right_cm = 0.0
    lower_lip_left_cm = 1_000_000.0
    for cross_section in cross_sections:
        if cross_section.articulator == int(Articulator.TONGUE):
            tongue_tip_right_cm = cross_section.pos
        if cross_section.articulator == int(Articulator.LOWER_LIP) and cross_section.pos < lower_lip_left_cm:
            lower_lip_left_cm = cross_section.pos
    tongue_tip_left_cm = tongue_tip_right_cm - params.TONGUE_TIP_REGION_LENGTH_CM

    for cross_section in cross_sections:
        if cross_section.pos <= tongue_tip_left_cm:
            if cross_section.area < min_area_back:
                cross_section.area = min_area_back
            if cross_section.circ < min_circ_back:
                cross_section.circ = min_circ_back

        if tongue_tip_left_cm <= cross_section.pos <= lower_lip_left_cm:
            if cross_section.area < min_area_tip:
                cross_section.area = min_area_tip
            if cross_section.circ < min_circ_tip:
                cross_section.circ = min_circ_tip


def calc_cross_sections(vt: "VocalTract") -> None:
    invalid = params.INVALID_PROFILE_SAMPLE
    delta_squared = params.PROFILE_SAMPLE_LENGTH * params.PROFILE_SAMPLE_LENGTH
    cross_sections: List[params.CrossSection] = []

    if not hasattr(vt, "center_line") or len(vt.center_line) == 0:
        vt.cross_sections = cross_sections
        return

    for idx, center_line_point in enumerate(vt.center_line):
        upper_profile, lower_profile, articulator = get_cross_profiles(
            vt=vt,
            point=center_line_point.point,
            normal=center_line_point.normal,
            consider_tongue=True,
        )

        area = 0.0
        circ = 0.0

        for sample_idx in range(params.NUM_PROFILE_SAMPLES - 1):
            upper0 = upper_profile[sample_idx]
            upper1 = upper_profile[sample_idx + 1]
            lower0 = lower_profile[sample_idx]
            lower1 = lower_profile[sample_idx + 1]

            if upper0 == invalid or upper1 == invalid or lower0 == invalid or lower1 == invalid:
                continue

            section_a = upper0 - lower0
            section_b = upper1 - lower1
            delta_area = 0.5 * (section_a + section_b) * params.PROFILE_SAMPLE_LENGTH
            area += delta_area

            upper_delta = upper1 - upper0
            lower_delta = lower1 - lower0
            circ += math.sqrt(upper_delta * upper_delta + delta_squared)
            circ += math.sqrt(lower_delta * lower_delta + delta_squared)

        cross_sections.append(
            params.CrossSection(
                area=max(0.0, area),
                circ=max(0.0, circ),
                pos=vt.center_line[idx].pos,
                articulator=int(articulator),
            )
        )

    vt.cross_sections = cross_sections
    _postprocess_cross_sections(vt)


def calculate_all(vt: "VocalTract") -> None:
    calc_surfaces(vt)
    calc_center_line(vt)
    calc_cross_sections(vt)
