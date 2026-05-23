# -*- coding: utf-8 -*-
"""
声道几何核心类 (Vocal Tract Geometry)

该模块保留历史兼容布局，但当前服务默认使用项目自有 Python 运行时。

核心功能:
- 存储声道解剖结构 (Anatomy)
- 存储声道形状 (Shape)
- 计算声道表面 (Surfaces)
- 计算中线 (CenterLine)
- 计算横截面积 (CrossSections)
"""

import math
import xml.etree.ElementTree as ET
from typing import List, Optional, Dict, Tuple

from . import params
from . import anatomy_geometry
from . import dynamic_geometry
from ..utils.geometry import Point2D, Point3D, Vector2D, Line2D
from ..utils.splines import Spline3D, LineStrip3D, LineStrip2D, MAX_SPLINE_POINTS

# =============================================================================
# Default Anatomy XML
# =============================================================================

DEFAULT_ANATOMY_XML = """
<anatomy>
  <!--****************************************************************************-->
  <palate>
    <p0 x="0.2" z="-2.3" teeth_height="0.5" top_teeth_width="1.05" bottom_teeth_width="1.05" palate_height="1.3" palate_angle_deg="39.5"/>
    <p1 x="0.9" z="-2.2" teeth_height="0.5" top_teeth_width="1.05" bottom_teeth_width="1.05" palate_height="1.15" palate_angle_deg="39.5"/>
    <p2 x="1.8" z="-2.0" teeth_height="0.5" top_teeth_width="1.0" bottom_teeth_width="1.0" palate_height="1.425" palate_angle_deg="60.8"/>
    <p3 x="2.8" z="-1.8" teeth_height="0.5" top_teeth_width="1.0" bottom_teeth_width="1.0" palate_height="1.6" palate_angle_deg="60.8"/>
    <p4 x="3.5" z="-1.6" teeth_height="0.6" top_teeth_width="0.8" bottom_teeth_width="0.8" palate_height="1.4" palate_angle_deg="60.8"/>
    <p5 x="4.15" z="-1.4" teeth_height="0.7" top_teeth_width="0.7" bottom_teeth_width="0.7" palate_height="0.7" palate_angle_deg="38.0"/>
    <p6 x="4.55" z="-1.1" teeth_height="0.8" top_teeth_width="0.65" bottom_teeth_width="0.3" palate_height="0.15" palate_angle_deg="23.4"/>
    <p7 x="4.7" z="-0.6" teeth_height="0.8" top_teeth_width="0.8" bottom_teeth_width="0.2" palate_height="0.0" palate_angle_deg="0.0"/>
    <p8 x="4.7" z="0.0" teeth_height="0.8" top_teeth_width="0.85" bottom_teeth_width="0.2" palate_height="0.0" palate_angle_deg="0.0"/>
  </palate>
  <!--****************************************************************************-->
  <jaw fulcrum_x="-6.5" fulcrum_y="2.0" rest_pos_x="0.0" rest_pos_y="-1.2" tooth_root_length="0.8">
    <p0 x="0.2" z="-2.3" teeth_height="0.5"  top_teeth_width="1.05" bottom_teeth_width="1.05" jaw_height="1.5" jaw_angle_deg="69.5"/>
    <p1 x="1.2" z="-2.2" teeth_height="0.5"  top_teeth_width="1.1" bottom_teeth_width="1.1" jaw_height="1.5" jaw_angle_deg="69.5"/>
    <p2 x="2.2" z="-1.9" teeth_height="0.5"  top_teeth_width="1.05" bottom_teeth_width="1.05" jaw_height="1.5" jaw_angle_deg="69.5"/>
    <p3 x="3.2" z="-1.6" teeth_height="0.5"  top_teeth_width="0.9" bottom_teeth_width="0.9" jaw_height="1.5" jaw_angle_deg="69.5"/>
    <p4 x="3.9" z="-1.4" teeth_height="0.5"  top_teeth_width="0.75" bottom_teeth_width="0.75" jaw_height="1.0" jaw_angle_deg="42.2"/>
    <p5 x="4.3" z="-1.1" teeth_height="0.55"  top_teeth_width="0.6" bottom_teeth_width="0.7" jaw_height="0.4" jaw_angle_deg="35.8"/>
    <p6 x="4.5" z="-0.7" teeth_height="0.6"  top_teeth_width="0.3" bottom_teeth_width="0.8" jaw_height="0.13" jaw_angle_deg="31.4"/>
    <p7 x="4.55" z="-0.5" teeth_height="0.7"  top_teeth_width="0.2" bottom_teeth_width="0.9" jaw_height="0.0" jaw_angle_deg="0.0"/>
    <p8 x="4.55" z="0.0" teeth_height="0.7"  top_teeth_width="0.2" bottom_teeth_width="0.9" jaw_height="0.0" jaw_angle_deg="0.0"/>
  </jaw>
  <!--****************************************************************************-->
  <tongue>
    <tip radius="0.2000"/>
    <body radius_x="1.8000" radius_y="1.8000"/>
    <root automatic_calc="1" trx_slope="0.938" trx_intercept="-5.11" try_slope="0.831" try_intercept="-3.03"/>
    <tonsil length="2.2" height="0.5"/>
  </tongue>
  <!--****************************************************************************-->
  <lips width="1.3"/>
  <!--****************************************************************************-->
  <velum velum_angle_deg="50.0" uvula_width="0.7" uvula_height="0.9" uvula_depth="0.7">
    <low points="-1.25 -0.7 -1 -0.3 -0.7 0 -0.3 0.35 0 0.55 "/>
    <mid points="-1.75 0 -1.55 0.3 -1.15 0.5 -0.6 0.7 0 0.87 "/>
    <high points="-1.75 0.5 -1.5 0.95 -1.1 1.2 -0.6 1.25 0 1.3 "/>
  </velum>
  <!--****************************************************************************-->
  <pharynx fulcrum_x="-2.372" fulcrum_y="1.214" rotation_angle_deg="-98.0" top_rib_y="-1.4" upper_depth="3.8" lower_depth="3.4" back_side_width="1.5"/>
  <!--****************************************************************************-->
  <larynx upper_depth="1.0" lower_depth="1.0" epiglottis_width="0.5" epiglottis_height="1.6" epiglottis_depth="1.4" epiglottis_angle_deg="100.0">
    <narrow points="1.8 0 1.05 -0.2 1.55 -1.2 2.68 -3.2 1.48 -3.2 1.1 -1.2 0 -1 0 0 "/>
    <wide points="3.3 0 2.13 -0.2 1.9 -1.2 2.68 -3.2 1.48 -3.2 1.45 -1.2 0 -1 0 0 "/>  
  </larynx>
  <vocal_folds default_f0 = "110"/>
  <piriform_fossa length = "2.5" volume = "1.5"/>
  <subglottal_cavity length="23.0"/>
  <nasal_cavity length="11.4"/>
  <!--****************************************************************************-->
  <param index="0"  name="HX"   min="0.0"   max="1.0"   neutral="1.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="1"  name="HY"   min="-6.0"  max="-3.5"  neutral="-4.75"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="2"  name="JX"   min="-0.5"  max="0.0"   neutral="0.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="3"  name="JA"   min="-7.0"  max="0.0"   neutral="-2.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="4"  name="LP"   min="-1.0"  max="1.0"   neutral="-0.07"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="5"  name="LD"   min="-2.0"  max="4.0"   neutral="0.95"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="6"  name="VS"  min="0.0"   max="1.0"   neutral="0.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="7"  name="VO"  min="-0.1"   max="1.5"   neutral="-0.1"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="8"  name="TCX"  min="-3.0"  max="4.0"   neutral="-0.4"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="9"  name="TCY"  min="-3.0"  max="1.0"   neutral="-1.46"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="10" name="TTX"   min="1.5"   max="5.5"   neutral="3.5"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="11" name="TTY"  min="-3.0"  max="2.5"   neutral="-1.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="12" name="TBX"  min="-3.0"  max="4.0"   neutral="2.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="13" name="TBY"  min="-3.0"  max="5.0"   neutral="0.5"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="14" name="TRX"  min="-4.0"  max="2.0"   neutral="0.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="15" name="TRY"  min="-6.0"  max="0.0"   neutral="0.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="16" name="TS1"  min="0.0" max="1.0"  neutral="0.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="17" name="TS2"  min="0.0" max="1.0"  neutral="0.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
  <param index="18" name="TS3"  min="-1.0" max="1.0"  neutral="0.0"   positive_velocity_factor="1.0"   negative_velocity_factor="1.0"/>
</anatomy>
"""

# =============================================================================
# Helper for Surface
# =============================================================================

class Surface:
    MAX_TILES_X = 15
    MAX_TILES_Y = 15
    STANDARD_TILE_WIDTH = 0.5
    STANDARD_TILE_HEIGHT = 0.5

    def __init__(self, num_ribs: int = 0, num_points: int = 0):
        self.num_ribs = num_ribs
        self.num_points = num_points
        self.vertices: List[List[Point3D]] = []
        self.left_border = 0.0
        self.right_border = 0.0
        self.bottom_border = 0.0
        self.top_border = 0.0
        self.num_tiles_x = 1
        self.num_tiles_y = 1
        self.tile_width = self.STANDARD_TILE_WIDTH
        self.tile_height = self.STANDARD_TILE_HEIGHT
        self._tile_triangles: List[List[List[int]]] = [[[]]]
        if num_ribs > 0 and num_points > 0:
            self.init(num_ribs, num_points)
            
    def init(self, num_ribs: int, num_points: int):
        self.num_ribs = num_ribs
        self.num_points = num_points
        # self.vertices[num_ribs][num_points]
        self.vertices = [[Point3D() for _ in range(num_points)] for _ in range(num_ribs)]
        self.crease_angle_deg = 0.0
        self._build_topology()

    def _build_topology(self) -> None:
        if self.num_ribs <= 0 or self.num_points <= 0:
            self._edge_vertices: List[Tuple[int, int]] = []
            self._triangle_vertices: List[Tuple[int, int, int]] = []
            self._triangle_edges: List[Tuple[int, int, int]] = []
            return

        edge_vertices: List[Tuple[int, int]] = []
        triangle_vertices: List[Tuple[int, int, int]] = []
        triangle_edges: List[Tuple[int, int, int]] = []

        # Horizontal edges.
        for point in range(self.num_points):
            for rib in range(self.num_ribs - 1):
                edge_vertices.append(
                    (
                        rib * self.num_points + point,
                        (rib + 1) * self.num_points + point,
                    )
                )

        # Vertical edges.
        for rib in range(self.num_ribs):
            for point in range(self.num_points - 1):
                edge_vertices.append(
                    (
                        rib * self.num_points + point,
                        rib * self.num_points + point + 1,
                    )
                )

        # Diagonal edges.
        for point in range(self.num_points - 1):
            for rib in range(self.num_ribs - 1):
                edge_vertices.append(
                    (
                        rib * self.num_points + point,
                        (rib + 1) * self.num_points + point + 1,
                    )
                )

        a_offset = self.num_points * (self.num_ribs - 1)
        b_offset = a_offset + (self.num_points - 1) * self.num_ribs

        for rib in range(self.num_ribs - 1):
            for point in range(self.num_points - 1):
                triangle_vertices.append(
                    (
                        rib * self.num_points + point,
                        (rib + 1) * self.num_points + point + 1,
                        rib * self.num_points + point + 1,
                    )
                )
                triangle_edges.append(
                    (
                        b_offset + point * (self.num_ribs - 1) + rib,
                        (point + 1) * (self.num_ribs - 1) + rib,
                        a_offset + rib * (self.num_points - 1) + point,
                    )
                )

                triangle_vertices.append(
                    (
                        rib * self.num_points + point,
                        (rib + 1) * self.num_points + point,
                        (rib + 1) * self.num_points + point + 1,
                    )
                )
                triangle_edges.append(
                    (
                        b_offset + point * (self.num_ribs - 1) + rib,
                        point * (self.num_ribs - 1) + rib,
                        a_offset + (rib + 1) * (self.num_points - 1) + point,
                    )
                )

        self._edge_vertices = edge_vertices
        self._triangle_vertices = triangle_vertices
        self._triangle_edges = triangle_edges
        self.prepare_intersections()

    def _vertex_by_linear_index(self, vertex_index: int) -> Point3D:
        rib = vertex_index // self.num_points
        point = vertex_index % self.num_points
        return self.get_vertex(rib, point)
        
    def set_vertex(self, rib: int, point: int, p: Point3D):
        if 0 <= rib < self.num_ribs and 0 <= point < self.num_points:
            self.vertices[rib][point] = Point3D(p.x, p.y, p.z)
            
    def get_vertex(self, rib: int, point: int) -> Point3D:
        if 0 <= rib < self.num_ribs and 0 <= point < self.num_points:
            vertex = self.vertices[rib][point]
            return Point3D(vertex.x, vertex.y, vertex.z)
        return Point3D()

    def get_triangle_count(self) -> int:
        if self.num_ribs < 2 or self.num_points < 2:
            return 0
        return (self.num_ribs - 1) * (self.num_points - 1) * 2

    def get_triangle(self, index: int) -> Tuple[Point3D, Point3D, Point3D]:
        if index < 0 or index >= len(self._triangle_vertices):
            return (Point3D(), Point3D(), Point3D())
        v0, v1, v2 = self._triangle_vertices[index]
        return (
            self._vertex_by_linear_index(v0),
            self._vertex_by_linear_index(v1),
            self._vertex_by_linear_index(v2),
        )

    def prepare_intersections(self):
        extreme = 1_000_000.0
        epsilon = 0.1
        triangle_count = len(getattr(self, "_triangle_vertices", []))

        if self.num_ribs <= 0 or self.num_points <= 0 or triangle_count <= 0:
            self._triangle_index_cache = []
            self._tile_triangles = [[[]]]
            self.num_tiles_x = 1
            self.num_tiles_y = 1
            self.tile_width = self.STANDARD_TILE_WIDTH
            self.tile_height = self.STANDARD_TILE_HEIGHT
            self.left_border = 0.0
            self.right_border = 0.0
            self.bottom_border = 0.0
            self.top_border = 0.0
            self._vertex_side_cache = {}
            self._edge_intersection_cache = {}
            self._line_vector = Point2D(0.0, 1.0)
            self._line_point = Point2D(0.0, 0.0)
            self._left_line_point = Point2D(0.0, 0.0)
            self._right_line_point = Point2D(0.0, 0.0)
            return

        # Surface bounding box
        left_border = extreme
        right_border = -extreme
        bottom_border = extreme
        top_border = -extreme
        for rib in range(self.num_ribs):
            for point in range(self.num_points):
                vertex = self.vertices[rib][point]
                if vertex.x < left_border:
                    left_border = vertex.x
                if vertex.x > right_border:
                    right_border = vertex.x
                if vertex.y < bottom_border:
                    bottom_border = vertex.y
                if vertex.y > top_border:
                    top_border = vertex.y

        left_border -= epsilon
        bottom_border -= epsilon
        right_border += epsilon
        top_border += epsilon

        num_tiles_x = int((right_border - left_border) / self.STANDARD_TILE_WIDTH) + 1
        num_tiles_y = int((top_border - bottom_border) / self.STANDARD_TILE_HEIGHT) + 1
        num_tiles_x = max(1, min(self.MAX_TILES_X, num_tiles_x))
        num_tiles_y = max(1, min(self.MAX_TILES_Y, num_tiles_y))
        tile_width = (right_border - left_border) / float(num_tiles_x)
        tile_height = (top_border - bottom_border) / float(num_tiles_y)

        tile_triangles: List[List[List[int]]] = [[[] for _ in range(num_tiles_y)] for _ in range(num_tiles_x)]

        # Assign triangles to all overlapping tiles (duplicates allowed like C++).
        for triangle_idx, (v0_idx, v1_idx, v2_idx) in enumerate(self._triangle_vertices):
            p0 = self._vertex_by_linear_index(v0_idx)
            p1 = self._vertex_by_linear_index(v1_idx)
            p2 = self._vertex_by_linear_index(v2_idx)

            min_x = min(p0.x, p1.x, p2.x) - epsilon
            max_x = max(p0.x, p1.x, p2.x) + epsilon
            min_y = min(p0.y, p1.y, p2.y) - epsilon
            max_y = max(p0.y, p1.y, p2.y) + epsilon

            left_tile = int((min_x - left_border) / tile_width)
            right_tile = int((max_x - left_border) / tile_width)
            bottom_tile = int((min_y - bottom_border) / tile_height)
            top_tile_idx = int((max_y - bottom_border) / tile_height)

            left_tile = max(0, min(num_tiles_x - 1, left_tile))
            right_tile = max(0, min(num_tiles_x - 1, right_tile))
            bottom_tile = max(0, min(num_tiles_y - 1, bottom_tile))
            top_tile_idx = max(0, min(num_tiles_y - 1, top_tile_idx))

            for tile_x in range(left_tile, right_tile + 1):
                for tile_y in range(bottom_tile, top_tile_idx + 1):
                    tile_triangles[tile_x][tile_y].append(triangle_idx)

        self.left_border = left_border
        self.right_border = right_border
        self.bottom_border = bottom_border
        self.top_border = top_border
        self.num_tiles_x = num_tiles_x
        self.num_tiles_y = num_tiles_y
        self.tile_width = tile_width
        self.tile_height = tile_height
        self._tile_triangles = tile_triangles
        self._triangle_index_cache = list(range(triangle_count))
        self._vertex_side_cache: Dict[int, int] = {}
        self._edge_intersection_cache: Dict[int, Tuple[bool, Point2D]] = {}
        self._line_vector = Point2D(0.0, 1.0)
        self._line_point = Point2D(0.0, 0.0)
        self._left_line_point = Point2D(0.0, 0.0)
        self._right_line_point = Point2D(0.0, 0.0)
        
    def prepare_intersection(self, p: Point2D, v: Point2D):
        epsilon = 1.0e-6
        direction = Point2D(v.x, v.y)
        if direction.magnitude() < 1.0e-12:
            direction = Point2D(0.0, 1.0)
        direction.normalize()

        self._line_vector = direction
        self._line_point = Point2D(p.x, p.y)

        n = Point2D(-direction.y, direction.x)
        self._left_line_point = self._line_point + epsilon * n
        self._right_line_point = self._line_point - epsilon * n

        self._vertex_side_cache = {}
        self._edge_intersection_cache = {}
        
    def get_triangle_list(self, max_entries: int) -> List[int]:
        if not hasattr(self, "_tile_triangles") or self.num_tiles_x <= 0 or self.num_tiles_y <= 0:
            return list(self._triangle_index_cache)
        if max_entries <= 0:
            return []

        entries: List[int] = []
        line_point = Point2D(self._line_point.x, self._line_point.y)
        line_vector = Point2D(self._line_vector.x, self._line_vector.y)

        def push_tile(tile_x: int, tile_y: int) -> bool:
            if tile_x < 0 or tile_x >= self.num_tiles_x or tile_y < 0 or tile_y >= self.num_tiles_y:
                return True
            triangles = self._tile_triangles[tile_x][tile_y]
            if len(entries) + len(triangles) > max_entries:
                return False
            entries.extend(triangles)
            return True

        # Rather horizontal line
        if abs(line_vector.y) <= abs(line_vector.x) * self.tile_height / self.tile_width:
            if line_vector.x < 0.0:
                line_vector.x = -line_vector.x
                line_vector.y = -line_vector.y

            if abs(line_vector.x) < 1.0e-12:
                return list(self._triangle_index_cache[:max_entries])

            y_value = line_point.y + (self.left_border - line_point.x) * line_vector.y / line_vector.x
            delta_y = line_vector.y * self.tile_width / line_vector.x
            tile_y = int((y_value - self.bottom_border) / self.tile_height)

            if line_vector.y > 0.0:
                next_border_y = self.bottom_border + float(tile_y + 1) * self.tile_height
                for tile_x in range(self.num_tiles_x):
                    if not push_tile(tile_x, tile_y):
                        return entries
                    if y_value + delta_y > next_border_y:
                        tile_y += 1
                        next_border_y += self.tile_height
                        if not push_tile(tile_x, tile_y):
                            return entries
                    y_value += delta_y
            else:
                next_border_y = self.bottom_border + float(tile_y) * self.tile_height
                for tile_x in range(self.num_tiles_x):
                    if not push_tile(tile_x, tile_y):
                        return entries
                    if y_value + delta_y < next_border_y:
                        tile_y -= 1
                        next_border_y -= self.tile_height
                        if not push_tile(tile_x, tile_y):
                            return entries
                    y_value += delta_y
        # Rather vertical line
        else:
            if line_vector.y < 0.0:
                line_vector.x = -line_vector.x
                line_vector.y = -line_vector.y

            if abs(line_vector.y) < 1.0e-12:
                return list(self._triangle_index_cache[:max_entries])

            x_value = line_point.x + (self.bottom_border - line_point.y) * line_vector.x / line_vector.y
            delta_x = line_vector.x * self.tile_height / line_vector.y
            tile_x = int((x_value - self.left_border) / self.tile_width)

            if line_vector.x > 0.0:
                next_border_x = self.left_border + float(tile_x + 1) * self.tile_width
                for tile_y in range(self.num_tiles_y):
                    if not push_tile(tile_x, tile_y):
                        return entries
                    if x_value + delta_x > next_border_x:
                        tile_x += 1
                        next_border_x += self.tile_width
                        if not push_tile(tile_x, tile_y):
                            return entries
                    x_value += delta_x
            else:
                next_border_x = self.left_border + float(tile_x) * self.tile_width
                for tile_y in range(self.num_tiles_y):
                    if not push_tile(tile_x, tile_y):
                        return entries
                    if x_value + delta_x < next_border_x:
                        tile_x -= 1
                        next_border_x -= self.tile_width
                        if not push_tile(tile_x, tile_y):
                            return entries
                    x_value += delta_x

        return entries

    def get_triangle_intersection(self, index: int, p0_out: Point2D, p1_out: Point2D, n_out: Point2D) -> bool:
        ok, p0, p1, n = self._get_triangle_intersection_prepared(index)
        if not ok:
            return False

        p0_out.x, p0_out.y = p0.x, p0.y
        p1_out.x, p1_out.y = p1.x, p1.y
        n_out.x, n_out.y = n.x, n.y
        return True

    def _get_vertex_side(self, vertex_index: int) -> int:
        if vertex_index in self._vertex_side_cache:
            return self._vertex_side_cache[vertex_index]

        vertex = self._vertex_by_linear_index(vertex_index)
        side = 0

        dx = vertex.x - self._left_line_point.x
        dy = vertex.y - self._left_line_point.y
        d = dx * self._line_vector.y - dy * self._line_vector.x
        if d < 0.0:
            side = -1

        dx = vertex.x - self._right_line_point.x
        dy = vertex.y - self._right_line_point.y
        d = dx * self._line_vector.y - dy * self._line_vector.x
        if d > 0.0:
            side = 1

        self._vertex_side_cache[vertex_index] = side
        return side

    def _get_edge_intersection(self, edge_index: int) -> Tuple[bool, Point2D]:
        if edge_index in self._edge_intersection_cache:
            return self._edge_intersection_cache[edge_index]

        eps = 1.0e-6
        v0, v1 = self._edge_vertices[edge_index]
        side0 = self._get_vertex_side(v0)
        side1 = self._get_vertex_side(v1)

        is_intersected = False
        intersection = Point2D()
        if ((side0 >= 0 and side1 <= 0) or (side0 <= 0 and side1 >= 0)):
            p = self._vertex_by_linear_index(v0)
            u = self._vertex_by_linear_index(v1) - p
            r = Point3D(p.x - self._line_point.x, p.y - self._line_point.y, p.z)
            denominator = -u.x * self._line_vector.y + u.y * self._line_vector.x

            if denominator != 0.0:
                d = (-self._line_vector.x * r.y + self._line_vector.y * r.x) / denominator
                if -eps <= d < 1.0 + eps:
                    is_intersected = True
                    intersection.x = (
                        self._line_vector.x * (u.y * r.z - u.z * r.y)
                        + self._line_vector.y * (u.z * r.x - u.x * r.z)
                    ) / denominator
                    intersection.y = (-u.x * r.y + u.y * r.x) / denominator

        self._edge_intersection_cache[edge_index] = (is_intersected, intersection)
        return is_intersected, intersection

    def _get_triangle_intersection_prepared(self, index: int) -> Tuple[bool, Point2D, Point2D, Point2D]:
        if index < 0 or index >= len(self._triangle_edges):
            return False, Point2D(), Point2D(), Point2D()

        e0, e1, e2 = self._triangle_edges[index]
        intersections: List[Point2D] = []

        ok, point = self._get_edge_intersection(e0)
        if ok:
            intersections.append(point)
        ok, point = self._get_edge_intersection(e1)
        if ok:
            intersections.append(point)
        ok, point = self._get_edge_intersection(e2)
        if ok:
            intersections.append(point)

        if len(intersections) < 2:
            return False, Point2D(), Point2D(), Point2D()

        v0_index, v1_index, v2_index = self._triangle_vertices[index]
        p0 = self._vertex_by_linear_index(v0_index)
        p1 = self._vertex_by_linear_index(v1_index)
        p2 = self._vertex_by_linear_index(v2_index)

        edge1 = p1 - p0
        edge2 = p2 - p0
        normal = Point3D(
            edge1.y * edge2.z - edge1.z * edge2.y,
            edge1.z * edge2.x - edge1.x * edge2.z,
            edge1.x * edge2.y - edge1.y * edge2.x,
        )
        n_2d = Point2D(
            normal.z,
            normal.x * self._line_vector.x + normal.y * self._line_vector.y,
        )

        if len(intersections) == 2:
            return True, intersections[0], intersections[1], n_2d

        l0 = (intersections[0] - intersections[1]).magnitude()
        l1 = (intersections[1] - intersections[2]).magnitude()
        l2 = (intersections[2] - intersections[0]).magnitude()
        if l0 >= l1 and l0 >= l2:
            return True, intersections[0], intersections[1], n_2d
        if l1 >= l2:
            return True, intersections[1], intersections[2], n_2d
        return True, intersections[2], intersections[0], n_2d

    def get_triangle_intersection(self, index: int, plane_p: Point2D, plane_v: Point2D) -> Tuple[bool, Point2D, Point2D, Point2D]:
        direction = Point2D(plane_v.x, plane_v.y)
        if direction.magnitude() < 1.0e-12:
            return False, Point2D(), Point2D(), Point2D()
        direction.normalize()

        need_prepare = (
            abs(self._line_point.x - plane_p.x) > 1.0e-12
            or abs(self._line_point.y - plane_p.y) > 1.0e-12
            or abs(self._line_vector.x - direction.x) > 1.0e-12
            or abs(self._line_vector.y - direction.y) > 1.0e-12
        )
        if need_prepare:
            self.prepare_intersection(plane_p, direction)

        return self._get_triangle_intersection_prepared(index)

    def get_triangle_intersection_vtl(self, index: int, plane_p: Point2D, plane_v: Point2D) -> Tuple[bool, Point2D, Point2D, Point2D]:
        """Legacy compatibility shim for older VTL-flavored call sites."""
        return self.get_triangle_intersection(index, plane_p, plane_v)

    def swap_triangle_orientation(self):
        if not hasattr(self, "_triangle_vertices"):
            return
        swapped: List[Tuple[int, int, int]] = []
        for v0, v1, v2 in self._triangle_vertices:
            swapped.append((v2, v1, v0))
        self._triangle_vertices = swapped



# =============================================================================
# VocalTract Class
# =============================================================================

class VocalTract:
    def __init__(self):
        self.anatomy = params.Anatomy()
        self.params: List[params.ParamDef] = [params.ParamDef() for _ in range(params.ParamIndex.NUM_PARAMS)]
        
        # Surfaces (stored in dict or list? C++ uses array with enum index)
        # Using dict for readability or list for performance/matching C++
        self.surfaces: Dict[int, Surface] = {} 
        # Using list is better to match enum
        self.surfaces_list: List['Surface'] = [None] * params.NUM_SURFACES
        
        # Geometry - Outlines (2D Midsagittal)
        self.upper_outline = LineStrip2D()
        self.lower_outline = LineStrip2D()
        self.tongue_outline = LineStrip2D()
        self.epiglottis_outline = LineStrip2D()

        # Added for compatibility with anatomy_geometry
        self.upper_gums_inner_edge = [Point3D() for _ in range(params.NUM_JAW_RIBS)]
        self.upper_gums_outer_edge = [Point3D() for _ in range(params.NUM_JAW_RIBS)]
        self.lower_gums_inner_edge_orig = [Point3D() for _ in range(params.NUM_JAW_RIBS)]
        self.lower_gums_outer_edge_orig = [Point3D() for _ in range(params.NUM_JAW_RIBS)]
        self.lower_gums_inner_edge = [Point3D() for _ in range(params.NUM_JAW_RIBS)]
        self.lower_gums_outer_edge = [Point3D() for _ in range(params.NUM_JAW_RIBS)]
        
        self.wide_lip_corner_path = LineStrip3D()
        self.narrow_lip_corner_path = LineStrip3D()
        self.lip_corner_path = LineStrip3D()

        self.center_line_length = 0.0
        self.nasal_port_area_cm2 = 0.0
        self.nasal_port_pos_cm = 0.0
        self.incisor_pos_cm = 0.0
        self.tongue_tip_side_elevation = 0.0
        
        self.init()

    def init(self):
        # 1. Init Surfaces
        self._init_surface_grids()
        
        # 2. Read Anatomy
        root = ET.fromstring(DEFAULT_ANATOMY_XML)
        self.read_anatomy_xml(root)
        
        # 3. Init Control Params Names
        self._init_param_names()
        
        # 4. Calculate All (Initial)
        self.calculate_all()
        
        # 5. Store control params (reset)
        self._neutral_param_values = [float(p.value) for p in self.params]



    def read_anatomy_xml(self, root: ET.Element):
        # Parse XML and populate self.anatomy and self.params
        
        # Palate
        palate_node = root.find("palate")
        if palate_node is not None:
            self.anatomy.palate_points = []
            self.anatomy.palate_height_cm = []
            self.anatomy.palate_angle_deg = []
            self.anatomy.upper_teeth_height_cm = []
            self.anatomy.upper_teeth_width_top_cm = []
            self.anatomy.upper_teeth_width_bottom_cm = []
            
            for i in range(params.NUM_PALATE_RIBS):
                node = palate_node.find(f"p{i}")
                if node is not None:
                    x = float(node.get("x", 0))
                    z = float(node.get("z", 0))
                    self.anatomy.palate_points.append(Point3D(x, 0, z))
                    self.anatomy.palate_height_cm.append(float(node.get("palate_height", 0)))
                    self.anatomy.palate_angle_deg.append(float(node.get("palate_angle_deg", 0)))
                    self.anatomy.upper_teeth_height_cm.append(float(node.get("teeth_height", 0)))
                    self.anatomy.upper_teeth_width_top_cm.append(float(node.get("top_teeth_width", 0)))
                    self.anatomy.upper_teeth_width_bottom_cm.append(float(node.get("bottom_teeth_width", 0)))

        # Jaw
        jaw_node = root.find("jaw")
        if jaw_node is not None:
            self.anatomy.jaw_fulcrum = Point2D(float(jaw_node.get("fulcrum_x", 0)), float(jaw_node.get("fulcrum_y", 0)))
            self.anatomy.jaw_rest_pos = Point2D(float(jaw_node.get("rest_pos_x", 0)), float(jaw_node.get("rest_pos_y", 0)))
            self.anatomy.tooth_root_length_cm = float(jaw_node.get("tooth_root_length", 0))
            
            self.anatomy.jaw_points = []
            self.anatomy.jaw_height_cm = []
            self.anatomy.jaw_angle_deg = []
            self.anatomy.lower_teeth_height_cm = []
            self.anatomy.lower_teeth_width_top_cm = []
            self.anatomy.lower_teeth_width_bottom_cm = []
            
            for i in range(params.NUM_JAW_RIBS):
                node = jaw_node.find(f"p{i}")
                if node is not None:
                    x = float(node.get("x", 0))
                    z = float(node.get("z", 0))
                    self.anatomy.jaw_points.append(Point3D(x, 0, z))
                    self.anatomy.jaw_height_cm.append(float(node.get("jaw_height", 0)))
                    self.anatomy.jaw_angle_deg.append(float(node.get("jaw_angle_deg", 0)))
                    self.anatomy.lower_teeth_height_cm.append(float(node.get("teeth_height", 0)))
                    self.anatomy.lower_teeth_width_top_cm.append(float(node.get("top_teeth_width", 0)))
                    self.anatomy.lower_teeth_width_bottom_cm.append(float(node.get("bottom_teeth_width", 0)))

        # Tongue
        tongue_node = root.find("tongue")
        if tongue_node is not None:
            tip = tongue_node.find("tip")
            if tip is not None:
                self.anatomy.tongue_tip_radius_cm = float(tip.get("radius", 0))
            body = tongue_node.find("body")
            if body is not None:
                self.anatomy.tongue_center_radius_x_cm = float(body.get("radius_x", 0))
                self.anatomy.tongue_center_radius_y_cm = float(body.get("radius_y", 0))
            root_node = tongue_node.find("root")
            if root_node is not None:
                self.anatomy.automatic_tongue_root_calc = (root_node.get("automatic_calc") == "1")
                self.anatomy.tongue_root_trx_slope = float(root_node.get("trx_slope", 0))
                self.anatomy.tongue_root_trx_intercept = float(root_node.get("trx_intercept", 0))
                self.anatomy.tongue_root_try_slope = float(root_node.get("try_slope", 0))
                self.anatomy.tongue_root_try_intercept = float(root_node.get("try_intercept", 0))
            tonsil = tongue_node.find("tonsil")
            if tonsil is not None:
                self.anatomy.tongue_tonsil_length_cm = float(tonsil.get("length", 0))
                self.anatomy.tongue_tonsil_height_cm = float(tonsil.get("height", 0))

        # Lips
        lips_node = root.find("lips")
        if lips_node is not None:
            self.anatomy.lips_width_cm = float(lips_node.get("width", 0))

        # Velum
        velum_node = root.find("velum")
        if velum_node is not None:
            self.anatomy.uvula_width_cm = float(velum_node.get("uvula_width", 0))
            self.anatomy.uvula_height_cm = float(velum_node.get("uvula_height", 0))
            self.anatomy.uvula_depth_cm = float(velum_node.get("uvula_depth", 0))
            
            def parse_points_str(s):
                vals = s.split()
                pts = []
                for i in range(0, len(vals), 2):
                    pts.append(Point2D(float(vals[i]), float(vals[i+1])))
                return pts

            low = velum_node.find("low")
            if low is not None:
                self.anatomy.velum_low_points = parse_points_str(low.get("points", ""))
            mid = velum_node.find("mid")
            if mid is not None:
                self.anatomy.velum_mid_points = parse_points_str(mid.get("points", ""))
            high = velum_node.find("high")
            if high is not None:
                self.anatomy.velum_high_points = parse_points_str(high.get("points", ""))

        # Pharynx
        pharynx_node = root.find("pharynx")
        if pharynx_node is not None:
            self.anatomy.pharynx_fulcrum = Point2D(float(pharynx_node.get("fulcrum_x", 0)), float(pharynx_node.get("fulcrum_y", 0)))
            self.anatomy.pharynx_rotation_angle_deg = float(pharynx_node.get("rotation_angle_deg", 0))
            self.anatomy.pharynx_top_rib_y_cm = float(pharynx_node.get("top_rib_y", 0))
            self.anatomy.pharynx_upper_depth_cm = float(pharynx_node.get("upper_depth", 0))
            self.anatomy.pharynx_lower_depth_cm = float(pharynx_node.get("lower_depth", 0))
            self.anatomy.pharynx_back_width_cm = float(pharynx_node.get("back_side_width", 0))

        # Larynx
        larynx_node = root.find("larynx")
        if larynx_node is not None:
            self.anatomy.larynx_upper_depth_cm = float(larynx_node.get("upper_depth", 0))
            self.anatomy.larynx_lower_depth_cm = float(larynx_node.get("lower_depth", 0))
            self.anatomy.epiglottis_width_cm = float(larynx_node.get("epiglottis_width", 0))
            self.anatomy.epiglottis_height_cm = float(larynx_node.get("epiglottis_height", 0))
            self.anatomy.epiglottis_depth_cm = float(larynx_node.get("epiglottis_depth", 0))
            self.anatomy.epiglottis_angle_deg = float(larynx_node.get("epiglottis_angle_deg", 0))
            
            def parse_points_str(s):
                vals = s.split()
                pts = []
                for i in range(0, len(vals), 2):
                    pts.append(Point2D(float(vals[i]), float(vals[i+1])))
                return pts

            narrow = larynx_node.find("narrow")
            if narrow is not None:
                self.anatomy.larynx_narrow_points = parse_points_str(narrow.get("points", ""))
            wide = larynx_node.find("wide")
            if wide is not None:
                self.anatomy.larynx_wide_points = parse_points_str(wide.get("points", ""))

        # Parameters
        for param_node in root.findall("param"):
            idx = int(param_node.get("index", -1))
            if 0 <= idx < params.ParamIndex.NUM_PARAMS:
                p = self.params[idx]
                p.name = param_node.get("name", "")
                p.min_val = float(param_node.get("min", 0))
                p.max_val = float(param_node.get("max", 0))
                p.neutral_val = float(param_node.get("neutral", 0))
                p.value = p.neutral_val # Init to neutral
                self.anatomy.positive_velocity_factor[idx] = float(param_node.get("positive_velocity_factor", 1.0))
                self.anatomy.negative_velocity_factor[idx] = float(param_node.get("negative_velocity_factor", 1.0))

        # Vocal folds and static cavities
        vocal_folds_node = root.find("vocal_folds")
        if vocal_folds_node is not None:
            self.anatomy.default_f0_hz = float(vocal_folds_node.get("default_f0", self.anatomy.default_f0_hz))

        piriform_node = root.find("piriform_fossa")
        if piriform_node is not None:
            self.anatomy.piriform_fossa_length_cm = float(piriform_node.get("length", self.anatomy.piriform_fossa_length_cm))
            self.anatomy.piriform_fossa_volume_cm3 = float(piriform_node.get("volume", self.anatomy.piriform_fossa_volume_cm3))

        subglottal_node = root.find("subglottal_cavity")
        if subglottal_node is not None:
            self.anatomy.subglottal_cavity_length_cm = float(subglottal_node.get("length", self.anatomy.subglottal_cavity_length_cm))

        nasal_node = root.find("nasal_cavity")
        if nasal_node is not None:
            self.anatomy.nasal_cavity_length_cm = float(nasal_node.get("length", self.anatomy.nasal_cavity_length_cm))

        # Init Reference Surfaces
        anatomy_geometry.init_reference_surfaces(self)

    def _init_param_names(self):
        for index, param_def in enumerate(self.params):
            if param_def.name == "":
                param_def.name = f"P{index}"
            if param_def.abbr == "":
                param_def.abbr = param_def.name

    def _init_surface_grids(self):
        # Initialize surface grids with correct dimensions
        self.surfaces_list[params.SurfaceIndex.UPPER_TEETH] = Surface(params.NUM_TEETH_RIBS, params.NUM_TEETH_POINTS)
        self.surfaces_list[params.SurfaceIndex.LOWER_TEETH] = Surface(params.NUM_TEETH_RIBS, params.NUM_TEETH_POINTS)
        self.surfaces_list[params.SurfaceIndex.UPPER_COVER] = Surface(params.NUM_UPPER_COVER_RIBS, params.NUM_UPPER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.LOWER_COVER] = Surface(params.NUM_LOWER_COVER_RIBS, params.NUM_LOWER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.UPPER_LIP] = Surface(params.NUM_LIP_RIBS, params.NUM_LIP_POINTS)
        self.surfaces_list[params.SurfaceIndex.LOWER_LIP] = Surface(params.NUM_LIP_RIBS, params.NUM_LIP_POINTS)
        self.surfaces_list[params.SurfaceIndex.PALATE] = Surface(params.NUM_PALATE_RIBS, params.NUM_UPPER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.MANDIBLE] = Surface(params.NUM_JAW_RIBS, params.NUM_LOWER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.LOWER_TEETH_ORIGINAL] = Surface(params.NUM_TEETH_RIBS, params.NUM_TEETH_POINTS)
        self.surfaces_list[params.SurfaceIndex.LOW_VELUM] = Surface(params.NUM_VELUM_RIBS, params.NUM_UPPER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.MID_VELUM] = Surface(params.NUM_VELUM_RIBS, params.NUM_UPPER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.HIGH_VELUM] = Surface(params.NUM_VELUM_RIBS, params.NUM_UPPER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.NARROW_LARYNX_FRONT] = Surface(params.NUM_LARYNX_RIBS, params.NUM_LOWER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.NARROW_LARYNX_BACK] = Surface(params.NUM_LARYNX_RIBS, params.NUM_UPPER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.WIDE_LARYNX_FRONT] = Surface(params.NUM_LARYNX_RIBS, params.NUM_LOWER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.WIDE_LARYNX_BACK] = Surface(params.NUM_LARYNX_RIBS, params.NUM_UPPER_COVER_POINTS)
        self.surfaces_list[params.SurfaceIndex.TONGUE] = Surface(params.NUM_TONGUE_RIBS, params.NUM_TONGUE_POINTS)
        
        # Two-side surfaces
        self.surfaces_list[params.SurfaceIndex.UPPER_COVER_TWOSIDE] = Surface(
            params.NUM_UPPER_COVER_RIBS, params.NUM_UPPER_COVER_POINTS * 2 - 1
        )
        self.surfaces_list[params.SurfaceIndex.LOWER_COVER_TWOSIDE] = Surface(
            params.NUM_LOWER_COVER_RIBS, params.NUM_LOWER_COVER_POINTS * 2 - 1
        )
        self.surfaces_list[params.SurfaceIndex.UPPER_TEETH_TWOSIDE] = Surface(
            params.NUM_TEETH_RIBS * 2 - 1, params.NUM_TEETH_POINTS
        )
        self.surfaces_list[params.SurfaceIndex.LOWER_TEETH_TWOSIDE] = Surface(
            params.NUM_TEETH_RIBS * 2 - 1, params.NUM_TEETH_POINTS
        )
        self.surfaces_list[params.SurfaceIndex.UPPER_LIP_TWOSIDE] = Surface(
            params.NUM_LIP_RIBS * 2 - 1, params.NUM_LIP_POINTS
        )
        self.surfaces_list[params.SurfaceIndex.LOWER_LIP_TWOSIDE] = Surface(
            params.NUM_LIP_RIBS * 2 - 1, params.NUM_LIP_POINTS
        )
        
        self.surfaces_list[params.SurfaceIndex.LEFT_COVER] = Surface(params.NUM_FILL_RIBS, params.NUM_FILL_POINTS)
        self.surfaces_list[params.SurfaceIndex.RIGHT_COVER] = Surface(params.NUM_FILL_RIBS, params.NUM_FILL_POINTS)
        
        self.surfaces_list[params.SurfaceIndex.UVULA_ORIGINAL] = Surface(params.NUM_UVULA_RIBS, params.NUM_UVULA_POINTS)
        self.surfaces_list[params.SurfaceIndex.UVULA] = Surface(params.NUM_UVULA_RIBS, params.NUM_UVULA_POINTS)
        self.surfaces_list[params.SurfaceIndex.UVULA_TWOSIDE] = Surface(
            params.NUM_UVULA_RIBS, params.NUM_UVULA_POINTS * 2 - 1
        )
        
        self.surfaces_list[params.SurfaceIndex.EPIGLOTTIS_ORIGINAL] = Surface(params.NUM_EPIGLOTTIS_RIBS, params.NUM_EPIGLOTTIS_POINTS)
        self.surfaces_list[params.SurfaceIndex.EPIGLOTTIS] = Surface(params.NUM_EPIGLOTTIS_RIBS, params.NUM_EPIGLOTTIS_POINTS)
        self.surfaces_list[params.SurfaceIndex.EPIGLOTTIS_TWOSIDE] = Surface(
            params.NUM_EPIGLOTTIS_RIBS, params.NUM_EPIGLOTTIS_POINTS * 2 - 1
        )
        
        self.surfaces_list[params.SurfaceIndex.RADIATION] = Surface(params.NUM_RADIATION_RIBS, params.NUM_RADIATION_POINTS)

        self.surfaces_list[params.SurfaceIndex.UPPER_COVER].crease_angle_deg = 170.0
        self.surfaces_list[params.SurfaceIndex.LOWER_COVER].crease_angle_deg = 80.0
        self.surfaces_list[params.SurfaceIndex.UPPER_TEETH].crease_angle_deg = 40.0
        self.surfaces_list[params.SurfaceIndex.LOWER_TEETH].crease_angle_deg = 40.0
        self.surfaces_list[params.SurfaceIndex.UPPER_LIP].crease_angle_deg = 90.0
        self.surfaces_list[params.SurfaceIndex.LOWER_LIP].crease_angle_deg = 90.0
        self.surfaces_list[params.SurfaceIndex.TONGUE].crease_angle_deg = 90.0
        self.surfaces_list[params.SurfaceIndex.LEFT_COVER].crease_angle_deg = 170.0
        self.surfaces_list[params.SurfaceIndex.RIGHT_COVER].crease_angle_deg = 170.0
        self.surfaces_list[params.SurfaceIndex.EPIGLOTTIS].crease_angle_deg = 170.0
        self.surfaces_list[params.SurfaceIndex.UVULA].crease_angle_deg = 170.0

        self.surfaces_list[params.SurfaceIndex.UPPER_COVER_TWOSIDE].crease_angle_deg = 170.0
        self.surfaces_list[params.SurfaceIndex.LOWER_COVER_TWOSIDE].crease_angle_deg = 80.0
        self.surfaces_list[params.SurfaceIndex.UPPER_TEETH_TWOSIDE].crease_angle_deg = 40.0
        self.surfaces_list[params.SurfaceIndex.LOWER_TEETH_TWOSIDE].crease_angle_deg = 40.0
        self.surfaces_list[params.SurfaceIndex.UPPER_LIP_TWOSIDE].crease_angle_deg = 90.0
        self.surfaces_list[params.SurfaceIndex.LOWER_LIP_TWOSIDE].crease_angle_deg = 90.0
        self.surfaces_list[params.SurfaceIndex.EPIGLOTTIS_TWOSIDE].crease_angle_deg = 170.0
        self.surfaces_list[params.SurfaceIndex.UVULA_TWOSIDE].crease_angle_deg = 170.0

        self.surfaces_list[params.SurfaceIndex.LOWER_COVER].swap_triangle_orientation()
        self.surfaces_list[params.SurfaceIndex.UPPER_TEETH].swap_triangle_orientation()
        self.surfaces_list[params.SurfaceIndex.UPPER_LIP].swap_triangle_orientation()
        self.surfaces_list[params.SurfaceIndex.TONGUE].swap_triangle_orientation()
        self.surfaces_list[params.SurfaceIndex.LEFT_COVER].swap_triangle_orientation()
        self.surfaces_list[params.SurfaceIndex.EPIGLOTTIS].swap_triangle_orientation()

        self.surfaces_list[params.SurfaceIndex.LOWER_COVER_TWOSIDE].swap_triangle_orientation()
        self.surfaces_list[params.SurfaceIndex.UPPER_TEETH_TWOSIDE].swap_triangle_orientation()
        self.surfaces_list[params.SurfaceIndex.UPPER_LIP_TWOSIDE].swap_triangle_orientation()
        self.surfaces_list[params.SurfaceIndex.EPIGLOTTIS_TWOSIDE].swap_triangle_orientation()

    def calculate_all(self):
        dynamic_geometry.calculate_all(self)
