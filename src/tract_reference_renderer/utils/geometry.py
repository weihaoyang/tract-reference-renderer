# -*- coding: utf-8 -*-
"""
几何工具库 (Geometry Utils)

保留历史几何原语接口，供当前自有运行时与兼容层共用。

包含:
- Point2D, Point3D
- Vector2D, Vector3D
- Line2D, Line3D
- Circle, Ellipse
"""

import math
from typing import Tuple, Optional

# =============================================================================
# Point2D
# =============================================================================

class Point2D:
    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"Point2D({self.x}, {self.y})"

    def set(self, x: float, y: float):
        self.x = x
        self.y = y

    def __add__(self, other):
        return Point2D(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Point2D(self.x - other.x, self.y - other.y)

    def __mul__(self, scalar: float):
        return Point2D(self.x * scalar, self.y * scalar)
        
    def __rmul__(self, scalar: float): # float * Point2D
        return Point2D(self.x * scalar, self.y * scalar)

    def __truediv__(self, scalar: float):
        return Point2D(self.x / scalar, self.y / scalar)
    
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def magnitude(self) -> float:
        return math.sqrt(self.x**2 + self.y**2)
        
    def normalize(self):
        m = self.magnitude()
        if m > 1e-9:
            self.x /= m
            self.y /= m
        return self

    def turn_right(self):
        # (x, y) -> (y, -x)
        return Point2D(self.y, -self.x)

    def turn_left(self):
        # (x, y) -> (-y, x)
        return Point2D(-self.y, self.x)

    def turn(self, angle: float):
        s = math.sin(angle)
        c = math.cos(angle)
        nx = self.x * c - self.y * s
        ny = self.y * c + self.x * s
        self.x = nx
        self.y = ny

    def is_right_from(self, vector: 'Vector2D') -> bool:
        # Cross product 2D: (v.x * (P.y - y) - v.y * (P.x - x)) >= 0
        dx = vector.p.x - self.x
        dy = vector.p.y - self.y
        return (vector.v.x * dy - vector.v.y * dx) >= 0

    def is_left_from(self, vector: 'Vector2D') -> bool:
        dx = vector.p.x - self.x
        dy = vector.p.y - self.y
        return (vector.v.x * dy - vector.v.y * dx) <= 0

    def distance_from(self, vector: 'Vector2D') -> float:
        w = Point2D(vector.p.x - self.x, vector.p.y - self.y)
        denom = vector.v.x**2 + vector.v.y**2
        if denom < 1e-9:
            denom = 1e-4 # Avoid div 0
        return (vector.v.x * w.y - vector.v.y * w.x) / denom


# =============================================================================
# Point3D
# =============================================================================

class Point3D:
    def __init__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"Point3D({self.x}, {self.y}, {self.z})"

    def set(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        return Point3D(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Point3D(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar: float):
        return Point3D(self.x * scalar, self.y * scalar, self.z * scalar)
    
    def __rmul__(self, scalar: float):
        return Point3D(self.x * scalar, self.y * scalar, self.z * scalar)

    def __truediv__(self, scalar: float):
        return Point3D(self.x / scalar, self.y / scalar, self.z / scalar)

    def magnitude(self) -> float:
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def normalize(self):
        m = self.magnitude()
        if m > 1e-9:
            self.x /= m
            self.y /= m
            self.z /= m
        return self
    
    def to_point2d(self) -> Point2D:
        return Point2D(self.x, self.y)


# =============================================================================
# Helper Functions
# =============================================================================

def scalar_product_2d(p: Point2D, q: Point2D) -> float:
    return p.x * q.x + p.y * q.y

def scalar_product_3d(p: Point3D, q: Point3D) -> float:
    return p.x * q.x + p.y * q.y + p.z * q.z

def cross_product(p: Point3D, q: Point3D) -> Point3D:
    # (P.y*Q.z - P.z*Q.y, P.z*Q.x - P.x*Q.z, P.x*Q.y - P.y*Q.x)
    return Point3D(
        p.y * q.z - p.z * q.y,
        p.z * q.x - p.x * q.z,
        p.x * q.y - p.y * q.x
    )


# =============================================================================
# Vector2D (Line defined by Point + Direction)
# =============================================================================

class Vector2D:
    def __init__(self, p: Point2D = None, v: Point2D = None):
        self.p = p if p else Point2D()
        self.v = v if v else Point2D()

    def set(self, p: Point2D, v: Point2D):
        self.p = p
        self.v = v

    def normalize(self):
        self.v.normalize()

    def get_point(self, t: float) -> Point2D:
        return self.p + self.v * t


# =============================================================================
# Line2D (Line defined by 2 Points)
# =============================================================================

class Line2D:
    def __init__(self, p0: Point2D = None, p1: Point2D = None):
        self.p = [p0 if p0 else Point2D(), p1 if p1 else Point2D()]

    def set(self, p0: Point2D, p1: Point2D):
        self.p[0] = p0
        self.p[1] = p1

    def get_point(self, t: float) -> Point2D:
        # P0 + t * (P1 - P0)
        return self.p[0] + (self.p[1] - self.p[0]) * t

    def get_length(self) -> float:
        return (self.p[1] - self.p[0]).magnitude()

    def get_intersection(self, vector: Vector2D) -> Tuple[Point2D, float, bool]:
        # Returns (Point, t, ok)
        # Note: arg is Vector2D V
        # Intersection of this line (segment?) with infinite line V(P, v)?
        # Or V intersects this line?
        # 历史兼容接口: Point2D Line2D::getIntersection(Vector2D V, double& t, bool& ok)
        
        a = self.p[0]
        v = self.p[1] - self.p[0]
        b = vector.p
        w = vector.v
        
        orig_denom = v.x * w.y - v.y * w.x
        denom = orig_denom
        if abs(denom) < 1e-9:
            denom = 1e-4
            
        # t = (w.x*(A.y - B.y) - w.y*(A.x - B.x)) / denominator;
        t = (w.x * (a.y - b.y) - w.y * (a.x - b.x)) / denom
        
        ok = (-0.01 < t < 1.01) and (orig_denom != 0.0)
        
        return self.get_point(t), t, ok


# =============================================================================
# Circle
# =============================================================================

class Circle:
    def __init__(self, center: Point2D = None, radius: float = 1.0):
        self.m = center if center else Point2D()
        self.r = radius
        self.arc_angle = [0.0, 0.0]

    def set_valid_arc(self, angle0: float, angle1: float):
        while angle0 < 0: angle0 += 2 * math.pi
        while angle0 > 2 * math.pi: angle0 -= 2 * math.pi
        self.arc_angle[0] = angle0
        
        while angle1 < 0: angle1 += 2 * math.pi
        while angle1 > 2 * math.pi: angle1 -= 2 * math.pi
        self.arc_angle[1] = angle1

    def get_point(self, angle: float) -> Point2D:
        return Point2D(
            self.m.x + self.r * math.cos(angle),
            self.m.y + self.r * math.sin(angle)
        )
