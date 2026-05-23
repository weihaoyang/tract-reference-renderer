# -*- coding: utf-8 -*-
"""
样条曲线工具库 (Splines Utils)

保留历史几何工具接口，供当前自有运行时复用。

包含:
- Spline3D
- LineStrip3D
- LineStrip2D (Core for Vocal Tract outlines)
- BezierCurve3D (Rational Bezier Curve)
"""

import math
from typing import List, Tuple, Optional
from .geometry import Point2D, Point3D, Vector2D, scalar_product_3d

MAX_SPLINE_POINTS = 256

# =============================================================================
# Spline3D (Base)
# =============================================================================

class Spline3D:
    def __init__(self, points: List[Point3D] = None, weights: List[float] = None):
        self.p: List[Point3D] = []
        self.w: List[float] = []
        self.num_points = 0
        self.points_changed = False
        
        if points:
            self.set_points(points, weights)

    def reset(self, new_num_points: int):
        self.num_points = new_num_points
        if self.num_points < 0: self.num_points = 0
        if self.num_points > MAX_SPLINE_POINTS: self.num_points = MAX_SPLINE_POINTS
        
        self.p = [Point3D() for _ in range(self.num_points)]
        self.w = [1.0] * self.num_points
        self.points_changed = True

    def add_point(self, point: Point3D, weight: float = 1.0):
        if self.num_points >= MAX_SPLINE_POINTS: return
        self.p.append(Point3D(point.x, point.y, point.z))
        self.w.append(float(weight))
        self.num_points += 1
        self.points_changed = True

    def set_points(self, points: List[Point3D], weights: Optional[List[float]] = None):
        self.p = [Point3D(pt.x, pt.y, pt.z) for pt in points[:MAX_SPLINE_POINTS]] # Copy
        self.num_points = len(self.p)
        
        if weights and len(weights) == self.num_points:
            self.w = weights[:MAX_SPLINE_POINTS] 
        else:
            self.w = [1.0] * self.num_points
        self.points_changed = True

    def get_control_point(self, index: int) -> Point3D:
        if 0 <= index < self.num_points:
            point = self.p[index]
            return Point3D(point.x, point.y, point.z)
        return Point3D()

    def get_point(self, t: float) -> Point3D:
        if self.num_points < 2:
            return Point3D(0, 0, 0)
            
        num_segments = self.num_points - 1
        segment = int(t * num_segments)
        
        if segment < 0: segment = 0
        if segment >= num_segments: segment = num_segments - 1
        
        delta_t = 1.0 / num_segments
        u = (t - segment * delta_t) / delta_t
        
        if u < 0: u = 0.0
        if u > 1: u = 1.0
        
        # Linear Interp Base class default?
        # VTL: return (P[segment]*u1 + P[segment+1]*u);
        u1 = 1.0 - u
        return self.p[segment] * u1 + self.p[segment+1] * u

    def get_uniform_param(self, t: float) -> float:
        EPSILON = 0.000001
        N = 100
        pts = []
        s = []
        pos = []
        
        for i in range(N):
            si = i / (N - 1)
            s.append(si)
            pts.append(self.get_point(si))
            if i == 0:
                pos.append(0.0)
            else:
                pos.append(pos[i-1] + (pts[i] - pts[i-1]).magnitude())
                
        l = pos[-1]
        if l < EPSILON: l = EPSILON
        
        pos = [p / l for p in pos]
        
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        
        k = -1
        ratio = 0.0
        
        for i in range(N - 1):
            if pos[i] <= t <= pos[i+1]:
                k = i
                seg_len = pos[i+1] - pos[i]
                if seg_len < EPSILON: seg_len = EPSILON
                ratio = (t - pos[i]) / seg_len
                
        if k == -1:
            if t < 0.5:
                k = 0
                ratio = 0.0
            else:
                k = N - 2
                ratio = 1.0
                
        u = s[k] + ratio * (s[k+1] - s[k])
        return u

# =============================================================================
# LineStrip3D
# =============================================================================

class LineStrip3D(Spline3D):
    def __init__(self, points: List[Point3D] = None):
        super().__init__(points)
        self.pos = [] # Parametrization by length
        if points:
            self._calculate_params()

    def _calculate_params(self):
        # Calculate pos[i] based on cumulative length
        self.pos = [0.0] * self.num_points
        total_len = 0.0
        
        for i in range(1, self.num_points):
            dist = (self.p[i] - self.p[i-1]).magnitude()
            total_len += dist
            self.pos[i] = total_len
            
        if total_len > 0:
            for i in range(self.num_points):
                self.pos[i] /= total_len # Normalize to 0..1
        self.points_changed = False

    def get_point(self, t: float) -> Point3D:
        if self.points_changed: self._calculate_params()
        if self.num_points < 1:
            return Point3D(0.0, 0.0, 0.0)
        if self.num_points == 1:
            point = self.p[0]
            return Point3D(point.x, point.y, point.z)

        if t < 0.0:
            t = 0.0
        if t > 1.0:
            t = 1.0

        epsilon = 1.0e-6
        ratio = 0.0
        k = -1

        for i in range(self.num_points - 1):
            if (t >= self.pos[i] - epsilon) and (t <= self.pos[i + 1] + epsilon):
                k = i
                segment_length = self.pos[i + 1] - self.pos[i]
                if segment_length < epsilon:
                    segment_length = epsilon
                ratio = (t - self.pos[i]) / segment_length

        if k == -1:
            return Point3D(0.0, 0.0, 0.0)
        return self.p[k] + (self.p[k + 1] - self.p[k]) * ratio

    def get_intersection(self, q: Point3D, v: Point3D) -> float:
        if self.points_changed:
            self._calculate_params()
        if self.num_points < 2:
            return 0.0

        epsilon = 1.0e-6
        result = 0.0
        min_dist = 1_000_000.0

        for i in range(self.num_points - 1):
            p0 = Point3D(self.p[i].x, self.p[i].y, self.p[i].z)
            p1 = Point3D(self.p[i + 1].x, self.p[i + 1].y, self.p[i + 1].z)

            if i == 0:
                segment = self.p[i + 1] - self.p[i]
                length = segment.magnitude()
                if length > 0.0:
                    segment /= length
                    p0 -= segment * epsilon

            if i + 1 == self.num_points - 1:
                segment = self.p[i + 1] - self.p[i]
                length = segment.magnitude()
                if length > 0.0:
                    segment /= length
                    p1 += segment * epsilon

            d0 = scalar_product_3d(p0 - q, v)
            d1 = scalar_product_3d(p1 - q, v)

            if ((d0 <= 0.0) and (d1 >= 0.0)) or ((d1 <= 0.0) and (d0 >= 0.0)):
                segment = p1 - p0
                denominator = scalar_product_3d(segment, v)
                if denominator != 0.0:
                    t = -scalar_product_3d(p0 - q, v) / denominator
                    intersection = p0 + t * segment
                    dist = (intersection - q).magnitude()
                    if dist < min_dist:
                        result = self.pos[i] + t * (self.pos[i + 1] - self.pos[i])
                        min_dist = dist

        return result

# =============================================================================
# BezierCurve3D
# =============================================================================

class BezierCurve3D(Spline3D):
    def __init__(self, points: List[Point3D] = None, weights: List[float] = None):
        super().__init__(points, weights)
        self.A: List[Point3D] = [Point3D() for _ in range(MAX_SPLINE_POINTS)]
        self.B: List[float] = [0.0] * MAX_SPLINE_POINTS
        
    def get_point(self, t: float) -> Point3D:
        if self.points_changed: self._calculate_coeff()
        
        if self.num_points < 2: return Point3D()
        
        n = self.num_points - 1
        return self._evaluate_polynomial(t, n)

    def _evaluate_polynomial(self, t: float, n: int) -> Point3D:
        Q = Point3D()
        f = 1.0 # t^i
        denominator = 0.0
        
        for i in range(n + 1):
            Q += self.A[i] * f
            denominator += self.B[i] * f
            f *= t
            
        if abs(denominator) > 1e-9:
            Q /= denominator
        return Q

    def _calculate_coeff(self):
        n = self.num_points - 1
        
        # Init coeffs
        for j in range(n + 1):
            self.A[j] = Point3D()
            self.B[j] = 0.0
            
        coeff = [0.0] * MAX_SPLINE_POINTS
        
        for i in range(n + 1):
            self._get_bernstein_coeff(i, n, coeff)
            
            for j in range(n + 1):
                # A[j] += w[i] * coeff[j] * P[i]
                term = self.p[i] * (self.w[i] * coeff[j])
                self.A[j] += term
                self.B[j] += self.w[i] * coeff[j]
                
        self.points_changed = False

    def _get_bernstein_coeff(self, k: int, n: int, coeff: List[float]):
        # Implementation from Splines.cpp
        # Returns power basis coefficients for Bernstein polynomial B_{k,n}(t)
        
        # Reset
        for i in range(n + 1): coeff[i] = 0.0
        
        # t^k makes coeff[k] = 1
        coeff[k] = 1.0
        
        # Multiply by (1-t)^(n-k)
        new_coeff = [0.0] * MAX_SPLINE_POINTS
        
        for i in range(1, n - k + 1):
            for j in range(n + 1): new_coeff[j] = 0.0
            
            # (1-t) * P(t) = P(t) - t*P(t)
            # coeff[j] contributes to new_coeff[j] (times 1) and new_coeff[j+1] (times -1)
            for j in range(n + 1):
                new_coeff[j] += coeff[j]
                if j + 1 < MAX_SPLINE_POINTS:
                    new_coeff[j+1] -= coeff[j]
                    
            for j in range(n + 1): coeff[j] = new_coeff[j]
            
        # Multiply by binomial coefficient (n choose k)
        numerator = 1.0
        denominator = 1.0
        
        for j in range(2, n + 1): numerator *= j
        for j in range(2, k + 1): denominator *= j
        for j in range(2, n - k + 1): denominator *= j
        
        factor = numerator / denominator
        for j in range(n + 1): coeff[j] *= factor


# =============================================================================
# LineStrip2D
# =============================================================================

class LineStrip2D:
    def __init__(self, points: List[Point2D] = None):
        self.p: List[Point2D] = []
        self.pos: List[float] = []
        self.num_points = 0
        
        if points:
            self.set_points(points)

    def set_points(self, points: List[Point2D]):
        self.p = [Point2D(pt.x, pt.y) for pt in points[:MAX_SPLINE_POINTS]]
        self.num_points = len(self.p)
        self._calculate_params()
    
    def add_point(self, point: Point2D):
        if self.num_points >= MAX_SPLINE_POINTS: return
        self.p.append(Point2D(point.x, point.y))
        self.num_points += 1
        self._calculate_params() # Re-calc immediately for simplicity

    def reset(self, num_points: int):
        self.p = [Point2D(0,0) for _ in range(num_points)]
        self.num_points = num_points
        self._calculate_params()

    def set_point(self, index: int, point: Point2D):
        if 0 <= index < self.num_points:
            self.p[index] = Point2D(point.x, point.y)
            self._calculate_params()

    def del_point(self):
        if self.num_points <= 0:
            return
        self.p.pop()
        self.num_points -= 1
        self._calculate_params()

    def get_control_point(self, index: int) -> Point2D:
        if 0 <= index < self.num_points:
            point = self.p[index]
            return Point2D(point.x, point.y)
        return Point2D()

    def get_num_points(self) -> int:
        return self.num_points

    def get_curve_param(self, index: int) -> float:
        if self.num_points <= 0:
            return 0.0
        safe_index = max(0, min(self.num_points - 1, index))
        if not self.pos:
            self._calculate_params()
        return self.pos[safe_index]

    def get_function_value(self, x: float) -> float:
        if self.num_points < 1:
            return 0.0
        if self.num_points == 1:
            return self.p[0].y

        epsilon = 1.0e-6
        result = 0.0

        for idx in range(self.num_points - 1):
            p0 = self.p[idx]
            p1 = self.p[idx + 1]
            if (x >= p0.x - epsilon) and (x <= p1.x + epsilon):
                length = p1.x - p0.x
                if length < epsilon:
                    length = epsilon
                result = p0.y + (p1.y - p0.y) * (x - p0.x) / length

        return result

    def _calculate_params(self):
        if self.num_points == 0:
            self.pos = []
            return
            
        self.pos = [0.0] * self.num_points
        total_len = 0.0
        
        for i in range(1, self.num_points):
            dist = (self.p[i] - self.p[i-1]).magnitude()
            total_len += dist
            self.pos[i] = total_len
            
        if total_len > 1e-9:
            for i in range(self.num_points):
                self.pos[i] /= total_len
        else:
            self.pos = [0.0] * self.num_points

    def get_point(self, t: float) -> Point2D:
        if self.num_points < 1:
            return Point2D(0.0, 0.0)
        if self.num_points == 1:
            point = self.p[0]
            return Point2D(point.x, point.y)

        if t < 0.0:
            t = 0.0
        if t > 1.0:
            t = 1.0

        epsilon = 1.0e-6
        ratio = 0.0
        k = -1

        for i in range(self.num_points - 1):
            if (t >= self.pos[i] - epsilon) and (t <= self.pos[i + 1] + epsilon):
                k = i
                segment_length = self.pos[i + 1] - self.pos[i]
                if segment_length < epsilon:
                    segment_length = epsilon
                ratio = (t - self.pos[i]) / segment_length

        if k == -1:
            return Point2D(0.0, 0.0)
        return self.p[k] + (self.p[k + 1] - self.p[k]) * ratio

    def get_closest_intersection(self, q: Point2D, v: Point2D) -> Tuple[bool, float, Point2D]:
        epsilon = 1.0e-6
        if self.num_points < 2:
            return False, 1_000_000.0, Point2D(0.0, 0.0)

        n = Point2D(-v.y, v.x).normalize()
        p_left = q + epsilon * n
        p_right = q - epsilon * n

        found = False
        best_t = 1_000_000.0
        intersection = Point2D(0.0, 0.0)

        old_section = 0
        for idx in range(self.num_points):
            section = 0

            w = self.p[idx] - p_left
            d = w.x * v.y - w.y * v.x
            if d < 0.0:
                section = -1

            w = self.p[idx] - p_right
            d = w.x * v.y - w.y * v.x
            if d > 0.0:
                section = 1

            if idx > 0 and (((section >= 0) and (old_section <= 0)) or ((section <= 0) and (old_section >= 0))):
                r = q - self.p[idx - 1]
                w = self.p[idx] - self.p[idx - 1]

                denominator = v.x * w.y - v.y * w.x
                if denominator != 0.0:
                    s = (v.x * r.y - r.x * v.y) / denominator
                    if -epsilon <= s <= 1.0 + epsilon:
                        t = (w.x * r.y - r.x * w.y) / denominator
                        if abs(t) < abs(best_t):
                            best_t = t
                            intersection = q + best_t * v
                            found = True

            old_section = section

        return found, best_t, intersection
                    
    def get_intersection_with_greatest_t(self, q: Point2D, v: Point2D) -> Tuple[bool, float, Point2D]:
        epsilon = 1.0e-6
        if self.num_points < 2:
            return False, -1_000_000.0, Point2D(0.0, 0.0)

        n = Point2D(-v.y, v.x).normalize()
        p_left = q + epsilon * n
        p_right = q - epsilon * n

        found = False
        best_t = -1_000_000.0
        intersection = Point2D(0.0, 0.0)

        old_section = 0
        for idx in range(self.num_points):
            section = 0

            w = self.p[idx] - p_left
            d = w.x * v.y - w.y * v.x
            if d < 0.0:
                section = -1

            w = self.p[idx] - p_right
            d = w.x * v.y - w.y * v.x
            if d > 0.0:
                section = 1

            if idx > 0 and (((section >= 0) and (old_section <= 0)) or ((section <= 0) and (old_section >= 0))):
                r = q - self.p[idx - 1]
                w = self.p[idx] - self.p[idx - 1]

                denominator = v.x * w.y - v.y * w.x
                if denominator != 0.0:
                    s = (v.x * r.y - r.x * v.y) / denominator
                    if (-epsilon <= s <= 1.0 + epsilon):
                        t = (w.x * r.y - r.x * w.y) / denominator
                        if t > best_t:
                            best_t = t
                            intersection = q + best_t * v
                            found = True

            old_section = section

        return found, best_t, intersection

    def get_special_intersection(self, q: Point2D, v: Point2D) -> Tuple[bool, float, Point2D]:
        epsilon = 1.0e-6
        if self.num_points < 2:
            return False, 1_000_000.0, Point2D(0.0, 0.0)

        n = Point2D(-v.y, v.x).normalize()
        p_left = q + epsilon * n
        p_right = q - epsilon * n

        found = False
        best_t = 1_000_000.0
        intersection = Point2D(0.0, 0.0)

        old_section = 0
        for idx in range(self.num_points):
            section = 0

            w = self.p[idx] - p_left
            d = w.x * v.y - w.y * v.x
            if d < 0.0:
                section = -1

            w = self.p[idx] - p_right
            d = w.x * v.y - w.y * v.x
            if d > 0.0:
                section = 1

            if idx > 0 and (section >= 0) and (old_section <= 0):
                r = q - self.p[idx - 1]
                w = self.p[idx] - self.p[idx - 1]

                denominator = v.x * w.y - v.y * w.x
                if denominator != 0.0:
                    s = (v.x * r.y - r.x * v.y) / denominator
                    if (-epsilon <= s <= 1.0 + epsilon):
                        t = (w.x * r.y - r.x * w.y) / denominator
                        accept = False
                        if not found:
                            accept = True
                        elif t > 0.0:
                            if best_t < 0.0 or t < best_t:
                                accept = True
                        else:
                            if best_t < 0.0 and t > best_t:
                                accept = True

                        if accept:
                            best_t = t
                            intersection = q + best_t * v
                            found = True

            old_section = section

        return found, best_t, intersection
