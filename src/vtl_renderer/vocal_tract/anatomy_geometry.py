# -*- coding: utf-8 -*-
"""
解剖结构几何生成 (Anatomy Geometry)

保留历史几何初始化拆分结构，供当前自有运行时与兼容路径共同复用。

将 VocalTract 的几何初始化逻辑分离到此模块
"""

import math
from typing import List, TYPE_CHECKING
from ..utils.geometry import Point3D, Point2D, Vector2D
from ..utils.splines import BezierCurve3D
from . import params

if TYPE_CHECKING:
    from .geometry import VocalTract, Surface

# =============================================================================
# Helpers
# =============================================================================

def get_pharynx_back_x(vt: 'VocalTract', y: float) -> float:
    min_angle_deg = -135.0
    max_angle_deg = -45.0

    fulcrum = vt.anatomy.pharynx_fulcrum
    angle_deg = float(vt.anatomy.pharynx_rotation_angle_deg)
    if angle_deg > 0.0:
        angle_deg -= 2.0 * math.pi
    if angle_deg < min_angle_deg:
        angle_deg = min_angle_deg
    if angle_deg > max_angle_deg:
        angle_deg = max_angle_deg
    angle_rad = angle_deg * math.pi / 180.0
    sine = math.sin(angle_rad)
    if abs(sine) < 1.0e-9:
        sine = 1.0e-9 if sine >= 0.0 else -1.0e-9
    return fulcrum.x + (y - fulcrum.y) * math.cos(angle_rad) / sine


# =============================================================================
# Init Larynx
# =============================================================================

def init_larynx(vt: 'VocalTract'):
    EPSILON = 0.000001
    anatomy = vt.anatomy
    
    w = [1.0, 0.7, 1.0] # Weights
    curve = BezierCurve3D()
    
    y = [0.0] * params.NUM_LARYNX_RIBS
    xl = [0.0] * params.NUM_LARYNX_RIBS
    xr = [0.0] * params.NUM_LARYNX_RIBS
    xm = [0.0] * params.NUM_LARYNX_RIBS
    zm = [0.0] * params.NUM_LARYNX_RIBS
    
    s_nb = vt.surfaces_list[params.SurfaceIndex.NARROW_LARYNX_BACK]
    s_nf = vt.surfaces_list[params.SurfaceIndex.NARROW_LARYNX_FRONT]
    s_wb = vt.surfaces_list[params.SurfaceIndex.WIDE_LARYNX_BACK]
    s_wf = vt.surfaces_list[params.SurfaceIndex.WIDE_LARYNX_FRONT]
    
    points_lists = [anatomy.larynx_narrow_points, anatomy.larynx_wide_points]
    surfaces_back = [s_nb, s_wb]
    surfaces_front = [s_nf, s_wf]
    
    for k in range(2):
        L = points_lists[k]
        
        # Calculate ribs params
        y[0]  = L[3].y
        xl[0] = L[4].x
        xr[0] = L[3].x
        xm[0] = 0.5 * (xl[0] + xr[0])
        zm[0] = -0.5 * anatomy.larynx_lower_depth_cm

        y[1]  = L[2].y
        xl[1] = L[5].x
        xr[1] = L[2].x
        xm[1] = 0.5 * (xl[1] + xr[1])
        zm[1] = -0.5 * anatomy.larynx_upper_depth_cm

        y[2]  = L[6].y
        xl[2] = L[6].x
        d = L[1].y - L[2].y
        if d < EPSILON: d = EPSILON
        xr[2] = L[2].x + (L[1].x - L[2].x) * (L[6].y - L[2].y) / d
        xm[2] = 0.25 * xl[2] + 0.75 * xr[2]
        zm[2] = -0.5 * anatomy.pharynx_lower_depth_cm

        y[3]  = L[1].y
        d = L[7].y - L[6].y
        if d < EPSILON: d = EPSILON
        xl[3] = L[6].x + (L[7].x - L[6].x) * (L[1].y - L[6].y) / d
        xr[3] = L[1].x
        xm[3] = 0.2 * xl[3] + 0.8 * xr[3]
        zm[3] = -0.5 * anatomy.pharynx_lower_depth_cm

        y[4]  = L[0].y
        xl[4] = L[7].x
        xr[4] = L[0].x
        xm[4] = 0.5 * (xl[4] + xr[4])
        zm[4] = -0.5 * anatomy.pharynx_lower_depth_cm

        # Generate Surfaces
        for i in range(params.NUM_LARYNX_RIBS):
            # Back Side
            s = surfaces_back[k]
            if s:
                C = [Point3D(xm[i], y[i], zm[i]), 
                     Point3D(xl[i], y[i], zm[i]),
                     Point3D(xl[i], y[i], 0.0)]
                curve.set_points(C, w)
                
                for j in range(params.NUM_UPPER_COVER_POINTS):
                    u = float(j) / (params.NUM_UPPER_COVER_POINTS - 1)
                    d_param = curve.get_uniform_param(u)
                    s.set_vertex(i, j, curve.get_point(d_param))
            
            # Front Side
            s = surfaces_front[k]
            if s:
                C = [Point3D(xm[i], y[i], zm[i]),
                     Point3D(xr[i], y[i], zm[i]),
                     Point3D(xr[i], y[i], 0.0)]
                curve.set_points(C, w)
                
                for j in range(params.NUM_LOWER_COVER_POINTS):
                    u = float(j) / (params.NUM_LOWER_COVER_POINTS - 1)
                    d_param = curve.get_uniform_param(u)
                    s.set_vertex(i, j, curve.get_point(d_param))

    # 3. Epiglottis
    width = anatomy.epiglottis_width_cm
    height = anatomy.epiglottis_height_cm
    depth = anatomy.epiglottis_depth_cm
    
    s_epi = vt.surfaces_list[params.SurfaceIndex.EPIGLOTTIS_ORIGINAL]
    if s_epi:
        # Ribs 0, 1, 2
        d = 0.0
        P = [Point3D(0.0, d, 0.0),
             Point3D(-0.25*width, d, -0.75*0.5*depth),
             Point3D(-0.5*width, d, -0.5*depth),
             Point3D(-width, d, -0.75*0.5*depth),
             Point3D(-width, d, 0.0)]
             
        for i in range(params.NUM_EPIGLOTTIS_POINTS):
             P[i].y = 0.0
             s_epi.set_vertex(0, i, P[i])
             P[i].y = 0.25 * height
             s_epi.set_vertex(1, i, P[i])
             P[i].y = 0.75 * height
             s_epi.set_vertex(2, i, P[i])
             
        # Rib 3
        d = height
        P = [Point3D(-0.5*width, d, 0.0),
             Point3D(-0.5*width, d, -0.75*0.5*depth),
             Point3D(-0.5*width, d, -0.75*0.5*depth),
             Point3D(-0.5*width, d, -0.75*0.5*depth),
             Point3D(-0.5*width, d, 0.0)]
             
        for i in range(params.NUM_EPIGLOTTIS_POINTS):
            s_epi.set_vertex(3, i, P[i])

# =============================================================================
# Init Jaws
# =============================================================================

def init_jaws(vt: 'VocalTract'):
    anatomy = vt.anatomy
    
    upper_jaw = vt.surfaces_list[params.SurfaceIndex.PALATE]
    lower_jaw = vt.surfaces_list[params.SurfaceIndex.MANDIBLE]
    upper_teeth = vt.surfaces_list[params.SurfaceIndex.UPPER_TEETH]
    lower_teeth = vt.surfaces_list[params.SurfaceIndex.LOWER_TEETH_ORIGINAL]
    
    MIN_ANGLE = 0.000001
    MAX_ANGLE = 89.999999
    
    curve = BezierCurve3D()
    w = [1.0, 1.0, 1.0]
    
    # Teeth geometry storage
    th = [0.0] * params.NUM_TEETH_RIBS
    ttw = [0.0] * params.NUM_TEETH_RIBS
    tbw = [0.0] * params.NUM_TEETH_RIBS
    tp = [Point3D() for _ in range(params.NUM_TEETH_RIBS)]
    tn = [Point3D() for _ in range(params.NUM_TEETH_RIBS)]

    # Upper/Lower Jaw Ribs
    for i in range(params.NUM_JAW_RIBS):
        # Upper Jaw
        C = [Point3D(), Point3D(), Point3D()]
        C[0] = Point3D(anatomy.palate_points[i].x, anatomy.palate_points[i].y, anatomy.palate_points[i].z)
        
        height = anatomy.palate_height_cm[i]
        angle = anatomy.palate_angle_deg[i]
        if angle < MIN_ANGLE: angle = MIN_ANGLE
        if angle > MAX_ANGLE: angle = MAX_ANGLE
        
        C[1].x = C[0].x
        C[1].y = height
        C[1].z = C[0].z + height / math.tan(angle * math.pi / 180.0)
        if C[1].z > 0.0: C[1].z = 0.0
        
        C[2].x = C[0].x
        C[2].y = height
        C[2].z = 0.0
        
        curve.set_points(C, w)
        
        if upper_jaw:
            for j in range(params.NUM_UPPER_COVER_POINTS):
                u = float(j) / (params.NUM_UPPER_COVER_POINTS - 1)
                u = curve.get_uniform_param(u)
                upper_jaw.set_vertex(i, j, curve.get_point(u))
                
        # Lower Jaw
        C[0] = Point3D(anatomy.jaw_points[i].x, anatomy.jaw_points[i].y, anatomy.jaw_points[i].z)
        
        height = anatomy.jaw_height_cm[i]
        angle = anatomy.jaw_angle_deg[i]
        if angle < MIN_ANGLE: angle = MIN_ANGLE
        if angle > MAX_ANGLE: angle = MAX_ANGLE
        
        # Oblique posterior ribs logic
        delta = (anatomy.jaw_points[8].x - anatomy.jaw_points[0].x) * 1.0 / 4.5
        if C[0].x - anatomy.jaw_points[0].x < delta:
            bottom_x = anatomy.jaw_points[0].x + delta
        else:
            bottom_x = C[0].x
            
        C[1].x = bottom_x
        C[1].y = -height
        C[1].z = C[0].z + height / math.tan(angle * math.pi / 180.0)
        if C[1].z > 0.0: C[1].z = 0.0
        
        C[2].x = bottom_x
        C[2].y = -height
        C[2].z = 0.0
        
        curve.set_points(C, w)
        
        if lower_jaw:
            for j in range(params.NUM_LOWER_COVER_POINTS):
                u = float(j) / (params.NUM_LOWER_COVER_POINTS - 1)
                u = curve.get_uniform_param(u)
                lower_jaw.set_vertex(i, j, curve.get_point(u))
    
    # Upper Teeth
    GROOVE_DEPTH = 0.15
    
    # Reference points and normals for teeth
    tp[0] = Point3D(anatomy.palate_points[0].x, 0, anatomy.palate_points[0].z)
    tn[0] = Point3D(anatomy.palate_points[1].z - anatomy.palate_points[0].z, 0, -(anatomy.palate_points[1].x - anatomy.palate_points[0].x))
    tn[0].normalize()
    
    tp[params.NUM_TEETH_RIBS-1] = Point3D(anatomy.palate_points[params.NUM_JAW_RIBS-1].x, 0, anatomy.palate_points[params.NUM_JAW_RIBS-1].z)
    tn[params.NUM_TEETH_RIBS-1] = Point3D(1.0, 0, 0)
    
    for i in range(1, params.NUM_JAW_RIBS - 1):
        j = i * 3
        tp[j] = Point3D(anatomy.palate_points[i].x, 0, anatomy.palate_points[i].z)
        tn[j] = Point3D(anatomy.palate_points[i+1].z - anatomy.palate_points[i-1].z, 0, -(anatomy.palate_points[i+1].x - anatomy.palate_points[i-1].x))
        tn[j].normalize()

    # Interpolate for teeth ribs between jaw points
    for i in range(params.NUM_JAW_RIBS - 1):
        t = 0.75
        
        th[i*3+1] = anatomy.upper_teeth_height_cm[i]
        ttw[i*3+1] = t * anatomy.upper_teeth_width_top_cm[i] + (1.0-t) * anatomy.upper_teeth_width_top_cm[i+1]
        tbw[i*3+1] = t * anatomy.upper_teeth_width_bottom_cm[i] + (1.0-t) * anatomy.upper_teeth_width_bottom_cm[i+1]
        tp[i*3+1] = tp[i*3] * t + tp[i*3+3] * (1.0-t)
        tn[i*3+1] = tn[i*3] * t + tn[i*3+3] * (1.0-t)
        tn[i*3+1].normalize()
        
        th[i*3+2] = anatomy.upper_teeth_height_cm[i]
        ttw[i*3+2] = (1.0-t) * anatomy.upper_teeth_width_top_cm[i] + t * anatomy.upper_teeth_width_top_cm[i+1]
        tbw[i*3+2] = (1.0-t) * anatomy.upper_teeth_width_bottom_cm[i] + t * anatomy.upper_teeth_width_bottom_cm[i+1]
        tp[i*3+2] = tp[i*3] * (1.0-t) + tp[i*3+3] * t
        tn[i*3+2] = tn[i*3] * (1.0-t) + tn[i*3+3] * t
        tn[i*3+2].normalize()
        
        if anatomy.upper_teeth_height_cm[i] == 0.0:
            ttw[i*3+1] = 0.0
            tbw[i*3+1] = 0.0
            ttw[i*3+2] = 0.0
            tbw[i*3+2] = 0.0

    # Ribs at teeth interspaces
    t = 0.8
    th[0] = th[1] - GROOVE_DEPTH
    
    th[params.NUM_TEETH_RIBS-1] = th[params.NUM_TEETH_RIBS-2] * t
    ttw[params.NUM_TEETH_RIBS-1] = ttw[params.NUM_TEETH_RIBS-2] * t
    tbw[params.NUM_TEETH_RIBS-1] = tbw[params.NUM_TEETH_RIBS-2] * t
    
    for i in range(1, params.NUM_JAW_RIBS - 1):
        j = i * 3
        if th[j-1] < th[j+1]:
            th[j] = th[j-1] - GROOVE_DEPTH
        else:
            th[j] = th[j+1] - GROOVE_DEPTH
            
        ttw[j] = 0.5 * (ttw[j-1] + ttw[j+1]) * t
        tbw[j] = 0.5 * (tbw[j-1] + tbw[j+1]) * t
        
    # Create Upper Teeth
    if upper_teeth:
        for i in range(params.NUM_TEETH_RIBS):
            if ttw[i] < 0: ttw[i] = 0
            if tbw[i] < 0: tbw[i] = 0
            if th[i] < 0: th[i] = 0
            
            Q = tp[i]
            upper_teeth.set_vertex(i, 0, Q)
            upper_teeth.set_vertex(i, 4, Q)
            
            Q = tp[i] + tn[i] * (ttw[i] - tbw[i])
            Q.y = -th[i]
            upper_teeth.set_vertex(i, 1, Q)
            
            Q = tp[i] + tn[i] * ttw[i]
            Q.y = -th[i]
            upper_teeth.set_vertex(i, 2, Q)
            Q.y = 0.0
            upper_teeth.set_vertex(i, 3, Q)
            
    # Upper gums edge construction
    for i in range(params.NUM_JAW_RIBS):
        vt.upper_gums_inner_edge[i] = Point3D(anatomy.palate_points[i].x, 0, anatomy.palate_points[i].z)
        vt.upper_gums_outer_edge[i] = vt.upper_gums_inner_edge[i] + tn[3*i] * anatomy.upper_teeth_width_top_cm[i] # tn is calculated above

    # Lower Teeth (Skip identical logic for now or implement?)
    # The logic is very similar to upper teeth but inverted Y mostly.
    # To be perfectly 1:1, I will implement it.
    
    # Reset for lower teeth
    tp[0] = Point3D(anatomy.jaw_points[0].x, 0, anatomy.jaw_points[0].z)
    tn[0] = Point3D(anatomy.jaw_points[1].z - anatomy.jaw_points[0].z, 0, -(anatomy.jaw_points[1].x - anatomy.jaw_points[0].x))
    tn[0].normalize()
    
    tp[params.NUM_TEETH_RIBS-1] = Point3D(anatomy.jaw_points[params.NUM_JAW_RIBS-1].x, 0, anatomy.jaw_points[params.NUM_JAW_RIBS-1].z)
    tn[params.NUM_TEETH_RIBS-1] = Point3D(1.0, 0, 0)
    
    for i in range(1, params.NUM_JAW_RIBS - 1):
        j = i * 3
        tp[j] = Point3D(anatomy.jaw_points[i].x, 0, anatomy.jaw_points[i].z)
        tn[j] = Point3D(anatomy.jaw_points[i+1].z - anatomy.jaw_points[i-1].z, 0, -(anatomy.jaw_points[i+1].x - anatomy.jaw_points[i-1].x))
        tn[j].normalize()

    for i in range(params.NUM_JAW_RIBS - 1):
        t = 0.75
        
        th[i*3+1] = anatomy.lower_teeth_height_cm[i]
        ttw[i*3+1] = t * anatomy.lower_teeth_width_top_cm[i] + (1.0-t) * anatomy.lower_teeth_width_top_cm[i+1]
        tbw[i*3+1] = t * anatomy.lower_teeth_width_bottom_cm[i] + (1.0-t) * anatomy.lower_teeth_width_bottom_cm[i+1]
        tp[i*3+1] = tp[i*3] * t + tp[i*3+3] * (1.0-t)
        tn[i*3+1] = tn[i*3] * t + tn[i*3+3] * (1.0-t)
        tn[i*3+1].normalize()
        
        th[i*3+2] = anatomy.lower_teeth_height_cm[i]
        ttw[i*3+2] = (1.0-t) * anatomy.lower_teeth_width_top_cm[i] + t * anatomy.lower_teeth_width_top_cm[i+1]
        tbw[i*3+2] = (1.0-t) * anatomy.lower_teeth_width_bottom_cm[i] + t * anatomy.lower_teeth_width_bottom_cm[i+1]
        tp[i*3+2] = tp[i*3] * (1.0-t) + tp[i*3+3] * t
        tn[i*3+2] = tn[i*3] * (1.0-t) + tn[i*3+3] * t
        tn[i*3+2].normalize()
        
        if anatomy.lower_teeth_height_cm[i] == 0.0:
            ttw[i*3+1] = 0.0
            tbw[i*3+1] = 0.0
            ttw[i*3+2] = 0.0
            tbw[i*3+2] = 0.0

    t = 0.8
    th[0] = th[1] - GROOVE_DEPTH
    
    th[params.NUM_TEETH_RIBS-1] = th[params.NUM_TEETH_RIBS-2] * t
    ttw[params.NUM_TEETH_RIBS-1] = ttw[params.NUM_TEETH_RIBS-2] * t
    tbw[params.NUM_TEETH_RIBS-1] = tbw[params.NUM_TEETH_RIBS-2] * t
    
    for i in range(1, params.NUM_JAW_RIBS - 1):
        j = i * 3
        if th[j-1] < th[j+1]:
            th[j] = th[j-1] - GROOVE_DEPTH
        else:
            th[j] = th[j+1] - GROOVE_DEPTH
            
        ttw[j] = 0.5 * (ttw[j-1] + ttw[j+1]) * t
        tbw[j] = 0.5 * (tbw[j-1] + tbw[j+1]) * t
    
    if lower_teeth:
        for i in range(params.NUM_TEETH_RIBS):
            if ttw[i] < 0: ttw[i] = 0
            if tbw[i] < 0: tbw[i] = 0
            if th[i] < 0: th[i] = 0
            
            Q = tp[i]
            lower_teeth.set_vertex(i, 0, Q)
            lower_teeth.set_vertex(i, 4, Q)
            
            Q = tp[i] + tn[i] * (tbw[i] - ttw[i])
            Q.y = th[i]
            lower_teeth.set_vertex(i, 1, Q)
            
            Q = tp[i] + tn[i] * tbw[i]
            Q.y = th[i]
            lower_teeth.set_vertex(i, 2, Q)
            Q.y = 0.0
            lower_teeth.set_vertex(i, 3, Q)
            
    # Lower gums edge construction
    for i in range(params.NUM_JAW_RIBS):
        vt.lower_gums_inner_edge_orig[i] = Point3D(anatomy.jaw_points[i].x, 0, anatomy.jaw_points[i].z)
        vt.lower_gums_outer_edge_orig[i] = vt.lower_gums_inner_edge_orig[i] + tn[3*i] * anatomy.lower_teeth_width_bottom_cm[i]
    
    # Lip Corner Path
    vt.wide_lip_corner_path.reset(0)
    vt.narrow_lip_corner_path.reset(0)
    
    for i in range(4, 7): # i=4,5,6
        Q = vt.upper_gums_outer_edge[i]
        Q.y = 0.0
        vt.wide_lip_corner_path.add_point(Q)
        vt.narrow_lip_corner_path.add_point(Q)
        
    palate_length = anatomy.palate_points[8].x - anatomy.palate_points[0].x
    palate_depth = -2.0 * anatomy.palate_points[0].z
    
    lip_protrusion = palate_length / 4.5
    lip_front_min_z = -1.0 * palate_depth / 4.6
    lip_front_max_z = -0.4 * palate_depth / 4.6
    
    max_lip_corner_x = vt.upper_gums_outer_edge[8].x + lip_protrusion
    
    vt.wide_lip_corner_path.add_point(Point3D(max_lip_corner_x, 0.0, lip_front_min_z))
    vt.narrow_lip_corner_path.add_point(Point3D(max_lip_corner_x, 0.0, lip_front_max_z))

# =============================================================================
# Init Velum
# =============================================================================

def init_velum(vt: 'VocalTract'):
    anatomy = vt.anatomy
    
    s_low = vt.surfaces_list[params.SurfaceIndex.LOW_VELUM]
    s_mid = vt.surfaces_list[params.SurfaceIndex.MID_VELUM]
    s_high = vt.surfaces_list[params.SurfaceIndex.HIGH_VELUM]
    
    w = [1.0, 1.0, 1.0]
    curve = BezierCurve3D()
    
    # Target Z calculation
    source_z = -0.5 * anatomy.pharynx_upper_depth_cm
    t = anatomy.palate_angle_deg[0]
    if t <= 0.00001: t = 0.00001
    if t > 89.99999: t = 89.99999
    target_z = anatomy.palate_points[0].z + anatomy.palate_height_cm[1] / math.tan(t * math.pi / 180.0)
    
    posterior_end_point = Point3D()
    posterior_end_point.x = get_pharynx_back_x(vt, anatomy.pharynx_top_rib_y_cm) + anatomy.pharynx_back_width_cm
    posterior_end_point.y = anatomy.pharynx_top_rib_y_cm
    posterior_end_point.z = -0.5 * anatomy.pharynx_upper_depth_cm
    
    anterior_end_point = Point3D(0.0, 0.0, anatomy.palate_points[0].z)
    
    C = [Point3D(), Point3D(), Point3D()]
    
    for i in range(params.NUM_VELUM_RIBS):
        if i == 0:
            C[0] = posterior_end_point
        else:
            t_interp = (i - 1) / (params.NUM_VELUM_RIBS - 2)
            C[0] = posterior_end_point * (1.0 - t_interp) + anterior_end_point * t_interp
            
        t_rib = i / (params.NUM_VELUM_RIBS - 1)
        
        # LOW
        if s_low:
            if i == 0:
                C[2].x = get_pharynx_back_x(vt, anatomy.velum_low_points[0].y)
                C[2].y = anatomy.velum_low_points[0].y
            else:
                pt = anatomy.velum_low_points[i-1]
                C[2] = Point3D(pt.x, pt.y, 0.0)
            C[2].z = 0.0
            
            C[1].x = C[2].x
            C[1].y = C[2].y
            C[1].z = source_z * (1.0 - t_rib) + target_z * t_rib
            
            curve.set_points(C, w)
            for j in range(params.NUM_UPPER_COVER_POINTS):
                u = j / (params.NUM_UPPER_COVER_POINTS - 1)
                u = curve.get_uniform_param(u)
                s_low.set_vertex(i, j, curve.get_point(u))

        # MID
        if s_mid:
            if i == 0:
                C[2].x = get_pharynx_back_x(vt, anatomy.velum_mid_points[0].y)
                C[2].y = anatomy.velum_mid_points[0].y
            else:
                pt = anatomy.velum_mid_points[i-1]
                C[2] = Point3D(pt.x, pt.y, 0.0)
            C[2].z = 0.0
            
            C[1].x = C[2].x
            C[1].y = C[2].y
            C[1].z = source_z * (1.0 - t_rib) + target_z * t_rib
            
            curve.set_points(C, w)
            for j in range(params.NUM_UPPER_COVER_POINTS):
                u = j / (params.NUM_UPPER_COVER_POINTS - 1)
                u = curve.get_uniform_param(u)
                s_mid.set_vertex(i, j, curve.get_point(u))

        # HIGH
        if s_high:
            if i == 0:
                C[2].x = get_pharynx_back_x(vt, anatomy.velum_high_points[0].y)
                C[2].y = anatomy.velum_high_points[0].y
            else:
                pt = anatomy.velum_high_points[i-1]
                C[2] = Point3D(pt.x, pt.y, 0.0)
            C[2].z = 0.0
            
            C[1].x = C[2].x
            C[1].y = C[2].y
            C[1].z = source_z * (1.0 - t_rib) + target_z * t_rib
            
            curve.set_points(C, w)
            for j in range(params.NUM_UPPER_COVER_POINTS):
                u = j / (params.NUM_UPPER_COVER_POINTS - 1)
                u = curve.get_uniform_param(u)
                s_high.set_vertex(i, j, curve.get_point(u))
                
    # Uvula Init (uvula_original)
    width = anatomy.uvula_width_cm
    height = anatomy.uvula_height_cm
    depth = anatomy.uvula_depth_cm
    
    s_uvula = vt.surfaces_list[params.SurfaceIndex.UVULA_ORIGINAL]
    if s_uvula:
        for i in range(params.NUM_UVULA_POINTS):
            t = math.pi * i / (params.NUM_UVULA_POINTS - 1)
            # Rib 0
            P = Point3D(0.5*width*math.cos(t) - 0.5*width, 0.0, -0.5*depth*math.sin(t))
            s_uvula.set_vertex(0, i, P)
            # Rib 1
            P.y = -0.5*height
            s_uvula.set_vertex(1, i, P)
            # Rib 2
            P = Point3D(0.75*0.5*width*math.cos(t) - 0.5*width, -0.75*height, -0.75*0.5*depth*math.sin(t))
            s_uvula.set_vertex(2, i, P)
            # Rib 3
            P = Point3D(-0.5*width, -height, 0.0)
            s_uvula.set_vertex(3, i, P)

    
def init_reference_surfaces(vt: 'VocalTract'):
    init_larynx(vt)
    init_jaws(vt)
    init_velum(vt)
