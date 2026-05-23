# -*- coding: utf-8 -*-
"""
声道参数定义 (Vocal Tract Parameters)

保留历史参数布局与索引约定，供当前自有运行时与兼容路径共用。

包含:
- 物理常量 (Ribs count, Points count)
- ParamIndex 枚举
- Anatomy 数据结构
- TongueRib, CenterLinePoint 等辅助结构
"""

from enum import IntEnum
from dataclasses import dataclass, field
from typing import List, Tuple
from ..utils.geometry import Point2D, Point3D
from ..types import Articulator

# =============================================================================
# Constants
# =============================================================================

NUM_EPIGLOTTIS_RIBS = 4
NUM_EPIGLOTTIS_POINTS = 5
NUM_UVULA_RIBS = 4
NUM_UVULA_POINTS = 5

NUM_LARYNX_RIBS  = 5
NUM_PHARYNX_RIBS = 3
NUM_VELUM_RIBS   = 6
NUM_PALATE_RIBS  = 9    # for both upper and lower jaw
NUM_JAW_RIBS     = NUM_PALATE_RIBS
NUM_THROAT_RIBS  = NUM_PHARYNX_RIBS

NUM_UPPER_COVER_RIBS = NUM_LARYNX_RIBS + NUM_PHARYNX_RIBS + NUM_VELUM_RIBS + NUM_JAW_RIBS
NUM_UPPER_COVER_POINTS = 6

NUM_LOWER_COVER_RIBS = NUM_LARYNX_RIBS + NUM_THROAT_RIBS + NUM_JAW_RIBS
NUM_LOWER_COVER_POINTS = 5

NUM_DYNAMIC_TONGUE_RIBS = 33
NUM_STATIC_TONGUE_RIBS  = 4
NUM_TONGUE_RIBS = NUM_DYNAMIC_TONGUE_RIBS + NUM_STATIC_TONGUE_RIBS
NUM_TONGUE_POINTS = 11
MAX_TONGUE_RIBS_GLOBAL = 128

NUM_TEETH_RIBS = 3 * (NUM_JAW_RIBS - 1) + 1
NUM_TEETH_POINTS = 5
    
NUM_LIP_RIBS = NUM_JAW_RIBS
NUM_INNER_LIP_POINTS = 5
NUM_OUTER_LIP_POINTS = 5
NUM_LIP_POINTS = NUM_INNER_LIP_POINTS + NUM_OUTER_LIP_POINTS

NUM_FILL_RIBS = NUM_JAW_RIBS + 3
NUM_FILL_POINTS = 4

NUM_RADIATION_RIBS = NUM_LIP_RIBS + 4
NUM_RADIATION_POINTS = 6

# =============================================================================
# SurfaceIndex
# =============================================================================

class SurfaceIndex(IntEnum):
    UPPER_TEETH = 0
    LOWER_TEETH = 1
    UPPER_COVER = 2
    LOWER_COVER = 3
    UPPER_LIP = 4
    LOWER_LIP = 5
    PALATE = 6
    MANDIBLE = 7
    LOWER_TEETH_ORIGINAL = 8
    LOW_VELUM = 9
    MID_VELUM = 10
    HIGH_VELUM = 11
    NARROW_LARYNX_FRONT = 12
    NARROW_LARYNX_BACK = 13
    WIDE_LARYNX_FRONT = 14
    WIDE_LARYNX_BACK = 15
    TONGUE = 16
    UPPER_COVER_TWOSIDE = 17
    LOWER_COVER_TWOSIDE = 18
    UPPER_TEETH_TWOSIDE = 19
    LOWER_TEETH_TWOSIDE = 20
    UPPER_LIP_TWOSIDE = 21
    LOWER_LIP_TWOSIDE = 22
    LEFT_COVER = 23
    RIGHT_COVER = 24
    UVULA_ORIGINAL = 25
    UVULA = 26
    UVULA_TWOSIDE = 27
    EPIGLOTTIS_ORIGINAL = 28
    EPIGLOTTIS = 29
    EPIGLOTTIS_TWOSIDE = 30
    RADIATION = 31
    NUM_SURFACES = 32


NUM_SURFACES = int(SurfaceIndex.NUM_SURFACES)

# =============================================================================
# Constants for cross-sections
# =============================================================================
NUM_PROFILE_SAMPLES = 96
PROFILE_LENGTH = 7.0
PROFILE_SAMPLE_LENGTH = PROFILE_LENGTH / NUM_PROFILE_SAMPLES
INVALID_PROFILE_SAMPLE = 1_000_000.0
EXTREME_PROFILE_VALUE = 1_000_000.0
MIN_PROFILE_VALUE = -3.0
MAX_PROFILE_VALUE = 10.0

CENTERLINE_EPSILON = 1.0e-6
LEFT_INCISOR_MARGIN_CM = 0.5
RIGHT_INCISOR_MARGIN_CM = 0.3
MIN_INCISOR_AREA_CM2 = 0.15
TONGUE_TIP_REGION_LENGTH_CM = 2.0

# Tube 主路径常量（自有运行时固定拓扑）
MIN_VELUM_OPENING_CM2 = 1.0e-4

# 迁移历史补偿值（仅用于离线对比脚本，不允许进入生产主路径）
LEGACY_TUBE_HEURISTIC_COMPENSATION = {
    "area_floor_cm2": 0.33,
    "area_log1p_scale": 1.05,
    "expanded_nasal_port_ratio": 0.42,
    "target_to_raw_piecewise_blend": 0.0,
}

NUM_CENTERLINE_POINTS_EXPONENT = 7
NUM_CENTERLINE_POINTS = (1 << NUM_CENTERLINE_POINTS_EXPONENT) + 1

# =============================================================================
# ParamIndex (19D Model)
# =============================================================================

class ParamIndex(IntEnum):
    HX = 0
    HY = 1
    JX = 2
    JA = 3
    LP = 4
    LD = 5
    VS = 6
    VO = 7
    TCX = 8
    TCY = 9
    TTX = 10
    TTY = 11
    TBX = 12
    TBY = 13
    TRX = 14
    TRY = 15
    TS1 = 16
    TS2 = 17
    TS3 = 18
    NUM_PARAMS = 19

# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class ParamDef:
    value: float = 0.0
    limited_value: float = 0.0
    min_val: float = 0.0
    max_val: float = 0.0
    neutral_val: float = 0.0
    abbr: str = ""
    name: str = ""

@dataclass
class Anatomy:
    # Palate
    palate_points: List[Point3D] = field(default_factory=list) # [NUM_PALATE_RIBS]
    palate_angle_deg: List[float] = field(default_factory=list)
    palate_height_cm: List[float] = field(default_factory=list)
    upper_teeth_height_cm: List[float] = field(default_factory=list)
    upper_teeth_width_top_cm: List[float] = field(default_factory=list)
    upper_teeth_width_bottom_cm: List[float] = field(default_factory=list)
    
    # Jaw
    jaw_fulcrum: Point2D = field(default_factory=Point2D)
    jaw_rest_pos: Point2D = field(default_factory=Point2D)
    tooth_root_length_cm: float = 0.0
    jaw_points: List[Point3D] = field(default_factory=list)
    jaw_angle_deg: List[float] = field(default_factory=list)
    jaw_height_cm: List[float] = field(default_factory=list)
    lower_teeth_height_cm: List[float] = field(default_factory=list)
    lower_teeth_width_top_cm: List[float] = field(default_factory=list)
    lower_teeth_width_bottom_cm: List[float] = field(default_factory=list)

    # Tongue
    tongue_tip_radius_cm: float = 0.0
    tongue_center_radius_x_cm: float = 0.0
    tongue_center_radius_y_cm: float = 0.0
    automatic_tongue_root_calc: bool = False
    tongue_root_trx_slope: float = 0.0
    tongue_root_trx_intercept: float = 0.0
    tongue_root_try_slope: float = 0.0
    tongue_root_try_intercept: float = 0.0
    tongue_tonsil_length_cm: float = 0.0
    tongue_tonsil_height_cm: float = 0.0

    # Lips
    lips_width_cm: float = 0.0

    # Velum & Uvula
    uvula_width_cm: float = 0.0
    uvula_height_cm: float = 0.0
    uvula_depth_cm: float = 0.0
    velum_low_points: List[Point2D] = field(default_factory=list)
    velum_mid_points: List[Point2D] = field(default_factory=list)
    velum_high_points: List[Point2D] = field(default_factory=list)

    # Pharynx
    pharynx_fulcrum: Point2D = field(default_factory=Point2D)
    pharynx_rotation_angle_deg: float = 0.0
    pharynx_top_rib_y_cm: float = 0.0
    pharynx_upper_depth_cm: float = 0.0
    pharynx_lower_depth_cm: float = 0.0
    pharynx_back_width_cm: float = 0.0

    # Epiglottis
    epiglottis_width_cm: float = 0.0
    epiglottis_height_cm: float = 0.0
    epiglottis_depth_cm: float = 0.0
    epiglottis_angle_deg: float = 0.0

    # Larynx
    larynx_upper_depth_cm: float = 0.0
    larynx_lower_depth_cm: float = 0.0
    larynx_wide_points: List[Point2D] = field(default_factory=list) # [8]
    larynx_narrow_points: List[Point2D] = field(default_factory=list) # [8]

    # Defaults
    default_f0_hz: float = 120.0

    # Cavities
    piriform_fossa_length_cm: float = 0.0
    piriform_fossa_volume_cm3: float = 0.0
    subglottal_cavity_length_cm: float = 0.0
    nasal_cavity_length_cm: float = 0.0

    # Velocity Factors
    positive_velocity_factor: List[float] = field(default_factory=lambda: [0.0]*19)
    negative_velocity_factor: List[float] = field(default_factory=lambda: [0.0]*19)

@dataclass
class TongueRib:
    point: Point2D = field(default_factory=Point2D)
    left_side_height: float = 0.0
    right_side_height: float = 0.0
    
    # Internal
    left: Point2D = field(default_factory=Point2D)
    right: Point2D = field(default_factory=Point2D)
    normal: Point2D = field(default_factory=Point2D)
    min_x: float = 0.0
    max_x: float = 0.0
    min_y: float = 0.0
    max_y: float = 0.0

@dataclass
class CenterLinePoint:
    point: Point2D = field(default_factory=Point2D)
    normal: Point2D = field(default_factory=Point2D)
    pos: float = 0.0
    min_v: float = 0.0 # renamed from min to avoid conflict
    max_v: float = 0.0

@dataclass
class CrossSection:
    area: float = 0.0
    circ: float = 0.0
    pos: float = 0.0
    articulator: int = 0 # Tube::Articulator enum value


def tongue_side_param_to_min_area_cm2(param_value: float) -> float:
    min_area_cm2 = 0.0
    if param_value > 0.2:
        max_area_cm2 = 0.15
        min_area_cm2 = max_area_cm2 * (param_value - 0.2) / 0.2
        if min_area_cm2 > max_area_cm2:
            min_area_cm2 = max_area_cm2
    elif param_value < -0.2:
        max_area_cm2 = 0.20
        min_area_cm2 = max_area_cm2 * (param_value + 0.2) / (-0.2)
        if min_area_cm2 > max_area_cm2:
            min_area_cm2 = max_area_cm2
    return min_area_cm2
