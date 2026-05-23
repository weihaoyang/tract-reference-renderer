# -*- coding: utf-8 -*-
"""
核心类型定义 (Core Type Definitions)
"""

from __future__ import annotations

from enum import IntEnum
from typing import Any, List, Tuple

from pydantic import BaseModel, Field

try:
    import torch
except ModuleNotFoundError:  # pragma: no cover - optional runtime dependency
    torch = None  # type: ignore[assignment]

# =============================================================================
# PyTorch 类型别名
# =============================================================================

Scalar = torch.Tensor if torch is not None else Any
Vector = torch.Tensor if torch is not None else Any
Matrix = torch.Tensor if torch is not None else Any
Point3D = Tuple[float, float, float]
Point2D = Tuple[float, float]


# =============================================================================
# 枚举类型
# =============================================================================

class Articulator(IntEnum):
    """发音器官枚举 (对应 Tube::Articulator)"""

    VOCAL_FOLDS = 0
    TONGUE = 1
    LOWER_INCISORS = 2
    LOWER_LIP = 3
    OTHER_ARTICULATOR = 4
    NUM_ARTICULATORS = 5


class GlottisModelType(IntEnum):
    """声门模型类型"""

    GEOMETRIC_2025 = 0
    GEOMETRIC_2019 = 1
    TRIANGULAR = 2
    LF_PULSE = 3


# =============================================================================
# Pydantic 数据模型
# =============================================================================

class GlottisParams(BaseModel):
    """
    通用声门参数容器。
    """

    time_step: float = 1.0 / 48000.0
    sub_pressure: float = 800.0
    lower_pressure: float = 0.0
    upper_pressure: float = 0.0
    supra_pressure: float = 0.0
    f0: float = 120.0
    quality: float = 0.5
    aspiration_level: float = -40.0


class TubeSection(BaseModel):
    """单个管段参数。"""

    index: int
    length: float
    area: float
    articulator: Articulator
    pos_cm: float = 0.0
    circ_cm: float = 0.0
    pressure: float = 0.0
    flow: float = 0.0


class TubeState(BaseModel):
    """管道状态。"""

    sections: List[TubeSection] = Field(default_factory=list)
    incisor_pos_cm: float = 0.0
    tongue_tip_side_elevation: float = 0.0
    velum_opening_cm2: float = 0.0
    subglottal_cavity_length_cm: float = 23.0
    nasal_cavity_length_cm: float = 11.4
    piriform_fossa_length_cm: float = 3.0
    piriform_fossa_volume_cm3: float = 2.0


class CrossSectionData(BaseModel):
    """横截面数据。"""

    pos_cm: float
    area_cm2: float
    perimeter_cm: float
    articulator: Articulator


class TransferFunctionData(BaseModel):
    """传递函数数据容器。"""

    freqs_hz: List[float] = Field(default_factory=list)
    magnitude_db: List[float] = Field(default_factory=list)
    phase_rad: List[float] = Field(default_factory=list)


class SynthesisResult(BaseModel):
    """合成结果。"""

    audio: List[float] = Field(default_factory=list)
    sample_rate_hz: int = 48000
    num_samples: int = 0


class VocalTractShape(BaseModel):
    """声道形状参数 (19D 模型)。"""

    hx: float = Field(0.0, description="Hyoid X")
    hy: float = Field(0.0, description="Hyoid Y")
    jx: float = Field(0.0, description="Jaw X")
    ja: float = Field(0.0, description="Jaw Angle")
    lp: float = Field(0.0, description="Lip Protrusion")
    ld: float = Field(0.0, description="Lip Distance")
    vs: float = Field(0.0, description="Velum Shape")
    vo: float = Field(0.0, description="Velum Opening")
    tcx: float = Field(0.0, description="Tongue Body Center X")
    tcy: float = Field(0.0, description="Tongue Body Center Y")
    ttx: float = Field(0.0, description="Tongue Tip X")
    tty: float = Field(0.0, description="Tongue Tip Y")
    tbx: float = Field(0.0, description="Tongue Back X")
    tby: float = Field(0.0, description="Tongue Back Y")
    trx: float = Field(0.0, description="Tongue Root X")
    try_: float = Field(0.0, alias="try", description="Tongue Root Y")
    ts1: float = Field(0.0, description="Tongue Side 1")
    ts2: float = Field(0.0, description="Tongue Side 2")
    ts3: float = Field(0.0, description="Tongue Side 3")
