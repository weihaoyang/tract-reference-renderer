from __future__ import annotations

HELPER_NAME = "tract-reference-renderer"
HELPER_VERSION = "0.1.0"
PROTOCOL_VERSION = "coach.reference_renderer.v1"
LICENSE_STATUS = "external_gpl_helper"
DEFAULT_IPC_BOUNDARY = "localhost_websocket"


def build_health_payload() -> dict[str, object]:
    return {
        "status": "ok",
        "helper_name": HELPER_NAME,
        "helper_version": HELPER_VERSION,
        "protocol_version": PROTOCOL_VERSION,
        "license_status": LICENSE_STATUS,
        "ipc_boundary": DEFAULT_IPC_BOUNDARY,
        "truth_tier": "reference_visualization_not_patient_truth",
        "clinical_truth_claim_allowed": False,
    }


__all__ = [
    "DEFAULT_IPC_BOUNDARY",
    "HELPER_NAME",
    "HELPER_VERSION",
    "LICENSE_STATUS",
    "PROTOCOL_VERSION",
    "build_health_payload",
]
