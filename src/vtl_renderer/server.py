from __future__ import annotations

import asyncio
import json
import os
from typing import Any, Dict

from websockets.asyncio.server import serve

from .renderer import DEFAULT_HEIGHT_PX, DEFAULT_WIDTH_PX, render_svg_pair

HOST = os.environ.get("TRACT_RENDERER_HOST", os.environ.get("VTL_GPL_TRACT_HOST", "127.0.0.1"))
PORT = int(os.environ.get("TRACT_RENDERER_PORT", os.environ.get("VTL_GPL_TRACT_PORT", "8076")))


def _error_response(request_id: str, error: str) -> Dict[str, Any]:
    return {
        "request_id": request_id,
        "status": "error",
        "width_px": DEFAULT_WIDTH_PX,
        "height_px": DEFAULT_HEIGHT_PX,
        "error": error,
    }


def _handle_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    request_id = payload.get("request_id")
    if not isinstance(request_id, str) or request_id == "":
        return _error_response("", "request_id must be a non-empty string")
    current_tract_params = payload.get("current_tract_params")
    target_tract_params = payload.get("target_tract_params")
    render_target = payload.get("render_target", "current")
    if current_tract_params is None and target_tract_params is not None:
        current_tract_params = target_tract_params
    try:
        result = render_svg_pair(
            current_tract_params=current_tract_params,
            target_tract_params=target_tract_params,
            render_target=render_target,
            width_px=int(payload.get("width_px", DEFAULT_WIDTH_PX)),
            height_px=int(payload.get("height_px", DEFAULT_HEIGHT_PX)),
        )
    except KeyError as exc:
        return _error_response(request_id, f"missing required field: {exc.args[0]}")
    except Exception as exc:
        return _error_response(request_id, str(exc))
    return {
        "request_id": request_id,
        "status": "ok",
        "current_svg": result["current_svg"],
        "target_svg": result["target_svg"],
        "width_px": result["width_px"],
        "height_px": result["height_px"],
        "diagnostics": result["diagnostics"],
    }


async def _ws_handler(websocket) -> None:
    async for raw_message in websocket:
        request_id = ""
        try:
            payload = json.loads(raw_message)
            if isinstance(payload, dict):
                request_id = str(payload.get("request_id", ""))
                response = _handle_payload(payload)
            else:
                response = _error_response("", "payload must be a JSON object")
        except json.JSONDecodeError as exc:
            response = _error_response(request_id, f"invalid JSON: {exc.msg}")
        await websocket.send(json.dumps(response, ensure_ascii=False))


async def run_server(host: str = HOST, port: int = PORT) -> None:
    async with serve(_ws_handler, host, port):
        await asyncio.Future()


def main() -> None:
    asyncio.run(run_server())


if __name__ == "__main__":
    main()
