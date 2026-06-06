from __future__ import annotations

import asyncio
import json

import os

from websockets.asyncio.client import connect

from tract_reference_renderer.renderer import neutral_param_vector


async def main() -> None:
    port = os.environ.get("TRACT_RENDERER_PORT", "8076")
    payload = {
        "request_id": "smoke-ws-1",
        "current_tract_params": neutral_param_vector(),
        "render_target": "current",
    }
    drag_payload = {
        "type": "drag_solve",
        "request_id": "smoke-drag-1",
        "current_tract_params": [0.5] * 16 + [0.0, 0.0, 0.0],
        "control_id": "tongue_body_handle",
        "x_cm": 0.52,
        "y_cm": 0.62,
        "z_cm": 0.0,
    }
    async with connect(f"ws://127.0.0.1:{port}") as websocket:
        await websocket.send(json.dumps(payload))
        response = await websocket.recv()
        print(response)
        await websocket.send(json.dumps(drag_payload))
        drag_response = await websocket.recv()
        print(drag_response)


if __name__ == "__main__":
    asyncio.run(main())
