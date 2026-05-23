from __future__ import annotations

import asyncio
import json

import os

from websockets.asyncio.client import connect

from vtl_renderer.renderer import neutral_param_vector


async def main() -> None:
    port = os.environ.get("TRACT_RENDERER_PORT", os.environ.get("VTL_GPL_TRACT_PORT", "8076"))
    payload = {
        "request_id": "smoke-ws-1",
        "current_tract_params": neutral_param_vector(),
        "render_target": "current",
    }
    async with connect(f"ws://127.0.0.1:{port}") as websocket:
        await websocket.send(json.dumps(payload))
        response = await websocket.recv()
        print(response)


if __name__ == "__main__":
    asyncio.run(main())
