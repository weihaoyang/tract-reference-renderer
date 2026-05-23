from __future__ import annotations

import asyncio
import json

from websockets.asyncio.client import connect

from vtl_renderer.renderer import neutral_param_vector


async def main() -> None:
    payload = {
        "request_id": "smoke-ws-1",
        "current_tract_params": neutral_param_vector(),
        "render_target": "current",
    }
    async with connect("ws://127.0.0.1:8765") as websocket:
        await websocket.send(json.dumps(payload))
        response = await websocket.recv()
        print(response)


if __name__ == "__main__":
    asyncio.run(main())
