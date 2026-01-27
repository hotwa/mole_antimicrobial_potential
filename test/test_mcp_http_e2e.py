import os
import time

import httpx

try:
    import pytest
except ImportError:  # pragma: no cover
    class _PytestStub:
        @staticmethod
        def fail(message: str) -> None:
            raise AssertionError(message)

    pytest = _PytestStub()

MCP_URL = os.getenv("MCP_URL", "http://localhost:8001/mcp")
TIMEOUT_S = float(os.getenv("MCP_E2E_TIMEOUT", "60"))
_SESSION_ID: str | None = None
_INITIALIZED = False


def _post_jsonrpc(payload: dict) -> dict:
    global _SESSION_ID
    deadline = time.monotonic() + TIMEOUT_S
    while True:
        headers = {"Accept": "application/json"}
        if _SESSION_ID:
            headers["mcp-session-id"] = _SESSION_ID
        response = httpx.post(
            MCP_URL,
            json=payload,
            timeout=10.0,
            headers=headers,
        )
        session_id = response.headers.get("mcp-session-id")
        if session_id:
            _SESSION_ID = session_id
        if response.status_code == 400 and _SESSION_ID and "Missing session ID" in response.text:
            if time.monotonic() >= deadline:
                pytest.fail("Timed out retrying after session negotiation")
            continue
        if response.status_code == 202:
            if time.monotonic() >= deadline:
                pytest.fail("Timed out waiting for 202 Accepted to resolve")
            time.sleep(0.25)
            continue
        response.raise_for_status()
        return response.json()


def _send_notification(payload: dict) -> None:
    global _SESSION_ID
    deadline = time.monotonic() + TIMEOUT_S
    while True:
        headers = {"Accept": "application/json"}
        if _SESSION_ID:
            headers["mcp-session-id"] = _SESSION_ID
        response = httpx.post(
            MCP_URL,
            json=payload,
            timeout=10.0,
            headers=headers,
        )
        session_id = response.headers.get("mcp-session-id")
        if session_id:
            _SESSION_ID = session_id
        if response.status_code in (200, 202, 204):
            return
        if response.status_code == 400 and _SESSION_ID and "Missing session ID" in response.text:
            if time.monotonic() >= deadline:
                pytest.fail("Timed out retrying notification after session negotiation")
            continue
        response.raise_for_status()


def _ensure_initialized() -> None:
    global _INITIALIZED
    if _INITIALIZED:
        return
    init_payload = {
        "jsonrpc": "2.0",
        "id": "init",
        "method": "initialize",
        "params": {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {"name": "mcp-http-e2e", "version": "0.1.0"},
        },
    }
    _post_jsonrpc(init_payload)
    _send_notification(
        {"jsonrpc": "2.0", "method": "notifications/initialized", "params": {}}
    )
    _INITIALIZED = True


def _call(method: str, params: dict | None = None, request_id: str = "1") -> dict:
    _ensure_initialized()
    payload = {
        "jsonrpc": "2.0",
        "id": request_id,
        "method": method,
    }
    if params is not None:
        payload["params"] = params
    data = _post_jsonrpc(payload)
    assert data.get("jsonrpc") == "2.0"
    assert data.get("id") == request_id
    if "error" in data:
        pytest.fail(f"JSON-RPC error: {data['error']}")
    return data["result"]


def _tool_call(name: str, arguments: dict | None = None, request_id: str = "tool") -> dict:
    result = _call(
        "tools/call",
        params={"name": name, "arguments": arguments or {}},
        request_id=request_id,
    )
    structured = result.get("structuredContent")
    if structured is None:
        pytest.fail("Missing structuredContent in tool response")
    return structured


def test_tools_list():
    result = _call("tools/list", request_id="tools")
    tools = result.get("tools", [])
    names = {tool.get("name") for tool in tools}
    assert "status" in names
    assert "predict_antimicrobial_potential" in names


def test_resources_list():
    result = _call("resources/list", request_id="resources")
    resources = result.get("resources", [])
    uris = {resource.get("uri") for resource in resources}
    assert "resource://schema" in uris
    assert "resource://about" in uris
    assert "resource://strains" in uris


def test_resources_read_schema():
    result = _call("resources/read", params={"uri": "resource://schema"}, request_id="schema")
    contents = result.get("contents", [])
    assert contents
    text = contents[0].get("text")
    assert text


def test_resources_read_about():
    result = _call("resources/read", params={"uri": "resource://about"}, request_id="about")
    contents = result.get("contents", [])
    assert contents
    text = contents[0].get("text")
    assert text


def test_resources_read_strains():
    result = _call("resources/read", params={"uri": "resource://strains"}, request_id="strains")
    contents = result.get("contents", [])
    assert contents
    text = contents[0].get("text")
    assert text


def test_status_tool():
    status = _tool_call("status", request_id="status")
    required_keys = {"service", "version", "loaded", "device"}
    assert required_keys.issubset(status.keys())


def test_predict_tool_aggregate():
    result = _tool_call(
        "predict_antimicrobial_potential",
        arguments={"smiles": "CCO", "aggregate_scores": True},
        request_id="predict",
    )
    assert result["mode"] == "aggregate"
    items = result["items"]
    assert items
    sample = items[0]
    assert "chem_id" in sample
    assert "apscore_total" in sample
    assert "ginhib_total" in sample
    assert "broad_spectrum" in sample


if __name__ == "__main__":
    test_tools_list()
    test_resources_list()
    test_resources_read_schema()
    test_resources_read_about()
    test_resources_read_strains()
    test_status_tool()
    test_predict_tool_aggregate()
