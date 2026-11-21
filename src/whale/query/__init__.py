"""Read-only queries over White Whale's committed results, and the `whale mcp` server
(PLAN.md §3, §8). `tools` holds the questions as plain functions; `server` wraps them
for MCP clients and needs the optional `mcp` extra."""
from whale.query.store import Store
from whale.query.tools import QueryError

__all__ = ["Store", "QueryError"]
