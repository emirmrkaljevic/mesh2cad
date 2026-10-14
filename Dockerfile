# Stdio MCP server for Glama / local Docker use.
FROM python:3.12-slim

WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src ./src

RUN pip install --no-cache-dir ".[mcp]"

ENTRYPOINT ["mesh2cad", "mcp"]
