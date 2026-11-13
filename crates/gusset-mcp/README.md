# gusset-mcp

Gusset as an [MCP](https://modelcontextprotocol.io) server, so AI agents can model sheet
metal parts: the agent writes a part description (JSON, `schema/part.schema.json`), and
Gusset builds it and answers with the flat pattern, the folded part's mass properties, the
bend table and the shop checks, or with the kernel's message when it can't build the part.
Runs over stdio.

## Tools

| Tool | Input | What it returns |
|---|---|---|
| `build_part` | `part` (object or JSON text), optional `step_path`, `dxf_path`, `svg_path` | JSON: flat size, bounds and area; folded volume, area and centre of mass; bends (order, direction, angle, inner radius, K-factor, length, note); shop-check findings; files written. A part that can't be built (including one that folds into itself, with where) comes back as a tool error with the kernel's message. |
| `flat_svg` | `part` | The flat pattern as SVG text; areas that fold into each other marked with class `clash`. |
| `part_schema` | none | The part description's JSON Schema. |

The schema is also a resource, `gusset://schema/part`. Millimetres and degrees throughout.
File paths are resolved by the server process; pass absolute paths.

## Build

    cargo build -p gusset-mcp --release

The binary is `target/release/gusset-mcp` (`gusset-mcp.exe` on Windows).

## Use

Claude Code:

    claude mcp add gusset -- <path>/target/release/gusset-mcp.exe

Any client with an `mcpServers` config:

```json
{
  "mcpServers": {
    "gusset": {
      "command": "<path>/target/release/gusset-mcp.exe"
    }
  }
}
```

Then ask for a part, for example "a 100 x 60 x 2 mm plate with an 8 mm hole in the middle
and a 30 mm flange up along one long edge; give me the flat size and a STEP file".
