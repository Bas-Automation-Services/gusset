//! Gusset, an open sheet metal kernel, as an MCP server over stdio: an agent
//! describes a part as JSON and gets its flat pattern, folded mass, bend
//! table, shop checks and STEP, DXF or SVG files back.

/// The bracket from the README: a plate with a hole and one flange.
macro_rules! bracket {
    () => {
        r#"{"id": "bracket", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
 "base": {"outline": [[0, 0], [100, 0], [100, 60], [0, 60]],
          "holes": [{"center": [50, 30], "diameter": 8}]},
 "flanges": [{"id": "f", "on": {"base": 2}, "angle": 90, "length": 30,
              "datum": "outer", "position": "outside", "side": "up"}]}"#
    };
}

/// The part format in brief, for the server's instructions and build_part.
macro_rules! format_help {
    () => {
        concat!(
            "A part is JSON: a base plate (an outline of [x, y] points, holes, cut-outs) and the ",
            "features built on it (flanges on base edges, hems, contour flanges, folds, faces ",
            "joined by bends, punches, embosses). Millimetres and degrees. A minimal part, a ",
            "100 x 60 plate with an 8 mm hole and a 90 degree flange 30 mm high on base edge 2:\n",
            bracket!(),
            "\nCall part_schema for the full description."
        )
    };
}

mod part;

use part::{build, guarded, parse, svg, Outputs, SCHEMA};
use rmcp::handler::server::{router::tool::ToolRouter, wrapper::Parameters};
use rmcp::model::{
    Implementation, ListResourcesResult, PaginatedRequestParams, ReadResourceRequestParams,
    ReadResourceResponse, ReadResourceResult, Resource, ResourceContents, ServerCapabilities,
    ServerConfig,
};
use rmcp::schemars::{self, json_schema, JsonSchema, Schema, SchemaGenerator};
use rmcp::service::RequestContext;
use rmcp::{tool, tool_handler, tool_router, ErrorData, RoleServer, ServerHandler, ServiceExt};
use serde::Deserialize;
use serde_json::Value;

const SCHEMA_URI: &str = "gusset://schema/part";

fn part_field(_: &mut SchemaGenerator) -> Schema {
    json_schema!({
        "type": ["object", "string"],
        "description": "The part description (see part_schema), as a JSON object or as JSON text."
    })
}

#[derive(Deserialize, JsonSchema)]
struct PartArgs {
    #[schemars(schema_with = "part_field")]
    part: Value,
}

#[derive(Deserialize, JsonSchema)]
struct BuildArgs {
    #[schemars(schema_with = "part_field")]
    part: Value,
    /// Write the folded part as STEP (AP214, mm) to this path; absolute is safest.
    step_path: Option<String>,
    /// Write the flat pattern as DXF (R12, mm, bend lines and notes) to this path.
    dxf_path: Option<String>,
    /// Write the flat pattern as SVG to this path.
    svg_path: Option<String>,
}

#[derive(Clone)]
struct Gusset {
    tool_router: ToolRouter<Self>,
}

#[tool_router]
impl Gusset {
    #[tool(description = concat!(
        "Build a sheet metal part and sum it up: flat pattern size and area, folded volume, ",
        "surface area and centre of mass, the bend table (order, direction, angle, inner radius) ",
        "and shop-check findings (bend radius, leg length, hole distances and size). Optionally ",
        "writes STEP, flat DXF and flat SVG files. A part Gusset can't build (a mistake in the ",
        "description, something it doesn't build yet, a part that folds into itself) fails with ",
        "the kernel's message: fix the part and call again.\n\n",
        format_help!()
    ))]
    fn build_part(&self, Parameters(a): Parameters<BuildArgs>) -> Result<String, String> {
        guarded(|| {
            let spec = parse(&a.part)?;
            let out = Outputs {
                step: a.step_path,
                dxf: a.dxf_path,
                svg: a.svg_path,
            };
            Ok(serde_json::to_string_pretty(&build(&spec, &out)?).unwrap())
        })
    }

    #[tool(description = concat!(
        "The flat pattern as SVG text (mm): outline, holes, cut-outs, bend lines with their ",
        "notes, punch and emboss marks. Areas that fold into each other are marked with class ",
        "\"clash\". Takes the same part as build_part."
    ))]
    fn flat_svg(&self, Parameters(a): Parameters<PartArgs>) -> Result<String, String> {
        guarded(|| svg(&parse(&a.part)?))
    }

    #[tool(description = "The part description's full JSON Schema: every feature and option.")]
    fn part_schema(&self) -> String {
        SCHEMA.into()
    }
}

#[tool_handler(router = self.tool_router)]
impl ServerHandler for Gusset {
    fn get_info(&self) -> ServerConfig {
        ServerConfig::new(
            ServerCapabilities::builder()
                .enable_tools()
                .enable_resources()
                .build(),
        )
        .with_server_info(Implementation::new("gusset", env!("CARGO_PKG_VERSION")))
        .with_instructions(concat!(
            "Gusset models sheet metal parts exactly: folded solid, flat pattern, STEP/DXF/SVG. ",
            format_help!()
        ))
    }

    async fn list_resources(
        &self,
        _: Option<PaginatedRequestParams>,
        _: RequestContext<RoleServer>,
    ) -> Result<ListResourcesResult, ErrorData> {
        Ok(ListResourcesResult::with_all_items(vec![Resource::new(
            SCHEMA_URI,
            "part.schema.json",
        )
        .with_description("JSON Schema of a Gusset part description")
        .with_mime_type("application/schema+json")]))
    }

    async fn read_resource(
        &self,
        request: ReadResourceRequestParams,
        _: RequestContext<RoleServer>,
    ) -> Result<ReadResourceResponse, ErrorData> {
        if request.uri != SCHEMA_URI {
            return Err(ErrorData::resource_not_found(
                format!("no resource {}", request.uri),
                None,
            ));
        }
        Ok(ReadResourceResult::new(vec![
            ResourceContents::text(SCHEMA, SCHEMA_URI).with_mime_type("application/schema+json")
        ])
        .into())
    }
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let server = Gusset {
        tool_router: Gusset::tool_router(),
    }
    .serve(rmcp::transport::stdio())
    .await?;
    server.waiting().await?;
    Ok(())
}
