//! The server over its real transport: initialize, then list the tools.

use std::io::{BufRead, BufReader, Write};
use std::process::{Command, Stdio};

#[test]
fn the_server_answers_over_stdio() {
    let mut child = Command::new(env!("CARGO_BIN_EXE_gusset-mcp"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    let mut stdout = BufReader::new(child.stdout.take().unwrap());
    let mut ask = |line: &str| {
        writeln!(stdin, "{line}").unwrap();
        let mut answer = String::new();
        if line.contains(r#""id""#) {
            stdout.read_line(&mut answer).unwrap();
        }
        answer
    };

    let init = ask(
        r#"{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"test","version":"0"}}}"#,
    );
    assert!(init.contains(r#""name":"gusset""#), "{init}");
    ask(r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#);
    let tools = ask(r#"{"jsonrpc":"2.0","id":2,"method":"tools/list"}"#);
    for name in ["build_part", "flat_svg", "part_schema"] {
        assert!(tools.contains(&format!(r#""name":"{name}""#)), "{tools}");
    }

    drop(stdin);
    assert!(child.wait().unwrap().success());
}
