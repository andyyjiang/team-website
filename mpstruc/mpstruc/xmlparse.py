"""Lossless reader/writer for the mpstruc XML export.

The export is not well-formed XML: field text contains HTML named entities
(&Aring;, &nbsp;, ...), misspelled entities (&Ariing;, &alfa;) and broken
markup (<ahref='...'>, unclosed <a>, <super>). A standard XML parser rejects
it, and a "recover" parser silently mangles it.

The file is, however, perfectly line-regular: every line is exactly one of
  <tag attr="...">   start tag
  </tag>             end tag
  <tag/>             empty element
  <tag>RAW</tag>     leaf element whose RAW content is kept byte-for-byte
Lines are joined with "\n" and there is no trailing newline. This module
parses that grammar into a generic tree and writes it back byte-identically.

`normalize` accepts harmless variants of the same export (UTF-8 BOM, CRLF,
an XML declaration, a trailing newline) and reports what it changed; the
byte-identical round trip is then stated against the normalised bytes.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

_EMPTY = re.compile(r"<(\w+)((?: [^>]*)?)/>\Z")
_LEAF = re.compile(r"<(\w+)>(.*)</\1>\Z", re.S)
_END = re.compile(r"</(\w+)>\Z")
_START = re.compile(r"<(\w+)((?: [^>]*)?)>\Z")
_ATTR = re.compile(r'(\w+)="([^"]*)"')


class MpstrucParseError(ValueError):
    pass


def normalize(data: bytes) -> tuple[bytes, list[str]]:
    notes = []
    if data.startswith(b"\xef\xbb\xbf"):
        data, notes = data[3:], notes + ["removed UTF-8 BOM"]
    if b"\r\n" in data:
        data, notes = data.replace(b"\r\n", b"\n"), notes + ["converted CRLF line endings"]
    if data.startswith(b"<?xml"):
        data, notes = data.split(b"\n", 1)[1], notes + ["removed XML declaration"]
    if data.endswith(b"\n"):
        data, notes = data.rstrip(b"\n"), notes + ["removed trailing newline"]
    return data, notes


@dataclass
class Node:
    tag: str
    kind: str  # "elem" | "empty" | "leaf"
    attrs: str = ""  # raw attribute text including leading space, kept verbatim
    text: str | None = None  # raw inner content for leaf nodes
    children: list["Node"] = field(default_factory=list)

    def attr_dict(self) -> dict[str, str]:
        return dict(_ATTR.findall(self.attrs))


def parse(data: str | bytes) -> Node:
    if isinstance(data, bytes):
        data = data.decode("utf-8")
    stack: list[Node] = []
    root: Node | None = None
    for lineno, line in enumerate(data.split("\n"), 1):
        if m := _EMPTY.match(line):
            node = Node(m[1], "empty", attrs=m[2])
        elif m := _LEAF.match(line):
            if f"</{m[1]}>" in m[2] or re.search(rf"<{m[1]}[\s>/]", m[2]):
                raise MpstrucParseError(f"line {lineno}: more than one <{m[1]}> on one line")
            node = Node(m[1], "leaf", text=m[2])
        elif m := _END.match(line):
            if not stack or stack[-1].tag != m[1]:
                open_tag = stack[-1].tag if stack else None
                raise MpstrucParseError(f"line {lineno}: </{m[1]}> closes <{open_tag}>")
            done = stack.pop()
            if not stack:
                root = done
            continue
        elif m := _START.match(line):
            node = Node(m[1], "elem", attrs=m[2])
            if stack:
                stack[-1].children.append(node)
            elif root is not None:
                raise MpstrucParseError(f"line {lineno}: second root element")
            stack.append(node)
            continue
        else:
            raise MpstrucParseError(f"line {lineno}: unrecognised line {line[:80]!r}")
        if not stack:
            raise MpstrucParseError(f"line {lineno}: content outside root")
        stack[-1].children.append(node)
    if stack or root is None:
        raise MpstrucParseError("unexpected end of document")
    return root


def serialize(root: Node) -> str:
    out: list[str] = []

    def walk(n: Node) -> None:
        if n.kind == "empty":
            out.append(f"<{n.tag}{n.attrs}/>")
        elif n.kind == "leaf":
            out.append(f"<{n.tag}>{n.text}</{n.tag}>")
        else:
            out.append(f"<{n.tag}{n.attrs}>")
            for c in n.children:
                walk(c)
            out.append(f"</{n.tag}>")

    walk(root)
    return "\n".join(out)
