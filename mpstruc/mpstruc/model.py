"""Typed mpstruc model built from the lossless tree, and its exact inverse.

Every leaf field keeps its raw source text (None for an empty <tag/>), so
`to_xml(load(x)) == x` byte-for-byte. Display text is derived elsewhere
(mpstruc.text); nothing here decodes or cleans anything.

Record keys follow the reference xlsx: G## (group), S### (subgroup, numbered
within its group), P#### (parent within subgroup), D### (daughter within
parent), R### (related entry within its owner).
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

from .xmlparse import MpstrucParseError, Node, normalize, parse, serialize

BIB_FIELDS = ("pubMedId", "authors", "year", "title", "journal",
              "volume", "issue", "pages", "doi", "notes")
ENTRY_FIELDS = ("pdbCode", "name", "species", "taxonomicDomain",
                "expressedInSpecies", "resolution", "description")

TOPOLOGY = {
    "MONOTOPIC MEMBRANE PROTEINS": "Monotopic",
    "TRANSMEMBRANE PROTEINS: BETA-BARREL": "Beta-barrel transmembrane",
    "TRANSMEMBRANE PROTEINS: ALPHA-HELICAL": "Alpha-helical transmembrane",
}


@dataclass
class Bibliography:
    pubMedId: str | None = None
    authors: str | None = None
    year: str | None = None
    title: str | None = None
    journal: str | None = None
    volume: str | None = None
    issue: str | None = None
    pages: str | None = None
    doi: str | None = None
    notes: str | None = None


@dataclass
class Related:
    pdb: str
    key: str
    owner_key: str


@dataclass
class Entry:
    """A <protein> (parent) or <memberProtein> (daughter)."""
    kind: str  # "parent" | "daughter"
    key: str
    pdbCode: str
    name: str | None
    species: str | None
    taxonomicDomain: str | None
    expressedInSpecies: str | None
    resolution: str | None
    description: str | None
    bibliography: Bibliography
    secondary: list[Bibliography]
    related: list[Related]
    masterProteinPdbCode: str | None = None  # daughters only
    daughters: list["Entry"] = field(default_factory=list)  # parents only
    parent_key: str | None = None  # daughters only


@dataclass
class Subgroup:
    key: str
    name: str
    group_key: str
    parents: list[Entry]


@dataclass
class Group:
    key: str
    name: str
    subgroups: list[Subgroup]

    @property
    def topology(self) -> str:
        return TOPOLOGY.get(self.name, self.name)


@dataclass
class Database:
    attrs: str  # raw root attribute text
    caption: str | None
    groups: list[Group]
    sha256: str = ""  # of the bytes as supplied
    normalizations: list[str] = field(default_factory=list)

    @property
    def meta(self) -> dict[str, str]:
        return Node("mpstruc", "elem", attrs=self.attrs).attr_dict()

    def subgroups(self):
        for g in self.groups:
            yield from g.subgroups

    def entries(self):
        """Parents and daughters in document order."""
        for sg in self.subgroups():
            for p in sg.parents:
                yield p
                yield from p.daughters


# ---------------------------------------------------------------- reading

def _no_attrs(n: Node) -> None:
    if n.attrs.strip():
        raise MpstrucParseError(f"<{n.tag}> has unexpected attributes {n.attrs.strip()!r}")
    for c in n.children:
        _no_attrs(c)


def _expect(node: Node, tags: tuple[str, ...]) -> list[Node]:
    got = tuple(c.tag for c in node.children)
    if got != tags:
        raise MpstrucParseError(f"<{node.tag}> has children {got}, expected {tags}")
    return node.children


def _leaf(n: Node) -> str | None:
    if n.kind == "empty":
        return None
    if n.kind != "leaf":
        raise MpstrucParseError(f"<{n.tag}> should be a text field")
    return n.text


def _list(n: Node, child_tag: str) -> list[Node]:
    if n.kind == "empty":
        return []
    if n.kind != "elem" or not n.children:
        raise MpstrucParseError(f"<{n.tag}> should be a non-empty list or <{n.tag}/>")
    for c in n.children:
        if c.tag != child_tag:
            raise MpstrucParseError(f"<{n.tag}> contains <{c.tag}>, expected <{child_tag}>")
    return n.children


def _bib(n: Node) -> Bibliography:
    kids = _expect(n, BIB_FIELDS)
    return Bibliography(*(_leaf(k) for k in kids))


def _entry(n: Node, kind: str, key: str) -> Entry:
    if kind == "parent":
        tags = ENTRY_FIELDS + ("bibliography", "secondaryBibliographies",
                               "relatedPdbEntries", "memberProteins")
    else:
        tags = (("pdbCode", "masterProteinPdbCode") + ENTRY_FIELDS[1:]
                + ("bibliography", "secondaryBibliographies", "relatedPdbEntries"))
    k = dict(zip(tags, _expect(n, tags)))
    pdb = _leaf(k["pdbCode"])
    if not pdb or not pdb.strip():
        raise MpstrucParseError(f"{key}: missing pdbCode")
    related = [_leaf(r) for r in _list(k["relatedPdbEntries"], "pdbCode")]
    if any(not r or not r.strip() for r in related):
        raise MpstrucParseError(f"{key}: empty related pdbCode")
    e = Entry(
        kind=kind, key=key, pdbCode=pdb,
        **{f: _leaf(k[f]) for f in ENTRY_FIELDS[1:]},
        bibliography=_bib(k["bibliography"]),
        secondary=[_bib(b) for b in _list(k["secondaryBibliographies"], "bibliography")],
        related=[Related(r, f"{key}-R{i:03d}", key) for i, r in enumerate(related, 1)],
    )
    if kind == "daughter":
        e.masterProteinPdbCode = _leaf(k["masterProteinPdbCode"])
    else:
        e.daughters = [_entry(d, "daughter", f"{key}-D{i:03d}")
                       for i, d in enumerate(_list(k["memberProteins"], "memberProtein"), 1)]
        for d in e.daughters:
            d.parent_key = key
    return e


def from_tree(root: Node) -> Database:
    if root.tag != "mpstruc":
        raise MpstrucParseError(f"root is <{root.tag}>")
    for c in root.children:
        _no_attrs(c)
    caption, groups = _expect(root, ("caption", "groups"))
    db = Database(attrs=root.attrs, caption=_leaf(caption), groups=[])
    for gi, gn in enumerate(_list(groups, "group"), 1):
        name, gprot, subs = _expect(gn, ("name", "proteins", "subgroups"))
        if _list(gprot, "protein"):
            raise MpstrucParseError("group-level <proteins> is expected to be empty")
        g = Group(f"G{gi:02d}", _leaf(name) or "", [])
        for s_no, sn in enumerate(_list(subs, "subgroup"), 1):
            sname, prots = _expect(sn, ("name", "proteins"))
            skey = f"{g.key}-S{s_no:03d}"
            g.subgroups.append(Subgroup(
                skey, _leaf(sname) or "", g.key,
                [_entry(p, "parent", f"{skey}-P{pi:04d}")
                 for pi, p in enumerate(_list(prots, "protein"), 1)]))
        db.groups.append(g)
    return db


def load(path_or_bytes: str | Path | bytes) -> Database:
    data = path_or_bytes if isinstance(path_or_bytes, bytes) else Path(path_or_bytes).read_bytes()
    clean, notes = normalize(data)
    db = from_tree(parse(clean))
    db.sha256 = hashlib.sha256(data).hexdigest()
    db.normalizations = notes
    return db


# ---------------------------------------------------------------- writing

def _l(tag: str, v: str | None) -> Node:
    return Node(tag, "empty") if v is None else Node(tag, "leaf", text=v)


def _lst(tag: str, children: list[Node]) -> Node:
    return Node(tag, "elem", children=children) if children else Node(tag, "empty")


def _bib_node(b: Bibliography) -> Node:
    return Node("bibliography", "elem", children=[_l(f, getattr(b, f)) for f in BIB_FIELDS])


def _entry_node(e: Entry) -> Node:
    kids = [_l("pdbCode", e.pdbCode)]
    if e.kind == "daughter":
        kids.append(_l("masterProteinPdbCode", e.masterProteinPdbCode))
    kids += [_l(f, getattr(e, f)) for f in ENTRY_FIELDS[1:]]
    kids += [_bib_node(e.bibliography),
             _lst("secondaryBibliographies", [_bib_node(b) for b in e.secondary]),
             _lst("relatedPdbEntries", [_l("pdbCode", r.pdb) for r in e.related])]
    if e.kind == "parent":
        kids.append(_lst("memberProteins", [_entry_node(d) for d in e.daughters]))
    return Node("memberProtein" if e.kind == "daughter" else "protein", "elem", children=kids)


def to_tree(db: Database) -> Node:
    groups = [Node("group", "elem", children=[
        _l("name", g.name), Node("proteins", "empty"),
        _lst("subgroups", [Node("subgroup", "elem", children=[
            _l("name", sg.name),
            _lst("proteins", [_entry_node(p) for p in sg.parents])]) for sg in g.subgroups]),
    ]) for g in db.groups]
    return Node("mpstruc", "elem", attrs=db.attrs,
                children=[_l("caption", db.caption), _lst("groups", groups)])


def to_xml(db: Database) -> str:
    return serialize(to_tree(db))
