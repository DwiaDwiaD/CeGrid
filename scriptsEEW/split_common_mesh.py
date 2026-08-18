#!/usr/bin/env python3
"""Split one Gmsh 2D MSH2 mesh into unfilled and filled views.

The input must be a single mesh generated once by Gmsh with the airfoil
surface (elementary surface 99) present.  The splitter deliberately rebuilds
both outputs so that:

* every element not on surface 99 is common to both outputs;
* common elements appear first and in exactly the same order in both files;
* common nodes appear first and in exactly the same order in both files;
* filled-only nodes/elements are appended to the filled output;
* node and element tags are rewritten sequentially from those stable orders;
* common physical groups are retained;
* the unfilled output gets an "airfoil" physical curve group when curves 1,2
  are present, while the filled output does not expose that group.

This avoids depending on the numbering chosen by two independent Gmsh runs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass
class Node:
    tag: int
    x: float
    y: float
    z: float


@dataclass
class Element:
    old_tag: int
    etype: int
    tags: list[int]
    nodes: list[int]
    original_index: int


def read_msh2(path: Path):
    lines = path.read_text().splitlines()
    try:
        i_nodes = lines.index("$Nodes")
        i_elems = lines.index("$Elements")
    except ValueError as exc:
        raise ValueError(f"{path}: expected MSH2 $Nodes/$Elements sections") from exc

    physical_names: list[tuple[int, int, str]] = []
    if "$PhysicalNames" in lines:
        i = lines.index("$PhysicalNames")
        n = int(lines[i + 1])
        for row in lines[i + 2:i + 2 + n]:
            dim, tag, name = row.split(maxsplit=2)
            physical_names.append((int(dim), int(tag), name))

    n_nodes = int(lines[i_nodes + 1])
    nodes: dict[int, Node] = {}
    for row in lines[i_nodes + 2:i_nodes + 2 + n_nodes]:
        p = row.split()
        if len(p) != 4:
            raise ValueError(f"Malformed node row in {path}: {row!r}")
        tag = int(p[0])
        nodes[tag] = Node(tag, float(p[1]), float(p[2]), float(p[3]))

    n_elems = int(lines[i_elems + 1])
    elems: list[Element] = []
    for idx, row in enumerate(lines[i_elems + 2:i_elems + 2 + n_elems]):
        p = row.split()
        if len(p) < 3:
            raise ValueError(f"Malformed element row in {path}: {row!r}")
        etype = int(p[1])
        ntags = int(p[2])
        if len(p) < 3 + ntags:
            raise ValueError(f"Malformed element tags in {path}: {row!r}")
        tags = [int(x) for x in p[3:3 + ntags]]
        node_tags = [int(x) for x in p[3 + ntags:]]
        elems.append(Element(int(p[0]), etype, tags, node_tags, idx))

    return physical_names, nodes, elems


def elementary_entity(e: Element) -> int | None:
    # In MSH2, the first two optional tags are normally physical and
    # elementary entity tags.  The current CeGrid outputs use this layout.
    if len(e.tags) >= 2:
        return e.tags[1]
    return None


def rewrite_element_tags(e: Element, physical_tags: list[int], entity_tag: int | None) -> list[int]:
    # Preserve the elementary entity as the second tag whenever it exists.
    tags = list(physical_tags)
    if entity_tag is not None:
        tags.append(entity_tag)
    return tags


def write_msh2(
    path: Path,
    physical_names: list[tuple[int, int, str]],
    nodes: list[Node],
    elements: list[Element],
    element_physical_override: dict[int, list[int]] | None = None,
):
    element_physical_override = element_physical_override or {}
    node_id = {n.tag: i + 1 for i, n in enumerate(nodes)}

    out: list[str] = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat"]

    out += ["$PhysicalNames", str(len(physical_names))]
    out += [f'{d} {t} {name}' for d, t, name in physical_names]
    out += ["$EndPhysicalNames"]

    out += ["$Nodes", str(len(nodes))]
    for i, n in enumerate(nodes, 1):
        out.append(f"{i} {n.x:.17g} {n.y:.17g} {n.z:.17g}")
    out += ["$EndNodes"]

    out += ["$Elements", str(len(elements))]
    for i, e in enumerate(elements, 1):
        phys = element_physical_override.get(e.original_index)
        if phys is None:
            old_phys = e.tags[0:1] if e.tags else []
            phys = old_phys
        entity = elementary_entity(e)
        tags = rewrite_element_tags(e, phys, entity)
        conn = [node_id[n] for n in e.nodes]
        out.append(f"{i} {e.etype} {len(tags)} {' '.join(map(str, tags))} {' '.join(map(str, conn))}")
    out += ["$EndElements", ""]
    path.write_text("\n".join(out))


def build_outputs(src: Path, unfilled: Path, filled: Path):
    physical_names, nodes, elements = read_msh2(src)

    # Surface 99 is the airfoil-filled 2D patch. All other mesh elements are
    # common between the two views.
    surface99 = [
        e for e in elements
        if e.etype in (2, 3) and elementary_entity(e) == 99
    ]
    common = [e for e in elements if e not in surface99]
    filled_extra = surface99

    # Stable element order: common elements first; filled-only surface 99
    # elements last. The relative original order is preserved within each set.
    common.sort(key=lambda e: e.original_index)
    filled_extra.sort(key=lambda e: e.original_index)

    filled_elements = common + filled_extra

    common_node_tags = set(n for e in common for n in e.nodes)
    filled_node_tags = set(n for e in filled_elements for n in e.nodes)

    common_nodes = [nodes[tag] for tag in sorted(common_node_tags)]
    filled_only_nodes = [nodes[tag] for tag in sorted(filled_node_tags - common_node_tags)]
    filled_nodes = common_nodes + filled_only_nodes

    # Keep the source physical names, but add the airfoil curve group to the
    # unfilled mesh. Physical group tag 8 is unused by the supplied source file.
    common_names = list(physical_names)
    if not any(dim == 1 and name.strip('"') == 'airfoil' for dim, _, name in common_names):
        common_names.append((1, 8, '"airfoil"'))

    # In the current full/filled output the airfoil curves are not in a
    # physical group. For the unfilled output, assign curves on elementary
    # entities 1 and 2 to physical group 8. For filled, leave them unassigned.
    unfilled_overrides: dict[int, list[int]] = {}
    for e in common:
        ent = elementary_entity(e)
        if e.etype == 1 and ent in (1, 2):
            unfilled_overrides[e.original_index] = [8]

    # Filled output removes any airfoil physical tag while retaining the same
    # common line elements themselves.
    filled_overrides: dict[int, list[int]] = {}
    for e in filled_elements:
        ent = elementary_entity(e)
        if e.etype == 1 and ent in (1, 2):
            filled_overrides[e.original_index] = []

    # The surface-99 elements must retain the fluid physical tag if present.
    # No override is needed: they already inherit their source physical tag.

    write_msh2(unfilled, common_names, common_nodes, common, unfilled_overrides)
    write_msh2(filled, physical_names, filled_nodes, filled_elements, filled_overrides)

    print(f"Source combined mesh : {src}")
    print(f"Common elements      : {len(common)}")
    print(f"Filled-only elements : {len(filled_extra)}")
    print(f"Common nodes         : {len(common_nodes)}")
    print(f"Filled-only nodes    : {len(filled_only_nodes)}")
    print(f"Unfilled output      : {unfilled}")
    print(f"Filled output        : {filled}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("unfilled", type=Path)
    parser.add_argument("filled", type=Path)
    args = parser.parse_args()
    build_outputs(args.source, args.unfilled, args.filled)


if __name__ == "__main__":
    main()
