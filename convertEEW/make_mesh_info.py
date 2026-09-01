#!/usr/bin/env python3

"""
make_mesh_info.py
==================

Generate mesh.info / mesh.data / mesh1.info / mesh1.data from a 2-D Gmsh mesh,
for a VAYU run that supplies its own pre-built 3-D mesh via mien/mxyz/mrng
(made by stackEEW.py)(i.e. the flow-only path, where driver.F's
`goto 1234` skips the internal 2D->3D stacking and jumps straight to hypoflow).

WHY THESE FILES ARE NEEDED AT ALL (read from parseinput.F)
------------------------------------------------------------

parseinput.F reads mesh.info/mesh1.info purely to compute nn/ne for sizing
MPI-partition arrays:

    nn = numnp2d*(wing_sli-1) + numnp2dl*(nslices-wing_sli+1)
    ne = numel2d*(wing_sli-1) + numel2dl*(nslices-wing_sli)

mesh.info also carries boundary-group node lists (nbndry2d groups), which
parseinput.F uses ONLY to build `bump_bn`/`node_bn` (airfoil shape-optimizer
control points) from GROUP INDEX 5 specifically (`nbump = npts2d(5)`).

For a plain flow run (no shape optimization - the driver.F you sent already
bypasses via `goto 1234`), this data is inert.

FILE FORMATS (read exactly as-is from parseinput.F)
------------------------------------------------------

mesh.info / mesh1.info are free-format (Fortran `read(unit,*)`):

    numnp2d numel2d
    nbndry2d
    npts2d(1)
    nbnd2d(1,1)
    ...
    nbnd2d(1,npts2d(1))
    npts2d(2)
    ...

mesh1.info is just:

    numnp2dl numel2dl

mesh.data / mesh1.data are FIXED-FORMAT text:

    read(11,'(16i5)') ((ien2d(i,j),i=1,3),j=1,numel2d)

    read(11,'(5e16.9)') dummy, ((x2d(i,j),i=1,2),j=1,numnp2d)
"""


import sys
import numpy as np


def read_gmsh(path):
    with open(path, "r") as f:
        lines = f.readlines()

    nodes = {}
    node_order = []
    elements = []  # list of (etype, phys_tag, [gmsh node ids])
    physical_names = {}  # (dim, tag) -> name

    i, n = 0, len(lines)

    while i < n:
        line = lines[i].strip()

        if line == "$PhysicalNames":
            i += 1
            cnt = int(lines[i].strip())
            i += 1

            for _ in range(cnt):
                parts = lines[i].split()
                i += 1

                dim = int(parts[0])
                tag = int(parts[1])
                name = parts[2].strip().strip('"')

                physical_names[(dim, tag)] = name

            i += 1  # $EndPhysicalNames

        elif line == "$Nodes":
            i += 1
            cnt = int(lines[i].strip())
            i += 1

            for _ in range(cnt):
                parts = lines[i].split()
                i += 1

                nid = int(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])

                nodes[nid] = (x, y, z)
                node_order.append(nid)

            i += 1  # $EndNodes

        elif line == "$Elements":
            i += 1
            cnt = int(lines[i].strip())
            i += 1

            for _ in range(cnt):
                parts = lines[i].split()
                i += 1

                etype = int(parts[1])
                ntags = int(parts[2])

                tags = [
                    int(t)
                    for t in parts[3:3 + ntags]
                ]

                nodelist = [
                    int(t)
                    for t in parts[3 + ntags:]
                ]

                # First Gmsh element tag is the physical-group tag.
                phys_tag = tags[0] if ntags > 0 else 0

                elements.append(
                    (etype, phys_tag, nodelist)
                )

            i += 1  # $EndElements

        else:
            i += 1

    return nodes, node_order, elements, physical_names


def write_mesh_info_data(
    msh_path,
    out_dir,
    boundary_group_order=(
        "outlet",
        "ground",
        "inlet",
        "topwall",
        "airfoil",
    ),
):
    import os

    os.makedirs(out_dir, exist_ok=True)

    nodes, node_order, elements, physical_names = read_gmsh(msh_path)

    # ------------------------------------------------------------------
    # Node numbering
    # ------------------------------------------------------------------

    nn = len(node_order)

    # Convert Gmsh node IDs -> VAYU local 1-based numbering.
    node_map = {
        old: k + 1
        for k, old in enumerate(node_order)
    }

    # ------------------------------------------------------------------
    # Physical groups
    #
    # IMPORTANT:
    #
    # A physical tag is NOT globally unique across dimensions.
    #
    # Therefore:
    #
    #     (1, 5) != (2, 5)
    #
    # A physical curve and a physical surface may legally have the same
    # numerical tag.
    #
    # We therefore ALWAYS keep (dimension, tag) together.
    # ------------------------------------------------------------------

    fluid_groups = {
        (dim, tag)
        for (dim, tag), name in physical_names.items()
        if name == "fluid"
    }

    boundary_groups = {
        (dim, tag): name
        for (dim, tag), name in physical_names.items()
        if name != "fluid" and dim > 0
    }

    name_to_group = {
        name: (dim, tag)
        for (dim, tag), name in boundary_groups.items()
    }

    # ------------------------------------------------------------------
    # Volume elements
    #
    # Gmsh element type 2 = 3-node first-order triangle.
    #
    # Fluid is a 2-D physical group, so only (2, phys_tag) is accepted.
    # ------------------------------------------------------------------

    tri_elements = [
        nodelist
        for (etype, phys_tag, nodelist) in elements
        if etype == 2 and (2, phys_tag) in fluid_groups
    ]

    numel2d = len(tri_elements)
    numnp2d = nn

    ien2d = np.array(
        [
            [node_map[g] for g in nodelist]
            for nodelist in tri_elements
        ],
        dtype=np.int64,
    )

    # ------------------------------------------------------------------
    # Coordinates
    # ------------------------------------------------------------------

    coords = np.zeros(
        (nn, 2),
        dtype=np.float64,
    )

    for k, old in enumerate(node_order):
        x, y, z = nodes[old]
        coords[k] = (x, y)

    # ------------------------------------------------------------------
    # Boundary node groups
    #
    # Boundary groups are expected to be 1-D physical curves.
    #
    # Gmsh element type 1 = 2-node first-order line.
    #
    # Again, dimension is explicitly checked:
    #
    #     dim == 1
    #     etype == 1
    #
    # ------------------------------------------------------------------

    group_nodes = {}

    for name in boundary_group_order:

        group = name_to_group.get(name)

        if group is None:
            group_nodes[name] = []
            continue

        dim, tag = group

        node_ids = set()

        for etype, phys_tag, nodelist in elements:

            if (
                dim == 1
                and etype == 1
                and phys_tag == tag
            ):
                for g in nodelist:
                    node_ids.add(node_map[g])

        group_nodes[name] = sorted(node_ids)

    # ------------------------------------------------------------------
    # surf_elem.info
    # ------------------------------------------------------------------

    airfoil_nodes = group_nodes.get(
        "airfoil",
        []
    )

    with open(
        f"{out_dir}/surf_elem.info",
        "w",
    ) as f:

        f.write(
            f"{len(airfoil_nodes)}\n"
        )

        for nid in airfoil_nodes:

            f.write(
                f"{nid}\n"
            )

    # ------------------------------------------------------------------
    # mesh.info
    # ------------------------------------------------------------------

    with open(
        f"{out_dir}/mesh.info",
        "w",
    ) as f:

        f.write(
            f"{numnp2d} {numel2d}\n"
        )

        f.write(
            f"{len(boundary_group_order)}\n"
        )

        for name in boundary_group_order:

            ids = group_nodes[name]

            f.write(
                f"{len(ids)}\n"
            )

            for nid in ids:

                f.write(
                    f"{nid}\n"
                )

    # ------------------------------------------------------------------
    # mesh1.info
    #
    # Sizing only - identical counts, no boundary groups needed.
    # ------------------------------------------------------------------

    with open(
        f"{out_dir}/mesh1.info",
        "w",
    ) as f:

        f.write(
            f"{numnp2d} {numel2d}\n"
        )

    # ------------------------------------------------------------------
    # mesh.data / mesh1.data
    #
    # Fixed Fortran formats:
    #
    #   (16i5)       element connectivity
    #   (5e16.9)     coordinates
    # ------------------------------------------------------------------

    def write_mesh_data(
        path,
        ien,
        xy,
    ):

        with open(path, "w") as f:

            # ----------------------------------------------------------
            # Element connectivity
            # ----------------------------------------------------------

            # We need:
            #
            #   ((ien2d(i,j), i=1,3), j=1,numel2d)
            #
            # i must vary fastest.
            #
            # ien is currently shaped:
            #
            #   (numel2d, 3)
            #
            # C-order flattening therefore gives exactly:
            #
            #   elem1_node1 elem1_node2 elem1_node3
            #   elem2_node1 elem2_node2 elem2_node3
            #   ...
            #
            flat_ien = ien.flatten()

            for i in range(
                0,
                len(flat_ien),
                16,
            ):

                chunk = flat_ien[
                    i:i + 16
                ]

                f.write(
                    " ".join(
                        str(int(v))
                        for v in chunk
                    )
                    + "\n"
                )

            # ----------------------------------------------------------
            # Coordinates
            # ----------------------------------------------------------

            flat_xy = np.empty(
                1 + 2 * xy.shape[0],
                dtype=np.float64,
            )

            # Dummy value read by Fortran but not used downstream.
            flat_xy[0] = 0.0

            flat_xy[1::2] = xy[:, 0]
            flat_xy[2::2] = xy[:, 1]

            for i in range(
                0,
                len(flat_xy),
                5,
            ):

                chunk = flat_xy[
                    i:i + 5
                ]

                f.write(
                    "".join(
                        f"{v:16.9E}"
                        for v in chunk
                    )
                    + "\n"
                )

    write_mesh_data(
        f"{out_dir}/mesh.data",
        ien2d,
        coords,
    )

    write_mesh_data(
        f"{out_dir}/mesh1.data",
        ien2d,
        coords,
    )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    print()
    print(f"--- {msh_path} ---")
    print(
        f"  numnp2d={numnp2d}  "
        f"numel2d={numel2d}"
    )

    print()
    print("  Physical groups:")

    for (dim, tag), name in sorted(
        physical_names.items()
    ):
        print(
            f"    dim={dim} tag={tag}: "
            f"'{name}'"
        )

    print()
    print("  Boundary groups:")

    for name in boundary_group_order:

        group = name_to_group.get(name)

        if group is None:

            print(
                f"    '{name}': 0 nodes "
                f"[NOT FOUND in physical groups!]"
            )

        else:

            dim, tag = group

            print(
                f"    '{name}': "
                f"{len(group_nodes[name])} nodes "
                f"(dim={dim}, tag={tag})"
            )

    print()
    print(
        f"  wrote mesh.info, mesh.data, "
        f"mesh1.info, mesh1.data, "
        f"surf_elem.info -> {out_dir}"
    )
    print()


if __name__ == "__main__":

    if len(sys.argv) != 3:

        print(
            "usage: "
            "python make_mesh_info.py "
            "input_2D.msh output_dir/"
        )

        sys.exit(1)

    write_mesh_info_data(
        sys.argv[1],
        sys.argv[2],
    )