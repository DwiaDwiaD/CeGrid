#!/usr/bin/env python3
"""
make_mesh_info.py
==================
Generate mesh.info / mesh.data / mesh1.info / mesh1.data from a 2-D Gmsh mesh,
for a VAYU run that supplies its own pre-built 3-D mesh via mien/mxyz/mrng
(i.e. the flow-only path, where driver.F's `goto 1234` skips the internal
2D->3D stacking and jumps straight to hypoflow).

WHY THESE FILES ARE NEEDED AT ALL (read from parseinput.F)
------------------------------------------------------------
parseinput.F reads mesh.info/mesh1.info purely to compute nn/ne for sizing
MPI-partition arrays:

    nn = numnp2d*(wing_sli-1) + numnp2dl*(nslices-wing_sli+1)
    ne = numel2d*(wing_sli-1) + numel2dl*(nslices-wing_sli)

These nn/ne MUST equal the true node/element counts in your actual mxyz/mien
(built by gmsh_to_vayu.py) or the partitioner will misallocate. This script
solves that by setting:

    wing_sli = 2, nslices = 2   (put these in fsi.in)
    numnp2d = numnp2dl = (your 3-D nn) / 2      <- the 2-D mesh's node count
    numel2d = numel2dl = your 3-D ne            <- the 2-D mesh's element count

which collapses the formula to nn = 2*numnp2d, ne = numel2d*1 - i.e. exactly
a single-layer spanwise extrusion of the 2-D mesh, matching how your 3-D
meshes were actually built (nn_3D = 2 * nn_2D was verified directly against
your uploaded files).

mesh.info also carries boundary-group node lists (nbndry2d groups), which
parseinput.F uses ONLY to build `bump_bn`/`node_bn` (airfoil shape-optimizer
control points) from GROUP INDEX 5 specifically (`nbump = npts2d(5)`). For a
plain flow run (no shape optimization - which the driver.F you sent already
bypasses via `goto 1234`), this data is inert. This script writes
nbndry2d = 4 (inlet, outlet, walls, airfoil, in that order) with REAL node
lists pulled from your (correctly tagged) 2-D mesh; group index 5 is left
absent, which relies on Fortran COMMON-block BSS zero-initialization making
npts2d(5) = 0 (so the `do i=1,nbump` loop just does nothing - safe, but you
should confirm this is acceptable if you ever turn shape optimization on).

FILE FORMATS (read exactly as-is from parseinput.F)
------------------------------------------------------
mesh.info / mesh1.info are free-format (Fortran `read(unit,*)`), so any
whitespace/newline separation works:
    numnp2d numel2d
    nbndry2d
    npts2d(1)
    nbnd2d(1,1)
    ...
    nbnd2d(1,npts2d(1))
    npts2d(2)
    ...
mesh1.info is just: `numnp2dl numel2dl` (no boundary groups needed - only
`read(10,*)numnp2dl,numel2dl` is ever executed against it).

mesh.data / mesh1.data are FIXED-FORMAT text (Fortran formatted read):
    read(11,'(16i5)') ((ien2d(i,j),i=1,3),j=1,numel2d)   ! elements, 16 per line, width 5
    read(11,'(5e16.9)') dummy, ((x2d(i,j),i=1,2),j=1,numnp2d)  ! 1 dummy + x,y pairs, 5/line, width 16
"""

import sys
import numpy as np
from gmsh_to_vayu import read_gmsh


def write_mesh_info_data(msh_path, out_dir, boundary_group_order=("inlet", "outlet", "topwall", "ground", "airfoil")):
    import os
    os.makedirs(out_dir, exist_ok=True)

    nodes, node_order, elements, physical_names = read_gmsh(msh_path)
    nn = len(node_order)
    node_map = {old: k + 1 for k, old in enumerate(node_order)}  # 1-indexed

    fluid_tags = {tag for (dim, tag), name in physical_names.items() if name == "fluid"}
    boundary_names = {tag: name for (dim, tag), name in physical_names.items()
                       if name != "fluid" and dim > 0}
    name_to_tag = {name: tag for tag, name in boundary_names.items()}

    # volume (triangle) elements -> ien2d, in 1..numnp2d local numbering
    tri_elements = [nodelist for (etype, phys_tag, nodelist) in elements if phys_tag in fluid_tags]
    numel2d = len(tri_elements)
    numnp2d = nn
    ien2d = np.array([[node_map[g] for g in nodelist] for nodelist in tri_elements], dtype=np.int64)

    coords = np.zeros((nn, 2), dtype=np.float64)
    for k, old in enumerate(node_order):
        x, y, z = nodes[old]
        coords[k] = (x, y)

    # boundary node groups: unique nodes touched by each named 1-D physical group
    group_nodes = {}
    for name in boundary_group_order:
        tag = name_to_tag.get(name)
        if tag is None:
            group_nodes[name] = []
            continue
        node_ids = set()
        for etype, phys_tag, nodelist in elements:
            if phys_tag == tag:
                for g in nodelist:
                    node_ids.add(node_map[g])
        group_nodes[name] = sorted(node_ids)

    # --- surf_elem.info ---
    airfoil_nodes = group_nodes.get(
        "airfoil",
        []
    )

    with open(
        f"{out_dir}/surf_elem.info",
        "w"
    ) as f:

        f.write(
            f"{len(airfoil_nodes)}\n"
        )

        for nid in airfoil_nodes:

            f.write(
                f"{nid}\n"
            )

    # --- mesh.info ---
    with open(f"{out_dir}/mesh.info", "w") as f:
        f.write(f"{numnp2d} {numel2d}\n")
        f.write(f"{len(boundary_group_order)}\n")
        for name in boundary_group_order:
            ids = group_nodes[name]
            f.write(f"{len(ids)}\n")
            for nid in ids:
                f.write(f"{nid}\n")

    # --- mesh1.info (sizing only - identical counts, no boundary groups needed) ---
    with open(f"{out_dir}/mesh1.info", "w") as f:
        f.write(f"{numnp2d} {numel2d}\n")

    # --- mesh.data (fixed format: 16i5 for elements, then dummy + 5e16.9 for x,y) ---
    def write_mesh_data(path, ien, xy):
        with open(path, "w") as f:
            flat_ien = ien.flatten()  # i=1..3 fastest, then j=1..numel2d
            for i in range(0, len(flat_ien), 16):
                chunk = flat_ien[i:i + 16]
                f.write("".join(f"{v:5d}" for v in chunk) + "\n")
            flat_xy = np.empty(1 + 2 * xy.shape[0], dtype=np.float64)
            flat_xy[0] = 0.0  # 'dummy' - read but never used downstream
            flat_xy[1::2] = xy[:, 0]
            flat_xy[2::2] = xy[:, 1]
            for i in range(0, len(flat_xy), 5):
                chunk = flat_xy[i:i + 5]
                f.write("".join(f"{v:16.9E}" for v in chunk) + "\n")

    write_mesh_data(f"{out_dir}/mesh.data", ien2d, coords)
    write_mesh_data(f"{out_dir}/mesh1.data", ien2d, coords)  # identical, per the wing_sli=nslices=2 trick

    print(f"--- {msh_path} ---")
    print(f"  numnp2d={numnp2d}  numel2d={numel2d}")
    for name in boundary_group_order:
        print(f"  boundary group '{name}': {len(group_nodes[name])} nodes"
              + ("" if name_to_tag.get(name) is not None else "  [NOT FOUND in physical groups!]"))
    print(f"  wrote mesh.info, mesh.data, mesh1.info, mesh1.data -> {out_dir}")
    print(f"  ADD TO fsi.in: wing_sli 2   and   nslices 2")
    print(f"  (this makes nn = 2*{numnp2d} = {2*numnp2d}, ne = {numel2d} match your actual 3-D mesh)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python make_mesh_info.py input_2D.msh output_dir/")
        sys.exit(1)
    write_mesh_info_data(sys.argv[1], sys.argv[2])
