#!/usr/bin/env python3
"""
make_mesh_info_pair.py
========================
Generates mesh.info/mesh.data from an UNFILLED 2-D Gmsh mesh (the actual
wing/root cross-section, airfoil hole present) and mesh1.info/mesh1.data
from a FILLED 2-D Gmsh mesh (the beyond-tip cap, airfoil hole filled in) -
in one step, no manual copying/renaming needed.

See make_mesh_info.py for the full explanation of why these files exist and
what parseinput.F does with them; this script is the same logic, just run
once per mesh and written straight to the target file names.

Usage:
    python make_mesh_info_pair.py Unfilled2D.msh Filled2D.msh output_dir/
"""

import sys
import os
from make_mesh_info import write_mesh_info_data


def make_pair(unfilled_msh, filled_msh, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    # --- mesh.info / mesh.data <- UNFILLED mesh ---
    tmp_u = os.path.join(out_dir, "_tmp_unfilled")
    write_mesh_info_data(unfilled_msh, tmp_u, boundary_group_order=("outlet", "ground", "inlet", "topwall", "airfoil"))
    os.replace(os.path.join(tmp_u, "mesh.info"), os.path.join(out_dir, "mesh.info"))
    os.replace(os.path.join(tmp_u, "mesh.data"), os.path.join(out_dir, "mesh.data"))
    os.replace(os.path.join(tmp_u, "surf_elem.info"), os.path.join(out_dir, "surf_elem.info"))
    os.rmdir(tmp_u) if not os.listdir(tmp_u) else None

    # --- mesh1.info / mesh1.data <- FILLED mesh ---
    tmp_f = os.path.join(out_dir, "_tmp_filled")
    write_mesh_info_data(filled_msh, tmp_f, boundary_group_order=("outlet", "ground", "inlet", "topwall"))
    os.replace(os.path.join(tmp_f, "mesh1.info"), os.path.join(out_dir, "mesh1.info"))
    os.replace(os.path.join(tmp_f, "mesh1.data"), os.path.join(out_dir, "mesh1.data"))
    os.rmdir(tmp_f) if not os.listdir(tmp_f) else None

    print()
    print(f"Wrote to {out_dir}/:")
    print(f"  mesh.info,  mesh.data   <- {unfilled_msh}  (root/wing cross-section)")
    print(f"  mesh1.info, mesh1.data  <- {filled_msh}  (beyond-tip cap)")
    print()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: python make_mesh_info_pair.py Unfilled2D.msh Filled2D.msh output_dir/")
        sys.exit(1)
    make_pair(sys.argv[1], sys.argv[2], sys.argv[3])
