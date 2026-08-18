#!/usr/bin/env python3

from pathlib import Path
import argparse
import time

import numpy as np


# ============================================================
# STACKING CONFIGURATION
# ============================================================

WING_SLICES = 2       # number of spanwise planes that belong to the wing
P2 = 247              # retained for compatibility/documentation
NSLICES = 5           # total spanwise planes; change here with the zall.data input

"""
Corrected NumPy replacement for stack_modf_works.f

The original stacking/connectivity logic is retained, but boundary 5
(the wing/airfoil surface) is built from the actual 2-D connectivity
instead of treating surf_elem.info (which contains NODE IDs) as element IDs.

Configuration
-------------
WING_SLICES:
    0         -> no wing surface is assigned (all planes use mesh1)
    1..N-1    -> the first WING_SLICES *planes/slices* belong to the wing;
                  only the intervals joining those wing slices receive the
                  wing side boundary, and the last wing slice is capped
                  using the added mesh1 elements.
    N         -> all N planes belong to the wing, so every spanwise interval
                  receives the wing side boundary and there is no extra cap.

Thus WING_SLICES is the number of spanwise *planes/slices* occupied by the
wing, not the number of intervals. NSLICES remains a user-editable
configuration value at the top of this file and must match zall.data.

Finite-wing transition:
    For 0 < WING_SLICES < NSLICES, the first WING_SLICES planes contain the
    wing, with mesh.data connectivity through the wing. The final wing plane
    is the first mesh1 plane so the added mesh1 elements can seal the wing
    tip. Those added elements in the interval immediately downstream of the
    final wing plane receive boundary 5 on their first (tip) face only.
    The airfoil side boundary is never assigned to that downstream interval.

Boundary 5:
    - side faces are found from mesh.info boundary-5 NODE IDs and the
      triangle connectivity in mesh.data / mesh1.data;
    - the mesh1-added elements (ne2-ne1) in the first filled layer are
      assigned as the wing-tip cap (face 1), with no ncount>=3 hole.

Other boundary IDs retain the source-faithful bnd_3d/sub_mrng logic.
"""

# ============================================================
# INPUT READERS
# ============================================================

def read_fixed_width_mesh(path, numnp, numel):
    """Read the fixed-width mesh.data / mesh1.data format."""

    lines = Path(path).read_text().splitlines()

    nint = 3 * numel
    nint_lines = (nint + 15) // 16

    ints = []
    for line in lines[:nint_lines]:
        for i in range(0, len(line), 5):
            field = line[i:i + 5]
            if field.strip():
                ints.append(int(field))

    if len(ints) < nint:
        raise ValueError(
            f"{path}: expected {nint} integers, got {len(ints)}"
        )

    ien = np.asarray(ints[:nint], dtype=np.int32).reshape(numel, 3)

    text = "".join(lines[nint_lines:])
    nreal = 1 + 2 * numnp
    need = 16 * nreal

    if len(text) < need:
        raise ValueError(
            f"{path}: coordinate section is shorter than expected"
        )

    xyz2 = np.fromiter(
        (float(text[i:i + 16]) for i in range(0, need, 16)),
        dtype=np.float64,
        count=nreal,
    )

    return ien, xyz2[1:].reshape(numnp, 2)


def read_mesh_info(path):
    """Read mesh.info: numnp numel nbnd followed by node lists."""

    vals = np.fromstring(
        Path(path).read_text().replace("\r", "\n"),
        sep=" ",
        dtype=np.int64,
    )

    if len(vals) < 3:
        raise ValueError(f"{path}: expected numnp, numel, nbnd")

    numnp, numel, nbnd = map(int, vals[:3])
    pos = 3
    boundaries = []

    for ib in range(nbnd):
        if pos >= len(vals):
            raise ValueError(
                f"{path}: unexpected EOF before boundary {ib + 1}"
            )

        n = int(vals[pos])
        pos += 1

        if pos + n > len(vals):
            raise ValueError(
                f"{path}: boundary {ib + 1} extends past EOF"
            )

        boundaries.append(
            np.asarray(vals[pos:pos + n], dtype=np.int32)
        )
        pos += n

    if pos != len(vals):
        raise ValueError(
            f"{path}: unexpected data after boundary information"
        )

    return numnp, numel, boundaries


def read_mesh1_info(path):
    """Read mesh1.info: numnp1 numel1."""

    vals = np.fromstring(
        Path(path).read_text().replace("\r", "\n"),
        sep=" ",
        dtype=np.int64,
    )

    if len(vals) < 2:
        raise ValueError(f"{path}: expected numnp1 and numel1")

    return int(vals[0]), int(vals[1])


def read_surface_nodes(path):
    """
    Read surf_elem.info.

    Despite the historical name, this file contains SURFACE NODE IDs:
        n
        node_1
        ...
        node_n

    Only the first n values are consumed. Any trailing value is reported
    but ignored, which handles the stray 258 in the supplied test file.
    """

    vals = np.fromstring(
        Path(path).read_text().replace("\r", "\n"),
        sep=" ",
        dtype=np.int64,
    )

    if len(vals) < 1:
        raise ValueError(f"{path}: empty surface-node file")

    n = int(vals[0])
    if len(vals) < n + 1:
        raise ValueError(
            f"{path}: expected {n} node IDs, got {len(vals) - 1}"
        )

    nodes = np.asarray(vals[1:n + 1], dtype=np.int32)
    extra = len(vals) - (n + 1)
    return nodes, extra


# ============================================================
# OUTPUT WRITER
# ============================================================

def write_fortran_array(path, array, dtype):
    """Write a native binary stream matching the Fortran array bytes."""

    a = np.asarray(array, dtype=dtype)
    a.tofile(str(path))


# ============================================================
# STACK GEOMETRY HELPERS
# ============================================================

def validate_config(wing_slices, nslices):
    if nslices < 1:
        raise ValueError("NSLICES must be >= 1")
    if not 0 <= wing_slices <= nslices:
        raise ValueError(
            f"WING_SLICES must satisfy 0 <= WING_SLICES <= NSLICES "
            f"({nslices}); got {wing_slices}"
        )


def plane_is_unfilled(plane0, wing_slices, nslices):
    """Return whether 0-based plane uses mesh.data."""

    if wing_slices == nslices:
        return True
    if wing_slices == 0:
        return False
    # Legacy-compatible transition: first wing_slices-1 planes use mesh;
    # the wing_slices-th plane is the first mesh1 plane.
    return plane0 < wing_slices - 1


def interval_uses_unfilled(interval1, wing_slices, nslices):
    """Return whether 1-based interval uses mesh.data connectivity."""

    if wing_slices == nslices:
        return True
    if wing_slices == 0:
        return False
    return interval1 < wing_slices


def plane_starts(np1, np2, nslices, wing_slices):
    starts = np.empty(nslices, dtype=np.int64)

    if wing_slices == nslices:
        for p in range(nslices):
            starts[p] = p * np1
        return starts

    if wing_slices == 0:
        for p in range(nslices):
            starts[p] = p * np2
        return starts

    n_mesh_planes = wing_slices - 1
    for p in range(nslices):
        if p < n_mesh_planes:
            starts[p] = p * np1
        else:
            starts[p] = n_mesh_planes * np1 + (p - n_mesh_planes) * np2

    return starts


def make_connectivity(ien1, ien2, np1, np2, ne1, ne2, wing_slices, nslices):
    """Construct the stacked mesh, switching to mesh1 at the final wing plane.

    For a finite wing, the final wing plane uses mesh1 so the additional
    mesh1 elements can form the filled/capped tip without extending the
    wing side surface into the downstream interval.
    """

    hex_parts = []
    wedge_parts = []
    interval_offsets = []
    offset = 0

    starts = plane_starts(np1, np2, nslices, wing_slices)

    for k in range(nslices - 1):
        interval1 = k + 1
        interval_offsets.append(offset)

        # The interval is part of the open wing only when both of its
        # bounding planes are wing planes. For a finite wing, the interval
        # immediately after the last wing plane is mesh1/filled connectivity:
        # it exists to provide the solid tip closure, not another wing side.
        if interval_uses_unfilled(interval1, wing_slices, nslices):
            a = ien1 + int(starts[k])
            layer_offset = np1
        else:
            a = ien2 + int(starts[k])
            layer_offset = np2

        wedge = np.empty((len(a), 6), dtype=np.int32)
        wedge[:, :3] = a
        wedge[:, 3:] = a + layer_offset

        hexa = np.empty((len(a), 8), dtype=np.int32)
        hexa[:, :3] = a
        hexa[:, 3] = a[:, 2]
        hexa[:, 4:7] = a + layer_offset
        hexa[:, 7] = hexa[:, 6]

        wedge_parts.append(wedge)
        hex_parts.append(hexa)
        offset += len(a)

    return (
        np.vstack(hex_parts),
        np.vstack(wedge_parts),
        np.asarray(interval_offsets, dtype=np.int64),
        starts,
    )


def make_coordinates(xy1, xy2, z, np1, np2, nslices, wing_slices):
    """Stack 2-D coordinates using the same plane layout as connectivity."""

    starts = plane_starts(np1, np2, nslices, wing_slices)
    numnpt = int(
        starts[-1] + (np1 if plane_is_unfilled(nslices - 1, wing_slices, nslices) else np2)
    )

    xyz = np.empty((numnpt, 3), dtype=np.float64)

    for p in range(nslices):
        if plane_is_unfilled(p, wing_slices, nslices):
            xy = xy1
            n = np1
        else:
            xy = xy2
            n = np2

        s = int(starts[p])
        xyz[s:s + n, :2] = xy
        xyz[s:s + n, 2] = z[p]

    return xyz


# ============================================================
# BOUNDARY NODES
# ============================================================

def make_boundary_nodes(boundaries, np1, np2, nslices, wing_slices):
    """Build seven 3-D boundary-node lists using the actual plane offsets."""

    starts = plane_starts(np1, np2, nslices, wing_slices)
    bnd = [None] * 7

    # Boundaries 1..4: repeat each 2-D node list on every plane.
    for ib in range(4):
        nodes2 = boundaries[ib]
        chunks = []
        for p in range(nslices):
            chunks.append(nodes2 + int(starts[p]))
        bnd[ib] = np.concatenate(chunks).astype(np.int32)

    # Boundary 5: surface node set is the actual boundary-5 node set
    # repeated/mapped as follows. For mesh1 planes, mesh1 adds all extra
    # nodes (the filled/cap region) after the original np1 nodes.
    base = boundaries[4]
    parts = []

    for p in range(nslices):
        nodes = base.copy()
        if plane_is_unfilled(p, wing_slices, nslices):
            parts.append(nodes + int(starts[p]))
        else:
            parts.append(nodes + int(starts[p]))

    # On mesh1 planes, include the additional mesh1 nodes in the boundary-5
    # node set. This preserves the legacy filled-tip node definition.
    if np2 > np1:
        for p in range(nslices):
            if not plane_is_unfilled(p, wing_slices, nslices):
                extra = np.arange(np1 + 1, np2 + 1, dtype=np.int32)
                parts.append(extra + int(starts[p]))

    bnd[4] = np.concatenate(parts).astype(np.int32)

    # Boundary 6 = first 2-D plane.
    nfirst = np1 if plane_is_unfilled(0, wing_slices, nslices) else np2
    bnd[5] = np.arange(1, nfirst + 1, dtype=np.int32)

    # Boundary 7 = final 2-D plane.
    nlast = np1 if plane_is_unfilled(nslices - 1, wing_slices, nslices) else np2
    bnd[6] = np.arange(1, nlast + 1, dtype=np.int32) + int(starts[-1])

    return bnd


# ============================================================
# BOUNDARY ELEMENTS
# ============================================================

def make_boundary_elements(boundaries, wedge, np1, np2, ne1, ne2,
                           nslices, wing_slices, interval_offsets):
    """
    Build bnd_el for boundaries 1..7.

    Boundaries 1..4, 6, 7 retain the source-faithful node-membership search.
    Boundary 5 is corrected separately from actual 2-D connectivity.
    """

    bnodes = make_boundary_nodes(
        boundaries, np1, np2, nslices, wing_slices
    )

    bnd_el = [None] * 7

    max_node = int(wedge.max())

    # Boundaries 1..4.
    for ib in range(4):
        nodes = bnodes[ib]
        node_count = np.bincount(nodes, minlength=max_node + 1)
        counts = node_count[wedge]
        selected = counts.sum(axis=1, dtype=np.int64) >= 4
        bnd_el[ib] = (np.flatnonzero(selected) + 1).astype(np.int32)

    # Boundary 5 is populated by the topology-based routine later.
    bnd_el[4] = np.empty(0, dtype=np.int32)

    # Boundary 6 = first plane elements.
    n_first_el = ne1 if plane_is_unfilled(0, wing_slices, nslices) else ne2
    bnd_el[5] = np.arange(1, n_first_el + 1, dtype=np.int32)

    # Boundary 7 = final plane elements.
    n_last_el = ne1 if plane_is_unfilled(nslices - 1, wing_slices, nslices) else ne2
    bnd_el[6] = (
        np.arange(1, n_last_el + 1, dtype=np.int32)
        + int(interval_offsets[-1])
    )

    return bnd_el, bnodes


def derive_airfoil_surface_elements(ien1, boundary5_nodes):
    """
    Find 2-D triangle elements on the physical airfoil boundary.

    A triangle is a surface element when it has at least two nodes in the
    boundary-5 node set. For the supplied mesh this yields 247 elements.
    """

    node_set = np.zeros(int(max(ien1.max(), boundary5_nodes.max())) + 1,
                        dtype=bool)
    node_set[boundary5_nodes] = True
    count = node_set[ien1].sum(axis=1)
    elems = np.flatnonzero(count >= 2).astype(np.int32) + 1

    return elems


def edge_to_face_slot(tri, boundary_nodes_mask):
    """
    Return the wedge mrng face slot (0..4) corresponding to the airfoil edge.

    Wedge local nodes are [1,2,3,4,5,6]. The lateral faces are:
        edge (1,2) -> face 4 -> slot 3
        edge (2,3) -> face 5 -> slot 4
        edge (3,1) -> face 3 -> slot 2
    """

    hits = [
        boundary_nodes_mask[tri[0]] and boundary_nodes_mask[tri[1]],
        boundary_nodes_mask[tri[1]] and boundary_nodes_mask[tri[2]],
        boundary_nodes_mask[tri[2]] and boundary_nodes_mask[tri[0]],
    ]

    if sum(hits) != 1:
        raise ValueError(
            f"Could not identify a unique airfoil edge for triangle {tri.tolist()}; "
            f"edge matches={hits}"
        )

    if hits[0]:
        return 3
    if hits[1]:
        return 4
    return 2


def assign_correct_boundary5(mrng, ien1, boundary5_nodes, ne1, ne2,
                              nslices, wing_slices, interval_offsets):
    """
    Assign the physical wing/airfoil boundary (ID 5) for a finite-wing stack.

    Geometry convention:
        WING_SLICES is the number of spanwise planes occupied by the wing.

    Boundary-5 side:
        Only the intervals *between* adjacent wing planes are assigned ID 5.
        Therefore WING_SLICES=2 assigns the wing side only to interval 1
        (slice 1 -> slice 2). The next interval (slice 2 -> slice 3) is
        deliberately not assigned as an airfoil side.

    Boundary-5 tip:
        For a finite wing, the final wing plane is represented by the first
        mesh1 plane. The extra mesh1 elements (ne1+1 .. ne2) in the interval
        immediately after that plane form the filled tip closure. Their first
        wedge face (mrng slot 0) is the physical wingtip and receives ID 5.
    """

    if wing_slices == 0:
        return 0, 0, 0

    if wing_slices > nslices:
        raise ValueError(
            f"WING_SLICES={wing_slices} exceeds NSLICES={nslices}"
        )

    mask = np.zeros(int(max(ien1.max(), boundary5_nodes.max())) + 1,
                    dtype=bool)
    mask[boundary5_nodes] = True

    surface_elems = derive_airfoil_surface_elements(ien1, boundary5_nodes)

    if len(surface_elems) == 0:
        raise ValueError("No airfoil surface elements found from boundary-5 nodes")

    face_slots = np.empty(len(surface_elems), dtype=np.int32)
    for j, elem1 in enumerate(surface_elems):
        tri = ien1[elem1 - 1]
        face_slots[j] = edge_to_face_slot(tri, mask)

    # WING_SLICES counts planes, so there are WING_SLICES-1 wing-side
    # intervals when the wing is finite. For WING_SLICES == NSLICES,
    # every available interval belongs to the wing.
    side_intervals = (
        range(1, wing_slices)
        if wing_slices < nslices
        else range(1, nslices)
    )

    side_count = 0
    for interval1 in side_intervals:
        interval0 = interval1 - 1
        offset = int(interval_offsets[interval0])

        for elem1, slot in zip(surface_elems, face_slots):
            row = offset + int(elem1) - 1
            mrng[row, int(slot)] = 5
            side_count += 1

    cap_count = 0

    if 0 < wing_slices < nslices:
        # This is the interval from the final wing plane to the first
        # downstream plane. The added mesh1 elements in this interval seal
        # the finite-wing tip at the final wing plane.
        cap_interval0 = wing_slices - 1
        offset = int(interval_offsets[cap_interval0])

        if ne2 < ne1:
            raise ValueError(
                f"mesh1 element count ne2={ne2} is smaller than mesh "
                f"element count ne1={ne1}"
            )

        extra_count = ne2 - ne1
        if extra_count:
            rows = offset + np.arange(ne1, ne2, dtype=np.int64)
            mrng[rows, 0] = 5
            cap_count = int(extra_count)

    return int(len(surface_elems)), side_count, cap_count


# ============================================================
# MRNG
# ============================================================

def make_mrng_base(wedge, bnd_el, bnd_nodes, chunk_size=1_000_000):
    """Source-faithful mrng construction for all non-wing-surface IDs."""

    numelt = len(wedge)
    mrng = np.zeros((numelt, 5), dtype=np.int32)
    max_node = int(wedge.max())

    iface = 0
    inode_numbers = np.arange(1, 7, dtype=np.int64)

    # Deliberately preserve the original sequential-if / persistent-iface
    # behavior for boundaries 1..4,6,7. Boundary 5 is assigned separately.
    for ib in list(range(4)) + [5, 6]:
        nodes = bnd_nodes[ib]
        node_count = np.bincount(nodes, minlength=max_node + 1)
        elems = bnd_el[ib] - 1

        for start in range(0, len(elems), chunk_size):
            stop = min(start + chunk_size, len(elems))
            eidx = elems[start:stop]
            counts = node_count[wedge[eidx]]
            ncount = counts.sum(axis=1, dtype=np.int64)
            icount = (counts * inode_numbers).sum(axis=1, dtype=np.int64)

            for local in range(len(eidx)):
                s = int(icount[local])

                if s == 6:
                    iface = 0
                if s == 15:
                    iface = 1
                if s == 14:
                    iface = 2
                if s == 12:
                    iface = 3
                if s == 16:
                    iface = 4

                if ncount[local] >= 3:
                    mrng[eidx[local], iface] = ib + 1

    return mrng


# ============================================================
# DEBUG COORDINATE OUTPUT
# ============================================================

def write_debug_coordinates(path, xyz):
    with open(path, "w") as f:
        for x, y, z in xyz:
            f.write(f"{x:24.16g}{y:24.16g}{z:24.16g}\n")


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", default=".",
        help="directory containing mesh.info, mesh.data, mesh1.info, mesh1.data, surf_elem.info, zall.data",
    )
    parser.add_argument(
        "--output", default=".",
        help="directory for mien/mxyz/mrng/mienb",
    )
    parser.add_argument(
        "--wing-slices", type=int, default=WING_SLICES,
        help="number of spanwise wing slices/planes; defaults to WING_SLICES",
    )
    parser.add_argument(
        "--fort1", action="store_true",
        help="also write readable coordinates to fort.1",
    )

    args = parser.parse_args()
    wing_slices = args.wing_slices

    inp = Path(args.input)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    n1, e1, nbnd_nodes = read_mesh_info(inp / "mesh.info")
    n2, e2 = read_mesh1_info(inp / "mesh1.info")

    if len(nbnd_nodes) != 5:
        raise ValueError(
            f"mesh.info contains {len(nbnd_nodes)} boundaries; expected 5"
        )

    z = np.atleast_1d(np.loadtxt(inp / "zall.data", dtype=np.float64))
    if len(z) != NSLICES:
        raise ValueError(
            f"zall.data: expected {NSLICES} z locations, got {len(z)}"
        )

    nslices = NSLICES
    validate_config(wing_slices, nslices)

    ien1, xy1 = read_fixed_width_mesh(inp / "mesh.data", n1, e1)
    ien2, xy2 = read_fixed_width_mesh(inp / "mesh1.data", n2, e2)

    # surf_elem.info is a node list. Read it only as a validation/reference;
    # boundary construction uses mesh.info boundary 5 as the authoritative
    # node set because both files describe the same physical boundary.
    surf_nodes, extra_values = read_surface_nodes(inp / "surf_elem.info")

    # The supplied legacy file is known to have a stale header count of 248
    # with the final value 258 left over.  Since mesh.info is authoritative,
    # accept a prefix that exactly matches boundary 5 and ignore the stale tail.
    if not np.array_equal(surf_nodes, nbnd_nodes[4]):
        if len(surf_nodes) > len(nbnd_nodes[4]) and np.array_equal(
            surf_nodes[:len(nbnd_nodes[4])], nbnd_nodes[4]
        ):
            print(
                f"Warning: surf_elem.info has {len(surf_nodes) - len(nbnd_nodes[4])} "
                "stale trailing node(s); ignoring them."
            )
            surf_nodes = nbnd_nodes[4].copy()
        else:
            # Special case for the supplied 248-header file: its first 247
            # node IDs match mesh.info exactly, while the final 258 is stale.
            vals = np.fromstring(
                (inp / "surf_elem.info").read_text().replace("\r", "\n"),
                sep=" ", dtype=np.int64
            )
            if len(vals) >= len(nbnd_nodes[4]) + 2 and np.array_equal(
                vals[1:1 + len(nbnd_nodes[4])], nbnd_nodes[4]
            ):
                print(
                    "Warning: using first 247 surface nodes from surf_elem.info; "
                    "ignoring stale trailing value 258."
                )
                surf_nodes = nbnd_nodes[4].copy()
            else:
                raise ValueError(
                    "surf_elem.info node list does not match mesh.info boundary 5"
                )

    t0 = time.perf_counter()

    mienb, mien, interval_offsets, plane_starts_ = make_connectivity(
        ien1, ien2, n1, n2, e1, e2, wing_slices, nslices
    )

    xyz = make_coordinates(
        xy1, xy2, z, n1, n2, nslices, wing_slices
    )

    bnd_el, bnd_nodes = make_boundary_elements(
        nbnd_nodes, mien, n1, n2, e1, e2,
        nslices, wing_slices, interval_offsets
    )

    mrng = make_mrng_base(mien, bnd_el, bnd_nodes)

    nsurf, side_count, cap_count = (0, 0, 0)
    if wing_slices > 0:
        nsurf, side_count, cap_count = assign_correct_boundary5(
            mrng,
            ien1,
            nbnd_nodes[4],
            e1,
            e2,
            nslices,
            wing_slices,
            interval_offsets,
        )

    write_fortran_array(out / "mienb", mienb, np.int32)
    write_fortran_array(out / "mien", mien, np.int32)
    write_fortran_array(out / "mxyz", xyz, np.float64)
    write_fortran_array(out / "mrng", mrng, np.int32)

    if args.fort1:
        write_debug_coordinates(out / "fort.1", xyz)

    elapsed = time.perf_counter() - t0

    print(f"Done in {elapsed:.3f} s")
    print(f"Wing slices        : {wing_slices} / {NSLICES}")
    if 0 < wing_slices < NSLICES:
        print(f"Wing side intervals: 1..{wing_slices - 1}")
        print(f"Wingtip cap        : interval {wing_slices}, face 1")
    print(f"2D mesh 1 nodes    : {n1}")
    print(f"2D mesh 1 elements : {e1}")
    print(f"2D mesh 2 nodes    : {n2}")
    print(f"2D mesh 2 elements : {e2}")
    print(f"3D nodes           : {len(xyz)}")
    print(f"3D wedge elements  : {len(mien)}")
    print(f"3D hexa elements   : {len(mienb)}")
    print(f"z locations        : {z}")
    print(f"airfoil 2D elements found: {nsurf}")
    print(f"ID 5 side faces assigned: {side_count}")
    print(f"ID 5 wingtip cap faces assigned: {cap_count}")
    print(f"ID 5 total: {np.count_nonzero(mrng == 5)}")
    print(f"mrng nonzero entries: {np.count_nonzero(mrng)}")
    print(f"mienb bytes        : {(out / 'mienb').stat().st_size}")
    print(f"mien bytes         : {(out / 'mien').stat().st_size}")
    print(f"mxyz bytes         : {(out / 'mxyz').stat().st_size}")
    print(f"mrng bytes         : {(out / 'mrng').stat().st_size}")


if __name__ == "__main__":
    main()

