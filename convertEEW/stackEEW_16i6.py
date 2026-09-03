#!/usr/bin/env python3

from pathlib import Path
import argparse
import time

import numpy as np


# ============================================================
# STACKING CONFIGURATION
# ============================================================

WING_SLICES = 186       # number of spanwise planes that belong to the wing
P2 = 248                # retained for compatibility/documentation
NSLICES = 237           # total spanwise planes; change here with the zall.data input

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

def read_fixed_width_coordinates(coord_text, numnp):
    """
    Read the coordinate section.

    Preferred format:
        1 + 2*numnp values in fixed-width 16-character fields.

    Fallback:
        whitespace-separated floating-point values.
    """

    nreal = 1 + 2 * numnp
    need = 16 * nreal

    if len(coord_text) >= need:
        try:
            xyz2 = np.fromiter(
                (
                    float(coord_text[i:i + 16])
                    for i in range(0, need, 16)
                ),
                dtype=np.float64,
                count=nreal,
            )

            if len(xyz2) == nreal:
                return xyz2[1:].reshape(numnp, 2)

        except ValueError:
            pass

    coord_vals = np.fromstring(
        coord_text,
        sep=" ",
        dtype=np.float64,
    )

    if len(coord_vals) != nreal:
        raise ValueError(
            f"coordinate section: expected {nreal} values, "
            f"got {len(coord_vals)}"
        )

    return coord_vals[1:].reshape(numnp, 2)


def read_fixed_width_mesh(path, numnp, numel):
    """
    Read mesh.data / mesh1.data.

    The connectivity is written in the same fixed-width format expected by
    the legacy Fortran reader::

        (16i6)  connectivity
        (5e16.9) coordinates

    Each connectivity record contains up to 16 six-character integer fields.
    The final record may contain fewer than 16 integers when 3*numel is not
    divisible by 16.
    """

    path = Path(path)
    raw = path.read_text()
    lines = raw.splitlines()

    nint = 3 * numel
    nint_lines = (nint + 15) // 16

    if len(lines) < nint_lines:
        raise ValueError(
            f"{path}: expected at least {nint_lines} connectivity records, "
            f"got {len(lines)}"
        )

    ints = []

    # ------------------------------------------------------------
    # Fixed-width 16i6 connectivity.
    # ------------------------------------------------------------
    try:
        for line_no, line in enumerate(lines[:nint_lines], start=1):
            # A 16i6 record is at most 96 characters.
            if len(line) > 96:
                raise ValueError(
                    f"{path}: connectivity record {line_no} exceeds 96 "
                    f"characters; expected 16i6"
                )

            for i in range(0, len(line), 6):
                field = line[i:i + 6]
                if field.strip():
                    ints.append(int(field))

        if len(ints) != nint:
            raise ValueError(
                f"{path}: expected {nint} connectivity integers from 16i6, "
                f"got {len(ints)}"
            )

    except ValueError as exc:
        raise ValueError(
            f"{path}: invalid 16i6 connectivity data: {exc}"
        ) from exc

    ien = np.asarray(ints, dtype=np.int32).reshape(numel, 3)

    # Coordinates begin immediately after the fixed-width connectivity
    # records and are stored as 5e16.9 values per record.
    coord_text = "".join(lines[nint_lines:])
    xy = read_fixed_width_coordinates(coord_text, numnp)

    return ien, xy

def read_mesh_info(path):
    """
    Read mesh.info:

        numnp numel nbnd
        boundary_1_count
        boundary_1_nodes...
        ...
    """

    vals = np.fromstring(
        Path(path).read_text().replace("\r", "\n"),
        sep=" ",
        dtype=np.int64,
    )

    if len(vals) < 3:
        raise ValueError(
            f"{path}: expected numnp, numel, nbnd"
        )

    numnp = int(vals[0])
    numel = int(vals[1])
    nbnd = int(vals[2])

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
            np.asarray(
                vals[pos:pos + n],
                dtype=np.int32
            )
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
        raise ValueError(
            f"{path}: expected numnp1 and numel1"
        )

    return int(vals[0]), int(vals[1])


def read_surface_nodes(path):
    """
    Read surf_elem.info.

    The file contains surface NODE IDs, not element IDs.
    """

    vals = np.fromstring(
        Path(path).read_text().replace("\r", "\n"),
        sep=" ",
        dtype=np.int64,
    )

    if len(vals) < 1:
        raise ValueError(
            f"{path}: empty surface-node file"
        )

    n = int(vals[0])

    if len(vals) < n + 1:
        raise ValueError(
            f"{path}: expected {n} node IDs, "
            f"got {len(vals) - 1}"
        )

    nodes = np.asarray(
        vals[1:n + 1],
        dtype=np.int32
    )

    extra = len(vals) - (n + 1)

    return nodes, extra


# ============================================================
# OUTPUT WRITER
# ============================================================

def write_fortran_array(path, array, dtype):
    """Write native binary data matching the Fortran arrays."""

    np.asarray(
        array,
        dtype=dtype
    ).tofile(str(path))


# ============================================================
# STACK GEOMETRY HELPERS
# ============================================================

def validate_config(wing_slices, nslices):

    if nslices < 1:
        raise ValueError(
            "NSLICES must be >= 1"
        )

    if not 0 <= wing_slices <= nslices:
        raise ValueError(
            f"WING_SLICES must satisfy "
            f"0 <= WING_SLICES <= NSLICES ({nslices}); "
            f"got {wing_slices}"
        )


def plane_is_unfilled(
    plane0,
    wing_slices,
    nslices
):
    """Return True if the plane uses mesh.data."""

    if wing_slices == nslices:
        return True

    if wing_slices == 0:
        return False

    return plane0 < wing_slices - 1


def interval_uses_unfilled(
    interval1,
    wing_slices,
    nslices
):
    """
    Return True if a 1-based interval uses mesh.data connectivity.
    """

    if wing_slices == nslices:
        return True

    if wing_slices == 0:
        return False

    return interval1 < wing_slices


def plane_starts(
    np1,
    np2,
    nslices,
    wing_slices
):
    """
    Starting node index (0-based) for every stacked plane.
    """

    p = np.arange(
        nslices,
        dtype=np.int64
    )

    if wing_slices == nslices:

        return p * np1

    if wing_slices == 0:

        return p * np2

    n_mesh_planes = wing_slices - 1

    starts = np.empty(
        nslices,
        dtype=np.int64
    )

    wing_mask = p < n_mesh_planes

    starts[wing_mask] = (
        p[wing_mask] * np1
    )

    starts[~wing_mask] = (
        n_mesh_planes * np1
        + (p[~wing_mask] - n_mesh_planes) * np2
    )

    return starts


# ============================================================
# CONNECTIVITY
# ============================================================

def make_connectivity(
    ien1,
    ien2,
    np1,
    np2,
    ne1,
    ne2,
    wing_slices,
    nslices
):
    """
    Construct stacked wedge and hexa connectivity.

    This version preallocates the final arrays and fills them directly,
    avoiding the large np.vstack() copy at the end.
    """

    starts = plane_starts(
        np1,
        np2,
        nslices,
        wing_slices
    )

    nintervals = nslices - 1

    interval_ne = np.empty(
        nintervals,
        dtype=np.int64
    )

    use_ien1 = np.empty(
        nintervals,
        dtype=bool
    )

    for k in range(nintervals):

        use_ien1[k] = interval_uses_unfilled(
            k + 1,
            wing_slices,
            nslices
        )

        interval_ne[k] = (
            ne1
            if use_ien1[k]
            else ne2
        )

    total_elements = int(
        interval_ne.sum()
    )

    mien = np.empty(
        (total_elements, 6),
        dtype=np.int32
    )

    mienb = np.empty(
        (total_elements, 8),
        dtype=np.int32
    )

    interval_offsets = np.empty(
        nintervals,
        dtype=np.int64
    )

    offset = 0

    for k in range(nintervals):

        interval_offsets[k] = offset

        if use_ien1[k]:

            a = (
                ien1
                + int(starts[k])
            )

            layer_offset = np1

        else:

            a = (
                ien2
                + int(starts[k])
            )

            layer_offset = np2

        n = len(a)

        w = mien[
            offset:offset + n
        ]

        h = mienb[
            offset:offset + n
        ]

        # ----------------------------
        # Wedge
        # ----------------------------

        w[:, :3] = a
        w[:, 3:] = a + layer_offset

        # ----------------------------
        # Hexa
        # ----------------------------

        h[:, :3] = a

        h[:, 3] = a[:, 2]

        h[:, 4:7] = (
            a + layer_offset
        )

        h[:, 7] = h[:, 6]

        offset += n

    return (
        mienb,
        mien,
        interval_offsets,
        starts,
    )


# ============================================================
# COORDINATES
# ============================================================

def make_coordinates(
    xy1,
    xy2,
    z,
    np1,
    np2,
    nslices,
    wing_slices
):
    """
    Stack 2-D coordinates into 3-D coordinates.
    """

    starts = plane_starts(
        np1,
        np2,
        nslices,
        wing_slices
    )

    last_unfilled = plane_is_unfilled(
        nslices - 1,
        wing_slices,
        nslices
    )

    last_n = (
        np1
        if last_unfilled
        else np2
    )

    numnpt = int(
        starts[-1] + last_n
    )

    xyz = np.empty(
        (numnpt, 3),
        dtype=np.float64
    )

    for p in range(nslices):

        if plane_is_unfilled(
            p,
            wing_slices,
            nslices
        ):

            xy = xy1
            n = np1

        else:

            xy = xy2
            n = np2

        s = int(starts[p])

        xyz[
            s:s + n,
            :2
        ] = xy

        xyz[
            s:s + n,
            2
        ] = z[p]

    return xyz


# ============================================================
# BOUNDARY NODES
# ============================================================

def make_boundary_nodes(
    boundaries,
    np1,
    np2,
    nslices,
    wing_slices
):
    """
    Build seven 3-D boundary-node lists.
    """

    starts = plane_starts(
        np1,
        np2,
        nslices,
        wing_slices
    )

    bnd = [None] * 7

    # --------------------------------------------------------
    # Boundaries 1..4
    # --------------------------------------------------------

    for ib in range(4):

        nodes2 = boundaries[ib]

        bnd[ib] = np.concatenate(
            [
                nodes2 + int(s)
                for s in starts
            ]
        ).astype(
            np.int32,
            copy=False
        )

    # --------------------------------------------------------
    # Boundary 5
    # --------------------------------------------------------

    base = boundaries[4]

    parts = [
        base + int(s)
        for s in starts
    ]

    if np2 > np1:

        extra = np.arange(
            np1 + 1,
            np2 + 1,
            dtype=np.int32
        )

        for p in range(nslices):

            if not plane_is_unfilled(
                p,
                wing_slices,
                nslices
            ):

                parts.append(
                    extra + int(starts[p])
                )

    bnd[4] = np.concatenate(
        parts
    ).astype(
        np.int32,
        copy=False
    )

    # --------------------------------------------------------
    # Boundary 6
    # --------------------------------------------------------

    nfirst = (
        np1
        if plane_is_unfilled(
            0,
            wing_slices,
            nslices
        )
        else np2
    )

    bnd[5] = np.arange(
        1,
        nfirst + 1,
        dtype=np.int32
    )

    # --------------------------------------------------------
    # Boundary 7
    # --------------------------------------------------------

    nlast = (
        np1
        if plane_is_unfilled(
            nslices - 1,
            wing_slices,
            nslices
        )
        else np2
    )

    bnd[6] = (
        np.arange(
            1,
            nlast + 1,
            dtype=np.int32
        )
        + int(starts[-1])
    )

    return bnd


# ============================================================
# BOUNDARY ELEMENTS
# ============================================================

def make_boundary_elements(
    boundaries,
    wedge,
    np1,
    np2,
    ne1,
    ne2,
    nslices,
    wing_slices,
    interval_offsets
):
    """
    Build bnd_el for boundaries 1..7.

    Boundaries 1..4, 6, 7 retain the existing
    node-membership logic.
    """

    bnodes = make_boundary_nodes(
        boundaries,
        np1,
        np2,
        nslices,
        wing_slices
    )

    bnd_el = [None] * 7

    max_node = int(
        wedge.max()
    )

    # --------------------------------------------------------
    # Boundaries 1..4
    # --------------------------------------------------------

    for ib in range(4):

        nodes = bnodes[ib]

        node_count = np.bincount(
            nodes,
            minlength=max_node + 1
        )

        counts = node_count[wedge]

        selected = (
            counts.sum(
                axis=1,
                dtype=np.int64
            ) >= 4
        )

        bnd_el[ib] = (
            np.flatnonzero(selected)
            + 1
        ).astype(
            np.int32,
            copy=False
        )

    # Boundary 5 is handled separately.
    bnd_el[4] = np.empty(
        0,
        dtype=np.int32
    )

    # --------------------------------------------------------
    # Boundary 6
    # --------------------------------------------------------

    n_first_el = (
        ne1
        if plane_is_unfilled(
            0,
            wing_slices,
            nslices
        )
        else ne2
    )

    bnd_el[5] = np.arange(
        1,
        n_first_el + 1,
        dtype=np.int32
    )

    # --------------------------------------------------------
    # Boundary 7
    # --------------------------------------------------------

    n_last_el = (
        ne1
        if plane_is_unfilled(
            nslices - 1,
            wing_slices,
            nslices
        )
        else ne2
    )

    bnd_el[6] = (
        np.arange(
            1,
            n_last_el + 1,
            dtype=np.int32
        )
        + int(interval_offsets[-1])
    )

    return (
        bnd_el,
        bnodes
    )


# ============================================================
# AIRFOIL SURFACE
# ============================================================

def derive_airfoil_surface_elements(
    ien1,
    boundary5_nodes
):
    """
    Find 2-D triangle elements that lie on the physical airfoil
    boundary.

    A triangle belongs to the airfoil surface if at least two of its
    nodes belong to the boundary-5 node set.
    """

    max_index = int(
        max(
            ien1.max(),
            boundary5_nodes.max()
        )
    )

    node_set = np.zeros(
        max_index + 1,
        dtype=bool
    )

    node_set[
        boundary5_nodes
    ] = True

    count = node_set[
        ien1
    ].sum(
        axis=1
    )

    elems = (
        np.flatnonzero(
            count >= 2
        )
        .astype(
            np.int32,
            copy=False
        )
        + 1
    )

    return elems


def edge_to_face_slots_vectorized(
    ien1,
    surface_elems,
    boundary_nodes_mask
):
    """
    Vectorized equivalent of edge_to_face_slot().

    Returns mrng slot numbers:

        edge (1,2) -> 3
        edge (2,3) -> 4
        edge (3,1) -> 2
    """

    tri = ien1[
        surface_elems - 1
    ]

    hit01 = (
        boundary_nodes_mask[tri[:, 0]]
        &
        boundary_nodes_mask[tri[:, 1]]
    )

    hit12 = (
        boundary_nodes_mask[tri[:, 1]]
        &
        boundary_nodes_mask[tri[:, 2]]
    )

    hit20 = (
        boundary_nodes_mask[tri[:, 2]]
        &
        boundary_nodes_mask[tri[:, 0]]
    )

    hits = (
        hit01.astype(np.int8)
        +
        hit12.astype(np.int8)
        +
        hit20.astype(np.int8)
    )

    if np.any(hits != 1):

        bad = np.flatnonzero(
            hits != 1
        )

        elem = surface_elems[
            bad[0]
        ]

        raise ValueError(
            f"Could not identify a unique airfoil edge "
            f"for triangle {elem}"
        )

    slots = np.empty(
        len(surface_elems),
        dtype=np.int32
    )

    slots[hit01] = 3
    slots[hit12] = 4
    slots[hit20] = 2

    return slots


# ============================================================
# BOUNDARY 5
# ============================================================

def assign_correct_boundary5(
    mrng,
    ien1,
    boundary5_nodes,
    ne1,
    ne2,
    nslices,
    wing_slices,
    interval_offsets
):
    """
    Assign the physical airfoil/wing boundary (ID 5).

    Side surface:
        assigned only between wing planes.

    Finite-wing tip:
        the extra mesh1 elements in the first downstream interval
        receive ID 5 on mrng face slot 0.
    """

    if wing_slices == 0:
        return 0, 0, 0

    if wing_slices > nslices:
        raise ValueError(
            f"WING_SLICES={wing_slices} exceeds "
            f"NSLICES={nslices}"
        )

    # --------------------------------------------------------
    # Boundary-5 node mask
    # --------------------------------------------------------

    max_index = int(
        max(
            ien1.max(),
            boundary5_nodes.max()
        )
    )

    mask = np.zeros(
        max_index + 1,
        dtype=bool
    )

    mask[
        boundary5_nodes
    ] = True

    # --------------------------------------------------------
    # Find physical airfoil triangles.
    # --------------------------------------------------------

    surface_elems = derive_airfoil_surface_elements(
        ien1,
        boundary5_nodes
    )

    if len(surface_elems) == 0:
        raise ValueError(
            "No airfoil surface elements found "
            "from boundary-5 nodes"
        )

    # --------------------------------------------------------
    # Determine corresponding wedge face slot.
    # --------------------------------------------------------

    face_slots = edge_to_face_slots_vectorized(
        ien1,
        surface_elems,
        mask
    )

    # --------------------------------------------------------
    # Wing-side intervals.
    #
    # All assignments here are NumPy operations over only the
    # ~247 airfoil elements, so this portion is already tiny.
    # --------------------------------------------------------

    if wing_slices < nslices:

        interval_numbers = np.arange(
            1,
            wing_slices,
            dtype=np.int64
        )

    else:

        interval_numbers = np.arange(
            1,
            nslices,
            dtype=np.int64
        )

    side_count = 0

    for interval1 in interval_numbers:

        interval0 = int(
            interval1 - 1
        )

        offset = int(
            interval_offsets[interval0]
        )

        rows = (
            offset
            + surface_elems.astype(
                np.int64
            )
            - 1
        )

        mrng[
            rows,
            face_slots
        ] = 5

        side_count += len(rows)

    # --------------------------------------------------------
    # Wingtip cap.
    # --------------------------------------------------------

    cap_count = 0

    if 0 < wing_slices < nslices:

        cap_interval0 = (
            wing_slices - 1
        )

        offset = int(
            interval_offsets[
                cap_interval0
            ]
        )

        if ne2 < ne1:
            raise ValueError(
                f"mesh1 element count ne2={ne2} "
                f"is smaller than mesh element count ne1={ne1}"
            )

        extra_count = (
            ne2 - ne1
        )

        if extra_count:

            rows = (
                offset
                + np.arange(
                    ne1,
                    ne2,
                    dtype=np.int64
                )
            )

            mrng[
                rows,
                0
            ] = 5

            cap_count = extra_count

    return (
        int(len(surface_elems)),
        int(side_count),
        int(cap_count)
    )


# ============================================================
# MRNG
# ============================================================

def make_mrng_base(
    wedge,
    bnd_el,
    bnd_nodes,
    chunk_size=1_000_000
):
    """
    Construct mrng for boundaries 1..4, 6, 7.

    This is the main performance-critical routine.

    The original code performed a Python loop over every selected
    element. This version performs the element calculations and
    assignment with NumPy.

    The persistent iface behavior of the original code is retained.
    """

    numelt = len(wedge)

    mrng = np.zeros(
        (numelt, 5),
        dtype=np.int32
    )

    max_node = int(
        wedge.max()
    )

    # Original behavior:
    #
    #   1 -> node weight 1
    #   2 -> node weight 2
    #   ...
    #   6 -> node weight 6
    #
    # icount therefore identifies the wedge face.
    #
    # iface persists if an icount value isn't one of:
    # 6, 15, 14, 12, 16.
    iface = 0

    # --------------------------------------------------------
    # Process each boundary independently.
    # --------------------------------------------------------

    for ib in (0, 1, 2, 3, 5, 6):

        nodes = bnd_nodes[ib]

        # Convert to int32 to reduce memory traffic.
        node_count = np.bincount(
            nodes,
            minlength=max_node + 1
        ).astype(
            np.int32,
            copy=False
        )

        elems = (
            bnd_el[ib] - 1
        )

        # ----------------------------------------------------
        # Chunk large boundary-element lists.
        # ----------------------------------------------------

        for start in range(
            0,
            len(elems),
            chunk_size
        ):

            eidx = elems[
                start:start + chunk_size
            ]

            w = wedge[eidx]

            # ------------------------------------------------
            # Instead of:
            #
            #   counts = node_count[w]
            #   icount = (counts * [1..6]).sum(...)
            #
            # directly gather the six node counts.
            #
            # This avoids constructing a giant (N,6) matrix
            # with a multiplication temporary.
            # ------------------------------------------------

            c0 = node_count[w[:, 0]]
            c1 = node_count[w[:, 1]]
            c2 = node_count[w[:, 2]]
            c3 = node_count[w[:, 3]]
            c4 = node_count[w[:, 4]]
            c5 = node_count[w[:, 5]]

            ncount = (
                c0.astype(np.int64)
                + c1
                + c2
                + c3
                + c4
                + c5
            )

            icount = (
                c0.astype(np.int64)
                + 2 * c1
                + 3 * c2
                + 4 * c3
                + 5 * c4
                + 6 * c5
            )

            # ------------------------------------------------
            # Determine iface values without iterating through
            # individual elements in Python.
            #
            # raw_slot == -1 means "keep previous iface".
            # ------------------------------------------------

            n = len(icount)

            raw_slot = np.full(
                n,
                -1,
                dtype=np.int8
            )

            raw_slot[
                icount == 6
            ] = 0

            raw_slot[
                icount == 15
            ] = 1

            raw_slot[
                icount == 14
            ] = 2

            raw_slot[
                icount == 12
            ] = 3

            raw_slot[
                icount == 16
            ] = 4

            valid = (
                raw_slot >= 0
            )

            if np.any(valid):

                valid_positions = np.flatnonzero(
                    valid
                )

                last_valid = int(
                    valid_positions[-1]
                )

                # Index of the most recent valid iface for each
                # element in this chunk.
                positions = np.arange(
                    n,
                    dtype=np.int64
                )

                last_seen = np.maximum.accumulate(
                    np.where(
                        valid,
                        positions,
                        0
                    )
                )

                slot = raw_slot[
                    last_seen
                ].astype(
                    np.int8,
                    copy=False
                )

                # Values before the first valid iface use the
                # persistent iface from the previous chunk/boundary.
                first_valid = int(
                    valid_positions[0]
                )

                if first_valid > 0:

                    slot[
                        :first_valid
                    ] = iface

                # Update persistent iface for the next chunk.
                iface = int(
                    raw_slot[
                        last_valid
                    ]
                )

            else:

                # No iface transition anywhere in this chunk.
                slot = np.full(
                    n,
                    iface,
                    dtype=np.int8
                )

            # ------------------------------------------------
            # Boundary assignment.
            # ------------------------------------------------

            selected = (
                ncount >= 3
            )

            if np.any(selected):

                rows = eidx[
                    selected
                ]

                cols = slot[
                    selected
                ]

                mrng[
                    rows,
                    cols
                ] = ib + 1

    return mrng


# ============================================================
# DEBUG COORDINATE OUTPUT
# ============================================================

def write_debug_coordinates(
    path,
    xyz
):
    """
    Write fort.1 in the same 24-character coordinate style as
    the original Python writer, but using NumPy instead of a
    Python loop over every node.
    """

    np.savetxt(
        path,
        xyz,
        fmt="%24.16g%24.16g%24.16g"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        default=".",
        help=(
            "directory containing mesh.info, mesh.data, "
            "mesh1.info, mesh1.data, surf_elem.info, zall.data"
        )
    )

    parser.add_argument(
        "--output",
        default=".",
        help=(
            "directory for mien/mxyz/mrng/mienb"
        )
    )

    parser.add_argument(
        "--wing-slices",
        type=int,
        default=WING_SLICES,
        help=(
            "number of spanwise wing slices/planes"
        )
    )

    parser.add_argument(
        "--fort1",
        action="store_true",
        help=(
            "do not write fort.1; useful for measuring "
            "mesh stacking time without text output"
        )
    )

    args = parser.parse_args()

    wing_slices = args.wing_slices

    inp = Path(
        args.input
    )

    out = Path(
        args.output
    )

    out.mkdir(
        parents=True,
        exist_ok=True
    )

    # ========================================================
    # TOTAL TIMER
    # ========================================================

    total_t0 = time.perf_counter()

    # ========================================================
    # READ METADATA
    # ========================================================

    n1, e1, nbnd_nodes = read_mesh_info(
        inp / "mesh.info"
    )

    n2, e2 = read_mesh1_info(
        inp / "mesh1.info"
    )

    if len(nbnd_nodes) != 5:
        raise ValueError(
            f"mesh.info contains {len(nbnd_nodes)} boundaries; "
            f"expected 5"
        )

    z = np.atleast_1d(
        np.loadtxt(
            inp / "zall.data",
            dtype=np.float64
        )
    )

    if len(z) != NSLICES:
        raise ValueError(
            f"zall.data: expected {NSLICES} z locations, "
            f"got {len(z)}"
        )

    nslices = NSLICES

    validate_config(
        wing_slices,
        nslices
    )

    # ========================================================
    # READ MESH DATA
    # ========================================================

    t_read0 = time.perf_counter()

    ien1, xy1 = read_fixed_width_mesh(
        inp / "mesh.data",
        n1,
        e1
    )

    ien2, xy2 = read_fixed_width_mesh(
        inp / "mesh1.data",
        n2,
        e2
    )

    surf_nodes, extra_values = read_surface_nodes(
        inp / "surf_elem.info"
    )

    # ========================================================
    # Validate surface-node file
    # ========================================================

    if not np.array_equal(
        surf_nodes,
        nbnd_nodes[4]
    ):

        if (
            len(surf_nodes)
            >
            len(nbnd_nodes[4])
            and
            np.array_equal(
                surf_nodes[
                    :len(nbnd_nodes[4])
                ],
                nbnd_nodes[4]
            )
        ):

            print(
                f"Warning: surf_elem.info has "
                f"{len(surf_nodes) - len(nbnd_nodes[4])} "
                f"stale trailing node(s); ignoring them."
            )

            surf_nodes = nbnd_nodes[4].copy()

        else:

            vals = np.fromstring(
                (
                    inp / "surf_elem.info"
                ).read_text().replace(
                    "\r",
                    "\n"
                ),
                sep=" ",
                dtype=np.int64
            )

            nbase = len(
                nbnd_nodes[4]
            )

            if (
                len(vals) >= nbase + 2
                and
                np.array_equal(
                    vals[
                        1:1 + nbase
                    ],
                    nbnd_nodes[4]
                )
            ):

                print(
                    f"Warning: using first {nbase} "
                    f"surface nodes from surf_elem.info; "
                    f"ignoring stale trailing value(s)."
                )

                surf_nodes = (
                    nbnd_nodes[4].copy()
                )

            else:

                raise ValueError(
                    "surf_elem.info node list does not "
                    "match mesh.info boundary 5"
                )

    read_elapsed = (
        time.perf_counter()
        - t_read0
    )

    # ========================================================
    # STACKING TIMER
    # ========================================================

    t0 = time.perf_counter()

    # ========================================================
    # CONNECTIVITY
    # ========================================================

    t_stage = time.perf_counter()

    (
        mienb,
        mien,
        interval_offsets,
        plane_starts_
    ) = make_connectivity(
        ien1,
        ien2,
        n1,
        n2,
        e1,
        e2,
        wing_slices,
        nslices
    )

    connectivity_elapsed = (
        time.perf_counter()
        - t_stage
    )

    # ========================================================
    # COORDINATES
    # ========================================================

    t_stage = time.perf_counter()

    xyz = make_coordinates(
        xy1,
        xy2,
        z,
        n1,
        n2,
        nslices,
        wing_slices
    )

    coordinates_elapsed = (
        time.perf_counter()
        - t_stage
    )

    # ========================================================
    # BOUNDARY ELEMENTS
    # ========================================================

    t_stage = time.perf_counter()

    (
        bnd_el,
        bnd_nodes
    ) = make_boundary_elements(
        nbnd_nodes,
        mien,
        n1,
        n2,
        e1,
        e2,
        nslices,
        wing_slices,
        interval_offsets
    )

    boundary_elapsed = (
        time.perf_counter()
        - t_stage
    )

    # ========================================================
    # MRNG BASE
    # ========================================================

    t_stage = time.perf_counter()

    mrng = make_mrng_base(
        mien,
        bnd_el,
        bnd_nodes
    )

    mrng_elapsed = (
        time.perf_counter()
        - t_stage
    )

    # ========================================================
    # BOUNDARY 5
    # ========================================================

    nsurf = 0
    side_count = 0
    cap_count = 0

    t_stage = time.perf_counter()

    if wing_slices > 0:

        (
            nsurf,
            side_count,
            cap_count
        ) = assign_correct_boundary5(
            mrng,
            ien1,
            nbnd_nodes[4],
            e1,
            e2,
            nslices,
            wing_slices,
            interval_offsets
        )

    boundary5_elapsed = (
        time.perf_counter()
        - t_stage
    )

    # ========================================================
    # BINARY OUTPUT
    # ========================================================

    t_stage = time.perf_counter()

    write_fortran_array(
        out / "mienb",
        mienb,
        np.int32
    )

    write_fortran_array(
        out / "mien",
        mien,
        np.int32
    )

    write_fortran_array(
        out / "mxyz",
        xyz,
        np.float64
    )

    write_fortran_array(
        out / "mrng",
        mrng,
        np.int32
    )

    binary_elapsed = (
        time.perf_counter()
        - t_stage
    )

    # ========================================================
    # FORT.1
    # ========================================================

    fort1_elapsed = 0.0

    if args.fort1:

        t_stage = time.perf_counter()

        write_debug_coordinates(
            out / "fort.1",
            xyz
        )

        fort1_elapsed = (
            time.perf_counter()
            - t_stage
        )

    # ========================================================
    # TIMING
    # ========================================================

    stacking_elapsed = (
        time.perf_counter()
        - t0
    )

    total_elapsed = (
        time.perf_counter()
        - total_t0
    )

    # ========================================================
    # OUTPUT
    # ========================================================

    print()
    print("=" * 60)
    print("STACK COMPLETE")
    print("=" * 60)

    print(
        f"Wing slices              : "
        f"{wing_slices} / {NSLICES}"
    )

    if (
        0 < wing_slices < NSLICES
    ):

        print(
            f"Wing side intervals      : "
            f"1..{wing_slices - 1}"
        )

        print(
            f"Wingtip cap              : "
            f"interval {wing_slices}, face 1"
        )

    print()

    print(
        f"2D mesh 1 nodes           : {n1}"
    )

    print(
        f"2D mesh 1 elements        : {e1}"
    )

    print(
        f"2D mesh 2 nodes           : {n2}"
    )

    print(
        f"2D mesh 2 elements        : {e2}"
    )

    print(
        f"3D nodes                  : {len(xyz)}"
    )

    print(
        f"3D wedge elements        : {len(mien)}"
    )

    print(
        f"3D hexa elements         : {len(mienb)}"
    )

    print(
        f"z locations              : {z}"
    )

    print(
        f"airfoil 2D elements found: {nsurf}"
    )

    print(
        f"ID 5 side faces assigned : {side_count}"
    )

    print(
        f"ID 5 wingtip cap faces   : {cap_count}"
    )

    print(
        f"ID 5 total               : "
        f"{np.count_nonzero(mrng == 5)}"
    )

    print(
        f"mrng nonzero entries     : "
        f"{np.count_nonzero(mrng)}"
    )

    print()
    print("TIMING")
    print("-" * 60)

    print(
        f"Mesh/input reading       : "
        f"{read_elapsed:.3f} s"
    )

    print(
        f"Connectivity             : "
        f"{connectivity_elapsed:.3f} s"
    )

    print(
        f"Coordinates              : "
        f"{coordinates_elapsed:.3f} s"
    )

    print(
        f"Boundary construction    : "
        f"{boundary_elapsed:.3f} s"
    )

    print(
        f"MRNG base                : "
        f"{mrng_elapsed:.3f} s"
    )

    print(
        f"Boundary 5               : "
        f"{boundary5_elapsed:.3f} s"
    )

    print(
        f"Binary output            : "
        f"{binary_elapsed:.3f} s"
    )

    if not args.fort1:

        print(
            f"fort.1                   : "
            f"SKIPPED"
        )

    else:

        print(
            f"fort.1                   : "
            f"{fort1_elapsed:.3f} s"
        )

    print()

    print(
        f"Stacking subtotal        : "
        f"{stacking_elapsed:.3f} s"
    )

    print(
        f"Total elapsed            : "
        f"{total_elapsed:.3f} s"
    )

    print()

    print(
        f"mienb bytes              : "
        f"{(out / 'mienb').stat().st_size}"
    )

    print(
        f"mien bytes               : "
        f"{(out / 'mien').stat().st_size}"
    )

    print(
        f"mxyz bytes               : "
        f"{(out / 'mxyz').stat().st_size}"
    )

    print(
        f"mrng bytes               : "
        f"{(out / 'mrng').stat().st_size}"
    )

    if args.fort1:

        print(
            f"fort.1 bytes             : "
            f"{(out / 'fort.1').stat().st_size}"
        )

    print("=" * 60)


if __name__ == "__main__":
    main()