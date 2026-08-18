#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: ./mesh.sh <filename_without_extension>"
    exit 1
fi

# -------- INPUTS --------
name=$1
CHORD=$2
AOA_DEG=$3
SPAN=$4
SCRIPT_DIR=${5:-essential}
FIRST_LAYER=$6
LAYERS=$7
FILLED=$8
DIM=$9
QUADS=${10}
WING_LAYERS=${11}
WINGSPAN=${12}
GROUND_RES=${13}
WAKE_RES=${14}
INLET_RES=${15}
OUTLET_RES=${16}
IFBOTH=${17}

echo
echo "Airfoil: $name"
echo "Chord: $CHORD"
echo "AoA (deg): $AOA_DEG"
echo "Span: $SPAN"
echo

mkdir -p "meshesEEW"

# For 2D we deliberately generate the mesh ONCE with the filled airfoil
# surface present, then derive the filled and unfilled views from that same
# in-memory/generated mesh. This makes their common node/element numbering
# identical and avoids the numbering changes produced by two independent Gmsh
# runs.
if [ "$IFBOTH" = true ]; then
    echo "NOTE: 2D mode -> one Gmsh mesh generation; creating both filled and unfilled outputs"
    DIM=2 # this splitting is only valid for 2D mesh

    combined="meshesEEW/${name}Combined2D.msh"
    unfilled="meshesEEW/${name}Unfilled2D.msh"
    filled="meshesEEW/${name}Filled2D.msh"

    gmsh "$SCRIPT_DIR/scriptsEEW/MasterAirfoilEEW.geo" -2 \
        -setnumber chord "$CHORD" \
        -setnumber AoA "$AOA_DEG" \
        -setnumber span "$SPAN" \
        -setnumber wingspan "$WINGSPAN" \
        -setnumber firstlayer "$FIRST_LAYER" \
        -setnumber layers "$LAYERS" \
        -setnumber wing_layers "$WING_LAYERS" \
        -setnumber filled 1 \
        -setnumber dim "$DIM" \
        -setnumber quads "$QUADS" \
        -setnumber groundres "$GROUND_RES" \
        -setnumber wakeres "$WAKE_RES" \
        -setnumber inletres "$INLET_RES" \
        -setnumber outletres "$OUTLET_RES" \
        -format msh2 \
        -nt 0 \
        -o "$combined"

    python3 "$SCRIPT_DIR/scriptsEEW/split_common_mesh.py" \
        "$combined" "$unfilled" "$filled"

    rm -f "$combined"

    echo "2D mesh generation complete."
    echo "  Unfilled: $unfilled"
    echo "  Filled  : $filled"
    exit 0
fi


if [ "$FILLED" -eq 1 ]; then
    echo "NOTE: filled mesh"
    name="${name}Filled${DIM}D"
else
    echo "NOTE: Unfilled mesh"
    name="${name}Unfilled${DIM}D"
fi

gmsh "$SCRIPT_DIR/scriptsEEW/MasterAirfoilEEW.geo" -3 \
    -setnumber chord "$CHORD" \
    -setnumber AoA "$AOA_DEG" \
    -setnumber span "$SPAN" \
    -setnumber wingspan "$WINGSPAN" \
    -setnumber firstlayer "$FIRST_LAYER" \
    -setnumber layers "$LAYERS" \
    -setnumber wing_layers "$WING_LAYERS" \
    -setnumber filled "$FILLED" \
    -setnumber dim "$DIM" \
    -setnumber quads "$QUADS" \
    -setnumber groundres "$GROUND_RES" \
    -setnumber wakeres "$WAKE_RES" \
    -setnumber inletres "$INLET_RES" \
    -setnumber outletres "$OUTLET_RES" \
    -format msh2 \
    -nt 0 \
    -o "meshesEEW/${name}.msh"
