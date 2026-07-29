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

echo
echo "Airfoil: $name"
echo "Chord: $CHORD"
echo "AoA (deg): $AOA_DEG"
echo "Span: $SPAN"
echo

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
    -setnumber firstlayer "$FIRST_LAYER" \
    -setnumber layers "$LAYERS" \
    -setnumber filled "$FILLED" \
    -setnumber dim "$DIM" \
    -setnumber quads "$QUADS" \
    -format msh2 \
    -nt 0 \
    -o "meshesEEW/${name}.msh"
