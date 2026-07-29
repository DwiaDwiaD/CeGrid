#!/bin/bash


line() {
    pattern="$1"
    width=$(tput cols)
    out=""

    for ((i=0; i<width; i++)); do
        out+=${pattern:i%${#pattern}:1}
    done

    echo "$out"
}

figlet_centered() {
    local font="$1"
    shift
    local text="$1"
    shift

    local colors_local=("$@")
    if [ ${#colors_local[@]} -eq 0 ]; then
        colors_local=(196 202 208 214 220 226)
    fi

    mapfile -t lines < <(figlet -f "$font" "$text" 2>/dev/null || figlet "$text")

    max=0
    for l in "${lines[@]}"; do
        (( ${#l} > max )) && max=${#l}
    done

    width=$(tput cols)
    offset=$(( (width - max) / 2 ))

    i=0
    for l in "${lines[@]}"; do
        colored=$'\033[38;5;'"${colors_local[$i]}"'m'"$l"$'\033[0m'
        printf "%*s%s\n" "$offset" "" "$colored"
        ((i=(i+1)%${#colors_local[@]}))
    done
}

line "-_"
echo
# figlet_centered standard "LE - FIN" 27 33 39 45 51
figlet_centered standard "LE - FIN" 51 39 27 220 208 196
echo
line "_-"
echo
