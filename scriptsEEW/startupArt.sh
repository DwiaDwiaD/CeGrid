#!/bin/bash

# ------------------------------------
# Just making some startup art for fun
# ------------------------------------

# the color format: \033[38;2;R;G;Bm

: << 'COLORS'
# Regular
# \033[0;30m  # black
# \033[0;31m  # red
# \033[0;32m  # green
# \033[0;33m  # yellow
# \033[0;34m  # blue
# \033[0;35m  # magenta
# \033[0;36m  # cyan
# \033[0;37m  # white

# Bold (brighter)
# \033[1;31m  # bright red
# \033[1;32m  # bright green
# \033[1;34m  # bright blue
COLORS

center() {
        width=$(tput cols)
        text="$1"

        # remove ANSI escape sequences for length calculation
        cln_txt=$(echo -e "$text" | sed 's/\x1b\[[0-9;]*m//g')
        padding=$(( (width - ${#cln_txt}) / 2 ))
        printf "%*s%s\n" $padding "" "$text"
}

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

printf "\n\n"
line "-_"

echo
figlet_centered big "CeGrid"
figlet_centered small " 2-D AIRFOIL MESHER" 27 33 39 45 51 57
# figlet_centered small "MESHER" 27 33 39 45 51 57
echo
center "YeAh, BaBy!!!"
echo

line "_-"
printf "\n\n"

center "[ Initializing CFD Pipeline ]"
center "[ $(date '+%Y-%m-%d %H:%M:%S') ]"
echo
echo
