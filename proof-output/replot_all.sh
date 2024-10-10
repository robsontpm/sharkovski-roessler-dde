#!/bin/bash

OUTDIR="."

# # gather data for pictures:
# find $OUTDIR/ -mindepth 1 -type d | xargs -I{} bash -c "echo "" > {}-all.dat"
# find $OUTDIR/ -mindepth 1 -type d | xargs -I{} bash -c "find {} -type f | grep '.dat' | xargs -I XX grep -v '^\s*#' XX >> {}-all.dat"

# #ls | grep "\-all.dat" | sed 's/\(.*\)\-all.dat\.*/\1/' | xargs -I{} echo "sed 's/ \+/ /g' {}-all.dat | sort -k6,6r -k8,8nr > {}-srt.dat"
# ls | grep "\-all.dat" | sed 's/\(.*\)\-all.dat\.*/\1/' | xargs -I{} echo "sed 's/ \+/ /g' {}-all.dat > {}-srt.dat"

gnuplot plot_C0_tail_coords.gp
gnuplot plot_C1_tail_coords.gp
gnuplot plot_C2_tail_coords.gp
gnuplot plot_main.gp
gnuplot plot_X0_tail.gp