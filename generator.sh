#!/bin/bash

# 0: pulse
# 1: 35-d modulation

mode=1
psi=$(echo "0.0 * 8 * a(1)" | bc -l)

time ./pulse $mode $psi > output.data
sed '/^$/q' output.data > x_ray_flux_as_superorbital_cycle.data
sed '1,/^$/d' output.data > arcs.data

#gnuplot arcs_map.plt
