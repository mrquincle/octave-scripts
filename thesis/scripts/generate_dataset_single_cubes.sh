#!/bin/sh

cd ..

mkdir -p data/single_cubes

i=0
N=100

while [ "$i" -lt "$N" ]; do
  fname="single_cubes_$(date '+%Y%m%d_%H%M%S_%3N').pnts.data"
  echo "Generate $i: $fname"
  ./dpmobjectrnd.m config/single_cube_config.m "data/single_cubes/$fname" skip_graphs
  i=$(expr $i + 1)
done
