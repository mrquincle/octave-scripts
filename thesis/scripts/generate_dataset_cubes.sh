#!/bin/sh

cd ..

mkdir -p data/cubes

i=0
while [ "$i" -lt 100 ]; do
  fname="cubes_$(date '+%Y%m%d_%H%M%S_%3N').pnts.data"
  echo "Generate $i: $fname"
  ./dpmobjectrnd.m config/cube_config.m "data/cubes/$fname" skip_graphs
  i=$(expr $i + 1)
done
