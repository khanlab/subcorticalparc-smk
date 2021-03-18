#!/bin/bash

in_seg=$1
out_seg=$2
k=$3


i=1; 
for centroid in `c3d $in_seg -split -foreach  -centroid -endfor | grep VOX | awk -F ',' '{print $2}' | tail -n $k `;
do 
    echo $centroid,$i
    i=$((i+1))
done  > centroids_$k.txt


i=1;
temp_segs=""
c3d $in_seg -scale 0 -o $out_seg
for line in `sort -r centroids_$k.txt`
do
    echo $line
    in_lbl=${line##*,}
    out_lbl=$i
    echo c3d $in_seg -threshold $in_lbl $in_lbl $out_lbl 0 $out_seg -add -o $out_seg
     c3d $in_seg -threshold $in_lbl $in_lbl $out_lbl 0 $out_seg -add -o $out_seg
    i=$((i+1))
done


