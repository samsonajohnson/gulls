#!/bin/bash

if [ $# -lt 1 ] || [ $# -gt 3 ]; then
    echo "./listcols <file> {<row>} {<separator>}"
fi

if [ $# -ge 2 ]; then
    row=$2
    echo "Indexed at 1"
    if [ $# -eq 3 ]; then
	sep=$3
	head -n$row $1 | tail -n1 | awk -v FS=$sep '{for(i=1;i<=NF;i++){if($1~"^#$") j=-1; else j=0; print i+j,$(i)}}'
    else
	head -n$row $1 | tail -n1 | awk '{for(i=1;i<=NF;i++){if($1~"^#$") j=-1; else j=0; print i+j,$(i)}}'
    fi
    echo "Indexed at 1"

else

echo "Indexed at 1"
head -n1 $1 | awk '{for(i=1;i<=NF;i++) print i,$(i)}'
echo "Indexed at 1"

fi
