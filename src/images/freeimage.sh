#!/bin/bash

if [ $# -ne 8 ]; then
  echo "Usage: $0 <detector> <texp> <nstack> <nstars> <magcol> <maglimit> <field> <output root>"
  exit
fi

~/src/images/freetest $1 $2 $3 $4 $5 $6 ~/mabuls/starfields/EUCLID-h-faint/out-$7 0.000020 ~/mabuls/starfields/EUCLID-h-moder1/out-$7 0.003 ~/mabuls/starfields/EUCLID-h-moder2/out-$7 0.003 ~/mabuls/starfields/EUCLID-h-bright/out-$7 0.001000 $8 
