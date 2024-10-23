#!/bin/bash

if [ $# -ne 3 ]; then
  echo "Usage: $0 <detectorlist> <field> <output root>"
  exit
fi

/home/emeade/Documents/Roman/gulls/src/images/freeColour $1 /home/emeade/Documents/Roman/input/faint/out-$2 0.000020 /home/emeade/Documents/Roman/input/moderate/out-$2 0.000020 /home/emeade/Documents/Roman/input/bright/out-$2 0.001000 $3
