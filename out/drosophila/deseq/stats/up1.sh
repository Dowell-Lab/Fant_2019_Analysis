#!/bin/zsh
awk '{if ($5 == "NA") {} else if ($5 > 0) {print $1, $5} else if ($5 < 0) {print $1, -$5}}' | tail -n+2 | sort -nk2,2
