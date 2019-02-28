## Get Standard  Chromosomes from a BED file
awk '/^chr[0-9]*\t/ {printf("%s\t0\t%s\n",$1,$2);}'
