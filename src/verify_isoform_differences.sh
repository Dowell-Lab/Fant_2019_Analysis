#!/bin/bash

paste \
		<(cut -f 4 C413_1_S3_R1_001.trim.sorted.isoform_max) \
		<(cut -f 4 C413_2_S4_R1_001.trim.sorted.isoform_max) \
		<(cut -f 4 PO_1_S1_R1_001.trim.sorted.isoform_max) \
		<(cut -f 4 PO_2_S2_R1_001.trim.sorted.isoform_max) \
		> isoforms.allmax

awkcmd2='{
# We filter the broadest categories first to simplify logic for later ones
# All 4 Match
if ($1 == $2 && $1 == $3 && $1 == $4 && $2 == $3 && $2 == $4 && $3 == $4)
	 print "AllMatch"
# No Matches
else if ($1 != $2 && $1 != $3 && $1 != $4 && $2 != $3 && $2 != $4 && $3 != $4)
	 print "NoneMatch"
# 3 Match, 1 MisMatch
else if (($1 != $2) && ($2 == $3 && $2 == $4 && $3 == $4))
	print "FirstTermMisMatch"
else if (($2 != $1) && ($1 == $3 && $1 == $4 && $3 == $4))
	print "SecondTermMisMatch"
else if (($3 != $1) && ($1 == $2 && $1 == $4 && $2 == $4))
	print "ThirdTermMisMatch"
else if (($4 != $1) && ($1 == $2 && $1 == $3 && $2 == $3))
	print "FourthTermMisMatch"
# Two Pairwise Matches
else if ($1 == $2 && $3 == $4)
	 print "ConditionPairsMatch"
else if ($1 == $3 && $2 == $4)
	 print "DayPairsMatch"
else if ($1 == $4 && $2 == $3)
	 print "DayPairsMisMatch"
# One Pairwise Match
else if ($1 == $2)
	 print "ControlPairMatch"
else if ($3 == $4)
	 print "TreatmentPairMatch"
else if ($1 == $3)
	 print "FirstRepPairMatch"
else if ($2 == $4)
	 print "SecondRepPairMatch"
else if ($1 == $4)
	 print "1/4PairMatch"
else if ($2 == $3)
	 print "2/3PairMatch"
# No Matches
else
	print "ERROR"
}'

awk "$awkcmd2" isoforms.allmax | sort | uniq -c | sort -rnk1,1 | awk '{print $2, ($1 / 19259) * 100}' > matchCounts
