InFile="~/dowell_lab/DMSO2_3.BedGraph"
OutFile="~/dowell_lab/DMSO2_3_500kb_windows.BedGraph"

bedtools makewindows -b $InFile -w 500 > $OutFile
