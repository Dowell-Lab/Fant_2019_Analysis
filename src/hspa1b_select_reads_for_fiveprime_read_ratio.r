<PO_... awk '{if ($1 == "chr6" && $2 > 31827737 && $2 < 31828000) print $0}' > reads_ct.bed
<C413_1_S3_R1_001.pos.sort.bedGraph awk '{if ($1 == "chr6" && $2 > 31827737 && $2 < 31828000) print $0}' > reads_ko.bed
