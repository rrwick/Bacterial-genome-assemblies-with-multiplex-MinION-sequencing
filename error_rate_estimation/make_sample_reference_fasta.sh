mkdir -p references

for BARCODE_NUMBER in 01 02 03 04 05 06 07 08 09 10 11 12; do
    BARCODE=barcode$BARCODE_NUMBER

    # MUMmer: find sequences where ABySS and Velvet are in exact agreement
    # Grep: extract just the sequences
    # Awk: trim 100 bp from the sequence ends and build a FASTA file
    mummer -mum -b -s -n -l 10000 abyss/$BARCODE/abyss_assembly-1.fa velvet/$BARCODE/contigs.fa | grep -i -P "[acgt]{10000,}" | awk '{name += 1; print ">"name"\n"substr(toupper($0), 100, length($0)-200)}' > references/"$BARCODE".fasta

done;

