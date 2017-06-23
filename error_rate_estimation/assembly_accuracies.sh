
# Pass the directory with the assemblies (in subdirectories barcode01, barcode02, etc.) and the assembly filename.
ASSEMBLY_DIR=$1
ASSEMBLY_NAME=$2

UNICYCLER_ILLUMINA_ASSEMBLIES=10_unicycler_illumina_only_assemblies

echo "Barcode\tAccuracy\tError rate\tDistance between errors"

for BARCODE_NUMBER in 01 02 03 04 05 06 07 08 09 10 11 12; do
    BARCODE=barcode$BARCODE_NUMBER

    printf $BARCODE_NUMBER"\t"

    # Prepare the BLAST queries from a pre-RR assembly (big contigs only, ends trimmed off).
    head -n 25 $UNICYCLER_ILLUMINA_ASSEMBLIES/$BARCODE/002_overlaps_removed.gfa | cut -f1-3 | sed 's|S\t|>|' | sed 's|\t|\n|' > blast_queries.fasta
    python3 trim_blast_queries.py
    mv blast_queries2.fasta blast_queries.fasta
    ASSEMBLY=$ASSEMBLY_DIR/$BARCODE/$ASSEMBLY_NAME

    # Build BLAST database.
    makeblastdb -dbtype nucl -in $ASSEMBLY > /dev/null

    # Get the best BLAST hit for each query.
    blastn -db $ASSEMBLY -query blast_queries.fasta -outfmt 6 | sort -nk1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -nk1,1 > blast_hits
   
    # Get the total error rates.
    python3 get_error_rate.py blast_hits
    printf "\n"

    rm blast_queries.fasta blast_hits

done
