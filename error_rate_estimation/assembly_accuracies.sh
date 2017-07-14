
# Pass the directory with the assemblies (in subdirectories barcode01, barcode02, etc.) and the assembly filename.
ASSEMBLY_DIR=$1
ASSEMBLY_NAME=$2

echo "Barcode\tAccuracy\tError rate\tDistance between errors"

for BARCODE_NUMBER in 01 02 03 04 05 06 07 08 09 10 11 12; do
    BARCODE=barcode$BARCODE_NUMBER

    printf $BARCODE_NUMBER"\t"

    REFERENCE=references/"$BARCODE".fasta
    ASSEMBLY=$ASSEMBLY_DIR/$BARCODE/$ASSEMBLY_NAME

    # Build BLAST database.
    makeblastdb -dbtype nucl -in $ASSEMBLY > /dev/null

    # Get the best BLAST hit for each query.
    blastn -db $ASSEMBLY -query $REFERENCE -outfmt 6 | sort -nk1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -nk1,1 > blast_hits
   
    # Get the total error rates.
    python3 get_error_rate.py blast_hits
    printf "\n"

    rm blast_hits

done
