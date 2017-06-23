

echo "Sample\tReplicon\tSize\tIllumina depth\tNanopore depth\tIllumina depth (normalised)\tNanopore depth (normalised)"

for BARCODE_NUMBER in 01 02 03 04 05 06 07 08 09 10 11 12; do
    BARCODE=barcode$BARCODE_NUMBER

    if [ $BARCODE_NUMBER = "01" ]; then REPS="1 2 3"; fi
    if [ $BARCODE_NUMBER = "02" ]; then REPS="1 2 3"; fi
    if [ $BARCODE_NUMBER = "03" ]; then REPS="1 2 3 4 5 6 7"; fi
    if [ $BARCODE_NUMBER = "04" ]; then REPS="1 2 3 4 5 6 7"; fi
    if [ $BARCODE_NUMBER = "05" ]; then REPS="1 2 3 4 5 6 7"; fi
    if [ $BARCODE_NUMBER = "06" ]; then REPS="1 2 3 4 5 6 7"; fi
    if [ $BARCODE_NUMBER = "07" ]; then REPS="1 2 3"; fi
    if [ $BARCODE_NUMBER = "08" ]; then REPS="1"; fi
    if [ $BARCODE_NUMBER = "09" ]; then REPS="1 2 3 4 5"; fi
    if [ $BARCODE_NUMBER = "10" ]; then REPS="1 2 3 4 5 6 7 8"; fi
    if [ $BARCODE_NUMBER = "11" ]; then REPS="1 2 3 4 5 6 7"; fi
    if [ $BARCODE_NUMBER = "12" ]; then REPS="1 2 3 4 5 6 7 8 9 10 11"; fi

    for REP in $REPS; do
        REP_SIZE=$(awk -v R="$REP" '$1 == R {print $3}' $BARCODE"_illumina_depths" | wc -l | tr -d '\n')
        ILLUMINA_DEPTH=$(awk -v R="$REP" '$1 == R {print $3}' $BARCODE"_illumina_depths" | awk '{total += $1; count++} END {print total/count}' | tr -d '\n')
        NANOPORE_DEPTH=$(awk -v R="$REP" '$1 == R {print $3}' $BARCODE"_nanopore_depths" | awk '{total += $1; count++} END {print total/count}' | tr -d '\n')
        
        if [ $REP = "1" ]; then
            REP_NAME="chromosome"
            ILLUMINA_CHROMOSOME_DEPTH=$ILLUMINA_DEPTH
            NANOPORE_CHROMOSOME_DEPTH=$NANOPORE_DEPTH
        else
            REP_NAME="plasmid"
        fi

        ILLUMINA_NORMALISED_DEPTH=$(echo "scale=3; $ILLUMINA_DEPTH/$ILLUMINA_CHROMOSOME_DEPTH" | bc)
        NANOPORE_NORMALISED_DEPTH=$(echo "scale=3; $NANOPORE_DEPTH/$NANOPORE_CHROMOSOME_DEPTH" | bc)

        printf $BARCODE_NUMBER"\t"$REP_NAME"\t"$REP_SIZE"\t"$ILLUMINA_DEPTH"\t"$NANOPORE_DEPTH"\t"$ILLUMINA_NORMALISED_DEPTH"\t"$NANOPORE_NORMALISED_DEPTH"\n"
        
    done

done;
