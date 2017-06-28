mkdir -p velvet

for BARCODE_NUMBER in 01 02 03 04 05 06 07 08 09 10 11 12; do
    BARCODE=barcode$BARCODE_NUMBER

    ILLUMINA_1=$(realpath ../03_trimmed_illumina_reads/$BARCODE/*_1.fq.gz)
    ILLUMINA_2=$(realpath ../03_trimmed_illumina_reads/$BARCODE/*_2.fq.gz)

    # ABySS assembly.
    mkdir -p abyss/$BARCODE
    cd abyss/$BARCODE
    abyss-pe -j 20 name=abyss_assembly k=64 graph=gfa in="$ILLUMINA_1 $ILLUMINA_2"
    cd ../..

    # Velvet assembly.
    cd velvet
    ~/velvet_1.2.10/velveth $BARCODE 63 -shortPaired -fastq -separate $ILLUMINA_1 $ILLUMINA_2
    ~/velvet_1.2.10/velvetg $BARCODE
    cd ..

done;
