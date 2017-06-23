TRIMMED_ILLUMINA=../03_trimmed_illumina_reads
TRIMMED_NANOPORE=../07_trimmed_nanopore_reads

THREADS=20

for BARCODE_NUMBER in 01 02 03 04 05 06 07 08 09 10 11 12; do
    BARCODE=barcode$BARCODE_NUMBER

    FASTA="NB"$BARCODE_NUMBER"*.fasta"

    # Align Illumina reads.
    bowtie2-build $FASTA $FASTA
    bowtie2 --very-sensitive-local --threads $THREADS -I 0 -X 2000 -x $FASTA -1 $TRIMMED_ILLUMINA/$BARCODE/*_1.fq.gz -2 $TRIMMED_ILLUMINA/$BARCODE/*_2.fq.gz | samtools sort -o $BARCODE"_illumina.bam" -
    samtools index $BARCODE"_illumina.bam"

    # Align Nanopore reads.
    bwa index $FASTA
    bwa mem -x ont2d -t $THREADS $FASTA $TRIMMED_NANOPORE/$BARCODE/BC$BARCODE_NUMBER.fastq.gz | samtools sort -o $BARCODE"_nanopore.bam" -
    samtools index $BARCODE"_nanopore.bam"

    # Put depths (per base) into files.
    samtools depth -a $BARCODE"_illumina.bam" > $BARCODE"_illumina_depths"
    samtools depth -a $BARCODE"_nanopore.bam" > $BARCODE"_nanopore_depths"

    # Clean up.
    rm *.bam
    rm *.bam.bai
    rm *.bt2
    rm *.amb
    rm *.ann
    rm *.bwt
    rm *.pac
    rm *.sa
done
