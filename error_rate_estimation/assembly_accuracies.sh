
UNICYCLER_ILLUMINA_ASSEMBLIES=08_unicycler_illumina_only_assemblies
CANU_NANOPORE_ASSEMBLIES=09_canu_long_read_only_assemblies
UNICYCLER_NANOPORE_ASSEMBLIES=10_unicycler_long_read_only_assemblies
UNICYCLER_HYBRID_ASSEMBLIES=11_unicycler_hybrid_assemblies

printf "Barcode"
printf "\tCanu accuracy\tCanu distance between errors"
printf "\tUnicycler long-read-only accuracy\tUnicycler long-read-only distance between errors"
printf "\tUnicycler hybrid accuracy\tUnicycler hybrid distance between errors"
printf "\n"

for BARCODE_NUMBER in 01 02 03 04 05 06 07 08 09 10 11 12; do
    BARCODE=barcode$BARCODE_NUMBER

    printf $BARCODE_NUMBER"\t"

    # Prepare the BLAST queries from a pre-RR assembly (big contigs only, ends trimmed off).
    head -n 25 $UNICYCLER_ILLUMINA_ASSEMBLIES/$BARCODE/002_overlaps_removed.gfa | cut -f1-3 | sed 's|S\t|>|' | sed 's|\t|\n|' > blast_queries.fasta
    python3 trim_blast_queries.py
    mv blast_queries2.fasta blast_queries.fasta

    # Build BLAST databases.
    makeblastdb -dbtype nucl -in $CANU_NANOPORE_ASSEMBLIES/$BARCODE/canu.contigs.fasta > /dev/null
    makeblastdb -dbtype nucl -in $UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/assembly.fasta > /dev/null
    makeblastdb -dbtype nucl -in $UNICYCLER_HYBRID_ASSEMBLIES/$BARCODE/assembly.fasta > /dev/null

    # Get the best BLAST hit for each query.
    blastn -db $CANU_NANOPORE_ASSEMBLIES/$BARCODE/canu.contigs.fasta -query blast_queries.fasta -outfmt 6 | sort -nk1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -nk1,1 > canu_blast_hits
    blastn -db $UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/assembly.fasta -query blast_queries.fasta -outfmt 6 | sort -nk1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -nk1,1 > unicycler_long_read_only_blast_hits
    blastn -db $UNICYCLER_HYBRID_ASSEMBLIES/$BARCODE/assembly.fasta -query blast_queries.fasta -outfmt 6 | sort -nk1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -nk1,1 > unicycler_hybrid_blast_hits

    # Get the total error rates.
    python3 get_error_rate.py canu_blast_hits
    printf "\t"
    python3 get_error_rate.py unicycler_long_read_only_blast_hits
    printf "\t"
    python3 get_error_rate.py unicycler_hybrid_blast_hits
    printf "\n"

    rm blast_queries.fasta *_blast_hits

done
