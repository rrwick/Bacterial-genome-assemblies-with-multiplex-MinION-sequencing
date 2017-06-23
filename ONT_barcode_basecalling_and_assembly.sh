# This script assumes these two directories already exist
RAW_ILLUMINA=illumina_reads/raw               # Should have subdirectories named barcode01, barcode02, etc, with *_1.fastq.gz and *_2.fastq.gz in each.
NANOPORE_FAST5S=nanopore_reads/raw_fast5      # Should have the raw (before basecalling) Nanopore *.fast5 files (can be nested in multipled directories).

# This step is necessary for Illumina-only and hybrid assemblies.
TRIM_ILLUMINA_READS=true                      # Runs Trim Galore to trim adapters from the Illumina reads and do a some conservative quality trimming.

# These steps are necessary for Nanopore and hybrid assemblies.
ALBACORE_BASECALLING=true                     # Runs albacore to make fastq files from fast5 files.
GATHER_UP_NANOPORE_FASTQS=true                # Groups Nanopore reads into a single fastq per barcode bin.
TRIM_NANOPORE_READS=true                      # Runs Porechop to remove adapters from Nanopore reads and split some chimeras.
SUBSAMPLE_NANOPORE_READS=true                 # Runs fastq_to_fastq.py to subsample Nanopore reads using length and (if there are too many reads) quality.

# These steps turn particular assemblies on and off.
ASSEMBLE_ILLUMINA_READS_WITH_UNICYCLER=true   # Do an Illumina-only assembly with Unicycler.
ASSEMBLE_ILLUMINA_READS_WITH_SPADES=true      # Do an Illumina-only assembly with SPAdes.
ASSEMBLE_NANOPORE_READS_WITH_UNICYCLER=true   # Do a Nanopore-only assembly with Unicycler.
ASSEMBLE_NANOPORE_READS_WITH_CANU=true        # Do a Nanopore-only assembly with Canu.
ASSEMBLE_HYBRID_READS_WITH_UNICYCLER=true     # Do a hybrid assembly with Unicycler.
ASSEMBLE_HYBRID_READS_WITH_SPADES=true        # Do a hybrid assembly with SPAdes.

# These steps are necessary for Nanopolish.
ALBACORE_BASECALLING_TO_FAST5=true            # Runs Albacore basecalling, this time saving basecalling into the fast5 files (necessary for Nanopolish).
PREPARE_NANOPOLISH_READS=true                 # Runs Nanopolish extract and gets stuff ready for Nanopolish.

# These steps do Nanopolish on the Nanopore-only assemblies.
NANOPOLISH_CANU_ASSEMBLY=true                 # Run Nanopolish on the Canu assembly.
NANOPOLISH_UNICYCLER_NANOPORE_ASSEMBLY=true   # Run Nanopolish on the Nanopore-only Unicycler assembly.

# This step Pilon-polishes the Canu assembly.
PILON_CANU_ASSEMBLY=true                      # Run Pilon on the Canu assembly.

# The script will make these directories, as necessary.
TRIMMED_ILLUMINA=illumina_reads/trimmed
BASECALLED_NANOPORE_FASTQ=nanopore_reads/basecalling/to_fastq
BASECALLED_NANOPORE_FAST5=nanopore_reads/basecalling/to_fast5
RAW_NANOPORE=nanopore_reads/fastq/1_raw
TRIMMED_NANOPORE=nanopore_reads/fastq/2_trimmed
SUBSAMPLED_NANOPORE=nanopore_reads/fastq/3_subsampled
READS_FOR_NANOPOLISH=nanopore_reads/for_nanopolish
UNICYCLER_ILLUMINA_ASSEMBLIES=assemblies/illumina_only/unicycler
SPADES_ILLUMINA_ASSEMBLIES=assemblies/illumina_only/spades
CANU_NANOPORE_ASSEMBLIES=assemblies/nanopore_only/canu
NANOPOLISHED_CANU_ASSEMBLIES=assemblies/nanopore_only/canu_nanopolish
UNICYCLER_NANOPORE_ASSEMBLIES=assemblies/nanopore_only/unicycler
NANOPOLISHED_UNICYCLER_NANOPORE_ASSEMBLIES=assemblies/nanopore_only/unicycler_nanopolish
UNICYCLER_HYBRID_ASSEMBLIES=assemblies/hybrid/unicycler
SPADES_HYBRID_ASSEMBLIES=assemblies/hybrid/spades
PILON_POLISHED_CANU_ASSEMBLIES=assemblies/hybrid/canu_pilon


# Adjust the thread count as appropriate for the hardware.
THREADS=40
NANOPOLISH_PROCESSES=10
NANOPOLISH_THREADS_PER_PROCESS=4

# Nanopore reads shorter than this will be excluded.
MIN_READ_LENGTH=2000

# If barcode bins have more bases than this, the reads will be subsampled down (using quality scores) to this number of bases.
BASE_LIMIT=500000000

# Canu needs to know the approximate genome size.
GENOME_SIZE=5.5m




# Basecall the Nanopore reads using Albacore (with direct-to-fastq basecalling).
if $ALBACORE_BASECALLING; then
    read_fast5_basecaller.py --input $NANOPORE_FAST5S --recursive --worker_threads $THREADS --save_path $BASECALLED_NANOPORE_FASTQ --barcoding --flowcell FLO-MIN106 --kit SQK-LSK108 --output_format fastq --reads_per_fastq_batch 100000000
fi

# Basecall the Nanopore reads using Albacore, saving them as fast5 files suitable for Nanopolish.
if $ALBACORE_BASECALLING_TO_FAST5; then
    read_fast5_basecaller.py --input $NANOPORE_FAST5S --recursive --worker_threads $THREADS --save_path $BASECALLED_NANOPORE_FAST5 --barcoding --flowcell FLO-MIN106 --kit SQK-LSK108 --output_format fast5
fi

for BARCODE_NUMBER in 01 02 03 04 05 06 07 08 09 10 11 12; do
    BARCODE=barcode$BARCODE_NUMBER

    if $TRIM_ILLUMINA_READS; then
        mkdir -p $TRIMMED_ILLUMINA/$BARCODE
        trim_galore --paired --quality 10 --output_dir $TRIMMED_ILLUMINA/$BARCODE $RAW_ILLUMINA/$BARCODE/*_1.fastq.gz $RAW_ILLUMINA/$BARCODE/*_2.fastq.gz
    fi

    if $GATHER_UP_NANOPORE_FASTQS; then
        mkdir -p $RAW_NANOPORE
        # Some versions of Albacore put newlines between fastq records, so we use grep to only get non-empty lines.
        # This step also serves to group each barcode bin into a single fastq if there are more than one (Albacore makes a separate file for each run ID).
        grep -h . $BASECALLED_NANOPORE_FASTQ/workspace/$BARCODE/*.fastq | gzip > $RAW_NANOPORE/$BARCODE.fastq.gz
    fi

    if $TRIM_NANOPORE_READS; then
        mkdir -p $TRIMMED_NANOPORE
        # Even though Albacore already sorted the reads into barcode bins, I'm running Porechop with barcode binning on.
        # This is so I can exclude reads where Porechop and Albacore disagree on the proper barcode bin.
        porechop -i $RAW_NANOPORE/$BARCODE.fastq.gz -b $TRIMMED_NANOPORE/$BARCODE --threads $THREADS
        zless $TRIMMED_NANOPORE/$BARCODE/BC$BARCODE_NUMBER.fastq.gz | paste - - - - | cut -f 2 | awk '{ print length($0); }' > $TRIMMED_NANOPORE/$BARCODE"_read_lengths"
    fi

    if $SUBSAMPLE_NANOPORE_READS; then
        mkdir -p $SUBSAMPLED_NANOPORE
        fastq_to_fastq.py --min_length $MIN_READ_LENGTH --target_bases $BASE_LIMIT $TRIMMED_NANOPORE/$BARCODE/BC$BARCODE_NUMBER.fastq.gz | gzip > $SUBSAMPLED_NANOPORE/$BARCODE.fastq.gz
    fi

    if $ASSEMBLE_ILLUMINA_READS_WITH_UNICYCLER; then
        mkdir -p $UNICYCLER_ILLUMINA_ASSEMBLIES
        unicycler -1 $TRIMMED_ILLUMINA/$BARCODE/*_1.fq.gz -2 $TRIMMED_ILLUMINA/$BARCODE/*_2.fq.gz -o $UNICYCLER_ILLUMINA_ASSEMBLIES/$BARCODE --threads $THREADS
    fi

    if $ASSEMBLE_ILLUMINA_READS_WITH_SPADES; then
        mkdir -p $SPADES_ILLUMINA_ASSEMBLIES
        spades.py -1 $TRIMMED_ILLUMINA/$BARCODE/*_1.fq.gz -2 $TRIMMED_ILLUMINA/$BARCODE/*_2.fq.gz -o $SPADES_ILLUMINA_ASSEMBLIES/$BARCODE --threads $THREADS --careful
    fi

    if $ASSEMBLE_NANOPORE_READS_WITH_UNICYCLER; then
        mkdir -p $UNICYCLER_NANOPORE_ASSEMBLIES
        unicycler -l $SUBSAMPLED_NANOPORE/$BARCODE.fastq.gz -o $UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE --threads $THREADS
    fi

    if $ASSEMBLE_NANOPORE_READS_WITH_CANU; then
        mkdir -p $CANU_NANOPORE_ASSEMBLIES
        canu -p canu -d $CANU_NANOPORE_ASSEMBLIES/$BARCODE genomeSize=$GENOME_SIZE -nanopore-raw $SUBSAMPLED_NANOPORE/$BARCODE.fastq.gz
    fi

    if $ASSEMBLE_HYBRID_READS_WITH_UNICYCLER; then
        mkdir -p $UNICYCLER_HYBRID_ASSEMBLIES
        unicycler -1 $TRIMMED_ILLUMINA/$BARCODE/*_1.fq.gz -2 $TRIMMED_ILLUMINA/$BARCODE/*_2.fq.gz -l $SUBSAMPLED_NANOPORE/$BARCODE.fastq.gz -o $UNICYCLER_HYBRID_ASSEMBLIES/$BARCODE --threads $THREADS
    fi

    if $ASSEMBLE_HYBRID_READS_WITH_SPADES; then
        mkdir -p $SPADES_HYBRID_ASSEMBLIES
        spades.py -1 $TRIMMED_ILLUMINA/$BARCODE/*_1.fq.gz -2 $TRIMMED_ILLUMINA/$BARCODE/*_2.fq.gz --nanopore $SUBSAMPLED_NANOPORE/$BARCODE.fastq.gz -o $SPADES_HYBRID_ASSEMBLIES/$BARCODE --threads $THREADS --careful
    fi

    if $PREPARE_NANOPOLISH_READS; then
        mkdir -p $READS_FOR_NANOPOLISH
        nanopolish extract --recurse --type template $BASECALLED_NANOPORE_FAST5/workspace/$BARCODE/ > $READS_FOR_NANOPOLISH/$BARCODE"_before_filter.fasta"
        python3 nanopolish_read_filter.py $READS_FOR_NANOPOLISH/$BARCODE"_before_filter.fasta" $TRIMMED_NANOPORE/$BARCODE/BC$BARCODE_NUMBER.fastq.gz > $READS_FOR_NANOPOLISH/$BARCODE.fasta
    fi

    if $NANOPOLISH_UNICYCLER_NANOPORE_ASSEMBLY; then
        mkdir -p $NANOPOLISHED_UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE
        bwa index $UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/assembly.fasta
        bwa mem -x ont2d -t $THREADS $UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/assembly.fasta $READS_FOR_NANOPOLISH/$BARCODE.fasta | samtools sort -o $NANOPOLISHED_UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/reads.sorted.bam -T reads.tmp -
        samtools index $NANOPOLISHED_UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/reads.sorted.bam
        python ~/nanopolish/scripts/nanopolish_makerange.py $UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/assembly.fasta | parallel --results $NANOPOLISHED_UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/nanopolish.results -P $NANOPOLISH_PROCESSES nanopolish variants --consensus $NANOPOLISHED_UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/polished.{1}.fa -w {1} -r $READS_FOR_NANOPOLISH/$BARCODE.fasta -b $NANOPOLISHED_UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/reads.sorted.bam -g $UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/assembly.fasta -t $NANOPOLISH_THREADS_PER_PROCESS --min-candidate-frequency 0.1
        python ~/nanopolish/scripts/nanopolish_merge.py $NANOPOLISHED_UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/polished.*.fa > $NANOPOLISHED_UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/polished_genome.fa
        rm $NANOPOLISHED_UNICYCLER_NANOPORE_ASSEMBLIES/$BARCODE/polished.*.fa
    fi

    if $NANOPOLISH_CANU_ASSEMBLY; then
        mkdir -p $NANOPOLISHED_CANU_ASSEMBLIES/$BARCODE
        bwa index $CANU_NANOPORE_ASSEMBLIES/$BARCODE/canu.contigs.fasta
        bwa mem -x ont2d -t $THREADS $CANU_NANOPORE_ASSEMBLIES/$BARCODE/canu.contigs.fasta $READS_FOR_NANOPOLISH/$BARCODE.fasta | samtools sort -o $NANOPOLISHED_CANU_ASSEMBLIES/$BARCODE/reads.sorted.bam -T reads.tmp -
        samtools index $NANOPOLISHED_CANU_ASSEMBLIES/$BARCODE/reads.sorted.bam
        python ~/nanopolish/scripts/nanopolish_makerange.py $CANU_NANOPORE_ASSEMBLIES/$BARCODE/canu.contigs.fasta | parallel --results $NANOPOLISHED_CANU_ASSEMBLIES/$BARCODE/nanopolish.results -P $NANOPOLISH_PROCESSES nanopolish variants --consensus $NANOPOLISHED_CANU_ASSEMBLIES/$BARCODE/polished.{1}.fa -w {1} -r $READS_FOR_NANOPOLISH/$BARCODE.fasta -b $NANOPOLISHED_CANU_ASSEMBLIES/$BARCODE/reads.sorted.bam -g $CANU_NANOPORE_ASSEMBLIES/$BARCODE/canu.contigs.fasta -t $NANOPOLISH_THREADS_PER_PROCESS --min-candidate-frequency 0.1
        python ~/nanopolish/scripts/nanopolish_merge.py $NANOPOLISHED_CANU_ASSEMBLIES/$BARCODE/polished.*.fa > $NANOPOLISHED_CANU_ASSEMBLIES/$BARCODE/polished_genome.fa
        rm $NANOPOLISHED_CANU_ASSEMBLIES/$BARCODE/polished.*.fa
    fi

    if $PILON_CANU_ASSEMBLY; then
        mkdir -p $PILON_POLISHED_CANU_ASSEMBLIES/$BARCODE
        INPUT_FASTA=$CANU_NANOPORE_ASSEMBLIES/$BARCODE/canu.contigs.fasta

        for PILON_ROUND in 1 2 3 4 5; do
            bowtie2-build $INPUT_FASTA $INPUT_FASTA
            BAM=$PILON_POLISHED_CANU_ASSEMBLIES/$BARCODE/alignments_"$PILON_ROUND".bam
            bowtie2 --local --very-sensitive-local --threads $THREADS -I 0 -X 2000 -x $INPUT_FASTA -1 $TRIMMED_ILLUMINA/$BARCODE/*_1.fq.gz -2 $TRIMMED_ILLUMINA/$BARCODE/*_2.fq.gz | samtools sort -o $BAM -T reads.tmp -
            samtools index $BAM
            java -jar ~/pilon-1.22.jar --genome $INPUT_FASTA --frags $BAM --changes --output pilon_round_"$PILON_ROUND" --outdir $PILON_POLISHED_CANU_ASSEMBLIES/$BARCODE --fix all
            INPUT_FASTA=$PILON_POLISHED_CANU_ASSEMBLIES/$BARCODE/pilon_round_"$PILON_ROUND".fasta
        done
    fi

done
