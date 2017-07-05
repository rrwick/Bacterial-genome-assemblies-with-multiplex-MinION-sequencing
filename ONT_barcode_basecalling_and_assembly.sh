# This script assumes these two directories already exist
raw_illumina=illumina_reads/raw                # Should have subdirectories named barcode01, barcode02, etc, with *_1.fastq.gz and *_2.fastq.gz in each.
nanopore_fast5s=nanopore_reads/raw_fast5       # Should have the raw (before basecalling) Nanopore *.fast5 files (can be nested in multipled directories).

# This step is necessary for Illumina-only and hybrid assemblies.
trim_illumina_reads=true                       # Runs Trim Galore to trim adapters from the Illumina reads and do a some conservative quality trimming.

# These steps are necessary for Nanopore and hybrid assemblies.
albacore_basecalling=true                      # Runs albacore to make fastq files from fast5 files.
gather_up_nanopore_fastqs=true                 # Groups Nanopore reads into a single fastq per barcode bin.
trim_nanopore_reads=true                       # Runs Porechop to remove adapters from Nanopore reads and split some chimeras.
subsample_nanopore_reads=true                  # Runs fastq_to_fastq.py to subsample Nanopore reads using length and (if there are too many reads) quality.

# These steps turn particular assemblies on and off.
assemble_illumina_reads_with_unicycler=true    # Do an Illumina-only assembly with Unicycler.
assemble_illumina_reads_with_spades=true       # Do an Illumina-only assembly with SPAdes.
assemble_nanopore_reads_with_unicycler=true    # Do a Nanopore-only assembly with Unicycler.
assemble_nanopore_reads_with_canu=true         # Do a Nanopore-only assembly with Canu.
assemble_hybrid_reads_with_unicycler=true      # Do a hybrid assembly with Unicycler.
assemble_hybrid_reads_with_spades=true         # Do a hybrid assembly with SPAdes.

# These steps are necessary for Nanopolish.
albacore_basecalling_to_fast5=true             # Runs Albacore basecalling, this time saving basecalling into the fast5 files (necessary for Nanopolish).
prepare_nanopolish_reads=true                  # Runs Nanopolish extract and gets stuff ready for Nanopolish.

# These steps do Nanopolish on the Nanopore-only assemblies.
nanopolish_canu_assembly=true                  # Run Nanopolish on the Canu assembly.
nanopolish_canu_assembly_2=true                # Run a second round of Nanopolish on the Canu assembly.
nanopolish_unicycler_nanopore_assembly=true    # Run Nanopolish on the Nanopore-only Unicycler assembly.
nanopolish_unicycler_nanopore_assembly_2=true  # Run a second round of Nanopolish on the Nanopore-only Unicycler assembly.

# This step Pilon-polishes the Canu assembly.
pilon_canu_assembly=true                       # Run Pilon on the Canu assembly.

# The script will make these directories, as necessary.
trimmed_illumina=illumina_reads/trimmed
basecalled_nanopore_fastq=nanopore_reads/basecalling/to_fastq
basecalled_nanopore_fast5=nanopore_reads/basecalling/to_fast5
raw_nanopore=nanopore_reads/fastq/1_raw
trimmed_nanopore=nanopore_reads/fastq/2_trimmed
subsampled_nanopore=nanopore_reads/fastq/3_subsampled
reads_for_nanopolish=nanopore_reads/for_nanopolish
unicycler_illumina_assemblies=assemblies/illumina-only/unicycler
spades_illumina_assemblies=assemblies/illumina-only/spades
canu_nanopore_assemblies=assemblies/nanopore-only/canu
nanopolished_canu_assemblies=assemblies/nanopore-only/canu_nanopolish
nanopolished_canu_assemblies_2=assemblies/nanopore-only/canu_nanopolish_2
unicycler_nanopore_assemblies=assemblies/nanopore-only/unicycler
nanopolished_unicycler_nanopore_assemblies=assemblies/nanopore-only/unicycler_nanopolish
nanopolished_unicycler_nanopore_assemblies_2=assemblies/nanopore-only/unicycler_nanopolish_2
unicycler_hybrid_assemblies=assemblies/hybrid/unicycler
spades_hybrid_assemblies=assemblies/hybrid/spades
pilon_polished_canu_assemblies=assemblies/hybrid/canu_pilon


# Adjust the thread count as appropriate for the hardware.
threads=40
nanopolish_processes=10
nanopolish_threads_per_process=4

# Nanopore reads shorter than this will be excluded.
min_read_length=2000

# If barcode bins have more bases than this, the reads will be subsampled down (using quality scores) to this number of bases.
base_limit=500000000

# Canu needs to know the approximate genome size.
genome_size=5.5m




# Basecall the Nanopore reads using Albacore (with direct-to-fastq basecalling).
if $albacore_basecalling; then
    read_fast5_basecaller.py --input $nanopore_fast5s --recursive --worker_threads $threads --save_path $basecalled_nanopore_fastq --barcoding --flowcell FLO-MIN106 --kit SQK-LSK108 --output_format fastq --reads_per_fastq_batch 100000000
fi

# Basecall the Nanopore reads using Albacore, saving them as fast5 files suitable for Nanopolish.
if $albacore_basecalling_to_fast5; then
    read_fast5_basecaller.py --input $nanopore_fast5s --recursive --worker_threads $threads --save_path $basecalled_nanopore_fast5 --barcoding --flowcell FLO-MIN106 --kit SQK-LSK108 --output_format fast5
fi

for barcode_number in 01 02 03 04 05 06 07 08 09 10 11 12; do
    barcode=barcode$barcode_number

    if $trim_illumina_reads; then
        mkdir -p $trimmed_illumina/$barcode
        trim_galore --paired --quality 10 --output_dir $trimmed_illumina/$barcode $raw_illumina/$barcode/*_1.fastq.gz $raw_illumina/$barcode/*_2.fastq.gz
    fi

    if $gather_up_nanopore_fastqs; then
        mkdir -p $raw_nanopore
        # Some versions of Albacore put newlines between fastq records, so we use grep to only get non-empty lines.
        # This step also serves to group each barcode bin into a single fastq if there are more than one (Albacore makes a separate file for each run ID).
        grep -h . $basecalled_nanopore_fastq/workspace/$barcode/*.fastq | gzip > $raw_nanopore/$barcode.fastq.gz
    fi

    if $trim_nanopore_reads; then
        mkdir -p $trimmed_nanopore
        # Even though Albacore already sorted the reads into barcode bins, I'm running Porechop with barcode binning on.
        # This is so I can exclude reads where Porechop and Albacore disagree on the proper barcode bin.
        porechop -i $raw_nanopore/$barcode.fastq.gz -b $trimmed_nanopore/$barcode --threads $threads
        zless $trimmed_nanopore/$barcode/BC$barcode_number.fastq.gz | paste - - - - | cut -f 2 | awk '{ print length($0); }' > $trimmed_nanopore/$barcode"_read_lengths"
    fi

    if $subsample_nanopore_reads; then
        mkdir -p $subsampled_nanopore
        fastq_to_fastq.py --min_length $min_read_length --target_bases $base_limit $trimmed_nanopore/$barcode/BC$barcode_number.fastq.gz | gzip > $subsampled_nanopore/$barcode.fastq.gz
    fi

    if $assemble_illumina_reads_with_unicycler; then
        mkdir -p $unicycler_illumina_assemblies
        unicycler -1 $trimmed_illumina/$barcode/*_1.fq.gz -2 $trimmed_illumina/$barcode/*_2.fq.gz -o $unicycler_illumina_assemblies/$barcode --threads $threads
    fi

    if $assemble_illumina_reads_with_spades; then
        mkdir -p $spades_illumina_assemblies
        spades.py -1 $trimmed_illumina/$barcode/*_1.fq.gz -2 $trimmed_illumina/$barcode/*_2.fq.gz -o $spades_illumina_assemblies/$barcode --threads $threads --careful
    fi

    if $assemble_nanopore_reads_with_unicycler; then
        mkdir -p $unicycler_nanopore_assemblies
        unicycler -l $subsampled_nanopore/$barcode.fastq.gz -o $unicycler_nanopore_assemblies/$barcode --threads $threads
    fi

    if $assemble_nanopore_reads_with_canu; then
        mkdir -p $canu_nanopore_assemblies
        canu -p canu -d $canu_nanopore_assemblies/$barcode genomeSize=$genome_size -nanopore-raw $subsampled_nanopore/$barcode.fastq.gz
    fi

    if $assemble_hybrid_reads_with_unicycler; then
        mkdir -p $unicycler_hybrid_assemblies
        unicycler -1 $trimmed_illumina/$barcode/*_1.fq.gz -2 $trimmed_illumina/$barcode/*_2.fq.gz -l $subsampled_nanopore/$barcode.fastq.gz -o $unicycler_hybrid_assemblies/$barcode --threads $threads
    fi

    if $assemble_hybrid_reads_with_spades; then
        mkdir -p $spades_hybrid_assemblies
        spades.py -1 $trimmed_illumina/$barcode/*_1.fq.gz -2 $trimmed_illumina/$barcode/*_2.fq.gz --nanopore $subsampled_nanopore/$barcode.fastq.gz -o $spades_hybrid_assemblies/$barcode --threads $threads --careful
    fi

    if $prepare_nanopolish_reads; then
        mkdir -p $reads_for_nanopolish
        nanopolish extract --recurse --type template $basecalled_nanopore_fast5/workspace/$barcode/ > $reads_for_nanopolish/$barcode"_before_filter.fasta"
        python3 nanopolish_read_filter.py $reads_for_nanopolish/$barcode"_before_filter.fasta" $trimmed_nanopore/$barcode/BC$barcode_number.fastq.gz > $reads_for_nanopolish/$barcode.fasta
    fi

    if $nanopolish_unicycler_nanopore_assembly; then
        mkdir -p $nanopolished_unicycler_nanopore_assemblies/$barcode
        bwa index $unicycler_nanopore_assemblies/$barcode/assembly.fasta
        bwa mem -x ont2d -t $threads $unicycler_nanopore_assemblies/$barcode/assembly.fasta $reads_for_nanopolish/$barcode.fasta | samtools sort -o $nanopolished_unicycler_nanopore_assemblies/$barcode/reads.sorted.bam -T reads.tmp -
        samtools index $nanopolished_unicycler_nanopore_assemblies/$barcode/reads.sorted.bam
        python ~/nanopolish/scripts/nanopolish_makerange.py $unicycler_nanopore_assemblies/$barcode/assembly.fasta | parallel --results $nanopolished_unicycler_nanopore_assemblies/$barcode/nanopolish.results -P $nanopolish_processes nanopolish variants --consensus $nanopolished_unicycler_nanopore_assemblies/$barcode/polished.{1}.fa -w {1} -r $reads_for_nanopolish/$barcode.fasta -b $nanopolished_unicycler_nanopore_assemblies/$barcode/reads.sorted.bam -g $unicycler_nanopore_assemblies/$barcode/assembly.fasta -t $nanopolish_threads_per_process --min-candidate-frequency 0.1
        python ~/nanopolish/scripts/nanopolish_merge.py $nanopolished_unicycler_nanopore_assemblies/$barcode/polished.*.fa > $nanopolished_unicycler_nanopore_assemblies/$barcode/polished_genome.fa
        rm $nanopolished_unicycler_nanopore_assemblies/$barcode/polished.*.fa
    fi

    if $nanopolish_unicycler_nanopore_assembly_2; then
        mkdir -p $nanopolished_unicycler_nanopore_assemblies_2/$barcode
        bwa index $nanopolished_unicycler_nanopore_assemblies/$barcode/polished_genome.fa
        bwa mem -x ont2d -t $threads $nanopolished_unicycler_nanopore_assemblies/$barcode/polished_genome.fa $reads_for_nanopolish/$barcode.fasta | samtools sort -o $nanopolished_unicycler_nanopore_assemblies_2/$barcode/reads.sorted.bam -T reads.tmp -
        samtools index $nanopolished_unicycler_nanopore_assemblies_2/$barcode/reads.sorted.bam
        python ~/nanopolish/scripts/nanopolish_makerange.py $nanopolished_unicycler_nanopore_assemblies/$barcode/polished_genome.fa | parallel --results $nanopolished_unicycler_nanopore_assemblies_2/$barcode/nanopolish.results -P $nanopolish_processes nanopolish variants --consensus $nanopolished_unicycler_nanopore_assemblies_2/$barcode/polished.{1}.fa -w {1} -r $reads_for_nanopolish/$barcode.fasta -b $nanopolished_unicycler_nanopore_assemblies_2/$barcode/reads.sorted.bam -g $nanopolished_unicycler_nanopore_assemblies/$barcode/polished_genome.fa -t $nanopolish_threads_per_process --min-candidate-frequency 0.1
        python ~/nanopolish/scripts/nanopolish_merge.py $nanopolished_unicycler_nanopore_assemblies_2/$barcode/polished.*.fa > $nanopolished_unicycler_nanopore_assemblies_2/$barcode/polished_genome.fa
        rm $nanopolished_unicycler_nanopore_assemblies_2/$barcode/polished.*.fa
    fi

    if $nanopolish_canu_assembly; then
        mkdir -p $nanopolished_canu_assemblies/$barcode
        bwa index $canu_nanopore_assemblies/$barcode/canu.contigs.fasta
        bwa mem -x ont2d -t $threads $canu_nanopore_assemblies/$barcode/canu.contigs.fasta $reads_for_nanopolish/$barcode.fasta | samtools sort -o $nanopolished_canu_assemblies/$barcode/reads.sorted.bam -T reads.tmp -
        samtools index $nanopolished_canu_assemblies/$barcode/reads.sorted.bam
        python ~/nanopolish/scripts/nanopolish_makerange.py $canu_nanopore_assemblies/$barcode/canu.contigs.fasta | parallel --results $nanopolished_canu_assemblies/$barcode/nanopolish.results -P $nanopolish_processes nanopolish variants --consensus $nanopolished_canu_assemblies/$barcode/polished.{1}.fa -w {1} -r $reads_for_nanopolish/$barcode.fasta -b $nanopolished_canu_assemblies/$barcode/reads.sorted.bam -g $canu_nanopore_assemblies/$barcode/canu.contigs.fasta -t $nanopolish_threads_per_process --min-candidate-frequency 0.1
        python ~/nanopolish/scripts/nanopolish_merge.py $nanopolished_canu_assemblies/$barcode/polished.*.fa > $nanopolished_canu_assemblies/$barcode/polished_genome.fa
        rm $nanopolished_canu_assemblies/$barcode/polished.*.fa
    fi

    if $nanopolish_canu_assembly_2; then
        mkdir -p $nanopolished_canu_assemblies_2/$barcode
        bwa index $nanopolished_canu_assemblies/$barcode/polished_genome.fa
        bwa mem -x ont2d -t $threads $nanopolished_canu_assemblies/$barcode/polished_genome.fa $reads_for_nanopolish/$barcode.fasta | samtools sort -o $nanopolished_canu_assemblies_2/$barcode/reads.sorted.bam -T reads.tmp -
        samtools index $nanopolished_canu_assemblies_2/$barcode/reads.sorted.bam
        python ~/nanopolish/scripts/nanopolish_makerange.py $nanopolished_canu_assemblies/$barcode/polished_genome.fa | parallel --results $nanopolished_canu_assemblies_2/$barcode/nanopolish.results -P $nanopolish_processes nanopolish variants --consensus $nanopolished_canu_assemblies_2/$barcode/polished.{1}.fa -w {1} -r $reads_for_nanopolish/$barcode.fasta -b $nanopolished_canu_assemblies_2/$barcode/reads.sorted.bam -g $nanopolished_canu_assemblies/$barcode/polished_genome.fa -t $nanopolish_threads_per_process --min-candidate-frequency 0.1
        python ~/nanopolish/scripts/nanopolish_merge.py $nanopolished_canu_assemblies_2/$barcode/polished.*.fa > $nanopolished_canu_assemblies_2/$barcode/polished_genome.fa
        rm $nanopolished_canu_assemblies_2/$barcode/polished.*.fa
    fi

    if $pilon_canu_assembly; then
        mkdir -p $pilon_polished_canu_assemblies/$barcode
        input_fasta=$canu_nanopore_assemblies/$barcode/canu.contigs.fasta

        for pilon_round in 1 2 3 4 5; do
            bowtie2-build $input_fasta $input_fasta
            BAM=$pilon_polished_canu_assemblies/$barcode/alignments_"$pilon_round".bam
            bowtie2 --local --very-sensitive-local --threads $threads -I 0 -X 2000 -x $input_fasta -1 $trimmed_illumina/$barcode/*_1.fq.gz -2 $trimmed_illumina/$barcode/*_2.fq.gz | samtools sort -o $BAM -T reads.tmp -
            samtools index $BAM
            java -jar ~/pilon-1.22.jar --genome $input_fasta --frags $BAM --changes --output pilon_round_"$pilon_round" --outdir $pilon_polished_canu_assemblies/$barcode --fix all
            input_fasta=$pilon_polished_canu_assemblies/$barcode/pilon_round_"$pilon_round".fasta
        done
    fi

done
