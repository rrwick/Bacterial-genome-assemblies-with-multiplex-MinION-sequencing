<p align="center"><img src="logo.png" alt="Bacterial genome assemblies with multiplex-MinION sequencing"></p>

# Completing bacterial genome assemblies with multiplex MinION sequencing

This repository contains supplemental data and code for our paper: [Completing bacterial genome assemblies with multiplex MinION sequencing](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000132).

Here you will find scripts used to generate our data and figures, links to the reads and assemblies, and summaries of our results. We have also included data from other assemblers not mentioned in the paper to facilitate comparison. If you have different assembly methods that you would like to share, we are happy to include them here. Please do a GitHub pull request with your results or create an [issue](https://github.com/rrwick/Bacterial-genome-assemblies-with-multiplex-MinION-sequencing/issues) on this repo.



## Links to read data

* [Raw Illumina reads (before trimming)](https://figshare.com/articles/Raw_Illumina_reads/5170816)
* [Trimmed Illumina reads](https://figshare.com/articles/Trimmed_Illumina_reads/5170831)
* [Basecalled ONT reads (straight out of Albacore)](https://figshare.com/articles/Basecalled_ONT_reads/5170843)
* [Trimmed ONT reads](https://figshare.com/articles/Trimmed_ONT_reads/5170852)
* [Subsampled ONT reads](https://figshare.com/articles/Subsampled_ONT_reads/5171491)

I did not put the ONT fast5 files on figshare due to their size (157 GB before basecalling and 1.5TB after basecalling). If you are interested in these, please contact me and we can try to work something out.



## Basecalling and assembly

The [ONT_barcode_basecalling_and_assembly.sh](ONT_barcode_basecalling_and_assembly.sh) script carries out the following steps:
* Trimming Illumina reads
* Basecalling ONT reads
* Trimming and subsampling ONT reads
* Assembling Illumina-only read sets (with SPAdes and Unicycler)
* Assembling ONT-only read sets (with Canu and Unicycler)
* Assembling hybrid read sets (with SPAdes, Canu+Pilon and Unicycler)
* Polishing ONT-only assemblies with Nanopolish

Each of these steps can be turned on/off using the variables at the top of the script. Details for some of the steps are described below.


#### Software versions used

* [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/): v0.4.4
* Albacore: v1.1.2
* [Porechop](https://github.com/rrwick/Porechop): v0.2.1
* [Unicycler](https://github.com/rrwick/Unicycler): v0.4.0
* [Canu](http://canu.readthedocs.io/en/latest/): [snapshot v1.5 +54 changes](https://github.com/marbl/canu/tree/f356c2c3f2eb37b53c4e7bf11e927e3fdff4d747)
* [SPAdes](http://cab.spbu.ru/software/spades/): v3.10.1
* [Pilon](https://github.com/broadinstitute/pilon): v1.22
* [Nanopolish](https://github.com/jts/nanopolish): v0.7.0


#### ONT read processing

When basecalling ONT reads using Albacore (ONT's command-line basecaller), we used the `--barcoding` option to sort the reads into barcode bins. We then ran [Porechop](https://github.com/rrwick/Porechop) on each bin to remove adapter sequences and discard chimeric reads.

When running Porechop we used its barcode binning as well so we could keep only reads where both Albacore and Porechop agreed on the barcode bin. For example, Albacore put 95064 reads in the bin for barcode 1. Of these, Porechop put 90919 in the barcode 1 bin, 118 reads into bins for other barcodes, 3941 reads into no bin and 86 reads were discarded as chimeras. By using only the 90919 reads where Albacore and Porechop agree, minimised barcode cross-contamination.

All reads shorter than 2 kbp were discarded for each sample – due to the long read N50s this was a very small proportion of the reads. For samples which still had more than 500 Mbp of reads, we subsampled the read set down to 500 Mbp. This was done using read quality – specifically the reads' minimum qscore over a sliding window. This means that the discarded reads were the one which had the lowest quality regions, as indicated by their qscores. This was done with the `fastq_to_fastq.py` script in [this repo](https://github.com/rrwick/Fast5-to-Fastq).


#### Illumina-only assembly

We used the [trimmed Illumina reads](https://figshare.com/articles/Trimmed_Illumina_reads/5170831) as input for SPAdes and Unicycler:
* SPAdes: `spades.py -1 short_1.fastq.gz -2 short_2.fastq.gz -o out_dir --careful`
* Unicycler: `unicycler -1 short_1.fastq.gz -2 short_2.fastq.gz -o out_dir`

For SPAdes, the `contigs.fasta` file was taken as the final assembly.


#### ONT-only assembly

The [subsampled ONT reads](https://figshare.com/articles/Subsampled_ONT_reads/5171491) were used as input for Canu and Unicycler:
* Canu: `canu -p canu -d out_dir genomeSize=5.5m -nanopore-raw long.fastq.gz`
* Unicycler: `unicycler -l long.fastq.gz -o out_dir`


#### Hybrid assembly

The [trimmed Illumina reads](https://figshare.com/articles/Trimmed_Illumina_reads/5170831) and [subsampled ONT reads](https://figshare.com/articles/Subsampled_ONT_reads/5171491) were used as input for hybrid assemblies:
* SPAdes: `spades.py -1 short_1.fastq.gz -2 short_2.fastq.gz --nanopore long.fastq.gz -o out_dir --careful`
* Canu+Pilon:
  * `canu -p canu -d out_dir genomeSize=5.5m -nanopore-raw long.fastq.gz`
  * Then fives rounds of Pilon:
    * `bowtie2 --local --very-sensitive-local -I 0 -X 2000 -x before_polish.fasta -1 short_1.fastq.gz -2 short_2.fastq.gz | samtools sort -o alignments.bam -T reads.tmp -; samtools index alignments.bam`
    * `java -jar ~/pilon-1.22.jar --genome before_polish.fasta --frags alignments.bam --changes --output after_polish --outdir out_dir --fix all`
* Unicycler: `unicycler -1 short_1.fastq.gz -2 short_2.fastq.gz -l long.fastq.gz -o out_dir`


#### Polishing with Nanopolish

We used Nanopolish to get the most accurate possible ONT-only assemblies:
* `python nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 nanopolish variants --consensus polished.{1}.fa -w {1} -r reads.fa -b reads.sorted.bam -g draft.fa -t 4 --min-candidate-frequency 0.1`
* `python nanopolish_merge.py polished.*.fa > polished_genome.fa`

For this step we used the full set of [trimmed ONT reads](https://figshare.com/articles/Trimmed_ONT_reads/5170852) (before the read sets were subsampled). After using `nanopolish extract` to produce a fasta file from Albacore's output directory, we used [this script](nanopolish_read_filter.py) to exclude reads where Porechop disagreed on the bin.

We tried a second round of Nanopolish but found that it did not significantly change the results, so here we only report results from a single round of Nanopolish.



## Error rate estimation

The files in the [error_rate_estimation](error_rate_estimation) directory were used to get error rate estimates for assemblies. We...
* assembled each sample separately using both [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) and [ABySS](https://github.com/bcgsc/abyss) ([assemblies_for_error_rate_estimation.sh](error_rate_estimation/assemblies_for_error_rate_estimation.sh)). These were chosen as independent assemblers from the ones we are assessing.
* used [MUMmer](http://mummer.sourceforge.net/) to extract only large (10+ kbp) contigs where Velvet and ABySS were in perfect agreement ([make_sample_reference_fasta.sh](error_rate_estimation/make_sample_reference_fasta.sh)). The reference contigs made in this process can be downloaded [here](https://figshare.com/articles/Reference_contigs_for_error_rate_estimation/5172487).
* BLASTed the assemblies using these contigs, keeping only the best hit per contig and averaged the error rate over the hits ([assembly_accuracies.sh](error_rate_estimation/assembly_accuracies.sh)).

By using only large (10+ kbp) contigs, this method only covers non-repetitive DNA. Error rates in repetitive regions will possibly be higher.



## Depth per replicon

The files in the [depth_per_replicon](depth_per_replicon) directory were used to generate Figure S4 which shows the read depth for each plasmid, relative to the chromosomal depth, for both Illumina and ONT reads. It shows that small plasmids are very underrepresented in ONT reads.



## ONT-only error rates

The files in the [nanopore_only_error_rates](nanopore_only_error_rates) were used to generate Figure S3 which shows Canu error rates (before and after Nanopolish) against ONT read depth.



## Result table

The [results.xlsx](results.xlsx) file contains statistics on each read set and assembly. The summaries below were taken from this table.



## Results: Illumina-only assemblies

| Assembler | Mean contigs | Mean N50 | Complete large plasmids | Complete small plasmids | Estimated error rate |
| :-------: | -----------: | -------: | ----------------------: | ----------------------: | -------------------: |
| SPAdes    |        379.1 |  218,479 |                     n/a |                     n/a |              0.0001% |
| Unicycler |        191.8 |  293,648 |                  2 / 28 |                 12 / 29 |              0.0000% |

#### Links to assemblies

* [SPAdes v3.10.1](https://figshare.com/articles/SPAdes_v3_10_1_assemblies_Illumina-only_/5165836)
* [Unicycler v0.4.0](https://figshare.com/articles/Unicycler_v0_4_0_assemblies_Illumina-only_/5170744)

#### Metrics

* Mean contigs: the number of contigs in the assembly, averaged over all 12 samples (fewer is better).
* Mean N50: the assembly N50, averaged over all 12 samples (bigger is better).
    * Illumina-only assemblies of bacterial genomes are almost never complete, so this value is much less than the chromosome size.
* Complete large/small plasmids: how many plasmids completely assembled into a single contig, totalled over all 12 samples.
    * For Unicycler, we defined "completely assembled" as the plasmid being circularised in the graph (i.e. has a link joining its start and end).
    * The complete plasmid counts aren't available for SPAdes, as it outputs its final assembly in contig form (not as a graph) so there's no simple _de novo_ way to tell if a contig is a complete replicon.
    * Large plasmids were defined as over 10 kbp (though there were no plasmids between 7 and 60 kbp).
    * The total number of plasmids in the 12 isolates (28 large, 29 small) was determined from the manually-completed assemblies.
* Estimated error rate: estimate of the assembly's base-level error rate – method described above

#### Conclusions

Overall, Unicycler and SPAdes perform similarly when assembling the Illumina reads – not surprising, since Unicycler uses SPAdes to assemble Illumina reads. Unicycler achieves slightly better values because it uses a wider k-mer range than SPAdes does by default. Experimenting with larger values for SPAdes' `-k` option would probably give results close to Unicycler's.

Both Unicycler and SPAdes had extremely low error rates. These means that their assemblies are in very good agreement with ABySS and Velvet for the non-repetitive sequences assessed. This agreement between different assemblers supports our assumption that Illumina-only assemblies have near-perfect base-level accuracy (at least for non-repetitive sequence).

The SPAdes mean contig count is greatly inflated by sample INF163 which has some low-depth contamination. The SPAdes assembly has many contigs from this contamination, but they are filtered out in the Unicycler assembly. Excluding that sample, the mean contig count for SPAdes is 213.7, much closer to Unicycler's value.

As expected for short reads, neither assembler was very good at completing large plasmids, as they usually contained shared sequence with other replicons. Even though exact completed-plasmid counts aren't available for SPAdes, it seemed to perform similarly to Unicycler on small plasmids – assembling them into single contigs when they only contain unique sequence, assembling them into incomplete contigs when they share sequence with each other.



## Results: ONT-only assemblies

| Assembler | Mean N50  | Complete chromosomes | Complete large plasmids | Complete small plasmids | Estimated error rate (pre-Nanopolish) | Estimated error rate (post-Nanopolish) |
| :-------: | --------: | -------------------: | ----------------------: | ----------------------: | ------------------------------------: | -------------------------------------: |
| Canu      | 4,784,356 |               4 / 12 |                 23 / 28 |                  0 / 29 |                               1.2219% |                                0.6681% |
| Unicycler | 4,965,584 |               7 / 12 |                 27 / 28 |                  5 / 29 |                               1.0164% |                                0.6164% |

#### Links to assemblies

* Canu v1.5 (f356c2c): [before Nanopolish](https://figshare.com/articles/Canu_v1_5_f356c2c_assemblies_ONT-only_/5165833) and [after Nanopolish](https://figshare.com/articles/Canu_v1_5_f356c2c_Nanopolish_v0_7_0_assemblies_ONT-only_/5165845)
* Unicycler v0.4.0 [before Nanopolish](https://figshare.com/articles/Unicycler_v0_4_0_assemblies_ONT-only_/5170747) and [after Nanopolish](https://figshare.com/articles/Unicycler_v0_4_0_Nanopolish_v0_7_0_assemblies_ONT-only_/5170753)

#### Metrics

* Mean N50: as described above.
    * When an assembly completes the chromosome sequence, that will be the assembly N50 (about 5.3 Mbp in these samples).
* Complete chromosomes: how many chromosomes completely assembled into a single contig, totalled over all 12 samples.
    * For both Unicycler and Canu, we defined "completely assembled" as the chromosome being circularised in the graph.
* Complete large/small plasmids: as described above
* Estimated error rate: as described above

#### Conclusions

Neither Canu nor Unicycler was particular good recovering small plasmids. This is probably because the small plasmids are very underrepresented in the ONT reads. Unicycler did manage to assemble a few small plasmids and Canu didn't get any. Altering Canu's settings as described [here](http://canu.readthedocs.io/en/latest/faq.html#why-is-my-assembly-is-missing-my-favorite-short-plasmid) may help.

The estimated error rates for both Canu and Unicycler are much higher than Illumina-only assemblies: near 1% (i.e. one error per ~100 bp in the assembly). Unicycler's error rates were slightly lower than Canu's, probably due to its repeated application of [Racon](https://github.com/isovic/racon) to the assembly. Running Racon on Canu's assembly would most likely result in a similar error rate to Unicycler's assemblies. Nanopolish was able to repair about half of the errors.



## Results: hybrid assemblies

| Assembler  | Mean N50  | Complete chromosomes | Complete large plasmids | Complete small plasmids | 100% complete | Estimated error rate |
| :--------: | --------: | -------------------: | ----------------------: | ----------------------: | ------------: | -------------------: |
| SPAdes     | 4,391,534 |                  n/a |                     n/a |                     n/a |           n/a |              0.0000% |
| Canu+Pilon | 4,831,660 |               4 / 12 |                 23 / 28 |                  0 / 29 |        0 / 12 |              0.0039% |
| Unicycler  | 5,334,509 |              12 / 12 |                 28 / 28 |                 18 / 29 |        7 / 12 |              0.0000% |

#### Links to assemblies

* [SPAdes v3.10.1](https://figshare.com/articles/SPAdes_v3_10_1_assemblies_hybrid_Illumina_and_ONT_/5165842)
* Canu v1.5 (f356c2c): [1x Pilon](https://figshare.com/articles/Canu_v1_5_f356c2c_1_x_Pilon_v1_22_assemblies_hybrid_Illumina_and_ONT_/5170756), [2x Pilon](https://figshare.com/articles/Canu_v1_5_f356c2c_2_x_Pilon_v1_22_assemblies_hybrid_Illumina_and_ONT_/5170759), [3x Pilon](https://figshare.com/articles/Canu_v1_5_f356c2c_3_x_Pilon_v1_22_assemblies_hybrid_Illumina_and_ONT_/5170762), [4x Pilon](https://figshare.com/articles/Canu_v1_5_f356c2c_4_x_Pilon_v1_22_assemblies_hybrid_Illumina_and_ONT_/5170765), [5x Pilon](https://figshare.com/articles/Canu_v1_5_f356c2c_5_x_Pilon_v1_22_assemblies_hybrid_Illumina_and_ONT_/5170771)
* [Unicycler v0.4.0](https://figshare.com/articles/Unicycler_v0_4_0_assemblies_hybrid_Illumina_and_ONT_/5170750)

#### Metrics

* Mean N50: as described above
* Complete chromosomes: as described above
    * As was the case for Illumina-only assemblies, SPAdes doesn't provide its final assembly in graph form and it was therefore excluded from the "complete" counts.
    * Since the Canu+Pilon assemblies are polished versions of the Canu ONT-only assemblies, their "complete" counts are identical to those for the Canu ONT-only assemblies.
* Complete large/small plasmids: as described above
* 100% complete: how many of the assemblies have a complete chromosome and all plasmids are complete.
* Estimated error rate: as described above

#### Conclusions

Unicycler does quite well here because hybrid assemblies are its primary focus. Of its five assemblies which were not 100% complete, four were due to incomplete small plasmids (underrepresented in the ONT reads). The remaining incomplete assembly (sample INF164) was due to a discrepancy between the Illumina and ONT reads – an 18 kbp sequence was present in the Illumina sample but absent in the ONT sample, causing an incomplete component in the assembly graph. This was most likely caused by a biological change between cultures used for DNA extraction.

Unicycler and SPAdes both produce their hybrid assemblies by scaffolding an Illumina-only assembly graph. This explains why their error rates are as low as the Illumina-only assemblies.

For Canu+Pilon, the error rates for each of the five rounds of Pilon polishing were: 0.0427, 0.0051, 0.0039, 0.0039 and 0.0039. It plateaued at 0.0039%, suggesting that three rounds of Pilon polishing is sufficient. The error rate never got as low as the SPAdes/Unicycler error rate, indicating that Pilon could correct most but not all errors in the ONT-only assembly. I'm not sure what's causing this and it deserves closer investigation!
