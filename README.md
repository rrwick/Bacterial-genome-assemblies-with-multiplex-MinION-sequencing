<p align="center"><img src="logo.png" alt="Bacterial genome assemblies with multiplex-MinION sequencing"></p>
# Completing bacterial genome assemblies with multiplex MinION sequencing

This repository contains supplemental data and code for our paper: [Completing bacterial genome assemblies with multiplex MinION sequencing
](https://sdfosidhfsidfjaosdjiodifjodifjsdof).

Here you will find the scripts used to generate our data, links to the reads and assemblies, scripts used to generate figures, and summaries of our results. We have also included results from other assemblers not mentioned in the paper to facilitate comparison. If you have different assembly methods that you would like to share, we are happy to include the results here! You can do a GitHub pull-request with your results or else create an [issue](https://github.com/rrwick/Bacterial-genome-assemblies-with-multiplex-MinION-sequencing/issues) on this repo.


## Links to read data

* [Raw Illumina reads (before trimming)](https://figshare.com/)
* [Trimmed Illumina reads](https://figshare.com/)
* [Basecalled ONT reads (straight out of Albacore)](https://figshare.com/)
* [Trimmed ONT reads](https://figshare.com/)
* [Subsampled ONT reads](https://figshare.com/)

I did not put the ONT fast5 files on figshare due to their size (157 GB before basecalling, 1.5TB after basecalling). If anybody is interested in these, please contact me and we can try to work something out.


## Basecalling and assembly

The [ONT_barcode_basecalling_and_assembly.sh](ONT_barcode_basecalling_and_assembly.sh) script carries out the following steps:
* Trimming Illumina reads
* Basecalling ONT reads
* Trimming and subsampling ONT reads
* Assembling Illumina-only read sets (with SPAdes and Unicycler)
* Assembling ONT-only read sets (with Canu and Unicycler)
* Assembling hybrid read sets (with SPAdes, Canu+Pilon and Unicycler)
* Polishing ONT-only assemblies with [Nanopolish](https://github.com/jts/nanopolish)

Each of these steps can be turned on/off using the variables at the top of the script. Details for some of the steps are described below.


#### Software versions used

* Albacore: v1.1.2
* [Porechop](https://github.com/rrwick/Porechop): v0.2.1
* [Unicycler](https://github.com/rrwick/Unicycler): v0.4.0
* [Canu](http://canu.readthedocs.io/en/latest/): [snapshot v1.5 +54 changes](https://github.com/marbl/canu/tree/f356c2c3f2eb37b53c4e7bf11e927e3fdff4d747)
* [SPades](http://cab.spbu.ru/software/spades/): v3.10.1
* [Pilon](https://github.com/broadinstitute/pilon): v1.22
* [Nanopolish](https://github.com/jts/nanopolish): v0.6.1


#### Illumina read trimming

We used [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) to trim adapter sequences from the Illumina reads and remove low-quality sequence. We used a conservative quality threshold of 10 to only remove particularly bad sequences.


#### ONT read processing

When basecalling ONT reads using Albacore (ONT's command-line basecaller), we used the `--barcoding` option to sort the reads into barcode bins. We then ran [Porechop](https://github.com/rrwick/Porechop) on each bin to remove adapter sequences and discard chimeric reads.

Notably, when running Porechop we used its barcode binning as well. This was so we would only keep reads where both Albacore and Porechop agreed on the barcode bin. For example, Albacore put 95064 reads in the bin for barcode 1. Of these reads, Porechop put 90919 in the barcode 1 bin, 118 reads into bins for other barcodes, 3941 reads into no bin and 86 reads were discarded as chimeras. By using only the 90919 reads where Albacore and Porechop agree, we can minimise barcode cross-contamination.

All reads shorter than 2 kbp were discarded for each sample - due to the long read N50s this was a very small proportion of the reads. For samples which still had more than 500 Mbp of reads, we subsampled the read set down to 500 Mbp. This was done using read quality - specifically the reads' minimum qscore over a sliding window. This means that the discarded reads were the one which had the lowest quality regions, as indicated by their qscores. This was done with the `fastq_to_fastq.py` script in [this repo](https://github.com/rrwick/Fast5-to-Fastq).


#### Illumina-only assembly commands
* SPAdes: `spades.py -1 short_1.fastq.gz -2 short_2.fastq.gz -o out_dir --careful`
* Unicycler: `unicycler -1 short_1.fastq.gz -2 short_2.fastq.gz -o out_dir`

For SPAdes, the `contigs.fasta` file was taken as the final assembly.


#### ONT-only assembly commands
* Canu: `canu -p canu -d out_dir genomeSize=5.5m -nanopore-raw long.fastq.gz`
* Unicycler: `unicycler -l long.fastq.gz -o out_dir`
* Nanopolish: `python nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 nanopolish variants --consensus polished.{1}.fa -w {1} -r reads.fa -b reads.sorted.bam -g draft.fa -t 4 --min-candidate-frequency 0.1; python nanopolish_merge.py polished.*.fa > polished_genome.fa`


#### Hybrid assembly commands
* SPAdes: `spades.py -1 short_1.fastq.gz -2 short_2.fastq.gz --nanopore long.fastq.gz -o out_dir --careful`
* Canu+Pilon:
  * `canu -p canu -d out_dir genomeSize=5.5m -nanopore-raw long.fastq.gz`
  * Then fives rounds of Pilon:
    * `bowtie2 --local --very-sensitive-local -I 0 -X 2000 -x before_polish.fasta -1 short_1.fastq.gz -2 short_2.fastq.gz | samtools sort -o alignments.bam -T reads.tmp -; samtools index alignments.bam`
    * `java -jar ~/pilon-1.22.jar --genome before_polish.fasta --frags alignments.bam --changes --output after_polish --outdir out_dir --fix all`
* Unicycler: `unicycler -1 short_1.fastq.gz -2 short_2.fastq.gz -l long.fastq.gz -o out_dir`


#### Polishing with Nanopolish

We used Nanopolish on the ONT-only assemblies to get their base-level accuracy as high as possible. For this step we used _all_ ONT reads for which Albacore and Porechop agreed on the barcode bin (before the read sets were subsampled to 500 Mbp). After using `nanopolish extract` to produce a fasta file from Albacore's output directory, we used [this script](nanopolish_read_filter.py) to exclude reads where Porechop disagreed on the bin.

We tried a second round of Nanopolish but found that it barely improved the results, so the values here are for a single round of Nanopolish.


## Depth per replicon

The files in the [depth_per_replicon](depth_per_replicon) directory were used to generate the supplementary figure showing the read depth for each replicon. This demonstrates that small plasmids were very underrepresented in the ONT reads.



## Error rate estimation

To estimate error rates, we:
* took the 25 largest contigs from Unicycler's Illumina-only assembly
* trimmed 1 kbp off each end
* BLASTed the assemblies using these contigs, keeping only the best hit per contig
* Averaged the error rate over the hits

The code to carry this out is in the [error_rate_estimation](error_rate_estimation) directory.

We used this method because we trust the base calls in an Illumina-only assembly for non-repetitive sequence. By taking only the largest Illumina-only contigs, we avoid repeat sequences. Since these contigs usually end in repeat sequences (repeats being the most common contig-length-limiting feature for Illumina-only assembly), we trim off a kilobase of sequence to ensure our test sequences are not too close to a repeat.

This method for error rate estimation therefore only covers non-repetitive DNA. Error rates in repetitive regions will possibly be higher.



## Result table





## Results: Illumina-only assemblies

Links:
* [SPAdes v3.10.1](https://figshare.com/)
* [Unicycler v0.4.0](https://figshare.com/)

Metrics:
* Mean contigs: the number of contigs in the assembly, averaged over all 12 samples (fewer is better).
* Mean N50: the assembly N50, averaged over all 12 samples (bigger is better). Illumina-only assemblies almost never complete on their own, so this value is much less than the chromosome size.
* Complete large/small plasmids: how many plasmids completely assembled into a single contig, totaled over all 12 samples. For Unicycler, 'completely assembled' means the plasmid is circularised in the graph (i.e. a link joining its start and end). Large plasmids were defined as over 10 kbp (though there were no plasmids between 7 and 60 kbp). The total number of plasmids in the 12 isolates (28 large, 29 small) was determined from the manually-completed assemblies. The completed plasmid counts aren't available for SPAdes, as it outputs its final assembly in contig form, not as a graph, so there's no easy _de novo_ way to tell if a contig is a complete replicon.

| Assembler | Mean contigs | Mean N50 | Complete large plasmids | Complete small plasmids |
| :-------: | -----------: | -------: | ----------------------: | ----------------------: |
| SPAdes    |        379.1 |  218,479 |                     n/a |                     n/a |
| Unicycler |        191.8 |  293,648 |                  2 / 28 |                 12 / 29 |

Overall, Unicycler and SPAdes perform similarly when assembling the Illumina reads. It's worth remembering here that Unicycler uses SPAdes to assemble Illumina reads.

The SPAdes mean contig count is greatly inflated by sample INF163 which has some low-depth contamination (which makes contigs in the SPAdes assembly but is filtered out in the Unicycler assembly). Excluding that sample, the mean contig count for SPAdes is 213.7, much closer to Unicycler's value.

Unicycler achieves somewhat better N50 values because it uses a wider k-mer range than SPAdes does by default. Experimenting with larger values for SPAdes' `-k` option would probably result in N50 values like Unicycler's.

As expected for short reads, neither assembler was very good at completing large plasmids, as they usually contained shared sequence with other replicons. Even though exact completed-plasmid counts aren't available for SPAdes, it seemed to perform similarly to Unicycler on small plasmids - assembling them into single contigs when they only contain unique sequence, assembling them into incomplete contigs when they share sequence with each other.


## Results: ONT-only assemblies

Links:
* [Canu v1.5 (f356c2c)](https://figshare.com/)
* [Canu v1.5 (f356c2c) + Nanopolish v0.7.0](https://figshare.com/)
* [Unicycler v0.4.0](https://figshare.com/)
* [Unicycler v0.4.0 + Nanopolish v0.7.0](https://figshare.com/)

Metrics:
* Mean N50: as described above. When an assembly completes the chromosome sequence, that will be the assembly N50 (about 5.3 Mbp in these samples).
* Complete chromosomes: how many chromosomes completely assembled into a single contig, totaled over all 12 samples. For both Unicycler and Canu, 'completely assembled' means the chromosome is circularised in the graph.
* Complete large/small plasmids: as described above
* Estimated error rate: estimate of the assembly's base-level error rate - method described above.

| Assembler | Mean N50  | Complete chromosomes | Complete large plasmids | Complete small plasmids | Estimated error rate (pre-Nanopolish) | Estimated error rate (post-Nanopolish) |
| :-------: | --------: | -------------------: | ----------------------: | ----------------------: | ------------------------------------: | -------------------------------------: |
| Canu      | 4,784,356 |               4 / 12 |                 23 / 28 |                  0 / 29 |                                1.224% |                                 0.670% |
| Unicycler | 4,965,584 |               7 / 12 |                 27 / 28 |                  5 / 29 |                                1.021% |                                 0.619% |

Neither Canu nor Unicycler was particular good recovering small plasmids. This may be because the small plasmids are very underrepresented in the ONT reads. Unicycler did manage to assemble a few small plasmids. Canu didn't get any, but altering Canu's settings as described [here](http://canu.readthedocs.io/en/latest/faq.html#why-is-my-assembly-is-missing-my-favorite-short-plasmid) may help.

The estimated error rates of Unicycler's assemblies were lower than Canu's, probably due to its repeated application of [Racon](https://github.com/isovic/racon) to the assembly. Running Racon on Canu's assembly would most likely result in a similar error rate to Unicycler's assemblies.


## Results: hybrid assemblies

Links:
* [SPAdes v3.10.1](https://figshare.com/)
* [Canu v1.5 (f356c2c) + 1x Pilon v1.22](https://figshare.com/)
* [Canu v1.5 (f356c2c) + 2x Pilon v1.22](https://figshare.com/)
* [Canu v1.5 (f356c2c) + 3x Pilon v1.22](https://figshare.com/)
* [Canu v1.5 (f356c2c) + 4x Pilon v1.22](https://figshare.com/)
* [Canu v1.5 (f356c2c) + 5x Pilon v1.22](https://figshare.com/)
* [Unicycler v0.4.0](https://figshare.com/)

Metrics:
* Mean N50: as described above
* Complete chromosomes: as described above
* Complete large/small plasmids: as described above
* 100% complete: how many of the assemblies have a complete chromosome and all plasmids are complete.
* Estimated error rate: as described above

| Assembler  | Mean N50  | Complete chromosomes | Complete large plasmids | Complete small plasmids | 100% complete | Estimated error rate |
| :--------: | --------: | -------------------: | ----------------------: | ----------------------: | ------------: | -------------------: |
| SPAdes     | 4,391,534 |                  n/a |                     n/a |                     n/a |           n/a |              0.0000% |
| Canu+Pilon | 4,831,660 |               4 / 12 |                 23 / 28 |                  0 / 29 |        0 / 12 |              0.0058% |
| Unicycler  | 5,334,509 |              12 / 12 |                 28 / 28 |                 18 / 29 |        7 / 12 |              0.0000% |

As was the case for Illumina-only assemblies, since SPAdes doesn't provide its final assembly in graph form, it is difficult to tell whether or not a contig represents a complete replicon. It was therefore excluded from the 'complete' counts. SPAdes' N50 values reflect that it often but not always assembled the chromosome into a single contig.

Since the Canu+Pilon assemblies are polished versions of the Canu ONT-only assemblies, the 'complete' counts are unchanged. The mean N50 for Canu+Pilon is slightly higher than it was for Canu, revealing that Pilon is inserting bases more often than removing them.

Unicycler does quite well here because hybrid assemblies are its primary focus. Of its five assemblies which did not complete 100%, four were due to incomplete small plasmids. Small plasmids were very underrepresented in the ONT reads, and Unicycler failed to separate small plasmids with shared sequence. The remaining incomplete assembly (sample INF164) was due to a discrepancy between the Illumina and ONT reads - an 18 kbp sequence was present in the Illumina sample but absent in the ONT sample, causing an incomplete component in the assembly graph.

The estimated error rates are low for all hybrid assembly methods. The Canu+Pilon error rate decreased with the first couple rounds of Pilon polishing but plateaued after 3-4 rounds. The values shown here are after 5 rounds of Pilon polishing. Interestingly, the Canu+Pilon error rate never got as low as the SPAdes/Unicycler error rate, implying that Pilon did not correct all ONT consensus errors.

