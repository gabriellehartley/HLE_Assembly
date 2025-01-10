      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!   Hoolock leuconedys Assembly   !!
      !!           LCL : Betty           !! 
      !!         Assembly  Steps         !!
      !!        Gabrielle Hartley        !!
      !!   gabrielle.hartley@uconn.edu   !!
      !!          R O'Neill Lab          !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## **HLE Assembly Steps**

### 1. Flye Assembly from Raw ONT Reads
Assembled raw ONT fastq reads to generate fasta formatted assembly. ONT reads were basecalled with the SUP algorithm using Guppy v5.0.16-GPU and Flye was subsequently run with the ```--nano-hq``` flag for SUP reads.

```
module load flye/2.9
flye --nano-hq /core/projects/EBP/Oneill/reads/nanopore/promethion/gibbon/Combined_Reads_Super_Accurate/HLE_Super_Accurate_ONT_combined.fastq \
        --genome-size 2.8g \
        --out-dir HLE_Flye_2021DEC08 \
        --asm-coverage 30 \
        --threads 32
```

| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 94.0% (8671)             | 92.5% (12749)             |
| *Single copy*    | 91.8% (234)             | 91.5% (8439)             | 90.6% (12485)             |
| *Multi copy*     | 7.8% (20)               | 2.5% (232)               | 1.9% (264)                |
| Fragmented       | 0.4% (1)                | 2.4% (222)               | 2.6% (364)                |
| Missing          | 0.0% (0)                | 3.6% (333)               | 4.9% (667)                |



| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    1040       | No. contigs                 |    1038       |
| No. contigs (>= 1000 bp)    |    999        | Largest contig              |    121643331  |
| No. contigs (>= 5000 bp)    |    895        | Total length                |    2795581140 |   
| No. contigs (>= 10000 bp)   |    785        | GC (%)                      |    41.46      |
| No. contigs (>= 25000 bp)   |    663        | N50                         |    31463617   |
| No. contigs (>= 50000 bp    |    560        | N75                         |    13778777   |
| Total length (>= 0 bp)      |    2795581655 | L50                         |    28         |
| Total length (>= 1000 bp)   |    2795551663 | L75                         |    63         |
| Total length (>= 5000 bp)   |    2795270428 | No. N's per 100 kbp         |    0.00       |
| Total length (>= 10000 bp)  |    2794473240 |
| Total length (>= 25000 bp)  |    2792502699 |
| Total length (>= 50000 bp)  |    2788686322 |

### 2. Long Read Polishing using Medaka
Polished assembly using long ONT read data basecalled with SUP algorithm. 

```
module load medaka/1.4.3

medaka_consensus -b 100 -i /core/projects/EBP/Oneill/reads/nanopore/promethion/gibbon/Combined_Reads_Super_Accurate/HLE_Super_Accurate_ONT_combined.fastq \
     -d /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Flye_2021DEC08/assembly.fasta \
     -o MEDAKA_HLE_2021DEC09 \
     -t 12 \ 
     -m r941_prom_sup_g507 
```

| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 95.3% (8791)             | 94.7% (13042)             |
| *Single copy*    | 91.4% (233)             | 92.8% (8561)             | 92.8% (12784)             |
| *Multi copy*     | 8.2% (21)               | 2.5% (230)               | 1.9% (258)                |
| Fragmented       | 0.0% (0)                | 1.6% (147)               | 1.5% (212)                |
| Missing          | 0.4% (1)                | 3.1% (288)               | 3.8% (526)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    1040       | No. contigs                 |    1035       |
| No. contigs (>= 1000 bp)    |    999        | Largest contig              |    121718504  |
| No. contigs (>= 5000 bp)    |    895        | Total length                |    2797566732 |
| No. contigs (>= 10000 bp)   |    785        | GC (%)                      |    41.45      |  
| No. contigs (>= 25000 bp)   |    663        | N50                         |    31472656   |
| No. contigs (>= 50000 bp    |    560        | N75                         |    13784686   |
| Total length (>= 0 bp)      |    2797568188 | L50                         |    28         |
| Total length (>= 1000 bp)   |    2797539248 | L75                         |    63         |
| Total length (>= 5000 bp)   |    2797256969 | No. N's per 100 kbp         |    0.00       |
| Total length (>= 10000 bp)  |    2796458597 |
| Total length (>= 25000 bp)  |    2794485364 |
| Total length (>= 50000 bp)  |    2790664566 |

### 3. Short Read Polishing Using Pilon
Polished assembly using Illumina ChIP Input, WGS075, and WGS076 sequencing runs for Betty.

#### a.) Align and Process ChIP Input, WGS075, and WGS076 Samples to Assembly

##### i.) Index Nanopolished Genome 

```
module load bwa/0.7.17

bwa index -p  HLE_Medaka_2021DEC10 /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Medaka_2021DEC09/MEDAKA_HLE_2021DEC09/consensus.fasta
```

##### ii.) Map and Process ChIP Input, WGS075, and WGS76 Samples Using BWA 

```
module load bwa/0.7.17

bwa mem -t 16 HLE_Medaka_2021DEC10 /core/labs/Oneill/ghartley/Reads/GibbonChIPSeqReads/Betty-input-R1-combined.fq /core/labs/Oneill/ghartley/Reads/GibbonChIPSeqReads/Betty-input-R1-combined.fq > BettyinputChIPReads_vs_HLE_Medaka_2021DEC10_bwa.sam
```

```
module load bwa/0.7.17

bwa mem -t 16 HLE_Medaka_2021DEC10 /core/labs/Oneill/ghartley/Reads/WGS-HLE/BETTY/Betty-HLE-SW075_R1.fastq.gz /core/labs/Oneill/ghartley/Reads/WGS-HLE/BETTY/Betty-HLE-SW075_R2.fastq.gz > BettyWGS075Reads_vs_HLE_Medaka_2021DEC10_bwa.sam
```

```
module load bwa/0.7.17

bwa mem -t 16 HLE_Medaka_2021DEC10 /core/labs/Oneill/ghartley/Reads/WGS-HLE/BETTY/Betty-HLE-SW076_R1.fastq.gz /core/labs/Oneill/ghartley/Reads/WGS-HLE/BETTY/Betty-HLE-SW075_R2.fastq.gz > BettyWGS076Reads_vs_HLE_Medaka_2021DEC10_bwa.sam
```

##### iii.) Convert SAM to BAM
```
module load samtools/1.7

samtools view -Sb BettyinputChIPReads_vs_HLE_Medaka_2021DEC10_bwa.sam > BettyinputChIPReads_vs_HLE_Medaka_2021DEC10_bwa.bam
```

```
module load samtools/1.7

samtools view -Sb BettyWGS075Reads_vs_HLE_Medaka_2021DEC10_bwa.sam > BettyWGS075Reads_vs_HLE_Medaka_2021DEC10_bwa.bam
```

```
module load samtools/1.7

samtools view -Sb BettyWGS076Reads_vs_HLE_Medaka_2021DEC10_bwa.sam > BettyWGS076Reads_vs_HLE_Medaka_2021DEC10_bwa.bam
```

##### iv.) Sort BAM Files

```
module load samtools/1.7

samtools sort BettyinputChIPReads_vs_HLE_Medaka_2021DEC10_bwa.bam > BettyinputChIPReads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam
```

```
module load samtools/1.7

samtools sort BettyWGS075Reads_vs_HLE_Medaka_2021DEC10_bwa.bam > BettyWGS075Reads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam
```

```
module load samtools/1.7

samtools sort BettyWGS076Reads_vs_HLE_Medaka_2021DEC10_bwa.bam > BettyWGS076Reads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam
```

##### v.) Index BAM Files
```
module load samtools/1.7

samtools index BettyinputChIPReads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam
```

```
module load samtools/1.7

samtools index BettyWGS075Reads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam
```

```
module load samtools/1.7

samtools index BettyWGS076Reads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam
```

##### vi.) Get Mapping Statistics

```
module load samtools/1.7
samtools flagstat BettyinputChIPReads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam > BettyinputChIPReads_vs_HLE_Medaka_2021DEC10_bwa_STATS.txt
```

```
module load samtools/1.7
samtools flagstat BettyWGS075Reads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam > BettyWGS075Reads_vs_HLE_Medaka_2021DEC10_bwa_STATS.txt
```

```
module load samtools/1.7
samtools flagstat BettyWGS076Reads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam > BettyWGS076Reads_vs_HLE_Medaka_2021DEC10_bwa_STATS.txt
```

| BWA Stats v. 0.7.17   | ChIP Input | Betty WGS075  | Betty WGS076  |
| -------------         | ---        | ---           | ---           |
| Total Reads           | 74549938   | 469534006     | 88170048      |
| Mapped                | 99.42%     | 99.66%        | 98.92%        |


#### b.) Split Assembly Into Fragments 
The assembly was split into smaller fragments broken at fasta record boundaries in order to improve Pilon runtime from several days to several hours. 
```
module load GenomeBrowser/20180626

faSplit sequence /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Medaka_2021DEC09/MEDAKA_HLE_2021DEC09/consensus.fasta 100 consensus_split_
```

#### c.) Polish Assembly Fragments Using Pilon
Note: This step was run on multiple assembly fragments generated in step 3b substituted for the --genome file split at fasta record boundaries and combined at the end to improve runtime. Command below represents example.

```
module load java
module load pilon/1.22

java -Xmx200g -jar /isg/shared/apps/pilon/1.22/pilon-1.22.jar --genome consensus_split_000.fa 
    --frags BettyinputChIPReads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam \
    --frags BettyWGS075Reads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam \
    --frags BettyWGS076Reads_vs_HLE_Medaka_2021DEC10_bwa_sorted.bam \
    --threads 24 --output polishedfrag_0 
```

#### d.) Combine Polished Genome Fragments to Generate Complete Polished Assembly
```
cat polishedfrag_*fa > HLE_Pilon_Assembly.fasta
```
| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 95.8% (8846)             | 95.5% (13163)             |
| *Single copy*    | 91.4% (233)             | 93.3% (8611)             | 93.6% (12902)             |
| *Multi copy*     | 8.2% (21)               | 2.5% (235)               | 1.9% (261)                |
| Fragmented       | 0.0% (0)                | 1.3% (122)               | 1.2% (163)                |
| Missing          | 0.4% (1)                | 2.9% (258)               | 3.3% (454)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    1040       | No. contigs                 |    1035       |
| No. contigs (>= 1000 bp)    |    999        | Largest contig              |    121688009  |
| No. contigs (>= 5000 bp)    |    895        | Total length                |    2796976253 |
| No. contigs (>= 10000 bp)   |    785        | GC (%)                      |    41.45      |  
| No. contigs (>= 25000 bp)   |    663        | N50                         |    31464756   |
| No. contigs (>= 50000 bp    |    560        | N75                         |    13783658   |
| Total length (>= 0 bp)      |    2796977709 | L50                         |    28         |
| Total length (>= 1000 bp)   |    2796948788 | L75                         |    63         |
| Total length (>= 5000 bp)   |    2796666656 | No. N's per 100 kbp         |    0.00       |
| Total length (>= 10000 bp)  |    2795868549 |
| Total length (>= 25000 bp)  |    2793898063 |
| Total length (>= 50000 bp)  |    2790085328 |

### 4. Purge Haplotigs and Remove Artefacts

#### a.) Map and Sort All ONT Reads to Assembly Using Minimap2
```
module load minimap2/2.15
module load samtools/1.7

minimap2 -ax map-ont \
    -t 8 \
    /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Pilon_2021DEC09/Combined/HLE_Pilon_Assembly.fasta /core/projects/EBP/Oneill/reads/nanopore/promethion/gibbon/Combined_Reads_Super_Accurate/HLE_Super_Accurate_ONT_combined.fasta | samtools sort -o HLECombinedONTReads_vs_HLE_Genome_Pilon_2021DEC11.sorted.bam -T HLE_allreads_combined.tmp
```
#### b.) Generate a Coverage Histogram with Purge Haplotigs
```
module load purge_haplotigs/1.0
module load R/3.5.1
module load samtools/1.7
module load bedtools/2.29.0

purge_haplotigs readhist \
    -b /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_PurgeHaplotigs_2021DEC11/HLECombinedONTReads_vs_HLE_Genome_Pilon_2021DEC11.sorted.bam \
    -g /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Pilon_2021DEC09/Combined/HLE_Pilon_Assembly.fasta \
    -t 16
```

![Coverage Map](/Hoolock_leuconedys_Assembly/HLECombinedONTReads_vs_HLE_Genome_Pilon_2021DEC11.sorted.bam.histogram.png)

#### c.) Set Cutoffs with Purge Haplotigs Following Manual Review of Peaks
Based on the above coverage plot, the following were chosen: ```-l 10 -m 50 -h 55```.
```
module load purge_haplotigs/1.0
module load R/3.5.1
module load samtools/1.7
module load bedtools/2.29.0

purge_haplotigs contigcov \
    -i /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_PurgeHaplotigs_2021DEC11/HLECombinedONTReads_vs_HLE_Genome_Pilon_2021DEC11.sorted.bam.gencov \
    -l 10 -m 50 -h 55
```
#### d.) Run the Purging Pipeline Using Purge Haplotigs
```
module load purge_haplotigs/1.0
module load R/3.5.1
module load samtools/1.7
module load bedtools/2.29.0

purge_haplotigs purge -g /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Pilon_2021DEC09/Combined/HLE_Pilon_Assembly.fasta \
-c /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_PurgeHaplotigs_2021DEC11/coverage_stats.csv \
-b /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_PurgeHaplotigs_2021DEC11/HLECombinedONTReads_vs_HLE_Genome_Pilon_2021DEC11.sorted.bam
```

| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 95.9% (8842)             | 95.4% (13155)             |
| *Single copy*    | 91.4% (233)             | 93.4% (8615)             | 93.6% (12904)             |
| *Multi copy*     | 8.2% (21)               | 2.5% (227)               | 1.8% (251)                |
| Fragmented       | 0.0% (0)                | 1.3% (121)               | 1.2% (161)                |
| Missing          | 0.4% (1)                | 2.8% (263)               | 3.4% (464)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    472        | No. contigs                 |    470        |
| No. contigs (>= 1000 bp)    |    462        | Largest contig              |    121688009  |
| No. contigs (>= 5000 bp)    |    434        | Total length                |    2760122363 |
| No. contigs (>= 10000 bp)   |    414        | GC (%)                      |    41.27      |  
| No. contigs (>= 25000 bp)   |    391        | N50                         |    31776345   |
| No. contigs (>= 50000 bp    |    376        | N75                         |    13953348   |
| Total length (>= 0 bp)      |    2760122655 | L50                         |    27         |
| Total length (>= 1000 bp)   |    2760116563 | L75                         |    61         |
| Total length (>= 5000 bp)   |    2760045211 | No. N's per 100 kbp         |    0.00       |
| Total length (>= 10000 bp)  |    2759898480 |
| Total length (>= 25000 bp)  |    2759532357 |
| Total length (>= 50000 bp)  |    2758982978 |

### 5. Impose 3kb Contig Threshold
Remove all contigs less than 3 kb in length prior to scaffolding. Moving forward, assembly has the following name designation, **FMPP_3kbmin** (**F**lye, **M**edaka, **P**ilon, **P**urgeHaplotigs).

```
module load seqtk/1.2

seqtk seq -L 3000 curated.fasta > HLE_FMPP_3kbmin.fasta
```
| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 95.9% (8842)             | 95.4% (13155)             |
| *Single copy*    | 91.4% (233)             | 93.4% (8615)             | 93.6% (12904)             |
| *Multi copy*     | 8.2% (21)               | 2.5% (227)               | 1.8% (251)                |
| Fragmented       | 0.0% (0)                | 1.3% (121)               | 1.2% (161)                |
| Missing          | 0.4% (1)                | 2.8% (263)               | 3.4% (464)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    443        | No. contigs                 |    443        |
| No. contigs (>= 1000 bp)    |    443        | Largest contig              |    121688009  |
| No. contigs (>= 5000 bp)    |    434        | Total length                |    2760080048 |
| No. contigs (>= 10000 bp)   |    414        | GC (%)                      |    41.27      |  
| No. contigs (>= 25000 bp)   |    391        | N50                         |    31776345   |
| No. contigs (>= 50000 bp    |    376        | N75                         |    13953348   |
| Total length (>= 0 bp)      |    2760080048 | L50                         |    27         |
| Total length (>= 1000 bp)   |    2760080048 | L75                         |    61         |
| Total length (>= 5000 bp)   |    2760045211 | No. N's per 100 kbp         |    0.00       |
| Total length (>= 10000 bp)  |    2759898480 |
| Total length (>= 25000 bp)  |    2759532357 |
| Total length (>= 50000 bp)  |    2758982978 |

### 6. Scaffold Assembly with Juicer, 3D-DNA, and JBAT
Dovetail Omni-C was performed using a Betty LCL cell pellet. A small, 1-2 million read library was first generated and processed according to Dovetail's QC pipeline to ensure library quality, then sequenced to a depth of 300 million reads.

##### a.) Create a .genome file
A ```.genome``` file was created of contig sizes. This is needed as an input for Juicer.

```
module load samtools/1.7

samtools faidx HLE_FMPP_3kbmin.fasta
cut -f1,2 HLE_FMPP_3kbmin.fasta.fai > HLE_FMPP_3kbmin.fasta.genome
```

##### b.) Pre-Process Omni-C Reads with Juicer
A Juicer directory (```Juicer_HLE_2021DEC13```) was set up with the following folders: ```fastq```, ```references```, ```scripts```. The assembly (```HLE_FMPP_3kbmin.fasta```) was placed in ```references```. The Omni-C reads (```Gabby-HLE-OmniC_S1_L002_R1_001.fastq.gz``` and ```Gabby-HLE-OmniC_S1_L002_R2_001.fastq.gz```) were placed in ```fastq```. Updated Juicer scripts to work with the SLURM Xanadu cluster were placed in the ```scripts``` folder, see [here](/Juicer_Scripts/) for folder contents. In order to run Juicer, the following code was run. Once all of the split fastq files populated in the ```splits``` folder, the run was terminated. The same code was submitted again, which runs the alignment process. This is due to an unknown bug that prevents the Juicer scripts from progressing from the splits to alignment phase. Only the ```merged_nodups.txt``` file is needed to proceed with scaffolding.

```
module load java-sdk/1.8.0_92
module load bwa/0.7.5a

bash /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Scaffolding_2021DEC13/Juicer_HLE_2021DEC13/scripts/juicer.sh \
    -D /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Scaffolding_2021DEC13/Juicer_HLE_2021DEC13/ \
    -z /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Scaffolding_2021DEC13/Juicer_HLE_2021DEC13/references/HLE_FMPP_3kbmin.fasta \
    -s none \
    -p ./HLE_FMPP_3kbmin.fasta.genome
```
##### c.) Scaffolding with 3D-DNA (Pre-JBAT Review)
Scaffolding was performed with 3D-DNA. The ```wrap-fasta-sequence.awk``` script was first used on the assembly to format the fasta correctly; if this script is not used, scaffolding will take weeks as opposed to less than one day. 

```
awk -f wrap-fasta-sequence.awk HLE_FMPP_3kbmin.fasta > HLE_FMPP_3kbmin.wrapped.fasta
```

Then, the wrapped fasta was piped into the 3D-DNA scaffolding program. Note that the temp directories need to be set correctly otherwise the run will have errors.

```
module load gnu-parallel/20160622

hostname
echo "\nStart time:"
date

export TMPDIR="/core/projects/EBP/Oneill/Gibbon_Working/hic_temp3"
export TMP_DIR="/core/projects/EBP/Oneill/Gibbon_Working/hic_temp3"
export _JAVA_OPTIONS=-Djava.io.tmpdir=/core/projects/EBP/Oneill/Gibbon_Working/hic_temp3

bash ../3d-dna-master/run-asm-pipeline.sh ./HLE_FMPP_3kbmin.wrapped.fasta /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Scaffolding_2021DEC13/Juicer_HLE_2021DEC13/aligned/merged_nodups.txt

echo "\nEnd time:"
date
```

##### d.) Manual review in JBAT
The assembly was reviewed manually in JBAT and errors were corrected. One contig was inverted and contigs were joined as scaffolds.

Pre-Correction:
![Pre-Correction](Hoolock_leuconedys_Assembly/HLE_Assembly_PRE-REVIEW.HiCImage.svg)
Post-Correction:
![Post-Correction](Hoolock_leuconedys_Assembly/HLE_Assembly_POST-REVIEW.HiCImage.svg)

##### e.) Scaffolding with 3D-DNA (Post-JBAT Review)

```
module load gnu-parallel/20160622

hostname
echo "\nStart time:"
date

export TMPDIR="/core/projects/EBP/Oneill/Gibbon_Working/hic_temp3"
export TMP_DIR="/core/projects/EBP/Oneill/Gibbon_Working/hic_temp3"
export _JAVA_OPTIONS=-Djava.io.tmpdir=/core/projects/EBP/Oneill/Gibbon_Working/hic_temp3

bash ../3d-dna-master/run-asm-pipeline-post-review.sh -r HLE_FMPP_3kbmin.wrapped.FINAL.review.assembly HLE_FMPP_3kbmin.wrapped.fasta /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Scaffolding_2021DEC13/Juicer_HLE_2021DEC13/aligned/merged_nodups.txt

echo "\nEnd time:"
date
```

| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 95.8% (8845)             | 95.5% (13151)             |
| *Single copy*    | 91.4% (233)             | 93.5% (8630)             | 93.7% (12906)             |
| *Multi copy*     | 8.2% (21)               | 2.3% (215)               | 1.8% (245)                |
| Fragmented       | 0.0% (0)                | 1.3% (124)               | 1.2% (164)                |
| Missing          | 0.4% (1)                | 2.9% (257)               | 3.3% (465)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    398        | No. contigs                 |    398        |
| No. contigs (>= 1000 bp)    |    398        | Largest contig              |    224464253  |
| No. contigs (>= 5000 bp)    |    290        | Total length                |    2760540048 |
| No. contigs (>= 10000 bp)   |    238        | GC (%)                      |    41.27      |  
| No. contigs (>= 25000 bp)   |    171        | N50                         |    159659651  |
| No. contigs (>= 50000 bp    |    95         | N75                         |    110860885  |
| Total length (>= 0 bp)      |    2760540048 | L50                         |    8          |
| Total length (>= 1000 bp)   |    2760540048 | L75                         |    13         |
| Total length (>= 5000 bp)   |    2760295299 | No. N's per 100 kbp         |    16.66      |
| Total length (>= 10000 bp)  |    2759919281 |
| Total length (>= 25000 bp)  |    2758881201 |
| Total length (>= 50000 bp)  |    2756378502 |

### 7. Set 3 kb Contig Threshold Limit Post Scaffolding

```
module load seqtk/1.2

seqtk seq -L 3000 HLE_FMPP_3kbmin.wrapped.FINAL.fasta > HLE_FMPPJ3_3KBLimit.fasta
```

| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 95.8% (8845)             | 95.5% (13151)             |
| *Single copy*    | 91.4% (233)             | 93.5% (8630)             | 93.7% (12906)             |
| *Multi copy*     | 8.2% (21)               | 2.3% (215)               | 1.8% (245)                |
| Fragmented       | 0.0% (0)                | 1.3% (124)               | 1.2% (164)                |
| Missing          | 0.4% (1)                | 2.9% (257)               | 3.3% (465)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    321        | No. contigs                 |    321        |
| No. contigs (>= 1000 bp)    |    321        | Largest contig              |    224464253  |
| No. contigs (>= 5000 bp)    |    290        | Total length                |    2760418395 |
| No. contigs (>= 10000 bp)   |    238        | GC (%)                      |    41.27      |  
| No. contigs (>= 25000 bp)   |    171        | N50                         |    159659651  |
| No. contigs (>= 50000 bp    |    95         | N75                         |    110860885  |
| Total length (>= 0 bp)      |    2760418395 | L50                         |    8          |
| Total length (>= 1000 bp)   |    2760418395 | L75                         |    13         |
| Total length (>= 5000 bp)   |    2760295299 | No. N's per 100 kbp         |    16.66      |
| Total length (>= 10000 bp)  |    2759919281 |
| Total length (>= 25000 bp)  |    2758881201 |
| Total length (>= 50000 bp)  |    2756378502 |

### 8. Manual Alignment Review 

Jake VanCampen (vancampe@ohsu.edu) from Lucia Carbone (carbone@ohsu.edu) lab performed alignment of hg38 and HLE assembly. 

Filtered alignment blocks can be found [here](Hoolock_leuconedys_Assembly/HLE-HG38-Synteny-2022JAN10-FMPPJ3.xlsx).

Dot Plot:
![HLE - HG38 Dotplot](Hoolock_leuconedys_Assembly/HLE_FMPP.hg38.net.filt.maf.png)

#### a.) Correction of Chromosome Segment Inversions with JBAT

Based on the alignment and syntenic blocks reported in Capozzi et. al (2012), the following chromosomes were assigned. Subsequent changes were made in JBAT according to the table below based on the Omni-C contact map and the corresponding alignment blocks.

| HLE Contig      | Synteny to HG38 | HLE Chr | Correct HLE-HG38 Chr Synteny | Changes Needed         | Action in JBAT Based on Contact Map       |
| ---             | ---             | ---     | ---              | ---                    | ---                                       |
| HiC_Scaffold_1  | 1               | 9p      | 1                | add to 9q              | none                                      |
| HiC_Scaffold_10 | 6-16-5-17-4     | 8       | 4-17-5-16-6      | reverse                | none                                      |
| HiC_Scaffold_11 | 14-20-2-17      | 14      | 20-2-17-14       | reverse 20-17, reverse | reversed 20-17, should now be 14-17-2-20 [video here](Hoolock_leuconedys_Assembly/Chr14_Scaffold11_ManualReview.mov)  |
| HiC_Scaffold_12 | X               | X       | X                | none                   | none                                      |
| HiC_Scaffold_13 | 12-19-1         | 15      | 1-19-12          | reverse                | none                                      |
| HiC_Scaffold_14 | 1-19-12-19      | 18      | 1-19-12-19       | none                   | none                                      |
| HiC_Scaffold_15 | 8-4-17-9        | 1       | 8-4-17-9         | none                   | none                                      |
| HiC_Scaffold_16 | 13-21-10-4      | 2       | 13-21-10-4       | none                   | none                                      |
| HiC_Scaffold_17 | 7-11-18         | 3       | 18-11-7          | reverse                | none                                      |
| HiC_Scaffold_18 | 2-20            | 17      | 20-2             | reverse                | none                                      |
| HiC_Scaffold_19 | 17-2-10-1       | 16      | 17-2-1           | none (10 has low score)| none                                      |
| HiC_Scaffold_2  | 11-8-3-12       | 9q      | 12-3-8-11        | reverse, add to 9p     | none                                      |
| HiC_Scaffold_20 | 11-5-11         | 12      | 5-8-11           | remove 11              | moved small 11 segment to debris, should now be 11-5
| HiC_Scaffold_3  |14-10-1-15-22-16-5-17-14|13|1-14-17-5-16-22-15| reverse 1-14, reverse 15-14  (10 has low score)  | reversed 1-14, 15-14; removed debris at end, should now be 1-10-14-14-17-5-16-22-15 [video here](Hoolock_leuconedys_Assembly/Chr13_Scaffold3_ManualReview.mov) |
| HiC_Scaffold_4  | 15-7-2          | 7       | 15-7-2           | none                   | none                                      |
| HiC_Scaffold_5  | 22-15-7-11-3-8  | 11      | 22-15-7-11-3-8   | none                   | none                                      |
| HiC_Scaffold_6  | 16-5-16-5-16    | 10      | 5-16-5-16        | none                   | none                                      |
| HiC_Scaffold_7  | 3-12-19-3-11    | 6       | 3-12-19-11       | none                   | none                                      |
| HiC_Scaffold_8  | 10-4            | 5       | 4-10             | reverse                | none                                      |
| HiC_Scaffold_9  | 12-2-7-12-6     | 4       | 7-2-12-6         | reverse 12-7               | reversed 12-7, should now be 7-2-12-12-6 [video here](Hoolock_leuconedys_Assembly/Chr4_Scaffold9_ManualReview.mov)                          |

3D-DNA was used to generate a fasta with the corrected inversions.

```
module load gnu-parallel/20160622

hostname
echo "\nStart time:"
date

export TMPDIR="/core/projects/EBP/Oneill/Gibbon_Working/hic_temp3"
export TMP_DIR="/core/projects/EBP/Oneill/Gibbon_Working/hic_temp3"
export _JAVA_OPTIONS=-Djava.io.tmpdir=/core/projects/EBP/Oneill/Gibbon_Working/hic_temp3

bash ../3d-dna-master/run-asm-pipeline-post-review.sh -r HLE_FMPP_3kbmin.wrapped.manualreview.FINAL.review.assembly HLE_FMPP_3kbmin.wrapped.FINAL_Step2.fasta /core/projects/EBP/Oneill/Gibbon_Working/HLE_Assembly/HLE_Scaffolding_2021DEC13/Juicer_HLE_2021DEC13/aligned/merged_nodups.txt

echo "\nEnd time:"
date

```

#### b.) Correction of Whole Chromosome Inversions 
Chromosomes were separated into individual fasta records.

```
module load GenomeBrowser/20180626

faSplit sequence HLE_FMPP_3kbmin.wrapped.FINAL_Step2.fasta 405 manual_review_
```
Then, chromosomes needing to be flipped were inverted using the following command. These included: HiC_Scaffold_10, HiC_Scaffold_11, HiC_Scaffold_13, HiC_Scaffold_17, HiC_Scaffold_18, Hi-C_Scaffold_2, and Hi-C_Scaffold_8.

```
module load emboss/6.6.0

head -1 manual_review_009.fa 
>HiC_scaffold_10
revseq manual_review_009.fa -reverse -complement -outseq manual_review_009.revcomp.fa

head -1 manual_review_010.fa 
>HiC_scaffold_11
revseq manual_review_010.fa -reverse -complement -outseq manual_review_010.revcomp.fa

head -1 manual_review_012.fa 
>HiC_scaffold_13
revseq manual_review_012.fa -reverse -complement -outseq manual_review_012.revcomp.fa

head -1 manual_review_016.fa 
>HiC_scaffold_17
revseq manual_review_016.fa -reverse -complement -outseq manual_review_016.revcomp.fa

head -1 manual_review_017.fa 
>HiC_scaffold_18
revseq manual_review_017.fa -reverse -complement -outseq manual_review_017.revcomp.fa

head -1 manual_review_001.fa
>HiC_scaffold_2
revseq manual_review_001.fa -reverse -complement -outseq manual_review_001.revcomp.fa

head -1 manual_review_007.fa 
>HiC_scaffold_8
revseq manual_review_007.fa -reverse -complement -outseq manual_review_007.revcomp.fa
```

The sequences were concatenated together to create the final manually reviewed assembly. The assembly was named HLE_FMPPJ3M.fasta, for **F**lye, **M**edaka, **P**ilon, **P**urgeHaplotigs, **J**uicer, **3**D-DNA, **M**anual Review. 

```
rm manual_review_009.fa 
rm manual_review_010.fa 
rm manual_review_012.fa 
rm manual_review_016.fa 
rm manual_review_017.fa 
rm manual_review_001.fa 
rm manual_review_007.fa 

cat manual_review_* > HLE_FMPPJ3M.fasta
```

### 9. Set 3 kb Contig Threshold Limit Post Manual Review

Since the 3D-DNA output was used for manual review, the 3kb limit was again imposed on the assembly. 

```
module load seqtk/1.2

seqtk seq -L 3000 HLE_FMPPJ3M.fasta > HLE_FMPPJ3M_3kblimit.fasta
```

| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 95.8% (8834)             | 95.5% (13150)             |
| *Single copy*    | 91.4% (233)             | 93.5% (8622)             | 93.7% (12906)             |
| *Multi copy*     | 8.2% (21)               | 2.3% (212)               | 1.8% (244)                |
| Fragmented       | 0.0% (0)                | 1.3% (123)               | 1.2% (164)                |
| Missing          | 0.4% (1)                | 2.9% (269)               | 3.3% (466)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    324        | No. contigs                 |    324        |
| No. contigs (>= 1000 bp)    |    324        | Largest contig              |    224502753  |
| No. contigs (>= 5000 bp)    |    293        | Total length                |    2760876895 |
| No. contigs (>= 10000 bp)   |    241        | GC (%)                      |    41.27      |  
| No. contigs (>= 25000 bp)   |    174        | N50                         |    159492152  |
| No. contigs (>= 50000 bp    |    98         | N75                         |    108774101  |
| Total length (>= 0 bp)      |    2760876895 | L50                         |    8          |
| Total length (>= 1000 bp)   |    2760876895 | L75                         |    13         |
| Total length (>= 5000 bp)   |    2760753799 | No. N's per 100 kbp         |    33.27      |
| Total length (>= 10000 bp)  |    2760377781 |
| Total length (>= 25000 bp)  |    2759339701 |
| Total length (>= 50000 bp)  |    2756837002 |

### 10. Gap Fill with TGS-GapCloser
Gaps were filled using TGS-GapCloser. Must use script that Patrick Grady (patrick.gs.grady@uconn.edu) updated that makes ```--tgstype pb``` actually work for ONT reads. (ie: ```-   TGS reads type is pb . MINIMAP2_PARAM is  -x ava-ont --score-N 3   MIN_IDY is 0.2 . MIN_MATCH is 200 .```)

```
tgsgapcloser --scaff ./HLE_FMPPJ3M_3kblimit.fasta \
        --tgstype pb \
        --reads /core/projects/EBP/Oneill/reads/nanopore/promethion/gibbon/Combined_Reads_Super_Accurate/HLE_Super_Accurate_ONT_combined.fasta \
        --output HLE_FMPPJ3M_3kblimit_TGS-GapCloser \
        --ne \
	--minmap-arg \'--score-N 3\' \
        --thread 32
```

| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 95.8% (8836)             | 95.5% (13151)             |
| *Single copy*    | 91.4% (233)             | 93.5% (8623)             | 93.7% (12906)             |
| *Multi copy*     | 8.2% (21)               | 2.3% (213)               | 1.8% (245)                |
| Fragmented       | 0.0% (0)                | 1.3% (122)               | 1.2% (164)                |
| Missing          | 0.4% (1)                | 2.9% (268)               | 3.3% (465)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    324        | No. contigs                 |    324        |
| No. contigs (>= 1000 bp)    |    324        | Largest contig              |    224502753  |
| No. contigs (>= 5000 bp)    |    293        | Total length                |    2760876895 |
| No. contigs (>= 10000 bp)   |    241        | GC (%)                      |    41.28      |  
| No. contigs (>= 25000 bp)   |    174        | N50                         |    159724090  |
| No. contigs (>= 50000 bp    |    98         | N75                         |    108836962  |
| Total length (>= 0 bp)      |    2760876895 | L50                         |    8          |
| Total length (>= 1000 bp)   |    2760876895 | L75                         |    13         |
| Total length (>= 5000 bp)   |    2760753799 | No. N's per 100 kbp         |    13.00      |
| Total length (>= 10000 bp)  |    2760377781 |
| Total length (>= 25000 bp)  |    2759339701 |
| Total length (>= 50000 bp)  |    2756837002 |

### 11. Polish Gap Fill with Pilon
Pilon was run to polish filled gaps as according to Step 3 using Illumina ChIP Input, WGS075, and WGS076 sequencing runs for Betty.

#### a.) Align and Process ChIP Input, WGS075, and WGS076 Samples to Assembly

##### i.) Index Nanopolished Genome 

```
module load bwa/0.7.17

bwa index -p  HLE_TGS_2022JAN14 ./HLE_FMPPJ3MT_3kblimit.fasta
```

##### ii.) Map and Process ChIP Input, WGS075, and WGS76 Samples Using BWA 

```
module load bwa/0.7.17

bwa mem -t 16 HLE_Medaka_2021DEC10 /core/labs/Oneill/ghartley/Reads/GibbonChIPSeqReads/Betty-input-R1-combined.fq /core/labs/Oneill/ghartley/Reads/GibbonChIPSeqReads/Betty-input-R1-combined.fq > BettyinputChIPReads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.sam
```

```
module load bwa/0.7.17

bwa mem -t 16 HLE_Medaka_2021DEC10 /core/labs/Oneill/ghartley/Reads/WGS-HLE/BETTY/Betty-HLE-SW075_R1.fastq.gz /core/labs/Oneill/ghartley/Reads/WGS-HLE/BETTY/Betty-HLE-SW075_R2.fastq.gz > BettyWGS075Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.sam
```

```
module load bwa/0.7.17

bwa mem -t 16 HLE_Medaka_2021DEC10 /core/labs/Oneill/ghartley/Reads/WGS-HLE/BETTY/Betty-HLE-SW076_R1.fastq.gz /core/labs/Oneill/ghartley/Reads/WGS-HLE/BETTY/Betty-HLE-SW075_R2.fastq.gz > BettyWGS076Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.sam
```

##### iii.) Convert SAM to BAM
```
module load samtools/1.7

samtools view -Sb BettyinputChIPReads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.sam > BettyinputChIPReads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.bam
```

```
module load samtools/1.7

samtools view -Sb BettyWGS075Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.sam > BettyWGS075Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.bam
```

```
module load samtools/1.7

samtools view -Sb BettyWGS076Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.sam > BettyWGS076Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.bam
```

##### iv.) Sort BAM Files

```
module load samtools/1.7

samtools sort BettyinputChIPReads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.bam > BettyinputChIPReads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam
```

```
module load samtools/1.7

samtools sort BettyWGS075Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.bam > BettyWGS075Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam
```

```
module load samtools/1.7

samtools sort BettyWGS076Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa.bam > BettyWGS076Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam
```

##### v.) Index BAM Files
```
module load samtools/1.7

samtools index BettyinputChIPReads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam
```

```
module load samtools/1.7

samtools index BettyWGS075Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam
```

```
module load samtools/1.7

samtools index BettyWGS076Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam
```

##### vi.) Get Mapping Statistics

```
module load samtools/1.7
samtools flagstat BettyinputChIPReads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam > BettyinputChIPReads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_STATS.txt
```

```
module load samtools/1.7
samtools flagstat BettyWGS075Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam > BettyWGS075Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_STATS.txt
```

```
module load samtools/1.7
samtools flagstat BettyWGS076Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam > BettyWGS076Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_STATS.txt
```

| BWA Stats v. 0.7.17   | ChIP Input | Betty WGS075  | Betty WGS076  |
| -------------         | ---        | ---           | ---           |
| Total Reads           | 74557499   | 469761140     | 88207829      |
| Mapped                | 99.13%     | 99.50%        | 98.69%        |

#### b.) Split Assembly Into Fragments 
The assembly was split into smaller fragments broken at fasta record boundaries in order to improve Pilon runtime from several days to several hours. Each fragment will
```
module load GenomeBrowser/20180626

faSplit sequence ./HLE_FMPPJ3MT_3kblimit.fasta 100 consensus_split_
```

#### c.) Polish Assembly Fragments Using Pilon
Note: This step was run on multiple assembly fragments generated in step 3b substituted for the --genome file split at fasta record boundaries and combined at the end to improve runtime. Command below represents example.

```
module load java
module load pilon/1.22

java -Xmx200g -jar /isg/shared/apps/pilon/1.22/pilon-1.22.jar --genome consensus_split_000.fa 
    --frags BettyinputChIPReads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam \
    --frags BettyWGS075Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam \
    --frags BettyWGS076Reads_vs_HLE_FMPPJ3MT_2022JAN14_bwa_sorted.bam \
    --threads 24 --output polishedfrag_0 
```

#### d.) Combine Polished Genome Fragments to Generate Polished Assembly
```
cat polishedfrag_*fa > HLE_FMPPJ3MTP.fasta
```

| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 95.8% (8839)             | 95.4% (13151)             |
| *Single copy*    | 91.4% (233)             | 93.5% (8627)             | 93.6% (12903)             |
| *Multi copy*     | 8.2% (21)               | 2.3% (212)               | 1.8% (248)                |
| Fragmented       | 0.0% (0)                | 1.3% (121)               | 1.2% (164)                |
| Missing          | 0.4% (1)                | 2.9% (266)               | 3.4% (465)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    324        | No. contigs                 |    324        |
| No. contigs (>= 1000 bp)    |    324        | Largest contig              |    224726311  |
| No. contigs (>= 5000 bp)    |    293        | Total length                |    2764447495 |
| No. contigs (>= 10000 bp)   |    241        | GC (%)                      |    41.28      |  
| No. contigs (>= 25000 bp)   |    157        | N50                         |    159711632  |
| No. contigs (>= 50000 bp    |    93         | N75                         |    108834272  |
| Total length (>= 0 bp)      |    2764447495 | L50                         |    8          |
| Total length (>= 1000 bp)   |    2764447495 | L75                         |    13         |
| Total length (>= 5000 bp)   |    2764324408 | No. N's per 100 kbp         |    12.80      |
| Total length (>= 10000 bp)  |    2763949195 |
| Total length (>= 25000 bp)  |    2762486006 |
| Total length (>= 50000 bp)  |    2760160693 |

### 12. TGS-GapCloser - Round 2
To reduce Ns, gaps filling was repeated with TGS-GapCloser. Must use script that Patrick Grady (patrick.gs.grady@uconn.edu) updated that makes ```--tgstype pb``` actually work for ONT reads. (ie: ```-   TGS reads type is pb . MINIMAP2_PARAM is  -x ava-ont --score-N 3   MIN_IDY is 0.2 . MIN_MATCH is 200 .```)

```
tgsgapcloser --scaff ./HLE_FMPPJ3MTP.fasta \
        --tgstype pb \
        --reads /core/projects/EBP/Oneill/reads/nanopore/promethion/gibbon/Combined_Reads_Super_Accurate/HLE_Super_Accurate_ONT_combined.fasta \
        --output HLE_FMPPJ3M_3kblimit_TGS-GapCloser \
        --ne \
	--minmap-arg \'--score-N 3\' \
        --thread 32
```

| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.6% (254)             | 95.8% (8836)             | 95.5% (13151)             |
| *Single copy*    | 91.4% (233)             | 93.5% (8623)             | 93.7% (12906)             |
| *Multi copy*     | 8.2% (21)               | 2.3% (213)               | 1.8% (245)                |
| Fragmented       | 0.0% (0)                | 1.3% (122)               | 1.2% (164)                |
| Missing          | 0.4% (1)                | 2.9% (268)               | 3.3% (465)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    324        | No. contigs                 |    324        |
| No. contigs (>= 1000 bp)    |    324        | Largest contig              |    224752805  |
| No. contigs (>= 5000 bp)    |    293        | Total length                |    2764783899 |
| No. contigs (>= 10000 bp)   |    241        | GC (%)                      |    41.28      |  
| No. contigs (>= 25000 bp)   |    157        | N50                         |    159711632  |
| No. contigs (>= 50000 bp    |    93         | N75                         |    108834272  |
| Total length (>= 0 bp)      |    2764783899 | L50                         |    8          |
| Total length (>= 1000 bp)   |    2764783899 | L75                         |    13         |
| Total length (>= 5000 bp)   |    2764660812 | No. N's per 100 kbp         |    11.65      |
| Total length (>= 10000 bp)  |    2764285599 |
| Total length (>= 25000 bp)  |    2762822410 |
| Total length (>= 50000 bp)  |    2760497097 |

### 12. Get Chromosome Scaffolds


```
module load seqtk/1.2

seqtk seq -L 30000000 HLE_FMPPJ3MTPT.fasta > HLE_FMPPJ3MTP_Chrs.fasta
```

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    19         | No. contigs                 |    19         |
| No. contigs (>= 1000 bp)    |    19         | Largest contig              |    224752805  |
| No. contigs (>= 5000 bp)    |    19         | Total length                |    2721536245 |
| No. contigs (>= 10000 bp)   |    19         | GC (%)                      |    41.13      |  
| No. contigs (>= 25000 bp)   |    19         | N50                         |    159711632  |
| No. contigs (>= 50000 bp    |    19         | N75                         |    101836356  |
| Total length (>= 0 bp)      |    2721536245 | L50                         |    8          |
| Total length (>= 1000 bp)   |    2721536245 | L75                         |    13         |
| Total length (>= 5000 bp)   |    2721536245 | No. N's per 100 kbp         |    11.12      |
| Total length (>= 10000 bp)  |    2721536245 |
| Total length (>= 25000 bp)  |    2721536245 |
| Total length (>= 50000 bp)  |    2721536245 |

### 13. RagTag and Final Polishing

To reduce any misassemblies associated with manual curation and scaffolding, pre-scaffolded contigs were used to fill the final assembly with RagTag. Only chromosome scaffolds were taken forward for polishing. File name continued to be updated, **T**GS-GapCloser, **P**ilon, **Chr**omosomes. The final assembly was then polished once more to reduce errors.

#### a.) RagTag Patching
```
ragtag.py patch HLE_FMPPJ3MTP_Chrs.fasta HLE_FMPPJ3MTP_Chrs_R.fasta
```

#### b.) Align and Process PCR Free Illumina Reads 

##### i.) Index Nanopolished Genome 

```
module load bwa/0.7.17

bwa index -p  HLE_FMPPJ3MTP_Chrs_R HLE_FMPPJ3MTP_Chrs_R.fasta
```

##### ii.) Map and Process PCR Free Illumina Reads

```
module load bwa/0.7.17

bwa mem -t 16 HLE_FMPPJ3MTP_Chrs_R HLE_PCRFree_R_1_CA-q20-m50_Filtered.fastq HLE_PCRFree_R_2_CA-q20-m50_Filtered.fastq > BettyPCRFree_v_HLE_FMPPJ3MTP_Chrs_R.sam
```
##### iii.) Convert SAM to BAM
```
module load samtools/1.7

samtools view -Sb BettyPCRFree_v_HLE_FMPPJ3MTP_Chrs_R.sam > BettyPCRFree_v_HLE_FMPPJ3MTP_Chrs_R.bam
```

##### v.) Index BAM Files
```
module load samtools/1.7

samtools index BettyPCRFree_v_HLE_FMPPJ3MTP_Chrs_R.bam
```

##### vi.) Get Mapping Statistics

```
module load samtools/1.7
samtools flagstat BettyPCRFree_v_HLE_FMPPJ3MTP_Chrs_R.bam > BettyPCRFree_v_HLE_FMPPJ3MTP_Chrs_R_STATS.txt
```

| BWA Stats v. 0.7.17   | PCR Free     | 
| -------------         | ---          | 
| Total Reads           | 1050836244   | 
| Mapped                | 99.45%       | 

#### c.) Split Assembly Into Fragments 
The assembly was split into smaller fragments broken at fasta record boundaries in order to improve Pilon runtime from several days to several hours. 
```
module load GenomeBrowser/20180626

faSplit sequence HLE_FMPPJ3MTP_Chrs_R.fasta 20 consensus_split_
```

#### d.) Polish Assembly Fragments Using Pilon
Note: This step was run on multiple assembly fragments generated in step 3b substituted for the --genome file split at fasta record boundaries and combined at the end to improve runtime. Command below represents example.

```
module load java
module load pilon/1.22

java -Xmx200g -jar /isg/shared/apps/pilon/1.22/pilon-1.22.jar --genome consensus_split_000.fa 
    --frags BettyPCRFree_v_HLE_FMPPJ3MTP_Chrs_R.bam \
    --threads 24 --output polishedfrag_0 
```

#### e.) Combine Polished Genome Fragments to Generate Complete Polished Assembly
```
cat polishedfrag_*fa > HLE_Ragtag_Patch.FINAL.fasta
```
| BUSCO v. 5.0.0   | eukaryota_odb10 (n=255) | mammalia_odb10 (n=9226)  | primates_odb10 (n=13780)  |
| -------------    | -------------           | -------------            | -------------             |
| Complete         | 99.7% (254)             | 95.9% (8843)             | 95.5% (13162)             |
| *Single copy*    | 92.2% (235)             | 93.5% (8623)             | 93.7% (12915)             |
| *Multi copy*     | 7.5% (19)               | 2.4% (220)               | 1.8% (247)                |
| Fragmented       | 0.0% (0)                | 1.3% (117)               | 1.2% (160)                |
| Missing          | 0.3% (1)                | 2.8% (266)               | 3.3% (458)                |

| Quast v. 5.0.2              |    Stats      | Stats for contigs > 500 bp  |               |
| ---                         | ---           | ---                         | ---           |
| No. contigs (>= 0 bp)       |    19         | No. contigs                 |    19.        |
| No. contigs (>= 1000 bp)    |    19         | Largest contig              |    226217378  |
| No. contigs (>= 5000 bp)    |    19         | Total length                |    2760609521 |
| No. contigs (>= 10000 bp)   |    19         | GC (%)                      |    41.25      |  
| No. contigs (>= 25000 bp)   |    19         | N50                         |    159996863  |
| No. contigs (>= 50000 bp    |    19         | N75                         |    115684913  |
| Total length (>= 0 bp)      |    2760609521 | L50                         |    8          |
| Total length (>= 1000 bp)   |    2760609521 | L75                         |    13         |
| Total length (>= 5000 bp)   |    2760609521 | No. N's per 100 kbp         |    0.22       |
| Total length (>= 10000 bp)  |    2760609521 |
| Total length (>= 25000 bp)  |    2760609521 |
| Total length (>= 50000 bp)  |    2760609521 |
