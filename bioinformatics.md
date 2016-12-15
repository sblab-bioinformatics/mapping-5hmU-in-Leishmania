<!-- MarkdownTOC -->

- [Software tools](#software-tools)
- [Preparing reference genomes and annotation files](#preparing-reference-genomes-and-annotation-files)
    - [Reference genomes](#reference-genomes)
    - [Annotation files](#annotation-files)
- [Read trimming and alignment](#read-trimming-and-alignment)
- [Detection of read enrichment \(peak calling\)](#detection-of-read-enrichment-peak-calling)
    - [Jaccard index between replicates](#jaccard-index-between-replicates)
    - [Consensus regions of read enrichment](#consensus-regions-of-read-enrichment)
- [Overlap between modifications and genomic regions](#overlap-between-modifications-and-genomic-regions)
    - [Enrichment of Base J and 5hmU in genomic regions](#enrichment-of-base-j-and-5hmu-in-genomic-regions)
- [RNA expression in Base J and 5hmU peaks](#rna-expression-in-base-j-and-5hmu-peaks)
    - [Profiles of expression](#profiles-of-expression)
- [Motif analysis](#motif-analysis)
- [Peak files](#peak-files)

<!-- /MarkdownTOC -->

# Software tools

Data processing was performed under Linux environment with GNU coreutils tools. The following software and data files are required for the analysis described here:

* [cutadapt](http://cutadapt.readthedocs.org/en/stable/guide.html) version 1.8
* [bwa](https://github.com/lh3/bwa) version 0.7
* [samtools](http://www.htslib.org/) version 1.1
* [MACS2](https://github.com/taoliu/MACS/) version 2.1.0
* [bedtools](http://bedtools.readthedocs.org/en/latest/) version 2.25
* [GAT](http://gat.readthedocs.io/en/latest/)
* [deepTools](https://github.com/fidelram/deepTools) version 2.1.0
* [bedGraphToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
* [dreme](http://meme-suite.org/index.html)
* [R-3.2.3](https://cran.r-project.org/) 

See also [scripts](scripts) directory.

# Preparing reference genomes and annotation files

## Reference genomes

Reference fasta sequences have been downloaded from the Sanger Institute repositories:

```
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Leishmania/major/Current/LmjF_v6.1_20131105/fasta/LmjF_v6.1_all_20131105.fa
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Leishmania/infantum/V5210211/Linfantum.fa
```

The control spike-in sequences were appended to the reference genomes and the combined 
fasta indexed by bwa (`bwa index` command with default arguments).

```
>fk_hmu_pos_01
TTCTTGGCTGTGGCTCTGCGTCCTTGTCCTGAGGCCAHCACAGCGCAHGAACGACGAGGCACAACAGAGAGCAACACCGCCGAGGA
>fk_hmu_pos_02
ATCGAGAATCCCGGTGCCGAGCTACACCHACTCTTTGHAGAATTAAGTCTCCAGGCACGTGTCAGATATATACATCCGAT
>fk_hmu_neg_01
GCTCGCTTTGTTGGTTTCCTTGTTCTCTGTGAGGCCATCACAGCGCATGAACGACGAAAGCAGCGCGAGCAAGCGAGACAGGACAC
```

## Annotation files

Annotation files in GFF format for *L. donovani* and *L. major* were downloaded 
from the Sanger Institute repositories and they were complemented to include strand switch
regions, telomeric regions, and *intergenic* regions, *i.e.* the portion of the genome not covered by any 
known feature. The script below refers to *L. donovani*. a similar procedure was applied 
to *L. major*:

```
## Download annotations
# L. major
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Leishmania/major/Current/LmjF_v6.1_20131105/gff/LmjF_v6.1_20131105_all.gff
# L. donovani
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Leishmania/infantum/current_gff3/Linfantum.gff3.gz &&
gunzip Linfantum.gff3.gz


gff='Linfantum.gff3'
genome=Linfantum.fa.fai ## Chromosome sizes from fasta index file

## Prepare telomeric regions
## ----------------------------
awk '{print $1"\t1\t5000"}' $genome > /tmp/start.bed
awk '{print $1"\t"$2-5000"\t"$2}' $genome > /tmp/end.bed
sort -k1,1 -k2,2n /tmp/start.bed /tmp/end.bed > tmp.telo_5kb.bed &&
rm /tmp/start.bed /tmp/end.bed

## Start building GFF
## ==================
## Remove chroms not in fasta ref
sortBed -i $gff | grep -v -P '^LinJ.\d\d_.+?\t' > annotation.gff 

## Append Switch regions
awk -v OFS='\t' '{print $1, "chado", "switch_region", $2, $3, ".", $4, ".", "ID=NA"}' \
Linfantum.mRNA.switch_reg.bed >> annotation.gff

## Append telomeric regions
awk -v OFS='\t' '{print $1, "chado", "telo", $2, $3, ".", ".", ".", "ID=NA"}' tmp.telo_5kb.bed \
>> annotation.gff && 
rm tmp.telo_5kb.bed

## Infer intergenic regions
complementBed -i <(sortBed -i annotation.gff) -g <(cut -f1,2 $genome | sort -k1,1 -k2,2n) > /tmp/unkn.bed
awk -v OFS='\t' '{print $1, "chado", "unknown", $2, $3, ".", ".", ".", "ID=NA"}' /tmp/unkn.bed >> annotation.gff && 
rm /tmp/unkn.bed
```

Strand switch regions were prepared as follows:

```
genome=Linfantum.fa.fai
gff=Linfantum.gff3

grep -P '\tmRNA\t' $gff \
| grep -v -P '^##' \
| sort -k1,1 -k4,4n -k5,5n \
| awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $7}' \
| groupBy -g 1,4 -c 2,3 -o min,max \
| awk 'BEGIN{OFS="\t"}{print $1, $3, $4, $2}' > Linfantum.mRNA.strand_switch.bed

complementBed -g $genome \
    -i Linfantum.mRNA.strand_switch.bed \
| slopBed -b 1 -g $genome > /tmp/interswitch.bed

intersectBed -a /tmp/interswitch.bed -b Linfantum.mRNA.strand_switch.bed -wa -wb \
| sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n \
| groupBy -g 1,2,3 -c 7 -o collapse \
| slopBed -b -1 -g $genome \
| awk 'BEGIN{OFS="\t"}{if ($2 == 1) {$2=0}; sub(/,./, "", $4); print $0, ".", $4}' > /tmp/Linfantum.mRNA.intergenic_switch.bed
rm /tmp/interswitch.bed

awk '$2!=0 {print $0}' /tmp/Linfantum.mRNA.intergenic_switch.bed > /tmp/Linfantum.mRNA.switch_reg.bed
rm /tmp/Linfantum.mRNA.intergenic_switch.bed

## Remove distal telomeric regions
groupBy -g 1 -c 2,3 -o max,max -i /tmp/Linfantum.mRNA.switch_reg.bed \
| intersectBed -v -a /tmp/Linfantum.mRNA.switch_reg.bed -b - > Linfantum.mRNA.switch_reg.bed
rm /tmp/Linfantum.mRNA.switch_reg.bed
```

# Read trimming and alignment

Prior to alignment, fastq reads have been trimmed to remove adapter contamination 
using cutadapt and aligned to the reference genome with bwa mem:

```
cutadapt -m 3 -a AGATCGGAAGAGC $fq \
| bwa mem -t 5 -M $ref - \
| samtools view -F4 -u - \
| samtools sort  -@ 10 - $bname &&
samtools index ${bname}.bam" 
```

`$fq` is a variable for the fastq file name, `$ref` is the reference genome prepared as 
above, and `$bname` a the basename for the output bam file.

Paired-end RNA-Seq reads were processed in a similar way. In addition to the Illumina
sequence adapters, the splice leader sequence (CAGTTTCTGTACTTTATTG) was also trimmed:

```
cutadapt -a AGATCGGAAGAGC -b CAGTTTCTGTACTTTATTG -A AGATCGGAAGAGC -B CAGTTTCTGTACTTTATTG \
-o fastq_trimmed/${bname}.r1.fq.gz \
-p fastq_trimmed/${bname}.r2.fq.gz \
$fq1 $fq2

bwa mem -M -t 4 $ref fastq_trimmed/${bname}.r1.fq.gz fastq_trimmed/${bname}.r2.fq.gz \
| samtools sort -@8 -o - -O bam -T bam_clean/${bname}.tmp.bam > bam/${bname}.bam &&
samtools index bam/${bname}.bam" 
```

# Detection of read enrichment (peak calling)

Input bam files (`$inbam`) were filtered to remove reads mapping to spike-in controls 
and reads with low mapping quality:

```
samtools view -h -q 10 $inbam \
| grep -P -v '\tfk_hmu_.*' \
| samtools view -Sb - > $bam
```

Peaks of read enrichment have been detected using MACS2 as follows:

```
macs2 callpeak --nomodel -t $bam -c $ctrl --keep-dup all -g 30e6 -n $bname 
```

where `$ctrl` is the appropriate control bam file.

## Jaccard index between replicates

The jaccard index between pairs of peak files was computed using the `jaccard`
command in bedtools. 

## Consensus regions of read enrichment

The peak files generated by MACS2 were merged within replicates and only regions
overlapped by 2 or more replicates were retained as *consensus* enrichment regions:

```
## 5hmU - L. major
mergePeaks.sh fk050_oxhyd1.oxNPD_peaks.narrowPeak \
              fk051_oxhyd4.oxNPD_peaks.narrowPeak \
              fk054_Chem5.oxNPD_peaks.narrowPeak \
| awk '$5 > 1' fk_chem_hmu_lmaj.mrg.bed > fk_chem_hmu_lmaj.cns.bed

## 5hmU - L. donovani
mergePeaks.sh fk066_Ldono_chem1_peaks.narrowPeak \
              fk067_Ldono_chem2_peaks.narrowPeak \
| awk '$5 > 1' > fk_chem_hmu_ldono.cns.bed

## Base J - L. major
mergePeaks.sh fk113_B2_J-BEADS.160218.lmaj_peaks.narrowPeak \
              fk116_B3_J_BEADS.160218.lmaj_peaks.narrowPeak \
              fk117_B3_J_BEADS_2.160218.lmaj_peaks.narrowPeak \
| awk '$5 > 1' > fk_baseJ.160218.cns.bed
```

`mergePeaks.sh` can be found in the [scripts](scripts) directory.

# Overlap between modifications and genomic regions

Consensus peak regions were annotated by intersecting them with the genomic features
prepared above. The significance of the association between peak regions and genomic 
features was assessed by simulation by means of the GAT software.

First, a *workspace* was defined as the genomic space were read coverage is possible. 
To this end all the regions covered by at least one read in the control input files
were extracted. Input bam files were filtered in the same way as for MACS2 peak calling. 
Files `*.chromSize.txt` list chromosome names and chromosome sizes:


```
## Workspace
# L./ major
genomeCoverageBed -ibam fk108_Leish_batch3_NPD.160115.lmaj.bam \
    -g LmjF_v6.1_all_20131105.chromSize.txt -bga \
| grep -P '\t0$' \
| awk -v OFS='\t' '{print $0, $3-$2}' \
| complementBed -i - -g LmjF_v6.1_all_20131105.chromSize.txt > lmaj.workspace.bed

# L. donovani
genomeCoverageBed -ibam fk070_Ldono_oxNPD.ldon.bam \
    -g Linfantum.chromSize.txt -bga \
| grep -P '\t0$' \
| awk -v OFS='\t' '{print $0, $3-$2}' \
| sortBed \
| complementBed -i - -g Linfantum.chromSize.txt > ldon.workspace.bed
```

GAT was run as follows:
 
```
awk -v OFS='\t' '$0 ~ /^[^#]/ {print $1, $4-1, $5, $3}' ./annotation.gff > ldon.annotation.bed
gat-run.py -s fk_chem_hmu_ldono.cns.bed -a ldon.annotation.bed \
    -w ldon.workspace.bed -n 10000 -t 4 | grep -v '#' > fk_chem_hmu_ldono.gat.tmp

awk -v OFS='\t' '$0 ~ /^[^#]/ {print $1, $4-1, $5, $3}' ../20160301_annotation_pqs/annotation.gff > lmaj.annotation.bed
gat-run.py -s fk_chem_hmu_lmaj.cns.bed -a lmaj.annotation.bed \
    -w lmaj.workspace.bed -n 10000 -t 4 | grep -v '#' > fk_chem_hmu_lmaj.gat.tmp
```

From the output of GAT the column *percent_overlap_size_track* was used to plot 
the percentage peak overlapped by genomic features. 

## Enrichment of Base J and 5hmU in genomic regions

Genomic regions were partitioned among those containing only Base J, only 5hmU, or 
both modifications as in the representation below:

```
---------       Base J peak
     ---------  5hmU peak
Resulting regions:
-----           Only Base J
     ----       Both
         -----  Only 5hmU
```

`partitionBed.py` can be found in the [scripts](scripts) directory.

```
partitionBed.py fk_baseJ.160218.cns.bed fk_chem_hmu_lmaj.cns.bed \
| awk -v OFS="\t" '{if($NF == "a"){x="baseJ"} 
                    else if($NF == "b"){x= "5hmU"} 
                    else if($NF == "ab"){x= "baseJ_5hmU"} 
                    else {exit 1} 
                    print $1, $2, $3, x}' | sort -k1,1 -k2,2n > fk_chem_lmaj_baseJ-hmu.prt.bed
```

Then the three partitions (Base J only, 5hmU only, both) were assessed for enrichment
in different genomic regions: 

```
for x in `cut -f 4 fk_chem_lmaj_baseJ-hmu.prt.bed | sort | uniq`
do
    grep -w $x fk_chem_lmaj_baseJ-hmu.prt.bed > gat.${x}.tmp
    gat-run.py -s gat.${x}.tmp -a lmaj.annotation.bed -w lmaj.workspace.bed -n 10000 -t 4 | grep -v '#' > gat.${x}.txt.tmp
done
```

From the resulting output files the column *l2fold* was used to plot the enrichment while the 
confidence intervals for log2 fold enrichment was calculated as 
*l2foldCI95high= log2((observed+1) / (CI95high+1))* and *l2foldCI95low := log2((observed+1) / (CI95low+1))*

# RNA expression in Base J and 5hmU peaks

RNA expression in Base J and 5hmU peaks, partitioned as above, was measured by counting 
reads falling in each region. Since the three RNA-Seq libraries were very similar
to each other in terms of profiles of coverage, they were merged in a single file: 

```
samtools merge - \
    fk094_F6-51_rep1.bam \
    fk095_F6-51_rep2.bam \
    fk096_F6-51_rep3.bam \
| samtools view -b -F4 -q 5 -@ 12 - > fk_F6-51_rna.bam

samtools view -u -F 128 fk_F6-51_rna.bam \
| coverageBed -counts -a fk_chem_lmaj_baseJ-hmu.prt.bed -b - > fk_F6-51_rna.fk_hmU_baseJ.prt.cnt.bed
```

The read count in a random set of genomic regions was used as background level of
expression:

```
bedtools random -g LmjF_v6.1_all_20131105.chromSize.txt -l 1000 -n 10000 > random.tmp.bed
samtools view -u -F 128 fk_F6-51_rna.bam \
| coverageBed -counts -a random.tmp.bed -b - > fk_F6-51_rna.random.tmp.bed
```

Finally, the violin plot were produced with the following code:

```
R
library(ggplot2)
library(data.table)
mod<- fread('fk_F6-51_rna.fk_hmU_baseJ.prt.cnt.bed')
setnames(mod, names(mod), c('chrom', 'start', 'end', 'mod', 'nreads'))

rnd<- fread('fk_F6-51_rna.random.tmp.bed')
setnames(rnd, names(rnd), c('chrom', 'start', 'end', 'name', 'len', 'strand', 'nreads'))

libsize<- as.numeric(system('samtools view -c -F128 fk_F6-51_rna.bam', intern= TRUE)) # 2522177

mod[, mod := ifelse(mod == 'baseJ_5hmU', 'Base J and 5hmU', mod)]
mod[, mod := ifelse(mod == 'baseJ', 'Base J only', mod)]
mod[, mod := ifelse(mod == '5hmU', '5hmU only', mod)]

mod[, rpkm := log10((nreads+1) / (libsize/1e6 * (end-start)))]
rnd[, rpkm := log10((nreads+1) / (libsize/1e6 * (end-start)))]
rndavg<- mean(rnd$rpkm)
rndsd<- sd(rnd$rpkm)

gg<- ggplot(data= mod, aes(x= factor(mod, levels= c('Base J and 5hmU', 'Base J only', '5hmU only')), y= rpkm)) +
    geom_rect(xmin= -10, xmax= 10, ymin= rndavg-rndsd, ymax= rndavg+rndsd, colour= NA, fill= 'orange') +
    geom_hline(yintercept= rndavg, colour= 'grey20', size= 0.5, linetype= 'dashed')+
    geom_violin(draw_quantiles= 0.5, fill= 'white', alpha= 1) +
    geom_violin(draw_quantiles= 0.5, fill= 'blue', alpha= 0.5) +
    geom_jitter(colour= 'red', alpha= 0.2, width= 0.1) + 
    xlab('') + ylab('log10(RPKM)') + ggtitle('RNA expression in modification peaks\n(shaded: mean background expr. +/- 1 sd)')
ggsave('viol_rna_mods.prt.pdf', w= 12/2.54, h= 12/2.54)


## Significance of difference between groups

dat<- rbind(mod[, list(mod, rpkm)], rnd[, list(mod= 'Background', rpkm)])
xlm<- aov(rpkm ~ mod, data= dat)
hsd<- TukeyHSD(xlm)

#                             diff        lwr        upr    p adj
# Background-5hmU       -0.1398003 -0.3223368  0.0427361 0.200277
# baseJ-5hmU            -0.4846197 -0.6816427 -0.2875968 0.000000
# baseJ_5hmU-5hmU       -0.9692186 -1.1840549 -0.7543823 0.000000
# baseJ-Background      -0.3448194 -0.4211308 -0.2685080 0.000000
# baseJ_5hmU-Background -0.8294183 -0.9441355 -0.7147010 0.000000
# baseJ_5hmU-baseJ      -0.4845989 -0.6211951 -0.3480026 0.000000

par(las= 1, mar= c(3, 15, 2, 1), mgp= c(2, 0.5, 0))
plot(hsd)
```

## Profiles of expression

The average profiles of expression in and around base J and 5hmU modifications were
visualized using deepTools as follows:

```
## Convert BAM to bigWig
bamCoverage -p 4 --bam fk_F6-51_rna.bam -o fk_F6-51_rna.bw --normalizeUsingRPKM  --binSize 10 -of bigwig

## Random regions
bedtools random -g LmjF_v6.1_all_20131105.chromSize.txt -l 250 -n 1000 > random.tmp.bed

## Extract regions: Base J only, 5hmU only, both
for x in `cut -f 4 fk_chem_lmaj_baseJ-hmu.prt.bed | sort | uniq`
do
    grep -w $x fk_chem_lmaj_baseJ-hmu.prt.bed > deep.${x}.tmp
done

computeMatrix scale-regions \
             -R deep.*.tmp random.tmp.bed \
             -S fk_F6-51_rna.bw \
             -out tmp.mat.gz \
             -b 2000 \
             -a 2000 \
             -bs 10 \
             --skipZeros
plotProfile --matrixFile tmp.mat.gz -out profiles-rnaseq.prt.pdf --startLabel 'Start' --endLabel 'End' -T '' 
```

Similarly for GC content profiles:


```
## Prepare GC content profile along the genome
windowMaker -g LmjF_v6.1_all_20131105.chromSize.txt -w 10 \
| nucBed -fi LmjF_v6.1_all_20131105.fa -bed - \
| tail -n+2 \
| cut -f 1-3,5 > LmjF_v6.1_all_20131105.gc.bedGraph

bedGraphToBigWig LmjF_v6.1_all_20131105.gc.bedGraph \
    /data/sblab-data/common/reference_data/genomes/leishmania_major/LmjF_v6.1_all_20131105.chromSize.txt \
    LmjF_v6.1_all_20131105.gc.bw
rm LmjF_v6.1_all_20131105.gc.bedGraph

computeMatrix scale-regions \
             -R deep.*.tmp random.tmp.bed \
             -S LmjF_v6.1_all_20131105.gc.bw \
             -out tmp.mat.gz \
             -b 2000 \
             -a 2000 \
             -bs 25 \
             --skipZeros
plotProfile --matrixFile tmp.mat.gz -out profiles-pct_gc.prt.pdf --startLabel 'Start' --endLabel 'End' -T 'GC content profile' 
```

# Motif analysis

Motifs enriched in Base J and 5hmU sites only regions or in regions occupied by both modifications
were detected using the dreme software.

Telomeric regions were removed from input regions before testing:

```
for x in `cut -f 4 fk_chem_lmaj_baseJ-hmu.prt.bed | sort | uniq`
do
    grep -w $x fk_chem_lmaj_baseJ-hmu.prt.bed > ${x}.tmp
    intersectBed -f 0.1 -v -a ${x}.tmp -b <(grep 'PQS' lmaj.annotation.bed) > ${x}.tmp.notelo.bed
    rm ${x}.tmp
done
```

Enrichment was tested against random genomic regions produced by shuffling the test
regions in the mappable portion of genome (file `lmaj.workspace.bed` as produced above).

```
j=1
for bed in 5hmU.tmp.notelo.bed \
           baseJ_5hmU.tmp.notelo.bed \
           baseJ.tmp.notelo.bed
do 
    bname=`basename $bed .tmp.notelo.bed`
    >$bname.shuffle.tmp.bed
    for i in {1..10}
    do
        echo $j
        bedtools shuffle -seed $j -incl lmaj.workspace.bed -noOverlapping -i $bed \
            -g LmjF_v6.1_all_20131105.chromSize.txt >> $bname.shuffle.tmp.bed
        j=`echo $j + 1 | bc` 
    done
    fastaFromBed -fi LmjF_v6.1_all_20131105.fa -bed $bed > $bed.fa
    fastaFromBed -fi LmjF_v6.1_all_20131105.fa -bed $bname.shuffle.tmp.bed > $bname.shuffle.tmp.fa

    echo "dreme -maxk 32 -m 5 -oc dreme_prt_$bname -p $bed.fa -n $bname.shuffle.tmp.fa && 
          rm $bed.fa $bname.shuffle.tmp.fa $bname.shuffle.tmp.bed" > $bname.dreme.tmp.sh
done
ls *.dreme.tmp.sh | xargs -P 0 -n 1 bash 
```

# Peak files

Summary table of peak files

filename | species | pull-down | method | input | n_peaks
-------- | ------- | --------- | ------ | ----- | --------
fk066_Ldono_chem1.oxNPD | L. donovani | 5hmU | Chem | fk070_Ldono_oxNPD.ldon.bam | 132
fk067_Ldono_chem2.oxNPD | L. donovani | 5hmU | Chem | fk070_Ldono_oxNPD.ldon.bam | 126
fk068_Ldono_fU1.oxNPD | L. donovani | control | Chem | fk070_Ldono_oxNPD.ldon.bam | 1
fk069_Ldono_fU2.oxNPD | L. donovani | control | Chem | fk070_Ldono_oxNPD.ldon.bam | 1
fk047_LdonDIP1.oxNPD | L. donovani | 5hmU | DIP | fk070_Ldono_oxNPD.ldon.bam | 2087
fk048_LdonDIP2.oxNPD | L. donovani | 5hmU | DIP | fk070_Ldono_oxNPD.ldon.bam | 2609
fk050_oxhyd1.oxNPD | L. major | 5hmU | Chem | fk059_ox_NPD_spiked.bam | 206
fk051_oxhyd4.oxNPD | L. major | 5hmU | Chem | fk059_ox_NPD_spiked.bam | 188
fk054_Chem5.oxNPD | L. major | 5hmU | Chem | fk059_ox_NPD_spiked.bam | 145
fk113_B2_J-BEADS.160218.lmaj | L. major | base J | Chem | fk108_Leish_batch3_NPD.160115.lmaj.bam | 188
fk116_B3_J_BEADS.160218.lmaj | L. major | base J | Chem | fk108_Leish_batch3_NPD.160115.lmaj.bam | 198
fk117_B3_J_BEADS_2.160218.lmaj | L. major | base J | Chem | fk108_Leish_batch3_NPD.160115.lmaj.bam | 191
fk045_hyd1.oxNPD | L. major | control | Chem | fk059_ox_NPD_spiked.bam | 6
fk046_hyd2.oxNPD | L. major | control | Chem | fk059_ox_NPD_spiked.bam | 5
fk056_fu3.oxNPD | L. major | control | Chem | fk059_ox_NPD_spiked.bam | 5
fk057_fu4.oxNPD | L. major | control | Chem | fk059_ox_NPD_spiked.bam | 5
fk041_F5_10_DIP1.oxNPD | L. major | 5hmU | DIP | fk059_ox_NPD_spiked.bam | 1175
fk043_F5_18_DIP2.oxNPD | L. major | 5hmU | DIP | fk059_ox_NPD_spiked.bam | 237
fk052_DIP4.oxNPD | L. major | 5hmU | DIP | fk059_ox_NPD_spiked.bam | 47
fk059_ox_NPD_spiked | L. major | input | na | fk058_NPD_spiked.clean.bam | 3

