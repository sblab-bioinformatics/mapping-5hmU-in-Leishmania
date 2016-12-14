<!-- MarkdownTOC -->

- [Comparison between MACS2 and SICER in 5hmU detection](#comparison-between-macs2-and-sicer-in-5hmu-detection)

<!-- /MarkdownTOC -->


Comparison between MACS2 and SICER in 5hmU detection
====================================================

The size and number of peaks detected in a pull down experiment can depend on
the  choice of peak caller. Consequently, downstream analyses and interpretation
may be dependant on the peak calling algorithm. Here, we compare our
peak caller of choice, MACS2, with an alternative peak calling program, SICER
(version 1.1, [PMID: 19505939](https://www.ncbi.nlm.nih.gov/pubmed/19505939)).
MACS2 ans SICER rely on substantially differ peak algorithms since SICER, in
contrast to MACS2, is designed for detecting broad and diffused peaks. 

We compare MACS2 and SICER in terms of peak number and peak size in the
detection of 5hmU in *Leishmania major* samples. Overall, we found that when the
output of SICER is filtered to remove peaks with low enrichment over the input,
it produces peaks similar to MACS2 in terms of size and number.

Here we execute SICER using th wrapper
[SICER.py](https://github.com/dariober/SICERpy) for convenience. The output was
filtered by the `awk` step to remove peaks supported by relatively few reads
(50) and with fold enrichment below 2.

```
for bam in fk050_oxhyd1.Lmaj.bam \
           fk051_oxhyd4.Lmaj.bam \
           fk054_Chem5.bam
do
     SICER.py --treatment $bam --control fk059_ox_NPD_spiked.bam --windowSize 50 \
     | awk '$4 > 50 && $7 > 2' > `basename $bam .bam`.sicer.bed &
done
```

The peaks detected by SICER were then merged across the three replicates and a
consensus peak set was created by retaining merged peaks supported by at least
two replicates. This step is identical to the one used for MACS2.

```
mergePeaks.sh fk050_oxhyd1.Lmaj.sicer.bed \
              fk051_oxhyd4.Lmaj.sicer.bed \
              fk054_Chem5.sicer.bed \
| awk '$5 > 1' > fk_chem_hmu_lmaj.sicer.cns.bed
```

Finally, we merged the consensus peak set from SICER with the consensus peak set 
from MACS2 previously prepared (`fk_chem_hmu_lmaj.cns.bed`):

```
mergePeaks.sh fk_chem_hmu_lmaj.sicer.cns.bed fk_chem_hmu_lmaj.cns.bed > sicer.macs.mrg.bed
```

The table below shows the number of consensus peaks present only in one peak set and in both peak
sets. Most of the peaks were detected by both methods:

| N. peaks | Peak set |
|----------|----------|
|  32      | MACS2    |
| 114      | Both     |
|   5      | SICER    | 

In order to compare the size of the peaks detected by the two methods, we intersected 
the MACS2 and SICER consensus peak sets and we extracted the width of the peaks:

```
intersectBed -wa -wb -a fk_chem_hmu_lmaj.cns.bed -b fk_chem_hmu_lmaj.sicer.cns.bed  \
| awk -v OFS='\t' '{print $3-$2,  $9-$8}' > sizes.txt

R
library(ggplot2)
library(data.table)
peaks<- fread('sizes.txt')
setnames(peaks, names(peaks), c('macs_width', 'sicer_width'))

gg<- ggplot(data= peaks, aes(x= macs_width, y= sicer_width)) + 
    geom_point(size= 0.5) +
    xlab('MACS2 width (bp)') +
    ylab('SICER width (bp)') +
    ggtitle('Peak size comparison between\nSICER and MACS2') +
    geom_abline(intercept= 0, slope= 1, color = 'blue', linetype= 'dotted')
ggsave('macs_vs_sicer_peak_size.png', width= 10, height= 10, units= 'cm')
```

<img src=figures/macs_vs_sicer_peak_size.png width= 500>

Th peak size is largely similar between MACS2 and SICER. SICER tends to give
wider  peaks as per its design. It should be noted that at least
some of the larger inconsistencies (data points on the top right of the plot)
are cases were MACS2 identified two separate peaks close to each other which SICER
merged together. It also useful to note that peak width tends to vary by approximately
one order of magnitude, from few hundred to few thousand kilo bases.
