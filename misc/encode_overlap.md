<!-- MarkdownTOC -->

- [Overlap between ChIP-Seq replicates from Encode](#overlap-between-chip-seq-replicates-from-encode)
    - [Dataset](#dataset)
    - [Overlap](#overlap)
    - [Scripts](#scripts)

<!-- /MarkdownTOC -->

Overlap between ChIP-Seq replicates from Encode
===============================================

We want to have an overview of how consistent typical pull down experiments are
in detecting peaks of read enrichment. To this end we compute the overlap between
pairs of ChIP-Seq replicates from the Encode project. _Overlap_ is computed as number of base pairs
in common between replicates over the union of the two replicates.

Dataset
-------

Regions (peaks) of enrichments have been downloaded from Encode:

<!--
cd /nas/sblab_data1/group_folders/berald01/projects/20140818_fumi_hmu_pull_down/20161124_encode_overlap/
-->

```
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/*.broadPeak.gz ./

for rep1 in *Rep1.broadPeak.gz
do
    rep2=${rep1/Rep1/Rep2}
    if [ -f $rep2  ]
    then
        # Merge regions overlapping in the two replicates
        mergePeaks.sh $rep1 $rep2 > ${rep1%Rep1.broadPeak.gz}.mrg.bed 2> /dev/null
    else 
        echo "$rep2 NOT FOUND"
    fi
done
```

Overlap
-------

Here we compute the overlap within each ChIP-Seq experiment and plot the result.

```
R
library(data.table)
library(ggplot2)

## Concatenate merged files from above
mrg<- fread('tableCat.py -i *.mrg.bed -r ".mrg.bed"')
setnames(mrg, names(mrg), c('chrom', 'start', 'end', 'reps', 'n', 'chip'))
mrg[, len := end - start]
mrg[, cns_len := ifelse(n > 1, len, 0)]

ovl<- mrg[, list(union= sum(len), cns= sum(cns_len)), by= chip]

## Mean and stdev fo the overlap
mu<- mean(ovl[, 100 * cns/union])
std<- sd(ovl[, 100 * cns/union])

gg<- ggplot(data= ovl, aes(x= 100 * cns/union)) + 
    geom_histogram(color= 'white') +
    ggtitle(sprintf('Overlap between pairs of replicates from ENCODE\n(n= %s pairs, mean= %.2f, sd= %.2f)', nrow(ovl), mu, std)) +
    geom_vline(xintercept= 66, linetype= 'dotted', color= 'red') +
    xlab('% Overlap as [bp common]/[bp union]')
ggsave('encode_overlap.png', width= 12/2.54, height= 12/2.54)
```

<!-- 
system('rsync --remove-source-files encode_overlap.png 10.20.12.18:~/git_sblab/mapping-5hmU-in-Leishmania/trunk/misc/figures') 
-->

<img src=figures/encode_overlap.png width=600>

Scripts
-------

* [tableCat.py](https://github.com/dariober/bioinformatics-cafe/tree/master/tableCat)

* [mergePeaks.sh](../scripts/mergePeaks.sh)