# testing_featurecounts

In the internal RNAseq pipeline, reads are mapped to the reference transcriptome by STAR/hisat/Salmon. Gene expression levels are then calculated by featureCounts. However, I have noted that, in the way the tool is implemented now, the reads are counted as they were SE, even when the input data is PE.

At the current stage, the pipeline does check whether the inut data is SE or PE, and when PE, adds ```-p``` to the featureCounts option - as per below:
```
featureCounts -B -C -s 0 -p -T $cpu -a $gtf -o $out $bam
```

However, as per featureCounts manual, when the reads are PE ```-p --countReadPairs``` should be added, and not only ```-p```!!

To define how badly this is influencing the final gene expression levels, I have here run two analyses with/without ```--countReadPairs``` and compared the results.

## Analysis

```
module load subread
module load r
bam=/cluster/work/nme/data/fansaloni/projs/emanuel/EMT_A549_MCF10a_cells/rnaseq/SEQ00221/01_genexpression/output/aligned/bam/S221_10_hsMCF10a_EMT_RNA-1_GRCh38.p14_star.Aligned.out.sorted.bam
gtf=/cluster/work/nme/data/fansaloni/projs/pipeline_issues/rnaseq/01.featurecounts/test_ENSG00000116151/star/Homo_sapiens_GRCh38_p14.gtf


# featureCounts SE
if [ ! -d featureCounts_SE ]; then mkdir featureCounts_SE; fi
cd featureCounts_SE
featureCounts -R BAM -B -C -s 0 -p -T $cpu -a $gtf -o res.txt $bam
cd ../


# featureCounts PE
if [ ! -d featureCounts_PE ]; then mkdir featureCounts_PE; fi
cd featureCounts_PE
featureCounts -R BAM -B -C -s 0 -p --countReadPairs -T $cpu -a $gtf -o res.txt $bam
cd ../


# test differences
cat featureCounts_SE/res.txt | grep -v "#" | cut -f 1,7 | grep ENS | sort -k1,1 > featureCounts_SE/short_SE.txt
cat featureCounts_PE/res.txt | grep -v "#" | cut -f 1,7 | grep ENS | sort -k1,1 > featureCounts_PE/short_PE.txt
Rscript merge_corr.R

```

## Results

The results are quite good as they highlight very high correlation between the two methods. Of course, the gene expression levels calculated by the SE method are ~2x of the PE as the two reads of the same pair are counted twice in the SE analysis and only once in the PE analysis. However, this does not seem to alter significantly the overall gene expression levels.

<img src="https://github.com/fansalon/testing_featurecounts/blob/main/rnaseq_cor_res.png" width="750" height="750"/>


As it appears from the plot above, there are some genes which are quantified differently by the two methods. The number of such genes is not incredibly high:
