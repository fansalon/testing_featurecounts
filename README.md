# testing_featurecounts

In the internal RNAseq pipeline, reads are mapped to the reference transcriptome by STAR/hisat/Salmon. Gene expression levels are then calculated by featureCounts. However, I have noted that, in the way the tool is implemented now, the reads are counted as they were SE, even when the input data is PE.

At the current stage, the pipeline does check whether the inut data is SE or PE, and when PE, adds ```-p``` to the featureCounts option - as per below:
```
featureCounts -B -C -s 0 -p -T $cpu -a $gtf -o $out $bam
```
