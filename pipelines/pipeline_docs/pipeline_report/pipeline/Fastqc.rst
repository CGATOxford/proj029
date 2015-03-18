.. _reads:


=============================
Quality control on raw reads
=============================


Number of reads per sample
===========================

Below is a summary of the number of raw reads that were obtained for each sample.
for each track. The number of reads is the actual number of reads. The number of pairs is
thus half this number.

.. report:: Fastqc.ReadSummary
   :render: table

   summary of read numbers



Per base quality scores
===========================

To get an idea of the quality of the sequencing we look at the distribution
of sequence qualities at each base in a sequencing cycle. Below is a representative
plot from this analysis - All samples looked very similar and so are not displayed here.


.. report:: Fastqc.PerBaseQuality
   :glob: /ifs/projects/proj029/full/DNA/fastqc1/export/fastqc/122260-WT-R1.fastq.*_fastqc/Images/per_base_quality.png
   :layout: columns-2
