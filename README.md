# Bd-RNA-seq-Analysis
This pipeline shows the RNA-seq analysis of the fungus Batrachochytrium dendrobatidis (Bd)

In UNIX you need to have Kallisto installed first
```
kallisto index -i transcripts.idx Batde5_best_transcripts.fasta 
```

For the analysis of multiple files, if you have the names of the fastq samples in a text filecalled "samples"
```
for s in 'cat samples';
do
kallisto quant -i transcripts.idx -o $s -b 100 $s_R1.fastq $s_R2.fastq;
done
```


If this works, now we move to R for statistical analysis

Loading required libraries
```
library('cowplot')
library('sleuth')

metada = read.csv("metadata.csv",header=T)

sample_ids = dir(file.path('KAL'))

sample_to_condition = dplyr::mutate(metada,path = file.path('KAL', sample_ids))

so <- sleuth_prep(sample_to_condition, ~Condition, extra_bootstrap_summary = TRUE)

plot_pca(so, color_by = 'Condition')

so <- sleuth_fit(so, ~Condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')

so <- sleuth_lrt(so, 'reduced', 'full')

full_results <- sleuth_results(so, 'reduced:full', 'wt', show_all = FALSE)

sleuth_significant <- dplyr::filter(full_results, qval <= 0.05)

so <- sleuth_wt(so, 'Condition')

results_table <- sleuth_results(so, 'Condition')

results_ordered <- results_table[order(results_table$qval), ]

write.table(results_ordered, file='results.txt', sep="\t", row.names=F, quote=F)
```

