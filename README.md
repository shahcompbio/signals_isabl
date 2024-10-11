## signals 

Script to run signals for use in signals isabl app. Command to run:

```
img="/data1/shahs3/users/william1/software/singularity/rstudio-docker_main.sif"

singularity exec --bind=/data1 --bind=/home \
    $img \
    Rscript run_signals.R  \
    --hmmcopyqc /data1/shahs3/isabl_data_lake/analyses/60/71/36071/results/SHAH_H002689_T03_01_DLP01_hmmcopy_metrics.csv.gz \
    --hmmcopyreads /data1/shahs3/isabl_data_lake/analyses/60/71/36071/results/SHAH_H002689_T03_01_DLP01_hmmcopy_reads.csv.gz \
    --allelecounts /data1/shahs3/isabl_data_lake/analyses/21/14/42114/results/haplotype_calling_count.csv.gz \
    --ncores 5 \
    --qcplot results/128701A/qcplot.png \
    --heatmap results/128701A/heatmap.png \
    --heatmapraw results/128701A/heatmapraw.png \
    --csvfile results/128701A/hscn.csv.gz \
    --qccsvfile results/128701A/qc.csv \
    --Rdatafile results/128701A/signals.Rdata \
    --sex MALE \
    --cell_list data/cell_list_128701A.txt
```

### v0.2.0 
For use with SIGNALS v0.11.0. Includes the following additions
- `sex` parameter for specifying patient sex, this determines what happens with chrX.
- `cell_list` only run for the set of cells listed in cell_list
- No additional filtering beyond HMMcopy
- Sets homozygous deletions to state 0|0

[Docker for v.11.0](https://hub.docker.com/layers/marcjwilliams1/signals/v0.11.0/images/sha256-f6f076ecd505a76e3fb13e4f0607647fa4da6e262b0e09a05854fabaddfacd20?context=explore). 
Path on iris: ``
