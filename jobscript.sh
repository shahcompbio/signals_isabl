#!/bin/bash
#SBATCH --job-name="signals_%j"
#SBATCH --partition=componc_cpu
#SBATCH --output=logs/signals_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00

singularity exec --bind=/data1 --bind=/home \
    /data1/shahs3/users/william1/software/singularity/rstudio-docker_main.sif \
    Rscript run_signals.R  \
    --hmmcopyqc /data1/shahs3/isabl_data_lake/analyses/18/69/31869/results/SHAH_H000708_T05_01_DLP01_hmmcopy_metrics.csv.gz /data1/shahs3/isabl_data_lake/analyses/19/36/31936/results/SHAH_H000708_T03_05_DLP01_hmmcopy_metrics.csv.gz \
    --hmmcopyreads /data1/shahs3/isabl_data_lake/analyses/18/69/31869/results/SHAH_H000708_T05_01_DLP01_hmmcopy_reads.csv.gz /data1/shahs3/isabl_data_lake/analyses/19/36/31936/results/SHAH_H000708_T03_05_DLP01_hmmcopy_reads.csv.gz \
    --allelecounts /data1/shahs3/isabl_data_lake/analyses/95/33/29533/results/haplotype_calling_count.csv.gz /data1/shahs3/isabl_data_lake/analyses/05/61/30561/results/haplotype_calling_count.csv.gz \
    --ncores 5 \
    --qcplot results/OV-107/qcplot.png \
    --heatmap results/OV-107/heatmap.png \
    --heatmapraw results/OV-107/heatmapraw.png \
    --csvfile results/OV-107/hscn.csv.gz \
    --qccsvfile results/OV-107/qc.csv \
    --Rdatafile results/OV-107/schnapps.Rdata

#singularity shell --bind=/juno --bind=/home /juno//work/shah/users/william1/singularity//schnapps_v0.5.4.sif
