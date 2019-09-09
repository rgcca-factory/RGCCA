#!/bin/bash


for i in $(seq 1 10); do srun --time=1:00:00 --cpus-per-task=4 --partition=normal --mem=2G --job-name=nucleipark --account=bioinfo --mail-user=etienne.camenen@icm-institute.org Rscript R/optspars.R --job-id out/${i} & done