#!/bin/bash

srun --time=1:00:00 --ntasks=50 --partition=bioinfo --mem=4G --job-name=nucleipark --mail-user=etienne.camenen@icm-institute.org conda Rscript R/optpars.R