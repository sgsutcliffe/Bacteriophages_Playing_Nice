#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4


module load sra-toolkit/2.9.6

fasterq-dump SRR828645
fasterq-dump SRR828660
fasterq-dump SRR828661
