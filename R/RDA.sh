#! /bin/bash
#SBATCH --time=4-0:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=gpu
#SBATCH --error=msg/rda.err
#SBATCH --output=msg/rda.out
#SBATCH --mail-user=emittman@iastate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load R
module load cuda

cd /home/emittman/computation-paper/R/
R CMD BATCH --vanilla --no-save RDA.R #run an R script using R
