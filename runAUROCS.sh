#!/bin/bash

#SBATCH -n 16
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --tmp=4000                        # per node!!
#SBATCH --job-name=netPropAUROC
#SBATCH --output=netPropAUROC.out
#SBATCH --error=netPropAUROC.err

env2lmod

module load gcc/11.4.0 r/4.3.2
module load gcc/11.4.0 glpk/5.0
module load gcc/11.4.0 gmp/6.2.1

Rscript /cluster/home/gmagnusson/netprop/Code/runNetProp.R
