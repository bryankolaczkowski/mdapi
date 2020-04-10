#!/bin/bash

#SBATCH --job-name=mdapi             # Job name
#SBATch --nodes=1                    # single node
#SBATCH --ntasks=1                   # single task
#SBATCH --cpus-per-task=8            # multi-threaded
#SBATCH --distribution=cyclic:block  # thread distribution
#SBATCH --mem=200MB                  # memory
#SBATCH --time=36:00:00              # time limit hrs:min:sec
#SBATCH --output=mdapi.out           # standard output log
#SBATCH --error=mdapi.err            # standard error log
#SBATCH --qos=bryankol-b             # queue
#SBATCH --partition=hpg2-compute     # make sure only intel processors used

date

module load intel/2018 gromacs/2019.2 pdb2pqr

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

prot=$1

mdapi.py --threads=$OMP_NUM_THREADS setup   ${prot}
mdapi.py --threads=$OMP_NUM_THREADS run     ${prot}
mdapi.py --threads=$OMP_NUM_THREADS analyze ${prot}

date
echo "done."
exit 0
