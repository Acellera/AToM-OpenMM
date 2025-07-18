#!/bin/bash
#
#SBATCH -J <JOBNAME>
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --account=<account>
#SBATCH --no-requeue
#SBATCH -t 10:00:00

jobname=<JOBNAME>

. $HOME/miniforge3/bin/activate
echo "Running on $(hostname)"

if [ ! -f ${jobname}_0.xml ]; then
   rbfe_structprep ${jobname}_asyncre.cntl || exit 1
fi

echo "localhost,0:0,1,CUDA,,/tmp" > nodefile
rbfe_production ${jobname}_asyncre.cntl
