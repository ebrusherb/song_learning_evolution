#!/bin/tcsh

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 20:00:00
#SBATCH -J maleboth_femaleboth
#SBATCH -o maleboth_femaleboth.out
# #SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=eleanor.brush@gmail.com

module load R

# mkdir -p /tmp/ebrush/$SLURM_JOB_ID

echo '-------------------------------'
date
echo 'Current directory:'
pwd
echo 'Running on'
hostname
echo '-------------------------------'

R --no-save < contsong_paramsweep_par_maleboth_femaleboth.R

# rm -rf /tmp/ebrush/$SLURM_JOB_ID
