#!/bin/bash -l
#SBATCH --time=00:30:00
#SBATCH --ntasks=2
#SBATCH --mem=14G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=$USER@umn.edu
#SBATCH --account lniedern
#SBATCH --array=0-6

module load rclone
source /users/6/${USER}/bin/mamba.init
datestamp=`date '+%y%m%d_permutation'`
samples=("H1-QYNNNRX_A1" "H1-QYNNNRX_D1" "H1-TMNBZ4N_A1" "H1-TMNBZ4N_D1" "H1-TQXC9QB_A1" "H1-TQXC9QB_D1" "H1-HDMGRMD_A1")

basedir=/scratch.global/$USER/${datestamp}
mkdir $basedir
mkdir $basedir/collated
mkdir $basedir/images
for sample in samples; do mkdir $basedir/${sample}; done

python3 /users/6/$USER/pylib/VisHD/make_genecat_table.py $basedir "${samples[@]}"

#python3 /users/6/${USER}/pylib/VisHD/Permutation_Runner.py $SLURM_ARRAY_TASK_ID

