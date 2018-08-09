#!/bin/sh
#SBATCH -S 24
#SBATCH -N 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=M.J.Tadema@student.rug.nl
#SBATCH -J "COBRA optimalization styrene"
#SBATCH -o ./cobra.log
#SBATCH --ntasks=24
#SBATCH --partition=gpu

module load MATLAB/2018a
module load GLPK/4.58-foss-2016a

matlab -nodisplay -nodesktop -r "run ./igemoptforcemain.m"
