#!/bin/bash
#SBATCH --job-name=Nextstrain_200420_4
#SBATCH --time=740:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=10240

module load beagle-lib/3.0.2-fosscuda-2018b

cd
cd Nextstrain_DTA

java -jar beast_1104_141118.jar -beagle_gpu -beagle_double -beagle_order 1 -overwrite Nextstrain_200420_4.xml

