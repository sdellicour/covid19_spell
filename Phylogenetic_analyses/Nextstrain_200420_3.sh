#!/bin/bash
#SBATCH --job-name=Nextstrain_200420_3
#SBATCH --time=740:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=10240
#SBATCH --qos=gpu_prio

module load beagle-lib/3.0.2-fosscuda-2018b

mkdir -p "$LOCALSCRATCH/$SLURM_JOB_ID"

cp -r "$SLURM_SUBMIT_DIR/*jar" "$LOCALSCRATCH/$SLURM_JOB_ID"
cp -r "$SLURM_SUBMIT_DIR/*xml" "$LOCALSCRATCH/$SLURM_JOB_ID"

cd  "$LOCALSCRATCH/$SLURM_JOB_ID"

java -jar beast_1104_141118.jar -beagle_gpu -beagle_double -beagle_order 1 -overwrite Nextstrain_200420_3.xml

cp -r "$LOCALSCRATCH/$SLURM_JOB_ID/*out" "$SLURM_SUBMIT_DIR/" &&\
cp -r "$LOCALSCRATCH/$SLURM_JOB_ID/*log" "$SLURM_SUBMIT_DIR/" &&\
cp -r "$LOCALSCRATCH/$SLURM_JOB_ID/*trees" "$SLURM_SUBMIT_DIR/" &&\

rm -rf "$LOCALSCRATCH/$SLURM_JOB_ID"
