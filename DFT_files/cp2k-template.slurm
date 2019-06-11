#!/bin/env bash
#SBATCH --job-name=SLURM_SYSNAME-energy
#SBATCH --output=SLURM_SYSNAME-energy.out
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=18
#SBATCH --cluster=mpi
#SBATCH --time=0-02:00:00

## ======================================================
## ===== Only make changes here!!!!  ====================
## ======================================================
sys_name="SLURM_SYSNAME"
xyz_file="EXTERNAL_XYZFILE"
template="EXTERNAL_INPUTFILE"
project_file="${sysname}-energy"
run_type="ENERGY"
print_level="MEDIUM"
box_size="EXTERNAL_BOX_SIZE"

## ======================================================
## ===== Prepare the Cores ==============================
## ======================================================
module purge
module load intel/2017.3.196 intel-mpi/2017.3.196 cp2k/6.1
ulimit -s unlimited
ulimit -l unlimited

## ======================================================
## ===== Copy all necessary files to all nodes ==========
## ======================================================
files=(${template} ${xyz_file} elements-PBE.inc)
for i in ${files[@]}; do
    sbcast $SLURM_SUBMIT_DIR/$i $SLURM_SCRATCH/$i
done

## ======================================================
## ===== Copy all files back at end of run ==============
## ======================================================
run_on_exit(){
    cp -r $SLURM_SCRATCH/* $SLURM_SUBMIT_DIR
    pkill --uid=$SLURM_JOB_USER cp2k.popt
}
trap run_on_exit EXIT

## ======================================================
## ===== Geometry Optimization Calculations =============
## ======================================================

cd $SLURM_SCRATCH
export I_MPI_FABRICS_LIST=ofa,tmi
export I_MPI_FALLBACK=0

## Fill in template
sed -e "s|EXTERNAL_PROJNAME|${project_file}|g" \
    -e "s/EXTERNAL_RUN_TYPE/${run_type}/g" \
    -e "s/EXTERNAL_PRINT_LEVEL/${print_level}/g" \
    -e "s/EXTERNAL_BOX_SIZE/${box_size}/g"
    -e "s/EXTERNAL_XYZ_FILE/${xyz_file}/g" \
    ${template} > "${project_file}-run.inp"

## Run calculation
echo "Beginning ${runtype} Calculation"
mpirun -np $SLURM_NTASKS cp2k.popt -i ${project_file}-run.inp -o ${project_file}-run.out
echo "Calculation finished successfully"
