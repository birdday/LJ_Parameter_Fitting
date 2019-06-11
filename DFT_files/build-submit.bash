#!/bin/env/bash

dft_input_files_dir="DFT_Input_Files"
dft_input_file="energy-nonperiodic"
xyz_dir="sotolon_graphite_xyzfiles"
run_dir="sotolon_graphite_dftcalcs"
box_size="19.688 19.688 15"

for i in *.xyz; do
  echo "${i%.*}" >> submissions.txt
  mkdir "../${run_dir}/${i%.*}"
  cp "${i}" ../${dft_input_files_dir}/* "../${run_dir}/${i%.*}"
  cd "../${run_dir}/${i%.*}"
  sed -e "s|SLURM_SYSNAME|${i%.*}|g" \
      -e "s/EXTERNAL_XYZFILE/${i}/g" \
      -e "s/EXTERNAL_INPUT_FILE/${dft_input_file}/g" \
      -e "s/EXTERNAL_BOX_SIZE/${box_size}/g"
      "cp2k-template.slurm" > "cp2k-submit.slurm"
  sbatch cp2k-submit.slurm
  cd "../${xyz_dir}"

done
