#!/bin/bash

# The purpose of this script and the extract_event_counts.r script is to go back to the grids
# and calculate new event times without having to re-generate the grids.

## The script was written by Sigurgeir Olafsson for the paper
## 'The landscape of somatic mutation in normal colorectal epithelial cells' by Henry Lee-Six et al.
## Please direct any questions or comments you may have to so11@sanger.ac.uk

## Submit to cluster:
# for i in {00..38}; do  bsub -o lt_0.05_event_counts_${i}.out -e lt_0.05_event_counts_${i}.err -M200 -R'span[hosts=1] select[mem>200] rusage[mem=200]' "bash ~/phd/somatic_ibd_p1/ABC/normal_fission_rate/run_extract_event_counts.sh ~/phd/somatic_ibd_p1/ABC/normal_fission_rate/biopsy_data_combined.txt lt_0.05 $PWD/fission_rate_files_tmp/lt_0.05${i}"; done

script_dir=/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/ABC/normal_fission_rate/
output_dir=/lustre/scratch114/projects/crohns/somatic_ibd_p1/ABC/normal_fission_rate/
biopsy_list=$1
batch=$2
crypt_fission_rate_list=$3
#crypt_fission_rate_list=/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/ABC/normal_fission_rate/crypt_fission_rate_${batch}.txt

event_times_dir=${output_dir}event_times/${batch}/
tree_dir=${output_dir}trees/${batch}/
grid_dir=${output_dir}grids/${batch}/
clusterOutput=${output_dir}clusterOutput/${batch}/

mkdir -p ${clusterOutput}
mkdir -p ${tree_dir}
mkdir -p ${grid_dir}
mkdir -p ${event_times_dir}

while read simulation_nr crypt_fission_rate; do
    /software/R-3.4.0/bin/Rscript ${script_dir}extract_event_counts.r ${crypt_fission_rate} ${simulation_nr} ${biopsy_list} ${tree_dir} ${grid_dir} ${event_times_dir}
done < ${crypt_fission_rate_list}



exit $?