#!/bin/bash

## This script takes care of submitting simulation jobs to the Sanger computing FARM.
## It submits a job-array with the same number of jobs as there are lines in the crypt_fission_rate_list file.

## The script was written by Sigurgeir Olafsson for the paper
## 'The landscape of somatic mutation in normal colorectal epithelial cells' by Henry Lee-Six et al.
## Please direct any questions or comments you may have to so11@sanger.ac.uk

script_dir=/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/ABC/normal_fission_rate/
#output_dir=/lustre/scratch119/humgen/teams/anderson/users/so11/somatic_ibd_p1/ABC/normal_fission_rate/
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


largest_CFR=$( awk '{print $2}' < ${crypt_fission_rate_list} | sort -rn | head -1 )
q=$( awk -v lar="$largest_CFR" 'BEGIN {if(lar<0.05) {print "normal"} else {print "long"}}' )

nrSimulations=$( wc -l ${crypt_fission_rate_list} | awk '{print $1}' )

if [[ "$q" == "long" ]]; then
    ## Only allows 7000 jobs running at a time.
    bsub -J "crypt_fission_abc[1-${nrSimulations}]%7000" -q long -M500 -R'span[hosts=1] select[mem>500] rusage[mem=500]' \
    -e ${clusterOutput}crypt_fission_abc_errors.%J.%I  -o ${clusterOutput}crypt_fission_abc_output.%J.%I  \
    bash ${script_dir}run_abc_simulation_worker.sh ${script_dir}abc_simulation_worker.r \
    ${crypt_fission_rate_list} ${biopsy_list} ${tree_dir} ${grid_dir} ${event_times_dir}
else
    ## Only allows 7000 jobs running at a time.
    bsub -J "crypt_fission_abc[1-${nrSimulations}]%7000" -q normal -M5000 -R'span[hosts=1] select[mem>5000] rusage[mem=5000]' \
    -e ${clusterOutput}crypt_fission_abc_errors.%J.%I  -o ${clusterOutput}crypt_fission_abc_output.%J.%I  \
    bash ${script_dir}run_abc_simulation_worker.sh ${script_dir}abc_simulation_worker.r \
    ${crypt_fission_rate_list} ${biopsy_list} ${tree_dir} ${grid_dir} ${event_times_dir}
fi



exit $?