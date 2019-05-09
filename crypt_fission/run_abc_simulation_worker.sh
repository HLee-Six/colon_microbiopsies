#!/bin/bash

## The script was written by Sigurgeir Olafsson for the paper
## 'The landscape of somatic mutation in normal colorectal epithelial cells' by Henry Lee-Six et al.
## Please direct any questions or comments you may have to so11@sanger.ac.uk

script=$1
crypt_fission_rate_list=$2
biopsy_list=$3
tree_dir=$4
grid_dir=$5
event_times_dir=$6
jobIndex=${LSB_JOBINDEX}


crypt_fission_rate=$( awk -v jI="$jobIndex" 'NR==jI {print $2; exit}' < ${crypt_fission_rate_list} )
simulation_nr=$( awk -v jI="$jobIndex" 'NR==jI {print $1; exit}' < ${crypt_fission_rate_list} )
rundir=$PWD/
cd ${grid_dir}

## This is for when I split the biopsy_data into many lists. Unzip a previously compressed directory
## so I don't make a new one and overwrite the old.
if [ -f ${grid_dir}sim${simulation_nr}.tar.gz ]; then
    gunzip < sim${simulation_nr}.tar.gz | tar -xv
    rm sim${simulation_nr}.tar.gz
fi

mkdir -p ${grid_dir}sim${simulation_nr}

echo "/software/R-3.4.0/bin/Rscript ${script} ${crypt_fission_rate} ${simulation_nr} ${biopsy_list} ${tree_dir} ${grid_dir} ${event_times_dir}"
/software/R-3.4.0/bin/Rscript ${script} ${crypt_fission_rate} ${simulation_nr} ${biopsy_list} ${tree_dir} ${grid_dir}sim${simulation_nr}/ ${event_times_dir}

#tar -czf ${grid_dir}sim${simulation_nr}.tar  --directory=${grid_dir}sim${simulation_nr}/ ${grid_dir}sim${simulation_nr}/
tar -cvf sim${simulation_nr}.tar sim${simulation_nr}/
gzip ${grid_dir}sim${simulation_nr}.tar
rm -r ${grid_dir}sim${simulation_nr}/
rm -f ${grid_dir}sim${simulation_nr}.tar

cd ${rundir}

# to unzip later:
#gunzip < sim1.tar.gz | tar -xv

exit $?