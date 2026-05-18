#!/bin/bash
# In order to change parameters, use --export=param1=val1,param2=val2,...
# Parameter multiD used for running several values of D sequentially
for H_dir in {100,110,001}; do
    # Submit the job
    echo H_dir=$H_dir
    if [ $H_dir -eq 001 ]; then
        hmax=4.0
    else
        hmax=1.0
    fi
    
    echo "sbatch -A ctmc --export=H_dir=$H_dir,D=2,multiD=2,hmax=$hmax,Delta1=0.2,Delta2=0.2,eta=1.5 SSM.run"
    sbatch -A ctmc --export=H_dir=$H_dir,D=2,multiD=2,hmax=$hmax,Delta1=0.2,Delta2=0.2,eta=1.5 SSM.run
done