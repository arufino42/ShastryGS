#!/bin/bash
# In order to change parameters, use --export=param1=val1,param2=val2,...
# Parameter multiD used for running several values of D sequentially
for H_dir in {100,110,001}; do
    # Submit the job
    echo H_dir=$H_dir
    sbatch -A ctmc --export=H_dir=$H_dir SSM.run
done