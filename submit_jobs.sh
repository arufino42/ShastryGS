#!/bin/bash
# In order to change parameters, use --export=param1=val1,param2=val2,...
# Parameter multiD used for running several values of D sequentially
for h in $(seq -f "%2.2f" 0.5 0.02 2.0); do
    # Submit the job
    sbatch -A ctmc --export=hz=$h,D=2,multiD=0,Jx=1.4,Jy=0.6,Jz=0. SSM.run
    echo "Submitted job for h=$h"
done