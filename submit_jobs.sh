#!/bin/bash
# In order to change parameters, use --export=param1=val1,param2=val2,...
# Parameter multiD used for running several values of D sequentially
for h in $(seq -f "%2.2f" 0.0 0.03 3.0); do
    # Submit the job
    # sbatch -A ctmc --export=hz=$h,D=2,multiD=2,Jx=1.2,Jy=0.8,Jz=0. SSM.run
    sbatch -A ctmc --export=hz=$h,D=2,multiD=4,model=XYZ_stagH,J1=3.2,J2=-10.6,Delta1=1.0,Delta2=0.5 SSM.run
    echo "Submitted job for h=$h"
done