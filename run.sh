#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output test.%j.out
#SBATCH --account=ucbclass2_summit1
#SBATCH --qos debug
#SBATCH --time=00:15:00

start=`date +%s`
./main.x
end=`date +%s`

runtime=$((end-start))
