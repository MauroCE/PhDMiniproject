#!/bin/bash
#
#
#PBS -l nodes=1:ppn=3,walltime=32:00:00



# Define working directory
export WORK_DIR=$HOME/newjob/job_contracting

# Define executable
export EXE="/usr/bin/env Rscript snep_contracting.R"

# Add R module
module add languages/R-3.0.2

# Change into working directory
cd $WORK_DIR

# Do some stuff
echo JOB ID: $PBS_JOBID
echo Working Directory: `pwd`
echo Start Time: `date`

# Execute code
$EXE

echo End Time: `date`



sleep 20
