#!/usr/bin/zsh

### Job name
#SBATCH --job-name=rect-aspect2-0w01-{JOB_PARAM}

### File / path where STDOUT will be written, %J is the job id
#SBATCH --output=rect-aspect2-0w01-{JOB_PARAM}-job-out.%J

### Request the time you need for execution: 1 hour
#SBATCH --time=59:00
###SBATCH --time=1:00:00

### Request memory you need for your job in MB
#SBATCH --mem-per-cpu=9600

### Request number of tasks/MPI ranks
#SBATCH --ntasks=256

### Define a variable for the common string
### Change to the work directory
JOB_NAME="rect_aspect2_0w01_{JOB_PARAM}"
cd $HOME/fenicsR13
#
#### Execute the container with commands directly
$MPIEXEC apptainer exec ./fenicsR13.sif /bin/bash -c "cd 3d_crooks && fenicsR13 input/rect/0w01/input_${JOB_NAME}.yml"
