#!/usr/bin/zsh

# List of parameters
#params=("0kn005")
params=("0kn01" "0kn02" "0kn04" "0kn08" "0kn16" "0kn32" "0kn64" "1kn28" "2kn56" "5kn12")

# Template and output file
template="rect_aspect2_0w01_template.sh"

for param in ${params[@]}; do
	# Generate a job file name
	job_file="rect_aspect2_0w01_${param}.slurm"

	# Replace placeholder in template and create a new job file
	sed "s/{JOB_PARAM}/$param/g" $template > $job_file

	# Submit the job
	sbatch $job_file

	#Optionally remove the job_file
	#rm $job_file
done
