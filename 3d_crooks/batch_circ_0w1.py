import subprocess

# Long list of input files
input_files = [
    #"input/circ/0w1/input_circ_0w1_0kn005.yml",
    #"input/circ/0w1/input_circ_0w1_0kn01.yml",
    "input/circ/0w1/input_circ_0w1_0kn02.yml",
    "input/circ/0w1/input_circ_0w1_0kn04.yml",
    "input/circ/0w1/input_circ_0w1_0kn08.yml",
    #"input/circ/0w1/input_circ_0w1_0kn16.yml",
    "input/circ/0w1/input_circ_0w1_0kn32.yml",
    "input/circ/0w1/input_circ_0w1_0kn64.yml",
    "input/circ/0w1/input_circ_0w1_1kn28.yml",
    "input/circ/0w1/input_circ_0w1_2kn56.yml",
    "input/circ/0w1/input_circ_0w1_5kn12.yml"
]

# Loop through each input file and run the Python program with mpirun -n 8
for input_file in input_files:
    print(f"Running program with '{input_file}' using mpirun -n 8...")
    # Construct the command to run the program with mpirun
    command = ["mpirun", "-n", "8", "fenicsR13", input_file]
    
    # Execute the command
    subprocess.run(command, check=True)  # check=True raises an error if the command fails
    
    print(f"Finished running '{input_file}'")
    print("---------------------------------")

