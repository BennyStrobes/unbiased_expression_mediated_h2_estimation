#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-30:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



simulation_results_dir="$1"
expression_data_sample_size="$2"




python3 run_one_instance_of_expression_mediated_h2_simulation.py $simulation_results_dir $expression_data_sample_size