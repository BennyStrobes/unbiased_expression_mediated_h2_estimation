#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



simulation_results_dir="$1"
visualize_simulation_dir="$2"


if false; then
expression_data_sample_size="200"
sbatch run_one_instance_of_expression_mediated_h2_simulation.sh $simulation_results_dir $expression_data_sample_size

expression_data_sample_size="400"
sbatch run_one_instance_of_expression_mediated_h2_simulation.sh $simulation_results_dir $expression_data_sample_size

expression_data_sample_size="600"
sbatch run_one_instance_of_expression_mediated_h2_simulation.sh $simulation_results_dir $expression_data_sample_size

expression_data_sample_size="800"
sbatch run_one_instance_of_expression_mediated_h2_simulation.sh $simulation_results_dir $expression_data_sample_size

expression_data_sample_size="1000"
sbatch run_one_instance_of_expression_mediated_h2_simulation.sh $simulation_results_dir $expression_data_sample_size
fi