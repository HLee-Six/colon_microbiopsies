
This file documents how to run the accompanying scripts for doing ABC to estimate the crypt fission rate of normal human colon.

1. Run preprocess_spatial_matrices.r locally to create a list of all necessary information for each biopsy.

2. Run find_observed_coalescence_times.r locally to count the nr of coalescence events in each time-bin in the observed data.

3. draw_priors.r (Can be run locally, no need to submit a job).
This will draw crypt fission rate values from a uniform distribution of the specified range. For me the range was 1/80 to 1
which was then divided into bins of size 0.05.

4. run_abc_simulation.sh. This should also be run locally. This will submit an array of jobs, where each job runs
run_abc_simulation_worker.sh, which in turn calls abc_simulation_worker.r, which sources abc_function_archive.r

run_abc_simulation_worker.sh, abc_simulation_worker.r and abc_function_archive.r need never be run on their own.

5. After the simulations finish running, you can extract event counts for the biopsies excluding events happening before
certain times (to exclude the effect of neonatal expansion) by running run_extract_event_counts.sh. I found it best to
break the crypt fission rate files into more manageable chunks and run the script on each chunk.
This script calls extract_event_counts.r

6. Finally, I run estimate_posterior.r locally (not as a job).









