# VNSRR
This repo includes the code to implement Variable Neighborhood Searching Rerandomization(VNSRR) and reproduce the numerical experiments.

## About VNSRR: 
Rerandomization discards undesired treatment assignments to ensure covariate balance in randomized experiments. However, rerandomization based on acceptance-rejection sampling is computationally inefficient, especially when numerous independent assignments are required to perform randomization-based statistical inference. Existing acceleration methods are rudimentary and are not applicable in structured experiments, including stratified experiments and experiments with clusters. Based on metaheuristics in combinatorial optimization, we propose a novel variable neighborhood searching rerandomization(VNSRR) method to efficiently draw balanced assignments in various experiments. We derive the unbiasedness and a lower bound for the variance reduction of the treatment effect estimator under VNSRR. Simulation studies and a real data example indicate that our method maintains the appealing statistical properties of rerandomization and can sample thousands of treatment assignments within seconds, even in cases where existing methods require an hour to complete the task.

## About the code:
The files with index from 1a to 1e generate the data in the simulation studies.

The files with index from 2a to 2g apply VNSRR and the competing methods to the simulated data.

The files with index from 3a to 3b analyze the real dataset.

The files whose file names start with 'rerand' implements VNSRR and the competing methods in different types of randomized experiments.

The files within the folder 'chart' summarize the simulation results and compare the performances of different methods.
