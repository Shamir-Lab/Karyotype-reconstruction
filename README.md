# Karyotype reconstruction
Reconstruct a rearranged karyotype from bridge and copy number data using ILP.
You can also use this code to generate simulations of rearranged karyotypes.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

* The ILP code is in java, so you need to set up JRE and have a compiler for java, [eclipse](https://www.eclipse.org/) is free and recommended.
* [IBM CPLEX](https://onboarding-oaas.docloud.ibmcloud.com/software/analytics/docloud) package for java (a jar file). License is free for academics.
* Python 2.7

### Setting up
* Make sure the CPLEX jar file (cplex.jar) is included in your projects' classpath
* Compile ex1.java with CPLEX (or you can run it straight from eclipse)

example (paths may differ):
```
javac -classpath workspace/LP_ex1/bin:/usr/local/stow/cplex125/lib/cplex125/lib/cplex.jar -d ./work/ ./workspace/LP_ex1/src/ex1.java
```

## Running

### Creating simulations

The main script is CreateGraph.py running which takes one argument. Typing

```
python CreateGraph.py
```

will print out a menu for all possible arguments:

* 1: create_simulation_batch - Create a batch of simulation (100 by default)
* 2: convert_karyotype_list_to_graph - Create files of the correct rearrangement for reference.
* 3: process_simulation_results - After the ILP created result files, this process the results and creates a tab delimited text file.
* 4: process_simulation_results_simul10 - Divides the reults into a hundred groups and analyse each seperately (to produce the results shown in Fig5)
* 5: process_simulation_results_for_cn_score - Calculates the CN score for each simulation.
* 6: clean_all - deletes all created files.


After creating the files for the simulations, you can run the ILP code using the file ex1.java with one parameter containing the list of samples (this is created by the python script).
Make sure the simulation files are accessible to the ILP code.

Example for a simulation list file named "simul_list". Note the paths may differ.
```
java -classpath workspace/LP_ex1/bin:/usr/local/stow/cplex125/lib/cplex125/lib/cplex.jar -Djava.library.path=/usr/local/lib/cplex125/bin/ ex1 "simul_list"
```

Once the ILP is done you can run any of the results processing scripts and get a tab delimited text file.

an example run may be:
```
python CreateGraph.py 1
python CreateGraph.py 2
java -classpath workspace/LP_ex1/bin:/usr/local/stow/cplex125/lib/cplex125/lib/cplex.jar -Djava.library.path=/usr/local/lib/cplex125/bin/ ex1 "simul_list"
python CreateGraph.py 3
python CreateGraph.py 6
```

### Changing the simulation parameteres

All of the parameters are defined in parameters.py.
For convinience the different set ups used in the paper are included and commented out.

The main parameters changed for each simulation are:

* alpha_values = A list of all alpha values to be checked. NOTE - this parameter has to be modified in ex1.java as well.
* eps_cnv - Mean noise level for the CN data.
* max_ops - Largest number of rearrangements a karytype goes through
* min_ops - Lowest number of rearrangements a karyotype goes through
* num_chrs - Max number of chromosomes involved in a rearrangement
* chromosomes - A list of strings representing the chromosomes names (default is ['1', '2', '3', ....])
* DEFAULT_CN = How many copies from each chromosome (2 for diploidity)
* p_miss_brdg - the probability of missing a bridge
* noise_kar_ratio - the percent of samples taken from a different karyotype
* noise_kar_options - A list detailing the rearrangement the "noise" karyotype goes through in relation to the main karyotype. the options are:
  * 'none' - The noisy karyotype doesn not go through a rearrangement (if this is the only option, the noisy karyotype will represent a healthy tissue)
  * 'same' - The noisy karyotype goes through the same rearrangement. 
  * 'dif' - The noisy karyotype goes through a different rearrangement.


In addition these parameters can also be set:

* SIMULATION - a boolean flag signfying if this is simulated or real data
* num_kars - the number of karyotypes in each simulation batch.
* simul_name - a name for the specific simulation that will be used as file name headers
* min_seg_size -The minimal size in bases of each segment (default - 300)
* total_size - The length (in bases) of the chromosome (default - 3,000,000)
* max_dist_relative - the maximum distance a bp can span in percentage of the total length (default - 0.1)
* supp_dist_lambda - The lambda for the distribution of support scores (default - 0.1866)
* ILP_ver - a string for the current version that will be concatentaed to the file names. NOTE: if changed, this must be identical to the value in ex1.java
* rearrangement - a list of the possible rearrangement types (default - ['INV', 'DUP', 'DEL', 'TRL'])

NOTE: the list of alpha values should also be set in ex1.java and be identical to the ones in params.py. the line to change is:
```
static double[] alpha_values = {0.1};
```

### Running on a real data.
To run the script on real data do the following steps:

* Create a tab delimited file called BP.txt containing all the bridge data. See below for details.
* Create a tab delimited file called CNV.txt conatining all the copy number data. See below for details.
* Change SIMULATION flag in params.py to False
* Run CreateGraph.py with argument 1, this will create the input files for the ILP
```
python CreateGraph.py 1
```
* create a text file samp_list.txt containing all the sample names you wish to run the ILP on.
* Run the ILP code
```
java -classpath workspace/LP_ex1/bin:/usr/local/stow/cplex125/lib/cplex125/lib/cplex.jar -Djava.library.path=/usr/local/lib/cplex125/bin/ ex1 "samp_list"
```

You can use the function CreateGraph.create_gv_input_list(samp_list) to convert the solution files into a dot language text file that can be visualised (using gvedit eg).



#### Create BP.txt and CNV.txt
In the BP.txt file each row represents a bridge between two breakpoints and has the following columns:
* chrom1 - Chromosome of the first breakpoint
* pos1 - Start position of the first breakpoint
* pos2 - End position of the first breakpoint
* chrom2 - Chromosome of the second breakpoint
* pos3 - Start position of the second breakpoint
* pos4 - End position of the second breakpoint
* support - Support score of the bridge
* strand1	- Strand of the first breakpoint (which "side" of the break is joined)
* strand2	- Strand of the second breakpoint (which "side" of the break is joined)
* variantClass - (INV, DUP, DEL, INTRA, INTER) (unused)
* sample - The sample the break was found in	

A sample BP.txt file can be found in the repository (based on data taken from Malhorta et al.)

In the CNV.txt file each row represents on copy number change, and has the following column headers:
* Chromosome - Chromosome of the change-point
* Start - Start coordinate of the change-point
* End	- End coordinate of the change-point
* Left Segment Zscore - Zscore of the flanking segment on the left of the change-point (unused)
* Right Segment Zscore - Zscore of the flanking segment on the right of the change-point (unused)
* Delta Zscore - Difference in Zscores (unused)
* Left Segment Copy Number - Estimated copy number of the flanking segment to the left of the change-point
* Right Segment Copy Number - Estimated copy number of the flanking segment to the right of the change-point
* Sample name - The sample the changepoint is found in.

A sample CNV.txt file can be found in the repository (based on data taken from Malhorta et al.)


## Authors

* **Rami Eitan** - [RamiEitan2337](https://github.com/RamiEitan2337)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

