simul_debug_fname = 'simulations_debug.txt'
OBSOLETE = ''
num_kars = 100
simul_name = 'no simulation selected'

# For CreateGraph:
out_fname = '_graph.txt'
graph_viz = '_gv.txt'


# For simulations
min_seg_size = 10000	# The minimal size in bases of each segment.
total_size = 3000000	# The length (in bases) of the chromosome)
max_dist_relative = 0.1		# the maximum distance a bp can span in percentage of the total length.
num_segs = total_size/min_seg_size						#The number of segments in a karyotype
max_dist = 	max_dist_relative*total_size/min_seg_size	#The max dist in segments
supp_dist_lambda = 0.1866

# For BP_graph
graph_fname = '_graph.txt'
MAX_BP_WINDOW = min_seg_size/2
TERMINTOR_EDGE_SIZE = min_seg_size*4
DEFAULT_END_POSITION = total_size


## new simulation sets

# General for all simulations
SIMULATION = True
ILP_ver = '4.2' #This is just used to easily distinguish between result files in case of multiple runs. It needs to be the same as the one in the java code.
ver = ILP_ver
eps_dist = OBSOLETE
rearrangement = ['INV', 'DUP', 'DEL', 'TRL']

# default values
alpha_values = [0.1]
eps_cnv = [0.28]    # The mean noise value to be added or subtracted from the real CN value
max_ops=5
min_ops=5
num_chrs = 5
chromosomes = [str(x) for x in range(num_chrs+1)[1:]]
DEFAULT_CN = 2	#This should be 2 for diploidity / 1 for single copy
p_miss_brdg = 0.05  # Probability to miss a bridge
noise_kar_ratio = 0
noise_kar_options = ['none', 'same', 'dif']
SIMULATION = True

#values for base scenario 100*100
#iter = 100
#simul_name = 'simul810'+'_'+str(iter)
#CNV_fname = simul_name+'_cnv.txt'
#BP_fname = simul_name+'_bp.txt'
#defaults_list_name = simul_name+'_list.txt'
#num_kars = 10000

#values for simulating different alpha values (please note that the java code also needs to be modified in this case):
simul_name = 'simul811'
CNV_fname = simul_name+'_cnv.txt'
BP_fname = simul_name+'_bp.txt'
defaults_list_name = simul_name+'_list.txt'
alpha_values = [0.0,0.01,0.02,0.03,0.04,0.05,0.1,0.25,0.5,1.0,2.0]

#values for simulating different noise values:
#simul_name = 'simul802'
#CNV_fname = simul_name+'_cnv.txt'
#BP_fname = simul_name+'_bp.txt'
#defaults_list_name = simul_name+'_list.txt'
#eps_cnv = [0,0.1,0.2,0.28,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

#values for simulating diofferent number of operations:
#simul_name = 'simul803'
#CNV_fname = simul_name+'_cnv.txt'
#BP_fname = simul_name+'_bp.txt'
#defaults_list_name = simul_name+'_list.txt'
#min_ops=1
#max_ops=30
#num_kars = 100

#values for simulating different number of chromosomes:
#simul_variant = '10'
#simul_name = 'simul804'+'_'+simul_variant
#CNV_fname = simul_name+'_cnv.txt'
#BP_fname = simul_name+'_bp.txt'
#defaults_list_name = simul_name+'_list.txt'
#num_chrs = int(simul_variant)
#chromosomes = [str(x) for x in range(num_chrs+1)[1:]]

#values for simulating different number of copies per chromosome:
#simul_variant = '1chr'
#simul_name = 'simul805'+'_'+simul_variant
#CNV_fname = simul_name+'_cnv.txt'
#BP_fname = simul_name+'_bp.txt'
#defaults_list_name = simul_name+'_list.txt'
#DEFAULT_CN = 2	#This should be 2 for diploidity / 1 for single copy

# change in the probability to miss a bridge
#values for simul 16: 0,0.025,0.05,0.075,0.1,0.125,0.15
#p_miss_brdg = 0.15
#simul_name = 'simul816'+'_'+str(p_miss_brdg)
#CNV_fname = simul_name+'_cnv.txt'
#BP_fname = simul_name+'_bp.txt'
#defaults_list_name = simul_name+'_list.txt'


#values for simulating heterogenity model with another rearranged karyotype
#noise_kar_options = ['none', 'same', 'dif']
#noise_kar_ratio = 0	#[0, 0.05, 0.1, 0.15, 0.2]
#simul_name = 'simul818a'+'_'+str(int(noise_kar_ratio*100))+'%'
#CNV_fname = simul_name+'_cnv.txt'
#BP_fname = simul_name+'_bp.txt'
#defaults_list_name = simul_name+'_list.txt'

#values for simulating contamination by a healthy tissue
#noise_kar_options = ['none']
#noise_kar_ratio = 0.3 #[0,0.05,0.1,0.15,0.2,0.3,0.45]
#simul_name = 'simul808b'+'_'+str(int(noise_kar_ratio*100))+'%'
#CNV_fname = simul_name+'_cnv.txt'
#BP_fname = simul_name+'_bp.txt'
#defaults_list_name = simul_name+'_list.txt'


#values for simulating different frequencies of operations
#simul_name = 'simul809_norm'
#CNV_fname = simul_name+'_cnv.txt'
#BP_fname = simul_name+'_bp.txt'
#defaults_list_name = simul_name+'_list.txt'
#rearrangement = ['DEL', 'DEL', 'DEL', 'DEL', 'DUP', 'DUP', 'DUP', 'DUP', 'INV', 'TRL']	#for malhorta dist uncomment this
#num_kars = 250


#For real data
#SIMULATION = False
#CNV_fname = 'CNV.txt'
#BP_fname = 'BP.txt'
#num_chrs = 23
#chromosomes = [str(x) for x in range(num_chrs)[1:]]+['X','Y']