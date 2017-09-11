import sys
import os
import shutil
from decimal import *
import simulations
import BP_graph
import params
import numpy
import time
import scores

from BP_graph import *
from params import *

		
def dir2int(dir_chr):
	if dir_chr=='+':
		return 0
	if dir_chr=='-':
		return 1
def int2dir(dir_int):
	if dir_int==0:
		return '+'
	if dir_int==1:
		return '-'

def neighbours(node, graph):
	neighbours = []
	candidates=[node]
	done = []
	
	while candidates!=[]:		
		neighbours += [x for x in graph[candidates[0]] if x not in neighbours]
		if candidates[0] not in done:
			done.append(candidates[0])
		candidates += [x for x in graph[candidates[0]] if x not in done and x not in candidates]
		candidates.remove(candidates[0])
	
	return list(set(neighbours))
	
# Create an ILP input file for every sample on the list without the summary
def create_ILP_input_list(samp_list=''):
	if samp_list == '':
		samp_list = params.defaults_list_name
	#TODO: change so samps of the form SAMP_T will also work
	s_list_file = open(samp_list)
	for line in s_list_file.readlines():
		g=None
		sys.stdout.write('\rdoing '+line[:-1])
		l = line.replace('\n','\t').split('\t')
		g = build_graph(l[0],l[1])[0]
		g.create_ILP_file()
	s_list_file.close()
	

# Return a list of all chrs connected to chr in samp (via translocations)
def get_chr_list(samp, chr):
	#sanity check of the input
	if samp[-2:]=='_T':
		samp=samp[:-2]
	
	cnv_f = open(CNV_fname)
	bp_f = open(BP_fname)
	
	chrs_graph={}	#This can be done once per sample
	for ch in chromosomes:
		chrs_graph[ch] = [ch]		
	
	lines = bp_f.readlines()
	for line in lines[1:]:
		l = line.replace('\n','\t').split('\t')
		if l[10] != samp+'_T':
			continue

		if l[3] not in chrs_graph[l[0]]:
			chrs_graph[l[0]].append(l[3])
		if l[0] not in chrs_graph[l[3]]:
			chrs_graph[l[3]].append(l[0])
	
	cnv_f.close()
	bp_f.close()
	return neighbours(chr, chrs_graph)

	
def convert_karyotype_list_to_graph(kar_list_file=''):
	if kar_list_file == '':
		kar_list_file = params.defaults_list_name
	done = []
	kar_list = open(kar_list_file)
	for kar in kar_list:
		samp_name = '_'.join(kar.split('\t')[0].split('_')[:-1])
		if samp_name not in done:		
			sys.stdout.write('\rdoing '+ samp_name)
			convert_karyotype_to_graph(samp_name)
			done.append(samp_name)
	kar_list.close()

# This gets a correct solution file produce by simulations and converts it to a correct graph file format that close_solurions can understand
# The inpur file is a karyotype created by simulations, it only has chr and pos for bps and cnvs
# The output file is a correct graph with nodes and edges
def convert_karyotype_to_graph(kar_fname):
	kar_graph = open('correct\\' + kar_fname + '_correct_graph.txt','w')
	bp_g = build_graph(kar_fname+'_1','1')[0]	# The build_graph takes a chr which is mostly a legacy from the data. with simulations I just build graphs including all chrs. 
	kar_graph.write(str(bp_g.num_nodes)+'\n')	# Also note the build graph needs a full name of a sample but since the simulations add a noise param to the name we just concatenate '_1' to the name
	
	#read the correct solution data into 3 dictionaries
	simul_correct = open('correct\\' + kar_fname + '_correct.txt')
	for line in simul_correct.readlines():
		l = line.split('\t')
		type = l[0]
		chr1 = l[1]
		pos1 = int(l[2])
		chr2 = l[3]
		pos2 = int(l[4])
		dir1 = l[5]
		dir2 = l[6]
		right = int(l[7])
		left = int(l[8])
		
		if type == 'v':
			dic = bp_g.var_edges
		elif type == 'i':
			dic = bp_g.int_edges
		elif type == 'r':
			dic = bp_g.ref_edges
		else:
			print 'bad kar file: ', line
			return
			
		# int edges are a special case since the3 kar file contains "segments" and the graph interval of consecutive segments with no bp's between them
		# So check for edges with the same starting position
		if type == 'i':
			for edge in dic:
				if chr1 == edge.u.chr and pos1 == edge.u.pos:
					if chr2!=chr1 or chr2!=edge.v.chr or int2dir(edge.u.side)!='-' or int2dir(edge.v.side)!='+' or dir1!='-' or dir2!='+': #sanity check
						print 'something is not alright, int edge should be - to + on the same chr'
						print edge
						print l
						print 'chr1: ',chr1, ' chr2: ', chr2, ' edge.v.chr: ', edge.v.chr, ' edge.u.chr: ', edge.u.chr, ' edge.u.side: ', edge.u.side, ' edge.v.side: ', edge.v.side, ' dir1: ', dir1, ' dir2: ', dir2
						return
					else:
						kar_graph.write('i\t' + str(edge.u.id) + '\t' + str(edge.v.id) + '\t' + str(right) + '\t' + str(left) + '\n')
		else:						
		#Go over the var/ref edges in the graph and use the nodes positions to find their correct weight in the kar file
			for edge in dic:
				if chr1==edge.u.chr and pos1==edge.u.pos and dir1==BP_graph.side_str(edge.u.side) and chr2==edge.v.chr and pos2 == edge.v.pos and dir2==BP_graph.side_str(edge.v.side):
					kar_graph.write(type + '\t' + str(edge.u.id) + '\t' + str(edge.v.id) + '\t' + str(right) + '\t' + str(left) + '\n')
				else:	# try the reverse edge
					if chr2==edge.u.chr and pos2==edge.u.pos and dir2==BP_graph.side_str(edge.u.side) and chr1==edge.v.chr and pos1 == edge.v.pos and dir1==BP_graph.side_str(edge.v.side):
						kar_graph.write(type + '\t' + str(edge.u.id) + '\t' + str(edge.v.id) + '\t' + str(left) + '\t' + str(right) + '\n')
		
		
	kar_graph.close()
	simul_correct.close()
	os.remove('correct\\' + kar_fname + '_correct.txt')
	
# build a bp graph from the raw data
def build_graph(sample, chr):	

	if sample[-2:]=='_T':
		sample=sample[:-2]
	
	bp_g = BPGraph(sample, chr)
	if SIMULATION:
		chr_cluster = chromosomes		
	else:
		chr_cluster = get_chr_list(sample, chr)

	
	cnv_f = open(CNV_fname)
	bp_f = open(BP_fname)
	
	lines = bp_f.readlines()
	
	# Get only relevant bps in a list
	raw_bp_list = []
	for line in lines[1:]:
		l = line.replace('\n','\t').split('\t')
		if l[10] != sample+'_T':
			continue
		if l[0] not in chr_cluster:
			continue
		pos1 = (int(l[1])+int(l[2]))/2
		pos2 = (int(l[4])+int(l[5]))/2
		raw_bp_list.append((l[0], pos1, l[3], pos2, dir2int(l[7]), dir2int(l[8]), int(l[6])))
	
	# get only relevant chr in a list
	cnv_list = []
	lines = cnv_f.readlines()
	for line in lines[1:]:
		l = line.replace('\n','\t').split('\t')
		if l[9] != sample+'_T':
			continue
		if l[0] not in chr_cluster:
			continue
		pos = (int(l[1])+int(l[2]))/2
		cnv_list.append((l[0],pos,Decimal(l[6]),Decimal(l[7])))
	
	bp_g.build_graph(raw_bp_list, cnv_list)
	bp_g.dir = get_dir(sample, chr)
	
	cnv_f.close()
	bp_f.close()
	
	return bp_g, raw_bp_list, cnv_list	
	
def process_simulation_results(samp_list=''):	
	if samp_list=='':
		samp_list = params.defaults_list_name
		
	print 'processing simulation result of batch: ', simul_name
	res = {}	# dictionary of results. res[ops, eps, alpha] = [total, identical to correct, equivalent (CN and var), equal CN profile, closer to data than to correct, 
				#	sum dist to data, sum dist to correct, sum % of var edges in relation to correct, sum % of var edges according to data

	list_file = open(samp_list)
	for line in list_file:
		sys.stdout.write('\rprocessing '+ line[:-1])
		samp = line.split('\t')[0]
		chr = line.split('\t')[1].strip()		# Currently this is not used. All chrs are part of the same graph.
		
		batch = '_'.join(samp.split('_')[0:-3])
		ops = samp.split('_')[-3]
		kar_num = samp.split('_')[-2]
		eps = samp.split('_')[-1]

		for alpha in params.alpha_values:
			correct_fname = 'correct\\'+batch+'_'+ops+'_'+kar_num+'_correct_graph.txt'
			input_graph_fname = 'Results\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'_graph.txt'
			result_graph_fname = 'Results\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'_ILP_v'+ver+'_alpha'+str(alpha)+'.txt'
			
			
			graph_scores = scores.calc_graph_dist(correct_fname,input_graph_fname,result_graph_fname)
			dist_from_data = graph_scores[0]
			supp_data = graph_scores[1]
			dist_from_correct = graph_scores[2]
			supp_correct = graph_scores[3]
			dist_cor_dat = graph_scores[4]
			cn_score = graph_scores[5]

			new_res = [0,0,0,0,0,0,0,0,0,0]
			if (ops,eps,alpha) in res.keys():
				new_res = res[ops,eps,alpha]
			new_res[0] += 1
			if scores.equal_solutions(result_graph_fname, correct_fname):
				new_res[1] += 1
			else:
				rev_circles_fname = scores.reverse_circles(result_graph_fname)
				if rev_circles_fname != '':
					if scores.equal_solutions(rev_circles_fname, correct_fname):
						new_res[1] += 1
			if dist_from_correct == 0 and scores.same_var_edge_usage(result_graph_fname, correct_fname): #equal profile
				new_res[2] += 1
			
			if dist_from_correct == 0:
				new_res[3] += 1
			if dist_from_correct == 0 or dist_from_data <= dist_from_correct:
				new_res[4] += 1
			new_res[5] += dist_from_data
			new_res[6] += dist_from_correct
			new_res[7] += supp_correct
			new_res[8] += supp_data
			
			if scores.correct_for_missing_bridges(correct_fname, result_graph_fname, input_graph_fname):
				new_res[9] += 1
				
			res[ops,eps,alpha] = new_res
	
	list_file.close()
	print_simulation_results(res, samp_list+'_res')
	print '\nDone. results are in: ', samp_list+'_res'


def process_simulation_results_for_cn_score(samp_list=''):
	
	if samp_list=='':
		samp_list = params.defaults_list_name
		
	print 'processing simulation result of batch: ', simul_name
	res = {}	# dictionary of results. res[ops, eps, alpha] = [total, identical to correct, equivalent (CN and var), equal CN profile, closer to data than to correct, 
				#	sum dist to data, sum dist to correct, sum % of var edges in relation to correct, sum % of var edges according to data

	res_file = open(samp_list+'_res_score.txt','w')
	res_file.write('samp\tnum ops\tnoise\talpha\tcn score\tsupp score\n')
	list_file = open(samp_list)
	for line in list_file:
		sys.stdout.write('\rdoing '+ line[:-1])
		samp = line.split('\t')[0]
		chr = line.split('\t')[1].strip()		# Currently this is not used. All chrs are part of the same graph.
		
		batch = '_'.join(samp.split('_')[0:-3])
		ops = samp.split('_')[-3]
		kar_num = samp.split('_')[-2]
		eps = samp.split('_')[-1]

		for alpha in params.alpha_values:
			correct_fname = 'correct\\'+batch+'_'+ops+'_'+kar_num+'_correct_graph.txt'
			input_graph_fname = 'Results\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'_graph.txt'
			result_graph_fname = 'Results\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'_ILP_v'+ver+'_alpha'+str(alpha)+'.txt'
			
			graph_scores = scores.calc_graph_dist(correct_fname,input_graph_fname,result_graph_fname)
			dist_from_data = graph_scores[0]
			supp_data = graph_scores[1]
			dist_from_correct = graph_scores[2]
			supp_correct = graph_scores[3]
			dist_cor_dat = graph_scores[4]
			cn_score = graph_scores[5]
			
			res_file.write(samp+'\t'+ops+'\t'+str(params.eps_cnv[int(eps)-1])+'\t'+str(alpha)+'\t'+str(cn_score)+'\t'+str(supp_correct)+'\n')

				
	
	list_file.close()
	res_file.close()
	print '\nDone. results are in: ', samp_list+'_res_score.txt'

# Process the simulation results for the base scenario where we want to split the big batch into smaller batches and calc the rate for each batch
def process_simulation_results_simul10(samp_list=''):
	
	if samp_list=='':
		samp_list = params.defaults_list_name
	
	res_f = open(samp_list+'_batch_res','w')
	res_cn = open(samp_list+'_batch_cn_res','w')	#holding all cn_scores regardless of batch
	res_cn.write('samp_name\tnum_ops\tepsilon\talpha\tsamp\tcn_score\n')
	print 'processing simulation result of batch: ', simul_name
	res = {}	# dictionary of results. res[batch_num] = [total, equivalent (CN and var), equal CN profile, EBS, cn_score] 
	
	samp_num = 0
	tot_cn_score = 0
	list_file = open(samp_list)
	for line in list_file:
		batch_num = samp_num/100
		samp_num+=1
		
		sys.stdout.write('\rdoing '+ line[:-1])
		samp = line.split('\t')[0]
		chr = line.split('\t')[1].strip()		# Currently this is not used. All chrs are part of the same graph.
		
		batch = '_'.join(samp.split('_')[0:-3])
		ops = samp.split('_')[-3]
		kar_num = samp.split('_')[-2]
		eps = samp.split('_')[-1]

		alpha = params.alpha_values[0]
		correct_fname = 'correct\\'+batch+'_'+ops+'_'+kar_num+'_correct_graph.txt'
		input_graph_fname = 'Results\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'_graph.txt'
		result_graph_fname = 'Results\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'_ILP_v'+ver+'_alpha'+str(alpha)+'.txt'
		
		
		graph_scores = scores.calc_graph_dist(correct_fname,input_graph_fname,result_graph_fname)
		dist_from_data = graph_scores[0]
		supp_data = graph_scores[1]
		dist_from_correct = graph_scores[2]
		supp_correct = graph_scores[3]
		dist_cor_dat = graph_scores[4]
		cn_score = graph_scores[5]
		
		if batch_num not in res.keys():
			res[batch_num] = [0,0,0,0,0,Decimal(0.0)]

		res[batch_num][0] += 1
		
		if dist_from_correct == 0 and scores.same_var_edge_usage(result_graph_fname, correct_fname): #equal profile
			res[batch_num][1] += 1
		
		if dist_from_correct == 0:
			res[batch_num][2] += 1
		if dist_from_correct == 0 or dist_from_data <= dist_from_correct:
			res[batch_num][3] += 1
		if scores.correct_for_missing_bridges(correct_fname, result_graph_fname, input_graph_fname):
				res[batch_num][4] += 1
		res[batch_num][5] += cn_score
		
		res_cn.write(batch+'\t'+ops+'\t'+str(params.eps_cnv[int(eps)-1])+'\t'+str(alpha)+'\t'+str(samp_num)+'\t'+str(cn_score)+'\n')
	
	res_f.write('batch\tnum_ops\tepsilon\talpha\ttotal\tcorrect\tECN\tEBS\tcorrect for missing bridges\taverage cn_score\n')
	for batch_num in res.keys():
		total = res[batch_num][0]
		correct = Decimal(res[batch_num][1])/total
		ECN = Decimal(res[batch_num][2])/total
		EBS = Decimal(res[batch_num][3])/total
		miss_brdg = Decimal(res[batch_num][4])/total
		cn_s = Decimal(res[batch_num][5])/total
		
		res_f.write(str(batch_num)+'\t'+ops+'\t'+str(params.eps_cnv[int(eps)-1])+'\t'+str(alpha)+'\t'+str(total)+'\t'+str(correct*100)+'%\t'+str(ECN*100)+'%\t'+str(EBS*100)+'%\t'+str(miss_brdg*100)+'%\t'+str(cn_s)+'\n')
	list_file.close()
	res_f.close()
	res_cn.close()
	print '\nDone. results are in: ', samp_list+'_batch_res and ', samp_list+'_batch_cn_res'
	
def print_simulation_results(res, fname):
# dictionary of results. res[ops, eps, alpha] = [total, identical to correct, equal CN profile, closer to data than to correct, 
				#	sum dist to data, sum dist to correct, sum % of var edges in relation to correct, sum % of var edges according to data
	out_f = open(fname,'w+')
	out_f.write('num_ops\tepsilon\talpha\t\ttotal\tcorrect\tequivalent(CN and var)\tequal CN\tcloser to data\tcorrect_for_missing_bridges\tavg dist from data\tavg dist from correct\tavg % of correct supp used\tavg % of data supp used\n')
	for ops,eps,alpha in res.keys():
		total, correct, equiv, equal, close, dist_data, dist_correct, supp_correct, supp_data, correct_for_missing_brdg = res[ops,eps,alpha]
		out_f.write(str(ops)+'\t'+str(params.eps_cnv[int(eps)-1])+'\t'+str(alpha)+'\t\t'+str(total)+'\t'+str(100*Decimal(correct)/total)+'%\t'+str(100*Decimal(equiv)/total)+'%\t'+str(100*Decimal(equal)/total)+'%\t'+str(100*Decimal(close)/total)+'%\t'+str(100*Decimal(correct_for_missing_brdg)/total)+'%\t'+str(dist_data/total)+'\t'+str(dist_correct/total)+'\t'+str(100*supp_correct/total)+'%\t'+str(100*supp_data/total)+'%\n')
	out_f.close()

	
# Return a working dir for a samp and chr - mainly to organize things
def get_dir(samp, chr):
	dir = os.getcwd() + '\\Results\\' + samp + '\\' + samp + '_' + chr
	if not os.path.exists(dir):
		os.makedirs(dir)
	return dir

def apply_bridge_filters(samp_list = ''):
	if samp_list == '':
		samp_list = params.defaults_list_name
	s_list_file = open(samp_list)
	for line in s_list_file.readlines():
		samp = line.split('\t')[0]
		chr = line.split('\t')[1].strip()		# Currently this is not used. There is one graph for all chrs even if they're not connected
		batch = '_'.join(samp.split('_')[0:-3])
		ops = samp.split('_')[-3]
		kar_num = samp.split('_')[-2]
		eps = samp.split('_')[-1]

		input_graph_fname = 'Results\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'\\'+batch+'_'+ops+'_'+kar_num+'_'+eps+'_'+chr+'_graph.txt'
		simulations.filter_bridges(input_graph_fname)
	
# Create a batch of simulated rearranged karyotypes
def create_simulation_batch(name_param=''):
	if (SIMULATION):
		if name_param=='':
			name = params.simul_name
		else:
			name = name_param
		print 'creating simulation batch ', name
		print 'operations: ', min_ops, '-',max_ops
		print num_kars, 'random karyotypes for each number of operations'
		print 'each with nosise values (epsilon): ', eps_cnv
		if name_param == '':
			simulations.create_simulated_karyotypes(max_ops, min_ops)
		else:
			simulations.create_simulated_karyotypes(max_ops, min_ops, name)
		print '\nCreated bp and cnv files: ', CNV_fname, ', ', BP_fname
		print 'Full list of batch is in: ', defaults_list_name, '. full kars are in: ', simul_name+'_debug.txt'
	
	print 'Preparing graph files for ILP run'
	if name_param=='':
		create_ILP_input_list()
	else:
		create_ILP_input_list(name+'_list.txt')
	print '\nILP input files created'
	
	
	print '\nApplying bridge filters'
	time.sleep(3)
	if name_param=='':
		apply_bridge_filters()
	else:
		apply_bridge_filters(name+'_list.txt')
		
	print '\nDone. ILP can be run'
	print 'Params for the ILP:'
	print 'smaples_file: ', defaults_list_name
	print 'ver: ', ver
	print 'alpha_values', alpha_values
	print 'When ILP is done, ready to run process_simulation_results'

def clean_all():
	try:
		os.remove(simul_name + '_bp.txt')
	except OSError:
		pass
	try:
		os.remove(simul_name + '_cnv.txt')
	except OSError:
		pass
	try:
		os.remove(simul_name + '_list.txt')
	except OSError:
		pass
	shutil.rmtree('correct/', ignore_errors=True)
	shutil.rmtree('Results/', ignore_errors=True)

def help():
	print 'please specify what to do as the first parameter: '
	print '1: create simulation batch'
	print '2: create correct karyotype files'
	print '3: process simulation results'
	print '4: process simulation results for base scenarion'
	print '5: process simulation results to produce cn score'
	print '6: clean all files'
	
if __name__ == '__main__':
	func = {'1': create_simulation_batch,
			'2': convert_karyotype_list_to_graph,
			'3': process_simulation_results,
			'4': process_simulation_results_simul10,
			'5': process_simulation_results_for_cn_score,
			'6': clean_all
			}
	
	print 'simulation name is : ', params.simul_name
	if not os.path.exists('correct/'):
		os.makedirs('correct/')
	if not os.path.exists('Results/'):
		os.makedirs('Results/')
		
	if len(sys.argv) < 2:
		help()
	else:
		try:
			func[sys.argv[1]]()
		except KeyError:
			help()
	
	
# Use the following functions to translate a solution file into a gv format
##create a gv for a specific result file
def build_gv_sample(samp, chr, type = 'solution', alpha=0.5, v=''):
	if v=='':
		v = params.ver
	g = build_graph(samp, chr)[0]
	if type == 'solution':
		res = os.getcwd() + '\\Results\\' + samp + '\\'+ samp + '_' + chr + '\\' + samp + '_' + chr + '_ILP_v' + v + '_alpha' + str(alpha) + '.txt'
	elif type == 'correct':
		res = os.getcwd() + '\\Results\\' + samp + '_' + chr + '\\' + samp + '_' + chr + '_correct.txt'
	elif type=='rev':	# for reversed circles
		res = os.getcwd() + '\\Results\\' + samp + '\\'+ samp + '_' + chr + '\\' + samp + '_' + chr + '_ILP_v' + v + '_alpha' + str(alpha) + '.txt_rev.txt'
	elif type == 'reversed':
		res = os.getcwd() + '\\Results\\' + samp + '\\'+ samp + '_' + chr + '\\' + samp + '_' + chr + '_ILP_v' + v + '_alpha' + str(alpha) + '.txt_reversed'
	g.build_results_gv(res)

 create a gv for input graph for every sample in list
def create_gv_input_list(samp_list):
	s_list_file = open(samp_list)
	for line in s_list_file.readlines():
		g=None
		print 'doing ',line
		l = line.replace('\n','\t').split('\t')
		g, bps, cnvs = build_graph(l[0],l[1])
		g.to_gv()
	s_list_file.close()
