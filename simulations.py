import sys
import os
import random
from decimal import *
from params import *
import CreateGraph
import numpy

BP_READ_LEN = 500
CNV_WINDOW_LEN = 5000

				

# Read a raw bp description of a simulation and create a BP.txt format from it
# NOTE - fields which are not used by CreateGraph are for now printed as 'TEST'
def create_simulation_BP(samp_file):
	sampf = open(samp_file)
	outf = open(samp_file+'.txt','w+')
	outf.write('chrom\tpos1\tpos2\tchrom2\tpos3\tpos4\tsupport\tstrand1\tstrand2\tvariantClass\tsample(s) the break was found in\n')
	for line in sampf.readlines()[1:]:
		l = line.split('\t')
		for i in range(len(eps_cnv)+1)[1:]:
			bp_supp = int(round(supp_noise()*Decimal(l[5])))
			if bp_supp>0:	# if bp_supp==0 that means it was too noisy to see it and so there should be no record of it at all #TODO: I don't think that happens at all anymore
				outf.write(l[1]+ '\t' + str(int(l[2])-BP_READ_LEN)+ '\t' + str(int(l[2])+BP_READ_LEN) + '\t' + l[3] + '\t' + str(int(l[4])-BP_READ_LEN)+ '\t' + str(int(l[4])+BP_READ_LEN) + '\t' + str(bp_supp) + '\t' + l[6] + '\t' + l[7].strip() + '\tTEST\t' + l[0] +'_' + str(i) + '_T\n')
	
	sampf.close()
	outf.close()

#Draw a support score from an exponential distribution
def supp_noise():
	return int(round(Decimal(0.5+random.expovariate(supp_dist_lambda))))
	
# simulate missing bridges by removing bridges from the specified graph	
def filter_bridges(fname):
	inp_file = open(fname)
	out_file = open('tmp','w+')
	out_file.write(inp_file.readline())
	for line in inp_file.readlines():
		if line.split('\t')[0] == 'v':
			if random.random() < p_miss_brdg:
				continue
		out_file.write(line)
	out_file.close()
	inp_file.close()
	os.remove(fname)
	os.rename('tmp',fname)
			

# Read a raw CNV description of a simulation and create a CNV.txt format from it for each epsilon value
def create_simulation_CNV(samp_file):
	sampf = open(samp_file)
	outf = open(samp_file+'.txt','w+')
	outf.write('Chromosome\tStart\tEnd\tLeft Segment Zscore\tRight Segment Zscore\tDelta Zscore\tLeft Segment Copy Number\tRight Segment Copy Number\tDatasetSpecID\tDatasetID\n')
	for line in sampf.readlines()[1:]:
		l = line.split('\t')
		for i,eps in enumerate(eps_cnv):
			outf.write(l[1] + '\t' + str(int(l[2])-CNV_WINDOW_LEN)+ '\t' + str(int(l[2])+CNV_WINDOW_LEN) + '\tTEST\tTEST\tTEST\t' + str(cnv_noise(Decimal(l[3]), eps)) + '\t' + str(cnv_noise(Decimal(l[4]), eps)) + '\tTEST\t' + l[0] + '_' + str(i+1) + '_T\n')
	
	sampf.close()
	outf.close()
	
def cnv_noise(cn, eps):
	return max(0, random.normalvariate(float(cn), eps))
	

test_file = 'Test'
def run():
	bp_file = test_file+'_BP.txt'
	cnv_file = test_file+'_CNV.txt'
	create_simulation_BP(bp_file)
	create_simulation_CNV(cnv_file)
	
# Return a rearranged karyotype. the inpout is the number of rearrangements. This also returns a secondary karyotype to serve as the "litter" of the simulation
def create_random_karyotype(operations):
	# randomly do operations on it
	kar = {}
	noise_kar = {}	#The secondary kar littering the sample
	for chr in chromosomes:
		kar[chr] = []
		noise_kar[chr] = []
		for copy in range(DEFAULT_CN):
			kar[chr].append([(chr,x) for x in range(num_segs+1)[1:]])
			noise_kar[chr].append([(chr,x) for x in range(num_segs+1)[1:]])
	for i in range(operations):
		op = random.choice(rearrangement)
		kar, pos = do_operation(kar, op)
		noise_kar_option = random.choice(noise_kar_options)
		if noise_kar_option == 'dif':
			op = random.choice(rearrangement)
			noise_kar = do_operation(noise_kar, op)[0]
		elif noise_kar_option == 'same':
			noise_kar = do_operation(noise_kar, op, pos)[0]
		
	return kar, noise_kar		# We don't add first and last 'technical' segments, but rather count them directly when creating the correct (count '[1' and 'n]')
	
# Perform the rearrangement op on kar
def do_operation(kar, op, op_md=[]):
	if op_md==[]:
		new_op = True
	else:
		new_op = False
	if not new_op:
		chr1 = op_md[0]
		chr1_copy = op_md[1]
		pos1 = op_md[2]
		pos2 = op_md[3]
		if op=='TRL':
			chr2=op_md[4]
			chr2_copy=op_md[5]
		
	else:
		chr1 = random.choice(kar.keys())
		chr1_copy = random.randint(0,len(kar[chr1])-1)	#Which copy of the chr to take
	
	kar1 = kar[chr1][chr1_copy]
	
	if len(kar1)<=2:
		return kar, []
	if new_op:
		pos1 = random.randint(1,len(kar1)-1)	# In accordance with breakpoint graphs we make sure simulations are not done on the edges
	
	if op=='TRL':	# Randomly pick a second chr
		if new_op:
			chr2 = random.choice(kar.keys())
			chr2_copy = random.randint(0,len(kar[chr2])-1)
		
		if chr1==chr2 and chr1_copy==chr2_copy:
			return kar, [chr1, chr1_copy, pos1, 0, chr2, chr2_copy]		# No point in doing a TRL on the same chr copy. This means that in some cases an operation will be moot. for a diploid 23 chrs kar this only happens for 1 in 46 cases though
		kar2 = kar[chr2][chr2_copy]
		if len(kar2)<=2:
			return kar, []
		if new_op:
			pos2 = random.randint(1,len(kar2)-1)
		ret_md = [chr1, chr1_copy, pos1, pos2, chr2, chr2_copy]
		kar1_a = kar1[0:pos1]
		kar1_b = kar1[pos1:]
		kar2_a = kar2[0:pos2]
		kar2_b = kar2[pos2:]
		if random.randint(0,1) == 0:	# Randomly pick an orientation for the trl
			kar1 = kar1_a + kar2_b
			kar2 = kar2_a + kar1_b
		else:
			kar1 = kar1_a + rev_kar(kar2_a)
			kar2 = rev_kar(kar1_b) + kar2_b
		
		kar[chr2][chr2_copy] = kar2
			
	else:
		if new_op:		# Make sure pos2 is also between 1 and len(kar)-1
			pos2 = random.randint(min(pos1+1,len(kar1)-1),min(pos1+max_dist, len(kar1)-1))		# This method skewes the randomness a bit as anything above the maximun segment number is truncated, so if the last seg is 30, pos2 has a higher probability to be 30
		ret_md = [chr1, chr1_copy, pos1, pos2]
		if op=='DEL':
			kar1 = kar1[0:pos1]+kar1[pos2:]
		elif op=='DUP':
			kar1 = kar1[0:pos2]+kar1[pos1:]
		elif op=='INV':
			kar1 = kar1[0:pos1]+rev_kar(kar1[pos1:pos2])+kar1[pos2:]
	
	kar[chr1][chr1_copy] = kar1
	
	return kar, ret_md

def rev_kar(kar):
	return [(x[0],-1*x[1]) for x in reversed(kar)]

# calls 'create_simulated_karyotypes' many times
def create_simulated_karyotypes(max_num_ops, min_num_ops=1, name='', debug=False):
	if name=='':
		simulation_filename = simul_name
	else:
		simulation_filename = name
	
	tmp_cnv_fname = simulation_filename + '_cnv'
	tmp_bp_fname = simulation_filename + '_bp'
	tmp_cnv = open(tmp_cnv_fname, 'w+')
	tmp_bp = open(tmp_bp_fname, 'w+')
	scheme_list = open(simulation_filename + eps_dist + '_list.txt', 'w+')
	if debug:
		dbg_f = open(simulation_filename + '_debug','w+')
		dbg_f.write('num_ops\tidx\tchr\tcopy\tkar\n')
		dbg_f_noise = open(simulation_filename + '_noise_debug','w+')
		dbg_f_noise.write('num_ops\tidx\tchr\tcopy\tkar\n')
	
	tmp_cnv.write('sample\tchr\tpos\tleft CN\tright CN\n')
	tmp_bp.write('sample\tchr1\tpos1\tchr2\tpos2\tused?\tstrand1\tstrand2\n')
	for num_ops in range(max_num_ops+1)[min_num_ops:]:
		sys.stdout.write('\rdoing ops '+ str(num_ops))
		for i in range(num_kars):
			samp_name = simulation_filename + '_RAND_' + str(num_ops) + '_' + str(i)
			correct_f = open('correct\\' + samp_name + '_correct.txt','w+')
			
			# Generate a radnom karyotype
			full_kar, full_kar_noise = create_random_karyotype(num_ops)
			if debug:
				print_kar(dbg_f, full_kar, num_ops, i)
				print_kar(dbg_f_noise, full_kar_noise, num_ops, i)
			
			# init the dics to hold the cnv and bp data
			cnv_dic_clean = {}	#cnv_dic[chr, seg] = [right cn, left cn]
			cnv_dic_noise = {}
			cnv_dic = {}
			for chr in full_kar.keys():
				for seg in range(num_segs+1)[1:]:
					cnv_dic_clean[chr, seg] = [0,0]
					cnv_dic_noise[chr, seg] = [0,0]
					cnv_dic[chr, seg] = [0,0]
			ref_dic = {}	# ref_dic[chr, prev seg] = [right cn, left cn]
			for chr in full_kar.keys():
				for seg in range(num_segs+1)[1:-1]:
					ref_dic[chr, seg] = [0,0]
			bp_dic = {}	# bp_dic[chr1,pos1,chr2,pos2,dir1,dir2] = [right cn, left cn]
			bp_dic_noise = {}
			
			for chr in full_kar.keys():
				for kar in full_kar[chr]:	# For each chr copy
					bp_dic = find_discordant_reads(kar, bp_dic)
					ref_dic = find_concordant_reads(kar, ref_dic)
					cnv_dic_clean = calc_cn(kar, cnv_dic_clean)
			for chr in full_kar_noise.keys():
				for kar in full_kar_noise[chr]:	# For each chr copy
					bp_dic_noise = find_discordant_reads(kar, bp_dic_noise)
					cnv_dic_noise = calc_cn(kar, cnv_dic_noise)
			
			# create the cnv for the hetrogenous sample
			for chr in full_kar.keys():
				for seg in range(num_segs+1)[1:]:
					cnv_dic[chr,seg] = (1-noise_kar_ratio)*sum(cnv_dic_clean[chr,seg]) + (noise_kar_ratio)*sum(cnv_dic_noise[chr,seg])
			
					
			# After going over all chrs in full_kar and updating the bp_dic, write discordant reads to file
			for chr1,pos1,chr2,pos2,dir1,dir2 in bp_dic.keys():
				if (chr1,pos1,chr2,pos2,dir1,dir2) in bp_dic_noise:
					noise_supp = (noise_kar_ratio)*sum(bp_dic_noise[chr1,pos1,chr2,pos2,dir1,dir2])
				else:
					noise_supp = 0
				supp = (1-noise_kar_ratio)*sum(bp_dic[chr1,pos1,chr2,pos2,dir1,dir2]) + noise_supp
				tmp_bp.write(samp_name + '\t' + str(chr1) + '\t' + str(pos1) + '\t' + str(chr2) + '\t' + str(pos2) + '\t' + str(supp) + '\t' + dir1 + '\t' + dir2 + '\n')
				correct_f.write('v\t' + str(chr1) + '\t' + str(pos1) + '\t' + str(chr2) + '\t' + str(pos2) + '\t' + dir1 + '\t' + dir2 + '\t' + str(bp_dic[chr1,pos1,chr2,pos2,dir1,dir2][0]) + '\t' + str(bp_dic[chr1,pos1,chr2,pos2,dir1,dir2][1]) + '\n')
			if noise_kar_ratio>0:	# Add discordant reads that exist in the noisy kar only. but not if its ratio is 0% because than it's just adding var edges (even if their weight is 0 the algorithm can still use them)
				for chr1,pos1,chr2,pos2,dir1,dir2 in bp_dic_noise.keys():
					if (chr1,pos1,chr2,pos2,dir1,dir2) not in bp_dic:	#A bp that exists only in the noise karyotype
						supp = (noise_kar_ratio)*sum(bp_dic_noise[chr1,pos1,chr2,pos2,dir1,dir2])
						tmp_bp.write(samp_name + '\t' + str(chr1) + '\t' + str(pos1) + '\t' + str(chr2) + '\t' + str(pos2) + '\t' + str(supp) + '\t' + dir1 + '\t' + dir2 + '\n')
					
			# Write CNV to file. This is done for each chr separately 
			for chr in full_kar.keys():
				left_cn = DEFAULT_CN
				for seg in range(num_segs+1)[1:]:
					if cnv_dic[chr, seg] != left_cn:
						right_cn = cnv_dic[chr, seg]
						tmp_cnv.write(samp_name + '\t' + str(chr) + '\t' + str(position(seg, 'sta')) + '\t' + str(left_cn) + '\t' + str(right_cn) + '\n')
						left_cn = right_cn
									
			# Write the correct solution to file
			for chr, seg in cnv_dic_clean:
				correct_f.write('i\t' + str(chr) + '\t' + str(position(seg, 'sta')) + '\t' + str(chr) + '\t' + str(position(seg, 'end')) + '\t-\t+\t' + str(cnv_dic_clean[chr, seg][0]) + '\t' + str(cnv_dic_clean[chr, seg][1]) + '\n')
			for chr, seg in ref_dic:
				correct_f.write('r\t' + str(chr) + '\t' + str(position(seg, 'end')) + '\t' + str(chr) + '\t' + str(position(seg+1, 'sta')) + '\t+\t-\t' + str(ref_dic[chr,seg][0]) + '\t' + str(ref_dic[chr,seg][1]) + '\n')
			
			# Compose a list of all the samples created (Easy to batch run the other scripts)
			for i in range(len(eps_cnv)+1)[1:]:
				scheme_list.write(samp_name + '_' + str(i) + '\t1\n')		# originally I composed a list of samp \t chr but that makes the correct files complicated so instead I just build graphs of the whole 23 chrs.
				
				#for chr in chrs_changed(bp_dic):
				#	scheme_list.write(samp_name + '\t' + str(i) + '\t' + str(chr) + '\n')
			correct_f.close()
			
			
	
	tmp_cnv.close()
	tmp_bp.close()
	scheme_list.close()
	if debug:
		dbg_f.close()
		dbg_f_noise.close()

	create_simulation_BP(tmp_bp_fname)
	create_simulation_CNV(tmp_cnv_fname)
	if not debug:
		os.remove(tmp_bp_fname)
		os.remove(tmp_cnv_fname)	

# Find discordant reads and produce a bp_dic out of them
# bp_dic[chr1,pos1,chr2,pos2,dir1,dir2] = [right cn, left cn]
def find_discordant_reads(kar, bp_dic):
	for ind, seg in enumerate(kar[:-1]):
		chr1 = seg[0]
		chr2 = kar[ind+1][0]
		seg1_num = seg[1]
		seg2_num = kar[ind+1][1]
		if chr1!=chr2 or seg1_num+1 != seg2_num:		# This is a discordant read
			# Determine which sides of the segments compose the discordant read
			if seg1_num > 0:	
				end1 = 'end'	# connect the end side of seg 1 (head)
				dir1 = '+'		# If the first segment is facing forward this is the plus side of the bp
			else:
				end1 = 'sta'
				dir1 = '-'
			pos1 = position(abs(seg1_num),end1)
			
			if seg2_num > 0:
				end2 = 'sta'	# connect the start side of seg 2 (tail)
				dir2 = '-'
			else:
				end2 = 'end'
				dir2 = '+'
			pos2 = position(abs(seg2_num),end2)
			
			# Add the discordant read to the bp_dic
			if (chr1,pos1,chr2,pos2,dir1,dir2) in bp_dic:
				bp_dic[chr1,pos1,chr2,pos2,dir1,dir2][0] += 1
			elif (chr2,pos2,chr1,pos1,dir2,dir1) in bp_dic:		#reverse traversal is still the same bp
				bp_dic[chr2,pos2,chr1,pos1,dir2,dir1][1] += 1
			else:
				bp_dic[chr1,pos1,chr2,pos2,dir1,dir2] = [1,0]
				
	return bp_dic
	
# ref_dic[chr, prev seg] = [right cn, left cn]
def find_concordant_reads(kar, ref_dic):
	for ind, seg in enumerate(kar[:-1]):
		chr1 = seg[0]
		chr2 = kar[ind+1][0]
		seg1_num = seg[1]
		seg2_num = kar[ind+1][1]
		if chr1==chr2 and seg1_num+1 == seg2_num:		# This is a concordant read
			# Determine which sides of the segments compose the read
			if seg1_num > 0:
				ref_dic[chr1, seg1_num][0] += 1
			else:	# we a connection of type -2 ==> -1
				ref_dic[chr1, abs(seg2_num)][1] += 1
	return ref_dic

def calc_cn(kar, cnv_dic):
	for seg in kar:
		chr = seg[0]
		seg_num = seg[1]
		if seg_num > 0:
			cnv_dic[chr, seg_num][0] += 1
		else:
			cnv_dic[chr, abs(seg_num)][1] += 1
	return cnv_dic

# get a seg number and a locator from ['sta','mid','end'] and return the position of its begining, middle or end
def position(seg, locator):
	if locator == 'sta':
		if seg == 1:		#special case
			return 0
		return int(min_seg_size*(seg-0.5))
	if locator == 'mid':
		return int(min_seg_size*seg)
	if locator == 'end':
		return int(min_seg_size*(seg+0.5))
	print 'bad segment locator: ',locator
		
	
def print_kar(dbg_f, kar, num_ops, idx):
	for chr in kar.keys():
		for copy_idx, copy in enumerate(kar[chr]):
			dbg_f.write(str(num_ops) + '\t' + str(idx) + '\t' + str(chr) + '\t' + str(copy_idx) +':\t' + str(copy) + '\n')

# get a bp dic of a simulation and return a list of chrs representing disconnected graph (ie only one chr for each connected component in the graph)
def chrs_changed(bp_dic):
	res = []
	used = []
	for bp in bp_dic:
		chr1 = bp[0]
		chr2 = bp[2]
		if chr1 not in used:
			used.append(chr1)
		if chr2 not in used:
			used.append(chr2)
		if chr1 not in res and chr2 not in res:
			res.append(chr1)
		elif chr1 in res and chr2 in res:
			res.remove(chr2)
	return res
	
#Asses how much of the karyotype is deleted for any number of rearrangements
def create_missing_karyotype_data(out_fname):
	max_op_num = 30
	num_iter = 100
	res = {}
	res_del = {}
	db={}
	for num_op in range(max_op_num+1)[1:]:
		res[num_op] = []
		res_del[num_op] = []
		db[num_op]=[]
		for i in range(num_iter):
			full_kar, full_kar_noise = create_random_karyotype(num_op)
			cnv_dic = {}	
				
			for chr in full_kar.keys():
				for seg in range(num_segs+1)[1:]:
					cnv_dic[chr, seg] = [0,0]		
			for chr in full_kar.keys():
				for kar in full_kar[chr]:	# For each chr copy
					cnv_dic = calc_cn(kar, cnv_dic)
			orig = 0
			cur = 0
			del_segs = 0
			for chr in full_kar.keys():
				for seg in range(num_segs+1)[1:]:
					orig += 2
					cur += min(2, sum(cnv_dic[chr, seg]))
					if sum(cnv_dic[chr, seg])<2:
						del_segs += 1
			res[num_op].append(float(float(cur)/float(orig)))
			res_del[num_op].append(del_segs)
			db[num_op].append((cur,orig, del_segs))
			
			
	
	out_f = open(out_fname, 'w+')
	out_f.write('num ops\tnum iter\taverage genomic content remained\tstandart deviation\tsegs deleted average\tsegs deleted std\n')
	for num_op in res:
		out_f.write(str(num_op) + '\t' + str(num_iter) + '\t' + str(numpy.mean(res[num_op])) + '\t' + str(numpy.std(res[num_op])) + '\t' + str(numpy.mean(res_del[num_op])) + '\t' + str(numpy.std(res_del[num_op])) + '\n')
	out_f.close()
	
	return res,db