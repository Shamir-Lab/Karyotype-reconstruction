import sys
import os
from decimal import *
import bisect
from params import *

'''
	This represents a bridge graph read from the data. Edges are tuples of nodes.
'''
class BPGraph:
	def __init__(self, samp, chr):
		self.dir = ''
		self.int_edges = []
		self.var_edges = []
		self.ref_edges = []
		self.num_nodes = 0
		self.num_chr = 0
		self.chrs = []
		self.chr_lengths = {}	# The length of each chromosome area of interest
		self.nodes = {}		#A dict of (ch,pos) -> (Node, Node) ; each pos has two sides
		self.samp = samp
		self.weights = {}	# A dict of double lists by chr. first list is of all positions and last one of average CN from that position onward
		self.chr = chr
		self.total_support = 0	# The total reads supporting all the var edges in this graph
		self.total_length = 0

	# calculate the lengths of chromsomes in the graph
	def calc_chr_lengths(self):
		for chr in self.chrs:
			self.chr_lengths[chr] = 0
		for e in self.int_edges:
			self.chr_lengths[e.u.chr] += e.get_len()
		for chr in self.chrs:
			self.total_length += self.chr_lengths[chr]
	
	#The major function - gets a list of BPs and CNV from the data files and builds it into a graph
	def build_graph(self, raw_bp_list, raw_cnv_list):
		self.build_nodes_dict(raw_bp_list)		# Build nodes dict - group adjacent breaks together as one node, this also creates the int and ref edges
		self.build_var_edges(raw_bp_list)		# Analyse the bp list and create all the var edges
		self.calc_int_weights(raw_cnv_list)		# Analyse the cnv list and infer the avg weights of the int edges
		self.calc_chr_lengths()

	# Write the graph into a file the ILP java code can read.
	def create_ILP_file(self):
		graph_f = open(self.dir + '\\' + self.samp + '_' + self.chr + graph_fname, 'w+')
		graph_f.write(str(self.num_nodes) + '\n')
		for e in self.int_edges:
			graph_f.write('i\t' + str(e.u.id) + '\t' + str(e.v.id) + '\t' + str(e.cn) + '\t' + str(Decimal(e.get_len())/self.total_length))
			if e.u.term:	#starting edge
				graph_f.write('\t->')
			elif e.v.term:	#end edge
				graph_f.write('\t->')
			graph_f.write('\n')
		for e in self.var_edges:
			if self.total_support == 0:
				normalized_supp = 0
			else:
				normalized_supp = Decimal(e.sup)/self.total_support
			graph_f.write('v\t' + str(e.u.id) + '\t' + str(e.v.id) + '\t' + str(normalized_supp) + '\n')
		for e in self.ref_edges:
			graph_f.write('r\t' + str(e.u.id) + '\t' + str(e.v.id) + '\n')
		graph_f.close()
	
	# Get a list of relevant bps (chr1, pos1, chr2, pos2, dir1, dir2, sup) and create the graph nodes, refrence and interval edges.
	def build_nodes_dict(self, bp_list):
		positions = {}		# a list of all breakpoints per chromosme
		for bp in bp_list:
			if bp[0] not in positions.keys():
				positions[bp[0]] = []
				self.chrs.append(bp[0])
				self.num_chr+=1
			if bp[2] not in positions.keys():
				positions[bp[2]] = []
				self.chrs.append(bp[2])
				self.num_chr+=1
			positions[bp[0]].append(bp[1])
			positions[bp[2]].append(bp[3])
		
		# Group adjacent breakpoints (AKA breakpoint filter)
		for ch in positions.keys():
			node2 = None
			while positions[ch] != []:
				adj_pos = [p for p in positions[ch] if p-min(positions[ch])<MAX_BP_WINDOW]	#get all bps that are close enough to the smallest one. this makes sure we do it in order
				avg_pos = sum(adj_pos)/len(adj_pos)
				node1 = Node(self.node_id_generator(), ch, avg_pos, 0)	#Each break has two nodes related to it
				# before giving node2 a new value we connet the older one with an interval edge. This is why it's important to do the positions in order
				if node2 != None:
					self.int_edges.append(IntEdge(node2, node1))
				else:	#Insert a start node for the ch
					self.int_edges.append(IntEdge(Node(self.node_id_generator(), ch, 0, 1, True), node1))
				node2 = Node(self.node_id_generator(), ch, avg_pos, 1)
				for pos in adj_pos:
					self.nodes[(ch,pos)] = node1,node2
				self.ref_edges.append(RefEdge(node1,node2))
				positions[ch] = [p for p in positions[ch] if p not in adj_pos]
			#Add end node for ch
			self.int_edges.append(IntEdge(node2, Node(self.node_id_generator(), ch, DEFAULT_END_POSITION, 0, True)))	#The default end position is used for simulations, I think the algorithm in general disregard it as it knows its a terminator edge
			
	# Calculate the weights of the interval edges from a list of CNV's
	def calc_int_weights(self, cnv_list):
		tmp_weights_dict = {}
		#read all CNV into a dict of sorted lists by chr
		for cnv in cnv_list:
			if cnv[0] not in tmp_weights_dict.keys():
				tmp_weights_dict[cnv[0]] = []
			tmp_weights_dict[cnv[0]].append(tuple(cnv[1:]))

		
		for ch, w_list in tmp_weights_dict.iteritems():
			w_list.sort()
			self.weights[ch] = ([],[])
			self.weights[ch][0].append(0)
			self.weights[ch][1].append(Decimal(w_list[0][1]))
			for i, cn in enumerate(w_list[:-1]):
				self.weights[ch][0].append(cn[0])
				self.weights[ch][1].append(Decimal(cn[2] + w_list[i+1][1])/2)
			self.weights[ch][0].append(w_list[-1][0])
			self.weights[ch][1].append(w_list[-1][2])
			
		
		for e in self.int_edges:
			e.set_cn(self.calc_edge_weight(e))
			
	
	# return the average cnv of a segment according to the cnv dict
	def calc_edge_weight(self, e):
		start_node = e.u
		end_node = e.v
		
		if start_node.chr != end_node.chr:
			print 'int edge must be on the same chr', start_node, end_node
			return
		
		ch = start_node.chr
		if ch not in self.weights.keys():
			return Decimal(DEFAULT_CN)
		weights = self.weights[ch]
		
		# left_cnv: latest position where the CN change that is left to the start of the segment
		# left_pos: the position where the segment starts
		# right_cnv: latest position where the CN change that is left to the end of the segment
		# right_pos: where the segment ends
		
		#Deal with start and end edges
		if start_node.term == True:
			left_pos = max(weights[0][1]-TERMINTOR_EDGE_SIZE,0)
		else:
			left_pos = start_node.pos
		left_cnv = bisect.bisect_left(weights[0],left_pos)
		
		if end_node.term == True:
			right_cnv = len(weights[0])
			right_pos = weights[0][-1]+TERMINTOR_EDGE_SIZE
		else:
			right_cnv = bisect.bisect_left(weights[0],end_node.pos)
			right_pos = end_node.pos
		
		if left_cnv==right_cnv:	#No spanning of CNV in this segment - we take the CN of where the entire segment lies which is the previous idx
			return weights[1][left_cnv-1]
			
		if right_pos == left_pos:	#This really only happens in simulations where we can have a bp on pos 0. In that case build_nodes_dict will simply add a "zero edge" at the beginning.
			sum = 0					# It should be noted that the algorithm in general can't handle rearranged karyotypes where the telomers are different than the reference
		else:
			sum = (Decimal(weights[0][left_cnv] - left_pos)/(right_pos-left_pos))*weights[1][left_cnv-1] + (Decimal(right_pos - weights[0][right_cnv-1])/(right_pos-left_pos))*weights[1][right_cnv-1]
			for i, p in enumerate(weights[0][left_cnv:right_cnv-1]):
				sum += (Decimal(weights[0][i+left_cnv+1]-p)/(right_pos-left_pos))*weights[1][i+left_cnv]
		
		if sum<0:
			print 'bad sum: ',sum
			print 'edge: ', e
			print 'weights: ', weights
		return sum
		

	#Build the graph's variant edges from a list of breakpoints
	def build_var_edges(self, bp_list):
		for bp in bp_list:
			u = self.nodes[(bp[0],bp[1])][bp[4]]
			v = self.nodes[(bp[2],bp[3])][bp[5]]
			if not SIMULATION or u!=v:							# With real data u and v are always different but after adjoining close bps they might be the same and than the var edge is a bit redundant,
				self.var_edges.append(VarEdge(u, v, bp[6]))		# but with simulated data we need to have these sort of edges.
				self.total_support += bp[6]						# One thing to consider is edges where u==v that represent an inverse duplicatiom. it might be a good idea to leave that out all together for real data as well
	
	def node_id_generator(self):
		self.num_nodes += 1
		return self.num_nodes-1
	
	def lower_bound_score(self):
		sum=0
		for e in self.int_edges:
			sum += abs(Decimal(round(e.cn))-Decimal(e.cn))
		return sum
	
	def avg_support(self):
		sum=0
		for e in self.var_edges:
			sum+=e.sup
		return Decimal(sum)/len(self.var_edges)
		
	def supp(self, res_file_name):
		sum = Decimal(0)
		count = 0
		
		res_file = open(res_file_name)
		for line in res_file.readlines()[1:]:
			l = line.split('\t')
			if l[0] != 'v':
				continue
			e = self.get_var_edge(int(l[1]),int(l[2]))
			if e!=None:
				sum += Decimal(l[3])*e.sup + Decimal(l[4])*e.sup
				count += Decimal(l[3]) + Decimal(l[4])
			else:
				print 'e not found. ', l[1], l[2]
		
		if count!=0:
			return sum/count
		return 0
			

	# return the edge connecting u,v (if there is more than one - choose the var)
	def get_var_edge(self, u, v):
		edge_list = [x for x in self.var_edges if x.u.id == u and x.v.id == v] + [x for x in self.var_edges if x.u.id == v and x.v.id == u]		# I should get to the bottom od this...
		if edge_list!=[]:
			return edge_list[0]
		return None
			
	
	# export to gv format - this is for the old version, it needs to be redone

	def to_gv(self):
		getcontext().prec=3
		gv_f = open(os.getcwd() + '\\Results\\' + self.samp + '\\' + self.samp + '_' + str(self.chr) + '\\' + self.samp + '_' + str(self.chr) + '_gv.txt', 'w+')
		
		gv_f.write('digraph G {\n')
		gv_f.write('\trankdir=LR\n')
		gv_f.write('\tranksep=0.1\n')
		gv_f.write('\tlabel="' + self.samp +  ', choromosme ' + str(self.chr) + '"\n')
		
		for e in self.int_edges:
			if e.u.term:
				gv_f.write('\t'+str(e.u.id)+' [shape=point, width=0.1, color=blue]\n')
			else:
				gv_f.write('\t'+str(e.u.id)+' [shape=point]\n')
			if e.v.term:
				gv_f.write('\t'+str(e.v.id)+' [shape=point, width=0.1, color=blue]\n')
			else:
				gv_f.write('\t'+str(e.v.id)+' [shape=point]\n')
			
			gv_f.write('\t' + str(e.u.id) + '->' + str(e.v.id) + ' [label=' + str(e.cn) + ', arrowhead=none, style=bold]\n')
			
		for e in self.ref_edges:
			gv_f.write('\t' + str(e.u.id) + '->' + str(e.v.id) + ' [arrowhead=none, style=dotted]\n')
		
		for e in self.var_edges:
			gv_f.write('\t' + str(e.u.id) + '->' + str(e.v.id) + ' [arrowhead=none, color=red, label=' + str(e.sup) +', constraint=false]\n')
		
		gv_f.write('}\n')
		gv_f.close()
	
	#TODO: this works but should be moved outside of this class.
	def build_results_gv(self, results_fname):
	
		if results_fname.split('\\')[-1].split('_')[0] != self.samp.split('_')[0]:
			print 'wrong solution file name, samp: ', self.samp.split('_')[0],  results_fname
			return
		if results_fname.split('\\')[-1].split('_')[1] != self.samp.split('_')[1]:
			print 'wrong solution file name, samp num: ', self.samp.split('_')[1],  results_fname
			return
		if results_fname.split('\\')[-1].split('_')[2] != str(self.chr):
			print 'wrong solution file name, char: ', str(self.chr), results_fname
			return
		
		getcontext().prec=1
		
		ilp_f = open(results_fname)
		gv_f = open(results_fname + '_gv.txt', 'w+')
			
		#use the solution file to build a 3 dictionaries for edges and their values)
		int_edge_dic = {}
		ref_edge_dic = {}
		var_edge_dic = {}
		
		lines = ilp_f.readlines()
		for line in lines[1:]:
			l = line.split('\t')
			if l[0] == 'v':
				var_edge_dic[int(l[1]),int(l[2])] = l[3]		#NOTE: I ignore l[5] which holds the original weight as it's stored in e.cn and was not used
				var_edge_dic[int(l[2]),int(l[1])] = l[4]
			elif l[0] == 'i':
				int_edge_dic[int(l[1]),int(l[2])] = l[3]
				int_edge_dic[int(l[2]),int(l[1])] = l[4]
			elif l[0] == 'r':
				ref_edge_dic[int(l[1]),int(l[2])] = l[3]
				ref_edge_dic[int(l[2]),int(l[1])] = l[4]
				continue
			else:
				print 'bad sol file ', results_fname
				return
				
		gv_f.write('digraph G {\n')
		gv_f.write('\trankdir=LR\n')
		gv_f.write('\tranksep=0.1\n')
		gv_f.write('\tlabel="' + self.samp +  ', choromosme ' + str(self.chr) + '"\n')
		
		for e in self.int_edges:
			if (e.u.id,e.v.id) in int_edge_dic.keys():
				rw = round(Decimal(int_edge_dic[e.u.id,e.v.id].strip()))	#right to left weight
				lw = round(Decimal(int_edge_dic[e.v.id,e.u.id].strip()))
			else:
				rw = 0
				lw = 0
				print '(',e.u.id,',',e.v.id,') not in int_edge_dic'
			
			if e.u.term:
				gv_f.write('\t'+str(e.u.id)+' [shape=point, width=0.1, color=blue]\n')
			else:
				gv_f.write('\t'+str(e.u.id)+' [shape=point]\n')
			if e.v.term:
				gv_f.write('\t'+str(e.v.id)+' [shape=point, width=0.1, color=blue]\n')
			else:
				gv_f.write('\t'+str(e.v.id)+' [shape=point]\n')
			
			if rw != 0:
				gv_f.write('\t' + str(e.u.id) + '->' + str(e.v.id) + ' [label="' + str(rw) + ' (' + str(Decimal(e.cn)) + ')" , style=bold]\n')
			else:
				gv_f.write('\t' + str(e.u.id) + '->' + str(e.v.id) + ' [label="0 (' + str(Decimal(e.cn)) + ')" , arrowhead=none, style=invis]\n')
			if lw != 0:
				gv_f.write('\t' + str(e.v.id) + '->' + str(e.u.id) + ' [label="' + str(lw) + ' (' + str(Decimal(e.cn)) + ')" , style=bold, constraint=false]\n')
				
		for e in self.ref_edges:
			if (e.u.id,e.v.id) in ref_edge_dic.keys():
				rw = round(Decimal(ref_edge_dic[e.u.id,e.v.id].strip()))	#right to left weight
				lw = round(Decimal(ref_edge_dic[e.v.id,e.u.id].strip()))
			else:
				rw = 0
				lw = 0
				print '(',e.u.id,',',e.v.id,') not in ref_edge_dic'
				
			if rw!=0:
				gv_f.write('\t' + str(e.u.id) + '->' + str(e.v.id) + ' [label="' + str(rw) + '", style=dotted]\n')
			else:
				gv_f.write('\t' + str(e.u.id) + '->' + str(e.v.id) + ' [style=invis]\n')		#invisble edges are used to make the graph looks nicer
			if lw!=0:
				gv_f.write('\t' + str(e.v.id) + '->' + str(e.u.id) + ' [label="' + str(lw) + '", style=dotted, constraint=false]\n')
		
		
		for e in self.var_edges:
			if (e.u.id,e.v.id) in var_edge_dic.keys():
				rw = round(Decimal(var_edge_dic[e.u.id,e.v.id].strip()))	#right to left weight
				lw = round(Decimal(var_edge_dic[e.v.id,e.u.id].strip()))
			else:
				lw = 0
				rw = 0
				print '(',e.u.id,',',e.v.id,') not in var_edge_dic'
				
			if rw!=0:
				gv_f.write('\t' + str(e.u.id) + '->' + str(e.v.id) + ' [label="' + str(rw) + '", color=red, constraint=false]\n')
			if lw!=0:
				gv_f.write('\t' + str(e.v.id) + '->' + str(e.u.id) + ' [label="' + str(lw) + '", color=red, constraint=false]\n')
			if rw==0 and lw==0:
				gv_f.write('\t' + str(e.u.id) + '->' + str(e.v.id) + ' [arrowhead=none, color=red, style=dotted, constraint=false]\n')
		
		gv_f.write('}\n')
		ilp_f.close()
		gv_f.close()
	
		
class RefEdge:
	def __init__(self, u, v):
		self.u = u
		self.v = v
	
	def __repr__(self):
		return 'Ref edge: %s ->\n%s\n' %(self.u,self.v)

class VarEdge:
	def __init__(self, u, v, sup):
		self.u = u
		self.v = v
		self.sup = sup
	
	def __repr__(self):
		return 'Var edge: %s ->\n%s\nsupp:%s' %(self.u,self.v,self.sup)
		

# represent an interval edge in the BPGraph with a weight
class IntEdge:
	def __init__(self, u, v, cn=0):
		self.u = u
		self.v = v
		self.cn = cn
	def set_cn(self, cn):
		self.cn = cn
		
	def get_len(self):
		if self.u.term:
			return TERMINTOR_EDGE_SIZE
		return abs(self.u.pos - self.v.pos)
	
	def __repr__(self):
		return 'int edge: %s ->\n%s\ncn:%s' %(self.u,self.v,self.cn)
	

# Represent a node in the BPGrpah. TODO: I should use a factory to make sure id's are unique
class Node:
	def __init__(self, id, chr, pos, side, terminator=False):
		self.id = id
		self.chr = chr
		self.side = side	# 0 for "+" (left side) or 1 for "-" (right side)
		if terminator:
			self.pos = 0
		else:	
			self.pos = pos
		self.term = terminator
	
	def __repr__(self):
		return 'node: %d. chr:%s, pos:%d, side:%s terminator:%s' %(self.id, self.chr, self.pos, side_str(self.side), self.term)
	
def side_str(side):
	if side == 0:
		return "+"
	elif side == 1:
		return "-"
	else:
		return str(side)
