from decimal import *

# The same as "equivalent" only disregarding affected bridges
def correct_for_missing_bridges(cor_sol, ILP_sol, ILP_input, dbg=None):
	
	cor_g_dic = read_sol_graph_to_dic(cor_sol)
	sol_g_dic = read_sol_graph_to_dic(ILP_sol)
	inp_g_dic = read_sol_graph_to_dic(ILP_input, struct_only = True)
	
	missing_bridges_list = find_missing_bridges(cor_g_dic, inp_g_dic)

	affected_edges = []
	for bridge in missing_bridges_list:
		u, v = bridge
		affected_edges.append(('v',u,v))
		affected_edges += find_affected_edges(cor_g_dic, bridge)

	for u,v in sol_g_dic['i'].keys():
		if ('i',u,v) in affected_edges or ('i',v,u) in affected_edges:
			continue
		if (u,v) not in cor_g_dic['i'].keys():
			return False
		if sol_g_dic['i'][u,v] + sol_g_dic['i'][v,u] != cor_g_dic['i'][u,v] + cor_g_dic['i'][v,u]:
			return False
			
	for u,v in cor_g_dic['v'].keys():
		if ('v',u,v) in affected_edges or ('v',v,u) in affected_edges:
			continue
		cor_uv = cor_g_dic['v'][u,v] + cor_g_dic['v'][v,u]		
		if (u,v) not in sol_g_dic['v'].keys():
			return False	#Note that the opposite is allowed, i.e. a bridge used by solution that is unused in correct (happens when there is a lesswer karyotype)		
		sol_uv = sol_g_dic['v'][u,v] + sol_g_dic['v'][v,u]
		if (sol_uv>0 and cor_uv == 0) or (sol_uv == 0 and cor_uv > 0):
			return False
	return True

# reads the file into a dictionary.
def read_sol_graph_to_dic(res_file, struct_only = False):
	res_sol = open(res_file)
	#use the solution file to build a 3 dictionaries for edges and their values)
	bp_dict = {}
	bp_dict['i'] = {}
	bp_dict['v'] = {}
	bp_dict['r'] = {}
	
	lines = res_sol.readlines()
	for line in lines[1:]:
		l = line.split('\t')
		bp_dict[l[0]][int(l[1]),int(l[2])] = 0
		bp_dict[l[0]][int(l[2]),int(l[1])] = 0
		if not struct_only:
			bp_dict[l[0]][int(l[1]),int(l[2])] = int(round(float(l[3])))
			bp_dict[l[0]][int(l[2]),int(l[1])] = int(round(float(l[4])))	
	res_sol.close()
	return bp_dict			

# Return a list of all bridges that were removed
def find_missing_bridges(cor_dic, inp_dic):
	missing_bridges = []
	for u,v in cor_dic['v'].keys():
		if (u,v) not in inp_dic['v'].keys():
			missing_bridges.append((u,v))
	return missing_bridges
	
# Return a list of edges affected by the removal of bridge
def find_affected_edges(cor_dic, bridge):
	affected_brdgs = []
	# This assumes nodes are numbered in an increasing order which is the case but there is no strict enforcement for now.
	u = min(bridge)
	v = max(bridge)
	
	u_dir = node_dir(u, cor_dic)
	v_dir = node_dir(v, cor_dic)

	#Translocation
	if not same_chr(u,v, cor_dic):
		affected_brdgs.append(ref_edge(u, cor_dic)) 
		affected_brdgs.append(ref_edge(v, cor_dic))
		affected_brdgs.append(('v',ref_node(u,cor_dic),ref_node(v,cor_dic)))
	#	print 'trl'
		return affected_brdgs
	
	#Deletion or duplication
	if u_dir != v_dir:
	#	print 'deldup', u, v
		return path_uv(u, v, cor_dic)
	
	#inversion
	if u_dir=='-':
		return path_uv(u, ref_node(v,cor_dic), cor_dic) + [('v', ref_node(u,cor_dic), ref_node(v,cor_dic))]
	else:
		return path_uv(ref_node(u,cor_dic), v, cor_dic) + [('v', ref_node(u,cor_dic), ref_node(v,cor_dic))]
	
# Return the reference edge connected to node node
def ref_edge(node, graph_dic):
	for (u,v) in graph_dic['r'].keys():
		if node in (u,v):
			return ('r',u,v)

# return the node connected to node via a reference edge
def ref_node(node, graph_dic):
	for (u,v) in graph_dic['r'].keys():
		if node == u:
			return v
		if node == v:
			return u

#return the side on the interval edge the node lies on
def node_dir(node, graph_dic):
	#find u interval edge (we use ref edges because of the 0 node problem (first node in every chrom seem to be the second, i.e first int edge is always 1->0 but the rest is ordered correctly)
	for edge in graph_dic['r'].keys():
		if edge[0] == node:
			if edge[1] > node:
				return '-'
			elif edge[1] < node:
				return '+'

# return True if both nodes are on the same chr
def same_chr(u, v, graph_dic):
	if u==v:
		return True
	if u>v:
		tmp=u;u=v;v=tmp
	
	next = u
	while next!=v:
	#	print next
		next = next_node(next, graph_dic)[0]
		if next<0:
			return False
		if next == v:
			return True

# return the i/r forward path from u to v. Note - this assumes both u,v are on the same chr and neither are the first telomere on the chr.
def path_uv(u, v, graph_dic):
	path = []
	if u==v:
		return []
	if u>v:
		tmp=u;u=v;v=tmp
		
	next = u
	while next!=v:
		next, next_e = next_node(next, graph_dic)
		path.append(next_e)
		if next<0:	#should not happen - they are not on the same chr
			return None
		if next == v:
			return path

# Return the next node on a forward moving chromosomal path and the edge leading to it
def next_node(node, graph_dic):
	telomere = True
	for edge in graph_dic['r'].keys():
		if min(edge) == node:
			return max(edge), ('r', edge[0], edge[1])
		if node in edge:
			telomere = False
	if not telomere:	#If we're here the node is on the + side
		for edge in graph_dic['i'].keys():
			if min(edge) == node:
				return max(edge), ('i', edge[0], edge[1])
	#If we're here that means the node is a telomere which is a special case since both sides are the min of their respective int edge
	for edge in graph_dic['i'].keys():
		if node in edge:
			other_node = min(edge)
			if node_dir(other_node, graph_dic) == '-':
				return other_node, ('i', edge[0], edge[1])
			return -1,None	#telomere
			
	
# return the distance of the solution from the correct file and the solution from the data. the distance is calculated on the int edges alone
def calc_graph_dist(correct_fname, noisy_fname, sol_fname):
	correct_g = open(correct_fname)
	noisy_g = open(noisy_fname)
	sol_g = open(sol_fname)
	
	int_edges = {}
	var_edges = {}
	v_counter = 0		# used for counting the var edges so i can normalize the support of the correct solution (TODO: in the future maybe a unified graph file format)
	
	#to improve readability I use 4 different dictionaries instead of one
	int_edges['sol'] = {}
	int_edges['cor'] = {}
	int_edges['dat'] = {}
	int_edges['len'] = {}
	
	dist_sol_dat = Decimal(0.0)
	supp_data = Decimal(0.0)
	dist_sol_cor = Decimal(0.0)
	supp_correct = Decimal(0.0)
	dist_dat_cor = Decimal(0.0)
	
	# read from sol file: for int - weight, for var - weight and normalized supp (in respect to data file)
	for line in sol_g.readlines()[1:]:
		l = line.split('\t')
		if l[0] == 'i':
			int_edges['sol'][int(l[1]),int(l[2])] = int(round((float(l[3])))) + int(round(float(l[4])))		# for some reason the ILP sometimes write the solution as 1.99999999999, this is probably due to the encoding in the java
			int_edges['sol'][int(l[2]),int(l[1])] = int(round((float(l[3])))) + int(round(float(l[4])))
		elif l[0] == 'v':
			v_counter += 1
			var_edges[int(l[1]),int(l[2])] = int(round((float(l[3])))) + int(round(float(l[4])))
			var_edges[int(l[2]),int(l[1])] = int(round((float(l[3])))) + int(round(float(l[4])))
			
	#read the noisy data and calc its distance from the solution
	for line in noisy_g.readlines()[1:]:
		l = line.split('\t')
		if l[0] == 'i':
			int_edges['dat'][int(l[1]),int(l[2])] = Decimal(l[3])	# noisy data value
			int_edges['dat'][int(l[2]),int(l[1])] = Decimal(l[3])
			int_edges['len'][int(l[1]),int(l[2])] = Decimal(l[4])	# normalized length
			int_edges['len'][int(l[2]),int(l[1])] = Decimal(l[4])	# normalized length

		elif l[0] == 'v':
			supp_data += min(1,var_edges[int(l[1]),int(l[2])]) * Decimal(l[3])		# Count the percentage of total support of edges used relative to the noisy data
	
			
	# Since the correct graph file does not contain a normalized support figures we need to sum it first so to normalize it
	sum_supp_correct = 0
	lines = correct_g.readlines()
	for line in lines[1:]:
		l = line.split('\t')
		if l[0] == 'v':
			sum_supp_correct += int(l[3]) + int(l[4])
			
	for line in lines[1:]:
		l = line.split('\t')
		if l[0] == 'i':
			sum = int(l[3]) + int(l[4])		# movement on the graph is in integers
			int_edges['cor'][int(l[1]),int(l[2])] = sum	# noisy data value
			int_edges['cor'][int(l[2]),int(l[1])] = sum
		elif l[0] == 'v':
			sum = int(l[3]) + int(l[4])
			if (int(l[1]),int(l[2])) in var_edges.keys():		# somtimes the var edge has no support in the simulated data
				supp_correct += min(1,var_edges[int(l[1]),int(l[2])])*((Decimal(l[3]) + Decimal(l[4]))/sum_supp_correct)
	
	for u,v in int_edges['sol']:
		# this edge case happens because my simulations are not limited to the edge of the karyotype, i'll fix it for better results
		if (u,v) not in int_edges['cor']:
			int_edges['cor'][u,v] = 0
		if (u,v) not in int_edges['len']:
			int_edges['len'][u,v] = 0
		if (u,v) not in int_edges['dat']:
			int_edges['dat'][u,v] = 0

		dist_sol_dat += abs(int_edges['sol'][u,v] - int_edges['dat'][u,v]) * int_edges['len'][u,v]
		dist_sol_cor += abs(int_edges['sol'][u,v] - int_edges['cor'][u,v]) * int_edges['len'][u,v]
		dist_dat_cor += abs(int_edges['dat'][u,v] - int_edges['cor'][u,v]) * int_edges['len'][u,v]
	
	# calculate the %cn score 
	cn_score = Decimal(0.0)
	for u,v in int_edges['sol']:
		if int_edges['sol'][u,v] == int_edges['cor'][u,v]:
			cn_score += int_edges['len'][u,v]
	cn_score /= 2	#since each edge is counted twice: u,v and v,u
	
	correct_g.close()
	noisy_g.close()
	sol_g.close()
	
	return dist_sol_dat, supp_data, dist_sol_cor, supp_correct, dist_dat_cor, cn_score

# In some cases a solution can contain a circle, which means it might look 'different',
#	this makes sure the solution is still the same even after reversing the circles
def reverse_circles(res_file):
	sol = open(res_file)
	
	edge_dic = {}
	lines = sol.readlines()
	num_nodes = int(lines[0])
	terminator_nodes = range(num_nodes)			# We determine the terminator nodes by elimination - nodes that dont have ref/var edges connected to them at all
	circles = []
	
	# read and compile an edge_dic: for each node the nodes conmnected to it via int and var/ref edges
	for line in lines[1:]:
		l = line.split('\t')
		e_type = l[0]
		u = int(l[1])
		v = int(l[2])
		uv = int(float(l[3]))
		vu = int(float(l[4]))
		
		
		if e_type == 'i':
			if uv>0:
				if (u,'i') not in edge_dic.keys():
					edge_dic[u,'i'] = []
				edge_dic[u,'i'] += [v]
				circles.append((u,v,'i'))
			if vu>0:
				if (v,'i') not in edge_dic.keys():
					edge_dic[v,'i'] = []
				edge_dic[v,'i'] += [u]
				circles.append((v,u,'i'))
		else:
			if u in terminator_nodes:
				terminator_nodes.remove(u)
			if v in terminator_nodes:
				terminator_nodes.remove(v)
			if uv>0:
				if (u,'rv') not in edge_dic.keys():
					edge_dic[u,'rv'] = []
				edge_dic[u,'rv'] += [v]
				circles.append((u,v,'rv'))
			if vu>0:
				if (v,'rv') not in edge_dic.keys():
					edge_dic[v,'rv'] = []
				edge_dic[v,'rv'] += [u]
				circles.append((v,u,'rv'))
	
	for e in spanning_tree(terminator_nodes, edge_dic):
		circles.remove(e)
	
	rev_graph_name = ''
	if len(circles) > 0:		# Print a new file with the circles reversed
		rev_graph_name = res_file+'_rev.txt'
		rev_graph = open(rev_graph_name, 'w+')
		rev_graph.write(str(num_nodes)+'\n')
		for line in lines[1:]:
			l = line.split('\t')
			e_type = l[0]
			u = int(l[1])
			v = int(l[2])
			if e_type == 'r' or e_type == 'v':
				e_type = 'rv'
			if (u,v,e_type) in circles:
				rev_graph.write(l[0]+'\t'+l[1]+'\t'+l[2]+'\t0.0\t'+str(int(float(l[3]))+int(float(l[4]))))
				if l[0] == 'i' or l[0] == 'v':
					rev_graph.write('\t'+l[5].strip())
				rev_graph.write('\n')
			elif (v,u,e_type) in circles:
				rev_graph.write(l[0]+'\t'+l[1]+'\t'+l[2]+'\t' + str(int(float(l[3]))+int(float(l[4]))) + '\t0.0')
				if l[0] == 'i' or l[0] == 'v':
					rev_graph.write('\t'+l[5].strip())
				rev_graph.write('\n')
			else:
				rev_graph.write(line)
		rev_graph.close()
	sol.close()
	
	return rev_graph_name
	
def spanning_tree(seeds, edge_dic):
	visited = []
	st_edges = []
	BFS_queue = []
	for node in seeds:
		if (node,'i') not in BFS_queue:
			BFS_queue.append((node,'i'))
	
	while BFS_queue != []:
		node, type = BFS_queue[0]
		visited.append((node,type))
		if type == 'i':
			if (node, 'i') in edge_dic.keys():
				for neighbour in edge_dic[node, 'i']:
					if (neighbour, 'rv') not in BFS_queue and (neighbour,'rv') not in visited:
						BFS_queue.append((neighbour,'rv'))
						st_edges.append((node, neighbour, 'i'))
		elif type == 'rv':
			if (node, 'rv') in edge_dic:
				for neighbour in edge_dic[node, 'rv']:
					if (neighbour, 'i') not in BFS_queue and (neighbour,'i') not in visited:
						BFS_queue.append((neighbour,'i'))
						st_edges.append((node, neighbour, 'rv'))
		BFS_queue.remove((node, type))
	
	return st_edges
		
# Return True if sol and cor use the same set of variant edges	
def same_var_edge_usage(sol_file, cor_file):
	cor_dic = read_sol_graph_to_dic(cor_file)
	sol_dic = read_sol_graph_to_dic(sol_file)
	
	for (u,v) in cor_dic['v']:
		if (u,v) not in sol_dic['v']:
			return False	#Solution is not using a bridge that correct is
		cor_var_weight = cor_dic['v'][u,v] + cor_dic['v'][v,u]
		sol_var_weight = sol_dic['v'][u,v] + sol_dic['v'][v,u]
		if cor_var_weight==0 and sol_var_weight != 0:
			return False
		if cor_var_weight!=0 and sol_var_weight == 0:
			return False
	for (u,v) in sol_dic['v']:
		if (u,v) not in cor_dic['v']:
			return False	#Solution is using a bridge that is not in correct - This happens when we introduce a contaminationg sample.
	return True
	
# checks if the graph solutions are identical (reverse order included)
def equal_solutions(res_file, correct_sol_file):
	correct_sol = open(correct_sol_file)		
	bp_dict = read_sol_graph_to_dic(res_file)
	
	res = True
	lines = correct_sol.readlines()
	for line in lines[1:]:
		l = line.split('\t')
		if (int(l[1]),int(l[2])) not in bp_dict[l[0]].keys():	#Happens when bridges are missed.
			return False
		if bp_dict[l[0]][int(l[1]),int(l[2])] != int(round(float(l[3]))) or bp_dict[l[0]][int(l[2]),int(l[1])] != int(round(float(l[4]))):
			res = False
			break
	
#	check reverse order
	if not res:
		res = True
		for line in lines[1:]:
			l = line.split('\t')
			if (int(l[1]),int(l[2])) not in bp_dict[l[0]].keys():	#Happens when bridges are missed.
				return False
			if bp_dict[l[0]][int(l[1]),int(l[2])] != int(round(float(l[4]))) or bp_dict[l[0]][int(l[2]),int(l[1])] != int(round(float(l[3]))):
				res = False
				break
				
	correct_sol.close()
	return res
	
