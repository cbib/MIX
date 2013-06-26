import re
import copy
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
import sys
import networkx as nx
#import pygraph, pydot
#import pygraphviz

from IPython.core.debugger import Tracer; debug_here = Tracer()


import logging
logger = logging.getLogger('mix_logger')

# if "logger" not in globals():
# 	logger = logging.getLogger('mix_graph')
# 	logger.setLevel(logging.DEBUG)

# 	# create console handler and set level to debug
# 	ch = logging.StreamHandler()
# 	ch.setLevel(logging.DEBUG)

# 	# create formatter
# 	formatter = logging.Formatter('%(asctime)s - %(filename)s - %(funcName)s - %(message)s',"%Y-%m-%d %H:%M:%S")
# 	# formatter = logging.Formatter('%(asctime)s - %(message)s')
# 	# add formatter to ch
# 	ch.setFormatter(formatter)

# 	# add ch to logger
# 	logger.addHandler(ch)


class graph ():
	##
	# @brief Graph representing the assembly 

	def __init__(self, alignments, contigs, included_contigs):
		##
		# @brief Graph constructor
		## the oriented graph representing the assembly
		self.GRAPH = nx.DiGraph()
		self.GRAPH.add_node("In", contig="In", bi="None", sens="None", coords="None")
		self.GRAPH.add_node("Out", contig="Out", bi="None", sens="None", coords="None")	
		self.append_alignments(alignments, contigs, included_contigs)

		# debug_here()
		#self.remove_cycles()

### Construction of the graph ####################################################################

	def append_contigs (self, contigs) : 
		for c in contigs.keys() : 
			c_already_in_GRAPH = False
			for n in self.GRAPH.nodes() : 
				if self.GRAPH.node[n]["contig"] == c :
					c_already_in_GRAPH = True
			if not c_already_in_GRAPH :
		 		n = len(self.GRAPH.nodes())-2
				self.GRAPH.add_node(n, contig=c, bi="b", sens="f", coords=1)
				self.GRAPH.add_node(n+1, contig=c, bi="b", sens="f", coords=contigs[c])	
				self.GRAPH.add_edge("In", n, length=0)
				self.GRAPH.add_edge(n, n+1, length=contigs[c])
				self.GRAPH.add_edge(n+1, "Out", length=0)

	def get_single_contigs(self, contigs):
		single_contigs = []
		for c in contigs.keys() : 
			c_already_in_GRAPH = False
			for n in self.GRAPH.nodes() : 
				if self.GRAPH.node[n]["contig"] == c :
					c_already_in_GRAPH = True
			if not c_already_in_GRAPH :
		 		single_contigs.append(c)
		return single_contigs

	def append_alignments(self, alignments, contigs, included_contigs) : 
		for a in alignments :
			#DEBUG SAM
			# Removed check for included contings
			if (a["TAGR"] not in included_contigs) and (a["TAGQ"] not in included_contigs) :
			# if True:
				#### DEBUG SAM  BOUNDARY SET OFF, NUMBER OF NT REQUIRED TO BE CONSIDERED EXTREMAL
				boudary_set_off = 52
				ref_boudary_set_off = boudary_set_off
				query_boudary_set_off = boudary_set_off
				#ref_boudary_set_off = (a["LENR"]/10000)+2
				#query_boudary_set_off = (a["LENQ"]/10000)+2
				ref_forward = a["S1"] < a["E1"]
				query_forward = a["S2"] < a["E2"]
				if ref_forward : 
					ref_aln_5 = a["S1"] < ref_boudary_set_off
					ref_aln_3 = a["E1"] > a["LENR"] - ref_boudary_set_off
				else : 
					ref_aln_5 = a["E1"] < ref_boudary_set_off
					ref_aln_3 = a["S1"] > a["LENR"] - ref_boudary_set_off
				if query_forward : 
					query_aln_5 = a["S2"] < query_boudary_set_off
					query_aln_3 = a["E2"] > a["LENQ"] - query_boudary_set_off
				else : 
					query_aln_5 = a["E2"] < query_boudary_set_off
					query_aln_3 = a["S2"] > a["LENQ"] - query_boudary_set_off
				# if ref_forward and query_forward : 
				# 	extends = (ref_aln_5 and query_aln_3) or (query_aln_5 and ref_aln_3)
				# else : 
				# 	extends = (ref_aln_5 and query_aln_5) or (query_aln_3 and ref_aln_3)
				# assert(extends)
				# if extends :
					# print a
				self.append_alignment(a, ref_forward, query_forward, ref_aln_5, query_aln_5, ref_aln_3, query_aln_3)
		self.link_5_aln_and_3_aln_of_the_same_contig()

	def append_alignment (self, a, ref_forward, query_forward, ref_aln_5, query_aln_5, ref_aln_3, query_aln_3) : 
		n = len(self.GRAPH.nodes())-2
		self.append_nodes_for_a_contig(n, a["TAGR"])#, a["S1"], a["E1"], a["LEN1"], ref_aln_5)
		self.append_nodes_for_a_contig(n+4, a["TAGQ"])#, a["S2"], a["E2"], a["LEN2"], query_aln_5)
		order = []
		if ref_forward and query_forward and ref_aln_3 and query_aln_5 : 
			order.append([n+0, n+1, n+5, n+4])
			order.append([n+6, n+7, n+3, n+2])
		elif ref_forward and (not query_forward) and ref_aln_3 and query_aln_3 : 
			order.append([n+0, n+1, n+7, n+6])
			order.append([n+4, n+5, n+3, n+2])
		elif ref_forward and query_forward and ref_aln_5 and query_aln_3 :
			order.append([n+4, n+5, n+1, n+0])
			order.append([n+2, n+3, n+7, n+6])
		elif ref_forward and (not query_forward) and ref_aln_5 and query_aln_5 :
			order.append([n+6, n+7, n+1, n+0])
			order.append([n+2, n+3, n+5, n+4])
		elif (not ref_forward) and query_forward and ref_aln_5 and query_aln_5 :
			order.append([n+2, n+3, n+5, n+4])
			order.append([n+6, n+7, n+1, n+0])
		elif (not ref_forward) and query_forward and ref_aln_3 and query_aln_3 :
			order.append([n+4, n+5, n+3, n+2])
			order.append([n+0, n+1, n+7, n+6])
		else : 
			sys.exit("Error : alignment configuration not supported.\n"+a)
		for i in range(len(order)) :
			self.GRAPH.add_edge(order[i][1], order[i][3], length=0 ) # A ajouter ?
			self.GRAPH.add_edge(order[i][0], order[i][2], length=0 )
			if order[i][0]%8 < 4 :
				self.GRAPH.add_edge("In", order[i][0], length=a["LENR"]-a["LEN1"])
				self.GRAPH.add_edge("In", order[i][2], length=0) # Inutile ?
				self.GRAPH.add_edge(order[i][0], order[i][1], length=a["LEN1"])
				self.GRAPH.add_edge(order[i][2], order[i][3], length=a["LEN2"])
				self.GRAPH.add_edge(order[i][3], "Out", length=a["LENQ"]-a["LEN2"])
				self.GRAPH.add_edge(order[i][1], "Out", length=0) # Inutile ?
				if i == 0 : 
					self.GRAPH.node[order[i][0]]["coords"], self.GRAPH.node[order[i][1]]["coords"] = a["S1"], a["E1"]
					self.GRAPH.node[order[i][2]]["coords"], self.GRAPH.node[order[i][3]]["coords"] = a["S2"], a["E2"]
				else : 
					self.GRAPH.node[order[i][0]]["coords"], self.GRAPH.node[order[i][1]]["coords"] = a["E1"], a["S1"]
					self.GRAPH.node[order[i][2]]["coords"], self.GRAPH.node[order[i][3]]["coords"] = a["E2"], a["S2"]
			else :
				self.GRAPH.add_edge("In", order[i][0], length=a["LENQ"]-a["LEN2"])
				self.GRAPH.add_edge("In", order[i][2], length=0) # Inutile ?
				self.GRAPH.add_edge(order[i][0], order[i][1], length=a["LEN2"])
				self.GRAPH.add_edge(order[i][2], order[i][3], length=a["LEN1"])
				self.GRAPH.add_edge(order[i][3], "Out", length=a["LENR"]-a["LEN1"])
				self.GRAPH.add_edge(order[i][1], "Out", length=0) # Inutile ?
				if i == 0 : 
					self.GRAPH.node[order[i][0]]["coords"], self.GRAPH.node[order[i][1]]["coords"] = a["S2"], a["E2"]
					self.GRAPH.node[order[i][2]]["coords"], self.GRAPH.node[order[i][3]]["coords"] = a["S1"], a["E1"]
				else : 
					self.GRAPH.node[order[i][0]]["coords"], self.GRAPH.node[order[i][1]]["coords"] = a["E2"], a["S2"]
					self.GRAPH.node[order[i][2]]["coords"], self.GRAPH.node[order[i][3]]["coords"] = a["E1"], a["S1"]

	def append_nodes_for_a_contig(self, n, TAG) :#, S, E, LEN, aln_5) :
		self.GRAPH.add_node(n, contig=TAG, bi="i", sens="f")#, coords=S)
		self.GRAPH.add_node(n+1, contig=TAG, bi="b", sens="f")#, coords=E)
		self.GRAPH.add_node(n+2, contig=TAG, bi="i", sens="r")#, coords=S)
		self.GRAPH.add_node(n+3, contig=TAG, bi="b", sens="r")#, coords=E)

	def link_5_aln_and_3_aln_of_the_same_contig (self) : 
		nodes_by_contig = {}
		for n in self.GRAPH.nodes() : 
			contig = self.GRAPH.node[n]["contig"]
			if contig in nodes_by_contig.keys():
				nodes_by_contig[contig].append(n)
			else : 
				nodes_by_contig[contig] = [n]
		for c in nodes_by_contig.keys() :
			nb_aln = len(nodes_by_contig[c])/4
			print "Nb aln tot pour le contig", c, ":", nb_aln
			if nb_aln > 1 : 
				#compter les forward
				forward = []
				for node in nodes_by_contig[c] :
					if self.GRAPH.node[node]["sens"] == "f" : 
						forward.append(node)
				#si plusieurs aln forward
				print "\tNb aln forward pour le contig", c, ":", len(forward)/2
				if len(forward)/2 > 1 : 
					#alors on lie les 5 aux 3
					IB = []
					BI = []
					for node in forward :
						if self.GRAPH.node[node]["bi"] == "i" : 
							ni = node
							nb = node+1
							if self.GRAPH.has_edge(ni, nb) :
								IB.append(node)
							else : 
								BI.append(node)
					print "\t\tNb ib :", len(IB), "\t\tNb bi :", len(BI)
					for ib in IB : 
						for bi in BI : 
							self.GRAPH.add_edge(bi, ib, length=self.GRAPH.node[ib]["coords"] - self.GRAPH.node[bi]["coords"])
							self.GRAPH.add_edge(ib+2, bi+2, length=0-(self.GRAPH.node[bi+2]["coords"] - self.GRAPH.node[ib+2]["coords"]))

	def link_all_aln_of_the_same_contig (self) : 
		nodes_by_contig = {}
		for n in self.GRAPH.nodes() : 
			contig = self.GRAPH.node[n]["contig"]
			if (contig in nodes_by_contig.keys()) :
				nodes_by_contig[contig].append(n)
			else : 
				nodes_by_contig[contig] = [n]
		for c in nodes_by_contig.keys() :
			nb_aln = len(nodes_by_contig[c])/4
			print "Nb aln tot pour le contig", c, ":", nb_aln
			if nb_aln > 1 : 
				for n1 in nodes_by_contig[c] :
					for n2 in nodes_by_contig[c] :
						if (n1 != n2) and  (self.GRAPH.node[n1]["bi"] == "i") and (self.GRAPH.node[n2]["bi"] == "i") and (self.GRAPH.node[n1]["sens"] == "f") and (self.GRAPH.node[n2]["sens"] == "f") :
							self.GRAPH.add_edge(n1, n2, length=self.GRAPH.node[n2]["coords"] - self.GRAPH.node[n1]["coords"])
							self.GRAPH.add_edge(n2, n1, length=self.GRAPH.node[n1]["coords"] - self.GRAPH.node[n2]["coords"])
							self.GRAPH.add_edge(n1+2, n2+2, length=self.GRAPH.node[n2+2]["coords"] - self.GRAPH.node[n1+2]["coords"])
							self.GRAPH.add_edge(n2+2, n1+2, length=self.GRAPH.node[n1+2]["coords"] - self.GRAPH.node[n2+2]["coords"])

	def remove_cycles(self) : 
		cycles =  nx.simple_cycles(self.GRAPH)
		while len(cycles) > 0 : 
			c = cycles[0]
			self.find_and_break_the_weakest_edge_in_cycle(c)
			cycles =  nx.simple_cycles(self.GRAPH)
	
	def find_and_break_the_weakest_edge_in_cycle(self, c) : 
		weakest_edge = [-1,-1,-1]
		for n in range (len(c)-1) : 
			n1, n2 = n, n+1
			contig1 = self.GRAPH.node[c[n1]]["contig"]
			contig2 = self.GRAPH.node[c[n2]]["contig"]
			not_in = (contig1 != "In") and (contig2 != "In")
			not_out = (contig1 != "Out") and (contig2 != "Out")
			if (contig1 != contig2) and not_in and not_out : 
				aln_len = self.get_aln_len_from_edge_between_two_different_contigs(c[n1], c[n2])
				if (weakest_edge[2] == -1) or (weakest_edge[2] > aln_len and weakest_edge[2] != -1) : 
					weakest_edge = [n1, n2, aln_len]
		self.GRAPH.remove_edge(c[weakest_edge[0]], c[weakest_edge[1]])

	def get_aln_len_from_edge_between_two_different_contigs(self, n1, n2) :
		if self.GRAPH.node[n1]["bi"]=="b" or self.GRAPH.node[n2]["bi"]=="b" : 
			for pre in self.GRAPH.predecessors(n1) :
				contig_pre = self.GRAPH.node[pre]["contig"]
				contig_n1 = self.GRAPH.node[n1]["contig"]
				if contig_pre == contig_n1 : 
					return self.GRAPH[pre][n1]["length"]
		return -1

### Path selection ##############################################################################################

	# def select_extensions (self, contigs) : 
	# 	self.remove_cycles()
	# 	graph_copy = self.GRAPH.copy()
	# 	longest_paths = []
	# 	while len(graph_copy.nodes()) > 2 : 
	# 		print "longest_path", len(graph_copy.nodes())
	# 		longest_path = self.longest_path_in_G (graph_copy)
	# 		longest_paths.append(longest_path)
	# 		self.remove_path_in_a_graph (graph_copy, longest_path)
	# 	print longest_paths
	# 	for p in longest_paths : 
	# 		for n in p :
	# 			self.GRAPH.node[n]["selected"]="True"
	# 	return self.simplify_paths (longest_paths, contigs)

	def select_extensions(self,contigs):
		if len(self.GRAPH)==2: # No alignments! just bail
			logger.info("alignment graph is empty,bailing out")
			return []
		### Test novel longest_path calculation 
		mix_graph,mapping,inv_mapping = prepare_mix_graph(copy.copy(self.GRAPH))
		self.acyclic_graph = remove_all_cycles(mix_graph)
		self.mapping_from_acyclic_to_self = inv_mapping
		self.longest_paths = maximal_independant_longest_path_for_acyclic_graph(self.acyclic_graph,self.mapping_from_acyclic_to_self)
#		debug_here()
		logger.debug("Longest paths:%s",self.longest_paths)
		for p in self.longest_paths : 
			for n in p :
				self.GRAPH.node[n]["selected"]="True"

		# DEBUG: Does she expect all paths to ber inversed ?
		longest_paths_inversed=[]
		for p in self.longest_paths:
			p.reverse()
			longest_paths_inversed.append(p)

		return self.simplify_paths(self.longest_paths, contigs)

	def simplify_paths (self, paths, contigs) : 
		new_paths = []
		for p in paths : 
			new_path = []
			for nb in p : 
				new_path.append(self.GRAPH.node[nb])
			new_paths.append(self.get_simplified_path(new_path, contigs))
		return new_paths

	def get_simplified_path (self, path, contigs) : 
		simplified_path = []
		nb = len(path) - 2
		while nb > 0 : 
			nb1 = nb
			node1 = path[nb1]
			contig = node1["contig"] 
			while path[nb-1]["contig"] == contig : 
				nb -= 1
			nb2 = nb
			node2 = path[nb2]
			forward = (node1["sens"] == "f")
			if forward : 
				if path[nb1+1]["contig"] == "In" : 
					start = 1
				else : 
					start = node1["coords"]
				if path[nb2-1]["contig"] == "Out" : 
					end = contigs[contig]
				else : 
					end = node2["coords"]
				simplified_path.append({"contig":contig,"start":start-1,"end":end, "sens":"f"})
			else :
				if path[nb1+1]["contig"] == "In" : 
					end = contigs[contig]
				else : 
					end = node1["coords"]
				if path[nb2-1]["contig"] == "Out" : 
					start = 1
				else : 
					start = node2["coords"]
				simplified_path.append({"contig":contig,"start":start-1,"end":end, "sens":"r"})
			nb -= 1
		return simplified_path

	def remove_path_in_a_graph (self, G, P) : 
		contigs_to_remove = []
		for np in P : 
			contig = G.node[np]["contig"]
			if contig != "In" and contig != "Out" :
				contigs_to_remove.append(contig)
				for ng in G.nodes() : 										#
					if G.node[ng]["contig"] == contig : 					#
						# pass
						for s in G.neighbors(ng) :							#
							contigs_to_remove.append(G.node[s]["contig"])	#
						for p in G.predecessors(ng) :						#
							contigs_to_remove.append(G.node[p]["contig"])	#
		contigs_to_remove = set(contigs_to_remove)
		contigs_to_remove.discard("In")										#
		contigs_to_remove.discard("Out")									#
		print "Dealing with path",P,[G.node[x]["contig"] for x in P]
		for ng in G.nodes() : 
			if G.node[ng]["contig"] in contigs_to_remove : 
				G.remove_node(ng)
		print "removed contigs "+str(contigs_to_remove)

	def longest_path_in_G (self, G) : 
		adjacency_matrix = self.get_adjacency_matrix_of_the_graph(G)
		#FW_predecessors, FW_weigth = self.Floyd_Warshall(adjacency_matrix)
		FW_predecessors = self.Floyd_Warshall(adjacency_matrix)
		del adjacency_matrix
		longest_path = self.extract_longest_path_from_FW(FW_predecessors)#, FW_weigth)
		del FW_predecessors
		return longest_path

	def extract_longest_path_from_FW(self, FW_predecessors):#, FW_weigth ) :
		n1, n2 = "In", "Out"
		path = []
		while (n2 != "In") : 
			#print "boucle? : ", n1, n2
			n = FW_predecessors[n1][n2]
			path.append(n2)	
			if n == n1 : 
				break
			n2 = n
		path.append("In")
		del FW_predecessors
		return path

	def get_adjacency_matrix_of_the_graph (self, G):
		"""
		Convert a directed graph to an adjacency matrix.
		"""
		total_edges_length = 0
		for e in G.edges() :
			total_edges_length += G.edge[e[0]][e[1]]["length"]
		vertices = G.nodes()
		dist = {}
		for i in vertices:
			dist[i] = {}
			for j in vertices:
				try:
					dist[i][j] = int(0 - (G[i][j]["length"]))
				except KeyError:
					# the distance from a node to itself is 0
					if i == j:
						dist[i][j] = int(0)
					# the distance from a node to an unconnected node is infinity
					else:
						dist[i][j] = int(total_edges_length*100) # +inf
		return dist
	 	 
	def Floyd_Warshall(self, g):
		"""
		Run the Floyd Warshall algorithm on an adjacency matrix.	 
		The Floyd Warshall algorithm computes the minimum cost of a simple path between each
		pair of vertices.
		"""
		vertices = g.keys()
		w = g.copy() # copy g
		p = {}
		for u in vertices :
			p[u] = {}
			for v in vertices : 
				p[u][v] = u
		for k in vertices:
			for i in vertices:
				for j in vertices:
					if int(w[i][j]) > int(w[i][k] + w[k][j]) : 
						w[i][j] = int(w[i][k] + w[k][j])
						p[i][j] = p[k][j]
		del w, vertices
		return p#, w

### drawing graph ####################################################################################

	def write_dot_graph(self, dot_adr):
		for n in self.GRAPH.nodes() : 
			self.GRAPH.node[n]["label"] = self.GRAPH.node[n]["contig"]+"_"+str(n)+"_"+self.GRAPH.node[n]["sens"]+"_"+self.GRAPH.node[n]["bi"]+"_"+str(self.GRAPH.node[n]["coords"])
			if self.GRAPH.node[n]["bi"] == "b" :
				self.GRAPH.node[n]["shape"] = "box"
			if self.GRAPH.node[n]["sens"] == "f" :
				self.GRAPH.node[n]["color"] = "green"
			elif self.GRAPH.node[n]["sens"] == "r" : 	
				self.GRAPH.node[n]["color"] = "red"
			if "selected" in self.GRAPH.node[n].keys() :
				self.GRAPH.node[n]["style"] = "filled"
		for e in self.GRAPH.edges() : 
			self.GRAPH.edge[e[0]][e[1]]["label"] = self.GRAPH.edge[e[0]][e[1]]["length"]
			if ((e[0] == "In" or e[0] == "Out") or (e[1] == "Out" or e[1] == "In"))  and (self.GRAPH.edge[e[0]][e[1]]["length"] == 0): 
				self.GRAPH.edge[e[0]][e[1]]["style"] = "dashed"
				self.GRAPH.edge[e[0]][e[1]]["color"] = "grey"
			elif self.GRAPH.node[e[0]]["contig"] == self.GRAPH.node[e[1]]["contig"] and self.GRAPH.node[e[0]]["sens"] == "f" : 
				self.GRAPH.edge[e[0]][e[1]]["color"] = "green"
			elif self.GRAPH.node[e[0]]["contig"] == self.GRAPH.node[e[1]]["contig"] and self.GRAPH.node[e[0]]["sens"] == "r" : 
				self.GRAPH.edge[e[0]][e[1]]["color"] = "red"
			elif ((e[0] == "In" or e[0] == "Out")  and self.GRAPH.node[e[1]]["sens"] == "f") or ((e[1] == "In" or e[1] == "Out")  and self.GRAPH.node[e[0]]["sens"] == "f") :
				self.GRAPH.edge[e[0]][e[1]]["color"] = "green"
			elif ((e[0] == "In" or e[0] == "Out")  and self.GRAPH.node[e[1]]["sens"] == "r") or ((e[1] == "In" or e[1] == "Out")  and self.GRAPH.node[e[0]]["sens"] == "r") :
				self.GRAPH.edge[e[0]][e[1]]["color"] = "red"
		nx.write_dot(self.GRAPH, dot_adr)

	def write_cytoscape_graph(self, prefix_path) : 
		nodes_file = open(prefix_path+"nodes_attributes.csv", "w")
		nodes_file.write("Number\tSens\tbi\tcoords\tselected")
		for n in self.GRAPH.nodes() : 
			contig = self.GRAPH.node[n]["contig"]
			sens = self.GRAPH.node[n]["sens"]
			bi = self.GRAPH.node[n]["bi"]
			coords = str(self.GRAPH.node[n]["coords"])
			selected = str("selected" in self.GRAPH.node[n].keys())
			nodes_file.write(str(n)+"\t"+sens+"\t"+bi+"\t"+coords+"\t"+selected+"\n")
		nodes_file.close()
		network_file = open(prefix_path+"network.csv", "w")
		for e in self.GRAPH.edges() : 
			n1 = e[0]
			n2 = e[1]
			length = str(self.GRAPH.edge[n1][n2]["length"])
			network_file.write(str(n1)+"\t"+str(n2)+"\t"+length+"\n")
		network_file.close()



### Longest path using SCC decomposition 

def bounded_bfs_paths(G,sources,max_depth=None):
	"""Produce edges in a breadth-first-search starting at source."""
	visited=set(sources)
	depth=0
	stack = []
	for source in sources:
		stack.append((source,iter(G[source]),depth,[source]))
	while stack:
		# print stack
		parent,children,last_depth,path = stack[0]
		depth=last_depth+1
		if (max_depth) and (depth > max_depth):
			break
		try:
			child = next(children)
			# if child not in visited:
#			yield parent,child
#			visited.add(child)
			if child not in path :
				yield path+[child]
				stack.append((child,iter(G[child]),depth,path+[child]))
		except StopIteration:
			stack.pop(0)


def prepare_mix_graph(g_dir):
	for src,tgt,mdata in g_dir.edges(data=True):
		g_dir[src][tgt]['weight']=mdata.get("length")


	# Build mapping to rename nodes
	mapping = {}
	for k in g_dir.node.keys():
		if type(k)==type((1,)): # If tuples are used to label nodes:
			mapping[k]="_".join(map(str,k))
		else:
			mapping[k]=str(k)

	inputs_nodes = [k for k,v in g_dir.in_degree().items() if v==0]
	outputs_nodes = [k for k,v in g_dir.out_degree().items() if v==0]
	assert(len(inputs_nodes) ==1)
	assert(len(outputs_nodes) ==1)

	for an_in in inputs_nodes:
		g_dir.node[an_in]['In']=True
	for an_out in outputs_nodes:
		g_dir.node[an_out]['Out']=True

	mapping[inputs_nodes[0]]="In"
	mapping[outputs_nodes[0]]="Out"

	g_dir=nx.relabel_nodes(g_dir,mapping)	
	return g_dir,mapping, {v:k for k, v in mapping.items()}


def remove_cycle(g):
	# Determine if there are cycles 
	n_scc = [x for x in nx.strongly_connected_component_subgraphs(g) if len(x)>1]
	if len(n_scc)==0:
		return g
	# scc=n_scc[0]
	# Mark nodes in the SCC 
	# Transform each SCC into a set of path 
	g_trans = nx.DiGraph(g)
	scc=n_scc[0] #take the largest 
	# Inputs are all the nodes from the scc having predecessors not in the scc, or nodes marked as inputs and part of the SCC 
	inputs=set()
	# Outputs are all the nodes from the scc having sucessors not in the scc, or nodes marked as outputs and part of the SCC 
	outputs=set()

	for n in scc:
		# Inputs 
		outside=set(g.predecessors(n)).difference(scc.nodes())
		for o in outside:
			inputs.add((o,n))

	for n in scc:
		# Inputs 
		outside=set(g.successors(n)).difference(scc.nodes())
		for o in outside:
			outputs.add((n,o))

	inputs=list(inputs)
	outputs=list(outputs)

	# Generate all paths between nodes in input and nodes in output 
	inputs_scc= list(set([x[1] for x in inputs]))
	outputs_scc= list(set([x[0] for x in outputs]))
	# We add nodes of the SCC that are labeled as I or O 

	inputs_scc=list(set(inputs_scc).union([k for k,v in scc.node.items() if ("In" in v) ]))
	outputs_scc=list(set(outputs_scc).union([k for k,v in scc.node.items() if ("Out" in v)]))
	# We only consider paths ending in outputs nodes that are not input nodes 
	all_paths = [path for path in bounded_bfs_paths(scc,inputs_scc,len(scc)+1) if path[-1] in outputs_scc]
	print "Found",len(all_paths),"paths"

	# Replace the SCC with all paths 
	g_trans.remove_nodes_from(scc)
	# Two path cannot share a node (if not, we can't guarantee we do not re-introduce shorter cycles)
	for idx in range(len(all_paths)):
	# for p in all_paths:
		if (len(all_paths)>5000) and ((idx%100)==0):
			print "\tProcessed",idx,"paths"
		p=all_paths[idx]

		p_renamed= [str(x)+"@"+str(idx) for x in p]
		g_trans.add_path(p_renamed)
		# Copy nodes and edges attributes from the original graph 
		for i in range(0,len(p)-1):
			src,tgt=p[i],p[i+1]
			src_r,tgt_r=p_renamed[i],p_renamed[i+1]
			g_trans[src_r][tgt_r].update(g[src][tgt])
			g_trans.node[src_r].update(g.node[src])
			g_trans.node[tgt_r].update(g.node[tgt])
		# Reconnect the path to the rest of the graph 

		for i in inputs:
			src,tgt=i
			if tgt in p:
				tgt_in_p=tgt+"@"+str(idx)
				g_trans.add_edge(src,tgt_in_p)
				g_trans[src][tgt_in_p].update(g[src][tgt])
		for o in outputs:
			src,tgt=o
			if src not in p:
				continue
			src_in_p= src+"@"+str(idx)
			g_trans.add_edge(src_in_p,tgt)
			g_trans[src_in_p][tgt].update(g[src][tgt])
	return g_trans

def remove_all_cycles(g):
	n_scc = [x for x in nx.strongly_connected_component_subgraphs(g) if len(x)>1]
	if len(n_scc)>0:
		print "Found",len(n_scc),"cycles, largest one of size",len(n_scc[0])
	else:
		print "No cycles"
		return g
	# Color them 
	for ccc_idx  in range(len(n_scc)): 
		for node in n_scc[ccc_idx].nodes():
			g.node[node]['in_scc']=ccc_idx

	#sys.exit(0)
	g_p=g
	for i in range(1,1000): #Up to a thousand iteration

		g_n= remove_cycle(g_p)
		if g_n.number_of_edges() == g_p.number_of_edges():
			break
	#	assert("10_4" in g_n)
		g_p=g_n
		n_scc_i = [x for x in nx.strongly_connected_component_subgraphs(g_p) if len(x)>1]
		print "Still",len(n_scc_i),"cycles"

	g_trans=g_p
	# Assert graph is acyclic 
	assert(max([len(x) for x in nx.strongly_connected_components(g_trans)])==1)

	return g_trans



def maximal_independant_longest_path_for_acyclic_graph(acyclic_graph,inv_mapping={}):
	# Compute longest-path-acyclic 
	g_dir_trans_io = nx.DiGraph(acyclic_graph)


	# After the cycle reduction, some nodes might have disappeared : 
	# This can happen in cases where we had a SCC but without inputs. Once this SCC is removed, some nodes might not have predecessors anymore 
	# We thus connect these nodes without pred back to the In node

	for inp in [x for x,v in g_dir_trans_io.node.items() if ("In" in v) and (x!="In")]+[x for x,v in g_dir_trans_io.in_degree().items() if (v==0) and (x!="In")]:
		g_dir_trans_io.add_edge("In",inp,attr_dict={"weight":0})

	for inp in [x for x,v in g_dir_trans_io.node.items() if ("Out" in v) and (x!="Out")]+[x for x,v in g_dir_trans_io.out_degree().items() if (v==0) and (x!="Out")]:
		g_dir_trans_io.add_edge(inp,"Out",attr_dict={"weight":0})

	nodes_in_order = nx.topological_sort(g_dir_trans_io)
	# The sort should start and end with resp. In and Out
	assert(nodes_in_order[0]=="In")
	assert(nodes_in_order[-1]=="Out")

	# Update each node attributes with the accumulated length from the beginning 
	# Q? Should we store each path separately or a greedy alg is ok? trying greedy 
	# Q? How to deal with mulitple equivalent paths ?

	for n in g_dir_trans_io:
		g_dir_trans_io.node[n]['in_longest']="False"
	for src,tgt in g_dir_trans_io.edges():
		g_dir_trans_io[src][tgt]["in_longest"]="False"

	all_longest_path_trans_io = []
	nodes_in_any_longest_path = set()

	inputs_nodes = [k for k,v in g_dir_trans_io.in_degree().items() if v==0]
	inputs_nodes.extend([x for x,v in g_dir_trans_io.node.items() if "In" in v])
	outputs_nodes = [k for k,v in g_dir_trans_io.out_degree().items() if v==0]
	outputs_nodes.extend([x for x,v in g_dir_trans_io.node.items() if "Out" in v])



	# Todo: Account for contig ID in longest path
	for long_path_index in range(0,2000):
		# Clear all longest path info 
		print "longest path iteration",long_path_index
		# print nodes_in_any_longest_path

		for n in g_dir_trans_io.node.keys():
			if "longest_path" in g_dir_trans_io[n]:
				del g_dir_trans_io[n]['longest_path']

		#Init
		g_dir_trans_io.node['In']['longest_path']=[0,[]]
		for n in nodes_in_order[1:]:
			contig = g_dir_trans_io.node[n].get("contig",n)
			if contig in nodes_in_any_longest_path:
				continue

			pred = [x for x in g_dir_trans_io.predecessors(n) if g_dir_trans_io.node[x].get("contig",x) not in nodes_in_any_longest_path]

			if len(pred)==0: 
				# This can happen in cases where we had a SCC but without inputs. Once this SCC is removed, some nodes might not have predecessors anymore 
				longest_path=[0,[]]
			else:
				longest_path= [(g_dir_trans_io.node[pred[0]]['longest_path'][0] + g_dir_trans_io[pred[0]][n]['weight']),\
							pred[0]]
			for p in pred:
				if (g_dir_trans_io.node[p]['longest_path'][0]+ g_dir_trans_io[p][n]['weight']) > longest_path[0]:
					longest_path= [(g_dir_trans_io.node[p]['longest_path'][0]+ g_dir_trans_io[p][n]['weight']),p]
			g_dir_trans_io.node[n]['longest_path']=longest_path


		# Rebuild the longest path, starting from "Out" up to "In" in predecessor order 
		current_node = "Out"
		longest_path_trans_io = []
		while (current_node != "In") and current_node!=[]: 
			longest_path_trans_io.append((current_node,g_dir_trans_io.node[current_node]["longest_path"][0]))
			current_node=g_dir_trans_io.node[current_node]['longest_path'][1]
		if current_node==[]:
			break
		longest_path_trans_io.append((current_node,g_dir_trans_io.node[current_node]["longest_path"][0]))
		longest_path_trans_io.reverse()

		print longest_path_trans_io
		print [g_dir_trans_io.node[x[0]].get("contig",x[0]) for x in longest_path_trans_io]
		if len(longest_path_trans_io)<=2: # only in and out?
			break
		all_longest_path_trans_io.append(longest_path_trans_io)

		# Mark this path in the graph

		g_dir_trans_io.node['Out']['in_longest']="True"


		for i in range(0, len(longest_path_trans_io)-1):
			n = longest_path_trans_io[i][0]
			g_dir_trans_io.node[n]['in_longest']=long_path_index
			g_dir_trans_io[longest_path_trans_io[i][0]][longest_path_trans_io[i+1][0]]['in_longest']=long_path_index
			# debug_here()
			if (n not in ["In","Out"]):
				nodes_in_any_longest_path.add(g_dir_trans_io.node[n]['contig'])
	# Trasnform the longest path into the original graph IDs
	original_longest_paths=[]
	for a_lp in all_longest_path_trans_io:
		lp = [inv_mapping.get(x[0].split("@")[0],x[0].split("@")[0]) for x in a_lp]
		original_longest_paths.append(lp)
	return original_longest_paths
