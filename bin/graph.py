import re
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
import sys
import networkx as nx
#import pygraph, pydot
#import pygraphviz

from IPython.core.debugger import Tracer; debug_here = Tracer()


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
		self.remove_cycles()

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

	def select_extensions (self, contigs) : 
		graph_copy = self.GRAPH.copy()
		longuest_paths = []
		while len(graph_copy.nodes()) > 2 : 
			print "longuest_path", len(graph_copy.nodes())
			longuest_path = self.longuest_path_in_G (graph_copy)
			longuest_paths.append(longuest_path)
			self.remove_path_in_a_graph (graph_copy, longuest_path)
		print longuest_paths
		for p in longuest_paths : 
			for n in p :
				self.GRAPH.node[n]["selected"]="True"
		return self.simplify_paths (longuest_paths, contigs)

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

	def longuest_path_in_G (self, G) : 
		adjacency_matrix = self.get_adjacency_matrix_of_the_graph(G)
		#FW_predecessors, FW_weigth = self.Floyd_Warshall(adjacency_matrix)
		FW_predecessors = self.Floyd_Warshall(adjacency_matrix)
		del adjacency_matrix
		longuest_path = self.extract_longuest_path_from_FW(FW_predecessors)#, FW_weigth)
		del FW_predecessors
		return longuest_path

	def extract_longuest_path_from_FW(self, FW_predecessors):#, FW_weigth ) :
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
