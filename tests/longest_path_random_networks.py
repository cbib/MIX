# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 00:05:57 2013

@author: hayssam
"""
import random
import networkx as nx 


# Make an accessible seed
import time
a = long(time.time() * 256) # use fractional seconds

# Bug for Seed is 351188051725 and grid 12,6
# a=351188051725

# Bug for see 351188626429, missing one node for the longest path 
# a=351188626429

random.seed(a)
print "Seed is",a
from IPython.core.debugger import Tracer; debug_here = Tracer() 


def find_all_paths(graph, start, end, path=[]):
	path = path + [start]
	if start == end:
		return [path]
	if not graph.has_key(start):
		return []
	paths = []
	for node in graph[start]:
		if node not in path:
			newpaths = find_all_paths(graph, node, end, path)
			for newpath in newpaths:
				paths.append(newpath)
	return paths


def path_to_str(path):
	if len(path)==0:
		return ""
	return str(path[0][0]) + "_" + "_".join([str(x[1]) for x in path])

def find_all_paths(graph, start,end):
	# Think about the case where start is also an end
	all_paths=[]
	current_path=[]
	def append_path(p):
		tot_weight = 0
		for src,tgt in current_path:
			tot_weight+=graph[src][tgt].get("weight",1)
		all_paths.append((path_to_str(current_path),tot_weight ))

	for e in nx.dfs_edges(nx.bfs_tree(graph,start)):
		if e[0]== start: #We start a new path
			if len(current_path):
				# Do we end in an out node ? 

				append_path(current_path)
				current_path=[]
		if e[1] in end: # We found a path
			append_path(current_path) 

		current_path.append(e)			
	if len(current_path):
		append_path(current_path)
	return all_paths


def bounded_bfs_paths(G,sources,max_depth=None):
	"""Produce edges in a breadth-first-search starting at source."""
	# Based on http://www.ics.uci.edu/~eppstein/PADS/BFS.py
	# by D. Eppstein, July 2004.
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
	
g= nx.grid_graph(dim=[15,15])
# Build mapping to rename nodes
mapping = {}
for k in g.node.keys():
	mapping[k]="_".join(map(str,k))
g=nx.relabel_nodes(g,mapping)

# Make a directed graph by selecting a random direction for each edge
g_dir=nx.DiGraph()
for e in g.edges():
	src,tgt = e
	if random.random()>=0.5:
		g_dir.add_edge(src,tgt)
	else:
		g_dir.add_edge(tgt,src)

# Randomly place an input node 
no_input = [k for k,v in g_dir.in_degree().items() if v==0] 
if len(no_input)>0:
	for n in no_input:
		g_dir.node[n]['In']=True
else:
	g_dir.node[random.choice(g_dir.node.keys())]['In']=True

no_output =[k for k,v in g_dir.out_degree().items() if v==0] 

if len(no_output)>0:
	for n in no_output:
		g_dir.node[n]['Out']=True

else:
	g_dir.node[random.choice(g_dir.node.keys())]['Out']=True


# Randomly place an output node 
g_dir.node[random.choice(g_dir.node.keys())]['Out']=True

# # Nested loops example 
# g_dir=nx.DiGraph()
# g_dir.add_path([1,2,3])
# g_dir.add_path([1,2,4,5,6,4])
# g_dir.add_edge(1,3)
# g_dir.add_edge(2,3)
# #g_dir.add_edge(4,3)
# g_dir.add_edge(5,3)
# g_dir.add_edge(6,4)
# g_dir.add_path([4,7,5])
# g_dir.add_path([4,7,8,5])

# g_dir.node[1]['In']=True
# g_dir.node[3]['Out']=True


# Add constant weights 
for src,tgt in g_dir.edges():
	g_dir[src][tgt]['weight']=1




def remove_cycle(g):
	# Determine if there are cycles (none for seed 1235)
	n_scc = [x for x in nx.strongly_connected_component_subgraphs(g) if len(x)>1]
	if len(n_scc)==0:
		return g
	# scc=n_scc[0]
	# Mark nodes in the SCC 
	# Transform each SCC into a set of path 
	g_trans = nx.DiGraph(g)
	scc=n_scc[0] #take the largest 
	# Inputs are all the nodes from the scc having predecessors not in the scc
	inputs=set()
	# Outputs are all the nodes from the scc having sucessors not in the scc
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
	# We only consider paths ending in outputs nodes that are not input nodes 
	all_paths = [path for path in bounded_bfs_paths(scc,inputs_scc,len(scc)+1) if path[-1] in outputs_scc]
	print "Found",len(all_paths),"paths"

	# Replace the SCC with all paths 
	g_trans.remove_nodes_from(scc)
	# Two path cannot share a node (if not, we can't guarantee we do not re-introduce shorter cycles)
	for idx in range(len(all_paths)):
	# for p in all_paths:
		p=all_paths[idx]

		p_renamed= [x+"@"+str(idx) for x in p]
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



n_scc = [x for x in nx.strongly_connected_component_subgraphs(g_dir) if len(x)>1]
print "Found",len(n_scc),"cycles"
# Color them 
for ccc_idx  in range(len(n_scc)): 
	for node in n_scc[ccc_idx].nodes():
		g_dir.node[node]['in_scc']=ccc_idx

g_p=g_dir
for i in range(1,1000): #Up to a thousand iteration

	g_n= remove_cycle(g_p)
	if g_n.number_of_edges() == g_p.number_of_edges():
		break
#	assert("10_4" in g_n)
	g_p=g_n
	n_scc_i = [x for x in nx.strongly_connected_component_subgraphs(g_p) if len(x)>1]
	print "Still",len(n_scc_i),"cycles"

g_dir_trans=g_p

# Assert we have all nodes and the graph is acyclic 
# assert(set(g_dir_trans.nodes()) == set(g_dir_trans.nodes()))
# Not true: Suppose that there's no single paths from input to output going through a node, then that node is removed 
# Cf node 6 in evernote:///view/141184/s2/8117ccaa-8d29-4b97-94f1-918acd0b14b9/8117ccaa-8d29-4b97-94f1-918acd0b14b9/

# Assert graph is acyclic 
assert(max([len(x) for x in nx.strongly_connected_components(g_dir_trans)])==1)

# Compute longest-path-acyclic 
g_dir_trans_io = nx.DiGraph(g_dir_trans)


# After the cycle reduction, some nodes might have disappeared : 
# This can happen in cases where we had a SCC but without inputs. Once this SCC is removed, some nodes might not have predecessors anymore 
# We thus connect these nodes without pred back to the In node

for inp in [x for x,v in g_dir_trans_io.node.items() if "In" in v]+[x for x,v in g_dir_trans_io.in_degree().items() if v==0]:
	g_dir_trans_io.add_edge("In",inp,attr_dict={"weight":0})

for inp in [x for x,v in g_dir_trans_io.node.items() if "Out" in v]+[x for x,v in g_dir_trans_io.out_degree().items() if v==0]:
	g_dir_trans_io.add_edge(inp,"Out",attr_dict={"weight":0})

nodes_in_order = nx.topological_sort(g_dir_trans_io)
# The sort should start and end with resp. In and Out
assert(nodes_in_order[0]=="In")
assert(nodes_in_order[-1]=="Out")

# Update each node attributes with the accumulated length from the beginning 
# Q? Should we store each path separately or a greedy alg is ok? trying greedy 
# Q? How to deal with mulitple equivalent paths ?

g_dir_trans_io.node['In']['longest_path']=[0,[]]
for n in nodes_in_order[1:]:
	pred = g_dir_trans_io.predecessors(n)
	if len(pred)==0: 
		# This can happen in cases where we had a SCC but without inputs. Once this SCC is removed, some nodes might not have predecessors anymore 
		longest_path=[0,[]]
	else:
		longest_path= [(g_dir_trans_io.node[pred[0]]['longest_path'][0] + g_dir_trans_io[pred[0]][n]['weight']),\
					pred[0]]
	for p in pred:
		if g_dir_trans_io.node[p]['longest_path'][0] > longest_path[0]:
			longest_path= [(g_dir_trans_io.node[p]['longest_path'][0]+ g_dir_trans_io[p][n]['weight']),p]
	g_dir_trans_io.node[n]['longest_path']=longest_path


# Rebuild the longest path, starting from "Out" up to "In" in predecessor order 
current_node = "Out"
longest_path = []
while current_node != "In":
	longest_path.append((current_node,g_dir_trans_io.node[current_node]["longest_path"][0]))
	current_node=g_dir_trans_io.node[current_node]['longest_path'][1]
longest_path.append((current_node,g_dir_trans_io.node[current_node]["longest_path"][0]))
longest_path.reverse()

# Mark this path in the graph
for n in g_dir_trans_io:
	g_dir_trans_io.node[n]['in_longest']="False"
for src,tgt in g_dir_trans_io.edges():
	g_dir_trans_io[src][tgt]["in_longest"]="False"

g_dir_trans_io.node['Out']['in_longest']="True"
for i in range(0, len(longest_path)-1):
	n = longest_path[i][0]
	g_dir_trans_io.node[n]['in_longest']="True"
	g_dir_trans_io[longest_path[i][0]][longest_path[i+1][0]]['in_longest']="True"






# Prepare for GML output 

for i in range(1, len(longest_path)-2):
	n = longest_path[i][0].split("@")[0]
	p = longest_path[i+1][0].split("@")[0]
	g_dir.node[n]['in_longest']="True"
	g_dir[n][p]['in_longest']="True"
# Add the last item
last_element = longest_path[-2][0].split("@")[0]
g_dir.node[last_element]['in_longest']="True"




for k,v in g_dir_trans_io.node.items():
	#g_dir_trans_io.node[k].update({"longest_path":"_".join(map(str,v['longest_path']))})
	del g_dir_trans_io.node[k]['longest_path']

	if "In" in g_dir_trans_io.node[k]:
		del g_dir_trans_io.node[k]['In']
	if "Out" in g_dir_trans_io.node[k]:
		del g_dir_trans_io.node[k]['Out']



# Add some weights to In and Out edges
for tgt in g_dir_trans_io['In']: 
	g_dir_trans_io['In'][tgt]['weight']=1

for src in g_dir_trans_io.predecessors('Out'):
	g_dir_trans_io[src]["Out"]['weight']=1

mapping = {}
for k in g_dir.node.keys():
	mapping[k]="N %s"%(k)
g_dir_trans_io=nx.relabel_nodes(g_dir_trans_io,mapping)
g_dir=nx.relabel_nodes(g_dir,mapping)


nx.write_gml(g_dir,"grid.gml")
nx.write_gml(g_dir_trans_io,"grid_trans_long_path.gml")
