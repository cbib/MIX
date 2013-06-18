MIX
===

**MIX: Combining multiple assemblies from NGS data**

Mix is a tool to combine multiple assemblies for NGS data. Its algorithm takes two assemblies and generates another one that mixes them in order to extend the length of resulting contigs. It builds an assembly graph in which all of the contigs are vertices and edges represent the best possible alignments between two contigs that have the potential of being used as basis for contig extension.
The resulting output assembly corresponds to a certain path in this assembly graph.