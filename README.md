# KappaUtilities
A simple collection of utilities for the purpose of handling and visualizing Kappa complexes, processing snapshots, determining whether two complexes are isomorphic, and identifying embeddings between a pattern and a host graph.

The typical usage is offline analysis of KaSim snapshots (https://github.com/Kappa-Dev/KappaTools).

For the official tools visit [KappaTools] (https://github.com/Kappa-Dev/KappaTools).

See also: https://github.com/hmedina/KaSaAn.

**kappathings.py**

The file defines the class *KappaComplex*, which takes as input the specification of a Kappa complex as a string in JSON format or as a Kappa expression.
It constructs an internal representation of the form

               self.agents[name] = { site_name: {'state': state, 'bond': bond_stub} }
            self.adjacency[name] = [ agent_names ]
                 self.info[name] = { 'id': identifier
                                     'type': agent_type
                                     'nbonds': {name: int e}  # e is the number of bonds to node 'name'
                                     'degree': int n }
                    self.nxbonds = [ (agent1, agent2, {'sites': site1-site2}) ]

 where everything, except e, is of type string, and 
  * the interface dictionary of an agent is sorted by site name (needs Python 3.7+)            
  * agent names are unique, consisting of type + identifier, eg Axin.42. (including the last dot), where the right and left separators are dots by default.
  * self.agents is a dictionary keyed by agent name with value a dictionary defining the interface of the agent
  * self.adjacency is a dictionary keyed by agent name and listing the agents bound to the key agent
  * self.info is a dictionary keyed by agent name providing additional information used by various utilities.
  * self.nxbonds is a list of unique tuples: [ (agent1, agent2, {'sites': 'site1-site2'}) ], lexicographically sorted on agent. It is used to make a networkx graph for the purpose of visualization.
  * bonds are stubs of the form name@site indicating the name of the bound agent and the site anchoring the bond.
            
A dictionary has no order by definition, but we can fake an order by iterating through it using an ordered list of its keys. This list is the purpose of self.name_list in the code.

Global properties:
* self.size: the number of agents in the complex
* self.composition: a dictionary keyed by agent type with value the abundance of that type in the complex.
* self.sum_formula: the sum formula of the complex.
* self.nxgraph the networkx representation of the complex for the purpose of visualization.
            
KappaComplex takes optional additional inputs: 
  * the number of complexes of this type (count=0 by default); usually determined by parsing the snapshot.
  * a flag asking to normalize the agent identifiers (normalize=True): this generates identifiers that form an ordered increasing sequence assigned to agents alphabetically ordered by type.
  * a flag asking to randomize the identifiers (randomize=False).
  
 Typical usage scenarios are provided in the top-level script environment \_\_main\_\_:
```Python
import kappathings as kt
import re

print("from Kappa:")
# input a simple Kappa string to create the KappaComplex representation
data = ' A(o[1], p[2] t{p}[3]), B(x[1] y[2] z[.]), C(w[3], y[z.B])'
c = kt.KappaComplex(data, count=175)
# print it
print(c)
# short version
print(c.show())

print(f'count: {c.count}')
print(f'size: {c.size}')
print(f'sum formula: {c.sum_formula}')
print('')

print("from JSON:")
# input a JSON string of the same complex as above
data = '[{"node_type":"A","node_sites":[{"site_name":"o","site_type":["port",{"port_links":[[[0,1],0]],"port_states":[]}]},{"site_name":"p","site_type":["port",{"port_links":[[[0,1],1]],"port_states":[]}]},{"site_name":"t","site_type":["port",{"port_links":[[[0,2],0]],"port_states":["p"]}]}]},{"node_type":"B","node_sites":[{"site_name":"x","site_type":["port",{"port_links":[[[0,0],0]],"port_states":[]}]},{"site_name":"y","site_type":["port",{"port_links":[[[0,0],1]],"port_states":[]}]},{"site_name":"z","site_type":["port",{"port_links":[],"port_states":[]}]}]},{"node_type":"C","node_sites":[{"site_name":"w","site_type":["port",{"port_links":[[[0,0],2]],"port_states":[]}]}]}]'
c = kt.KappaComplex(data)
print(c)

print("from file:")
# input a file containing one (large) Kappa string
line = open('TestData/bigly.ka', 'r').read()
# remove newlines that might occur in the file
line = re.sub(r'\n+', ' ', line)
# create a KappaComplex with whatever assignment of node identifiers arises
# (that's the normalize=False flag).
c1 = kt.KappaComplex(line, normalize=False)
print(c1)
print("list of nodes:")
# a list of node names of the complex
print(c1.nodes())
# the complex appears as an iterable
print('iterating through the nodes of the complex and showing adjacency views')
# testing __iter__ and __getitem__ methods
for n in c1:
    print(f'adjacency of {n}: {c1[n]}')
print('')

print("normalize identifiers:")
# "normalize" the identifiers: lexicographically sorted node types are assigned successively increasing identifiers.
c2 = kt.KappaComplex(line, normalize=True)
print(c2)

print("randomly permute identifiers:")
# assign identifiers randomly; this is useful for testing the isomorphism implementation
c3 = kt.KappaComplex(line, randomize=True)
print(c3)

print("list of nodes:")
# a list of node names of the complex
print(c3.nodes())
# the complex appears as an iterable
print('iterating through the nodes of the complex, showing adjacency views')
# testing __iter__ and __getitem__ methods
for n in c3:
    print(f'adjacency of {n}: {c3[n]}')
print('')
```

**kappasnap.py**

provides class *Snapshot*, which reads a snapshot file in JSON format or in Kappa (.ka) format and converts it into a list self.complex[], which contains the self.number_of_distinct_complexes as *KappaComplex* objects.  Typical usage scenarios are provided in the top-level script environment \_\_main\_\_.

**kappaviz.py**

provides a rendering of a *Kappacomplex* through render(), using networkx to do all the work. Again \_\_main\_\_ shows some usage scenarios.

**kappamorph.py**

provides two approaches to embedding a Kappa pattern graph G2 into a host graph G1:

* class *SiteGraphMatcher*(G1, G2) exploits the 'rigidity' of site graphs. It is simple and fast.
* class *GraphMatcher*(G1, G2) adapts the VF2 graph isomorphism implementation of networkx to Kappa site graphs. 

Both classes are initialized with two Kappa expressions G1 and G2. The mapping G2 -> G1 is the dictionary GM.mapping. Both classes provide an embed() method. In the case of GraphMatcher the embed() method has an optional keyword test; test='iso' (for isomorphism) and test='embed' (default) for embedding. The function all_embeddings(G1, G2, algo='sitegraph'|'graph') yields a list of all embeddings. If G1=G2 then these are, of course, the symmetries of G1. The algo keyword chooses the algorithm, with 'sitegraph' the default. Usage scenarios are in the top-level \_\_main\_\_. For example:

```Python
import kappathings as kt
import kappamorph as km

G1 = kt.KappaComplex('A(x[1],y[2]), A(x[1],y[3]), C(x[2]), C(x[3]{p})')
G2 = kt.KappaComplex('A(x[1],y[2]), A(x[1],y[3]), C(x[2]), C(x[3])')
GM = km.GraphMatcher(G1, G2)
print(f'G2 is embeddable in G1: {GM.embed(test="embed")}')
print(f'G1 and G2 are isomorphic: {GM.embed(test="iso")}')
GM = km.SiteGraphMatcher(G1, G2)
print(f'G2 is embeddable in G1: {GM.embed()}')
print(f'G1 and G2 are isomorphic: {GM.embed()}')

G1 = kt.KappaComplex('A(b[1] a[2]), A(b[3] a[2]), B(a[1] x{p}), B(a[3] x{u})')
G2 = kt.KappaComplex('A(b[2]), B(a[2] x)')
maps = km.all_embeddings(G1, G2, algo="graph")
print(f'number of embeddings from G2 into G1: {len(maps)}')
km.print_map(maps)
maps = km.all_embeddings(G1, G2, algo="sitegraph")
print(f'number of embeddings from G2 into G1: {len(maps)}')
km.print_map(maps)
```