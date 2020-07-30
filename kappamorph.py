"""
This is an adaptation of the undirected graph isomorphism code (VF2 algorithm) from networkx
(https://github.com/networkx/networkx/tree/master/networkx/algorithms/isomorphism)
to site graph patterns in the Kappa language.

VF2 references
--------------
[1]   Luigi P. Cordella, Pasquale Foggia, Carlo Sansone, Mario Vento,
      "A (Sub)Graph Isomorphism Algorithm for Matching Large Graphs",
      IEEE Transactions on Pattern Analysis and Machine Intelligence,
      vol. 26,  no. 10,  pp. 1367-1372,  Oct.,  2004.
      http://ieeexplore.ieee.org/iel5/34/29305/01323804.pdf
[2]   L. P. Cordella, P. Foggia, C. Sansone, M. Vento, "An Improved
      Algorithm for Matching Large Graphs", 3rd IAPR-TC15 Workshop
      on Graph-based Representations in Pattern Recognition, Cuen,
      pp. 149-159, 2001.
      http://amalfi.dis.unina.it/graph/db/papers/vf-algorithm.pdf
"""

# The original networkx code was based on work by Christopher Ellison

import sys


# The graphs are given by the internal representation in 'kappathings.py',
# which, for the purpose of this code, behaves like a networkx graph.

def print_map(maps):
    for i in range(0, len(maps)):
        print(f'map {i + 1}:')
        for k, v in maps[i].items():
            print(f'{k} --> {v}')


def all_embeddings(host, pattern):
    rarest_pattern_type = next(iter(pattern.composition))
    abundance = host.composition[rarest_pattern_type]

    mappings = []
    m = 0
    for i in range(0, abundance):
        GM = GraphMatcher(host, pattern)
        if GM.is_embeddable():
            # eliminate identical embeddings
            found = False
            for j in range(0, m):
                if GM.mapping == mappings[j]:
                    found = True
                    break
            if not found:
                mappings += [GM.mapping]
                m += 1
        else:
            break
        if i < abundance - 1:
            host.name_list = kt.shift(host.name_list)
    return mappings


class GraphMatcher:
    """
    Implementation of VF2 algorithm for matching undirected graphs.
    """

    def __init__(self, G1, G2):
        """
        Initialize GraphMatcher.
        Parameters
        ----------
        G1,G2: KappaComplex instances.
        The two graphs to check for isomorphism or for an embedding of G2 (pattern) into G1 (host).
        """
        self.G1 = G1
        self.G2 = G2
        self.G1_nodes = set(G1.nodes())
        self.G2_nodes = set(G2.nodes())
        self.G2_node_order = {n: i for i, n in enumerate(G2.name_list)}

        # Set recursion limit.
        self.old_recursion_limit = sys.getrecursionlimit()
        expected_max_recursion_level = len(self.G2)
        if self.old_recursion_limit < 1.5 * expected_max_recursion_level:
            # Give some breathing room.
            sys.setrecursionlimit(int(1.5 * expected_max_recursion_level))

        # default: Declare that we will be searching for a graph-graph isomorphism.
        self.test = 'iso'

        # Initialize state
        self.initialize()

    def reset_recursion_limit(self):
        """Restores the recursion limit."""
        # TODO:
        # Currently, we use recursion and set the recursion level higher.
        # It would be nice to restore the level, but because the
        # (Di)GraphMatcher classes make use of cyclic references, garbage
        # collection will never happen when we define __del__() to
        # restore the recursion level. The result is a memory leak.
        # So for now, we do not automatically restore the recursion level,
        # and instead provide a method to do this manually. Eventually,
        # we should turn this into a non-recursive implementation.
        sys.setrecursionlimit(self.old_recursion_limit)

    def initialize(self):
        """Reinitializes the state of the algorithm.
        This method should be redefined if using something other than GMState.
        If only subclassing GraphMatcher, a redefinition is not necessary.
        """

        # core_1[n] contains the index of the node paired with n, which is m,
        #           provided n is in the mapping.
        # core_2[m] contains the index of the node paired with m, which is n,
        #           provided m is in the mapping.
        self.core_1 = {}
        self.core_2 = {}

        # See the paper for definitions of M_x and T_x^{y}

        # inout_1[n]  is non-zero if n is in M_1 or in T_1^{inout}
        # inout_2[m]  is non-zero if m is in M_2 or in T_2^{inout}
        #
        # The value stored is the depth of the SSR tree when the node became
        # part of the corresponding set.
        self.inout_1 = {}
        self.inout_2 = {}
        # Practically, these sets simply store the nodes in the subgraph.

        self.state = GMState(self)

        # Provide a convenient way to access the isomorphism mapping.
        self.mapping = self.core_1.copy()

    def is_isomorphic(self):
        """
        Returns True if G1 and G2 are isomorphic graphs.
        The isomorphism is in GM.mapping
        """

        # Let's do two very quick checks!
        # QUESTION: Should we call faster_graph_could_be_isomorphic(G1,G2)?
        # For now, I just copy the code.

        # Check global properties
        if self.G1.order() != self.G2.order():
            return False
        if self.G1.sum_formula != self.G2.sum_formula:
            return False

        # Check local properties
        d1 = sorted(d for n, d in self.G1.degree())
        d2 = sorted(d for n, d in self.G2.degree())
        if d1 != d2:
            return False

        try:
            x = next(self.isomorphisms_iter())
            return True
        except StopIteration:
            return False

    def is_embeddable(self):
        """
        Returns True if G2 can be embedded in G1.
        The embedding is in GM.mapping
        """

        # Check size
        if self.G1.order() < self.G2.order():
            return False
        # Check composition
        for node_type in G2.composition:
            if G1.composition[node_type] < G2.composition[node_type]:
                return False

        try:
            x = next(self.embedding_iter())
            # invert the mapping for more natural reading
            m = {}
            for k, v in self.mapping.items():
                m[v] = k
            self.mapping = m
            return True
        except StopIteration:
            return False

    def isomorphisms_iter(self):
        """Generator over isomorphisms between G1 and G2."""
        # Declare that we are looking for a graph-graph isomorphism.
        self.test = 'iso'
        self.initialize()
        yield from self.match()

    def embedding_iter(self):
        """Generator over monomorphisms between a subgraph of G1 and G2."""
        # Declare that we are looking for graph-subgraph monomorphism.
        self.test = 'embed'
        self.initialize()
        yield from self.match()

    def match(self):
        """Extends the isomorphism mapping.
        This function is called recursively to determine if a complete
        isomorphism can be found between G1 and G2.  It cleans up the class
        variables after each recursive call. If an isomorphism is found,
        we yield the mapping.
        """
        if len(self.core_1) == len(self.G2):
            # Save the final mapping, otherwise garbage collection deletes it.
            self.mapping = self.core_1.copy()
            # The mapping is complete.
            yield self.mapping
        else:
            for G1_node, G2_node in self.candidate_pairs_iter():
                if self.syntactic_feasibility(G1_node, G2_node):
                    if self.semantic_feasibility(G1_node, G2_node):
                        # Recursive call, adding the feasible state.
                        newstate = self.state.__class__(self, G1_node, G2_node)
                        yield from self.match()

                        # restore data structures
                        newstate.restore()

    def candidate_pairs_iter(self):
        """Iterator over candidate pairs of nodes in G1 and G2."""

        # All computations are done using the current state!

        # G1_nodes = self.G1_nodes  # WF: not used
        G2_nodes = self.G2_nodes
        min_key = self.G2_node_order.__getitem__

        # First we compute the inout-terminal sets.
        T1_inout = [node for node in self.inout_1 if node not in self.core_1]
        T2_inout = [node for node in self.inout_2 if node not in self.core_2]

        # If T1_inout and T2_inout are both nonempty.
        # P(s) = T1_inout x {min T2_inout}
        if T1_inout and T2_inout:
            node_2 = min(T2_inout, key=min_key)
            for node_1 in T1_inout:
                yield node_1, node_2
        else:
            # If T1_inout and T2_inout were both empty....
            # P(s) = (N_1 - M_1) x {min (N_2 - M_2)}
            # if not (T1_inout or T2_inout):  # as suggested by  [2], incorrect
            if 1:  # as inferred from [1], correct
                # First we determine the candidate node for G2
                other_node = min(G2_nodes - set(self.core_2), key=min_key)
                for node in self.G1:
                    if node not in self.core_1:
                        yield node, other_node

        # For all other cases, we don't have any candidate pairs.

    def syntactic_feasibility(self, G1_node, G2_node):
        """Returns True if adding (G1_node, G2_node) is syntactically feasible.
        This function returns True if it is adding the candidate pair
        to the current partial isomorphism/embedding does not make it impossible
        for a final isomorphism/embedding to be found.
        """

        # edges here are explicit edges, not accounting for '_' or stubs or '#'

        ###
        # Test at each step to get a return value as soon as possible.
        ###

        # Look ahead 0

        # R_self

        # The number of self loops for G1_node must equal the number of
        # self-loops for G2_node. Without this check, we would fail on
        # R_neighbor at the next recursion level. But it is good to prune the
        # search tree now.

        if self.test == 'embed':
            if self.G1.number_of_edges(G1_node, G1_node) < self.G2.number_of_edges(G2_node, G2_node):
                return False
        else:
            if self.G1.number_of_edges(G1_node, G1_node) != self.G2.number_of_edges(G2_node, G2_node):
                return False

        # R_neighbor

        # For each neighbor n' of n in the partial mapping, the corresponding
        # node m' is a neighbor of m, and vice versa. Also, the number of
        # edges must be equal.
        if self.test == 'iso':
            for neighbor in self.G1[G1_node]:
                if neighbor in self.core_1:
                    if not (self.core_1[neighbor] in self.G2[G2_node]):
                        return False
                    elif self.G1.number_of_edges(neighbor, G1_node) != self.G2.number_of_edges(self.core_1[neighbor],
                                                                                               G2_node):
                        return False

        for neighbor in self.G2[G2_node]:
            if neighbor in self.core_2:
                if not (self.core_2[neighbor] in self.G1[G1_node]):
                    return False
                elif self.test == 'embed':
                    if self.G1.number_of_edges(self.core_2[neighbor], G1_node) < self.G2.number_of_edges(neighbor,
                                                                                                         G2_node):
                        return False
                else:
                    if self.G1.number_of_edges(self.core_2[neighbor], G1_node) != self.G2.number_of_edges(neighbor,
                                                                                                          G2_node):
                        return False

        if self.test == 'iso':
            # Look ahead 2

            # R_new
            # The number of neighbors of n that are neither in the core_1 nor
            # T_1^{inout} is equal to the number of neighbors of m
            # that are neither in core_2 nor T_2^{inout}.
            num1 = 0
            for neighbor in self.G1[G1_node]:
                if (neighbor in self.inout_1) and (neighbor not in self.core_1):
                    num1 += 1
            num2 = 0
            for neighbor in self.G2[G2_node]:
                if (neighbor in self.inout_2) and (neighbor not in self.core_2):
                    num2 += 1
            if not (num1 == num2):
                return False

            # Look ahead 1

            # R_terminout
            # The number of neighbors of n in T_1^{inout} is equal to the
            # number of neighbors of m that are in T_2^{inout}, and vice versa.
            num1 = 0
            for neighbor in self.G1[G1_node]:
                if neighbor not in self.inout_1:
                    num1 += 1
            num2 = 0
            for neighbor in self.G2[G2_node]:
                if neighbor not in self.inout_2:
                    num2 += 1
            if not (num1 == num2):
                return False

        # Otherwise, this node pair is syntactically feasible!
        return True

    def semantic_feasibility(self, G1_node, G2_node):
        """Returns True if adding (G1_node, G2_node) is semantically feasible.
        The semantic feasibility function should return True if it is
        acceptable to add the candidate pair (G1_node, G2_node) to the current
        partial isomorphism mapping. The logic should focus on semantic
        information contained in the edge data or a formalized node class.
        By acceptable, we mean that the subsequent mapping can still become a
        complete isomorphism mapping.
        """

        # type match
        G1_node_type = self.G1.info[G1_node]['type']
        G2_node_type = self.G2.info[G2_node]['type']
        if G1_node_type != G2_node_type:
            return False

        G1_node_iface = self.G1.agents[G1_node]
        G2_node_iface = self.G2.agents[G2_node]

        if self.test == 'embed':
            if len(G1_node_iface) < len(G2_node_iface):
                return False
        else:
            if len(G1_node_iface) != len(G2_node_iface):
                return False

        for site_name in G2_node_iface:
            if site_name not in G1_node_iface:
                return False
            else:
                if self.test == 'embed':
                    if G2_node_iface[site_name]['state'] != '#':
                        if G2_node_iface[site_name]['state'] != G1_node_iface[site_name]['state']:
                            return False
                elif self.test == 'iso':
                    if G2_node_iface[site_name]['state'] != G1_node_iface[site_name]['state']:
                        return False

                bond1 = G1_node_iface[site_name]['bond']
                bond2 = G2_node_iface[site_name]['bond']
                if '@' in bond1:
                    partner1, site1 = bond1.split(self.G1.bondsep)
                if '@' in bond2:
                    partner2, site2 = bond2.split(self.G2.bondsep)

                if self.test == 'iso':
                    if bond1 != '.' and bond2 != '.':
                        # both sites are bound
                        if partner2 in self.core_2:
                            if not (self.core_2[partner2] == partner1):
                                return False
                        if partner1 in self.core_1:
                            if not (self.core_1[partner1] == partner2):
                                return False
                        if site1 != site2:
                            return False
                    elif bond1 != bond2:
                        return False
                elif self.test == 'embed':  # _, don't care, stub, ., or bond
                    if bond2 == '.':
                        if bond1 != '.':
                            return False
                    elif '@' in bond2:  # specific bond
                        if not ('@' in bond1):
                            return False
                        else:
                            # both sites are bound
                            if partner2 in self.core_2:
                                if not (self.core_2[partner2] == partner1):
                                    return False
                                if site1 != site2:
                                    return False
                    elif bond2 == '_':  # unspecific bond
                        if bond1 == '.':
                            return False
                    elif '.' in bond2:  # stub ('.', as in free, is caught above)
                        if bond1 == '.':  # the site is free
                            return False
                        elif bond1 == '_':
                            return False  # is this True ?? (ask Pierre)
                        elif '@' in bond1:
                            ghost_site, ghost_type = bond2.split('.')
                            type1 = partner1.split('.')[0]
                            if ghost_type != type1 or ghost_site != site2:
                                return False
                        elif '.' in bond1:  # bond1 is also a stub
                            if bond2 != bond1:
                                return False
        return True


# ----------------------------------------------------------------------------------------

class GMState:
    """Internal representation of state for the GraphMatcher class.
    This class is used internally by the GraphMatcher class.  It is used
    only to store state specific data. There will be at most G2.order() of
    these objects in memory at a time, due to the depth-first search
    strategy employed by the VF2 algorithm.
    """

    def __init__(self, GM, G1_node=None, G2_node=None):
        """Initializes GMState object.
        Pass in the GraphMatcher to which this GMState belongs and the
        new node pair that will be added to the GraphMatcher's current
        isomorphism mapping.
        """
        self.GM = GM

        # Initialize the last stored node pair.
        self.G1_node = None
        self.G2_node = None
        self.depth = len(GM.core_1)

        if G1_node is None or G2_node is None:
            # Then we reset the class variables
            GM.core_1 = {}
            GM.core_2 = {}
            GM.inout_1 = {}
            GM.inout_2 = {}

        # Watch out! G1_node == 0 should evaluate to True.
        if G1_node is not None and G2_node is not None:
            # Add the node pair to the isomorphism mapping.
            GM.core_1[G1_node] = G2_node
            GM.core_2[G2_node] = G1_node

            # Store the node that was added last.
            self.G1_node = G1_node
            self.G2_node = G2_node

            # Now we must update the other two vectors.
            # We will add only if it is not in there already!
            self.depth = len(GM.core_1)

            # First we add the new nodes...
            if G1_node not in GM.inout_1:
                GM.inout_1[G1_node] = self.depth
            if G2_node not in GM.inout_2:
                GM.inout_2[G2_node] = self.depth

            # Now we add every other node...

            # Updates for T_1^{inout}
            new_nodes = set()
            for node in GM.core_1:
                new_nodes.update([neighbor for neighbor in GM.G1[node] if neighbor not in GM.core_1])
            for node in new_nodes:
                if node not in GM.inout_1:
                    GM.inout_1[node] = self.depth

            # Updates for T_2^{inout}
            new_nodes = set()
            for node in GM.core_2:
                new_nodes.update([neighbor for neighbor in GM.G2[node] if neighbor not in GM.core_2])
            for node in new_nodes:
                if node not in GM.inout_2:
                    GM.inout_2[node] = self.depth

    def restore(self):
        """Deletes the GMState object and restores the class variables."""
        # First we remove the node that was added from the core vectors.
        # Watch out! G1_node == 0 should evaluate to True.
        if self.G1_node is not None and self.G2_node is not None:
            del self.GM.core_1[self.G1_node]
            del self.GM.core_2[self.G2_node]

        # Now we revert the other two vectors.
        # Thus, we delete all entries which have this depth level.
        for vector in (self.GM.inout_1, self.GM.inout_2):
            for node in list(vector.keys()):
                if vector[node] == self.depth:
                    del vector[node]


if __name__ == '__main__':
    import kappathings as kt
    import re

    # usage scenarios

    # input a file containing one (large) Kappa string
    line = open('TestData/bigly.ka', 'r').read()
    # remove newlines that might occur in the file
    line = re.sub(r'\n+', ' ', line)
    # create a KappaComplex with whatever assignment of node identifiers arises
    # (that's the normalize=False flag).
    G1 = kt.KappaComplex(line, normalize=False)
    G2 = kt.KappaComplex(line, normalize=True)
    GM = GraphMatcher(G1, G2)
    print(f'G1 and G2 are isomorphic: {GM.is_isomorphic()}')
    print_map([GM.mapping])

    G1 = kt.KappaComplex('A(x[1],y[2]), B(x[2],y[3]), C(x[3],y[1])', normalize=False)
    for i in range(0, 10):
        # randomly permute identifiers:
        G2 = kt.KappaComplex('A(x[1],y[2]), B(x[2],y[3]), C(x[3],y[1])', randomize=True)
        GM = GraphMatcher(G1, G2)
        print(f'G1 and G2 are isomorphic: {GM.is_isomorphic()}')
    print('')

    G1 = kt.KappaComplex('A(x[1],y[2]), B(x[2],y[3]), C(x[3],y[1])', randomize=True)
    G2 = kt.KappaComplex('A(x[.],y[2]), B(x[2],y[3]), C(x[3],y[.])', normalize=True)
    GM = GraphMatcher(G1, G2)
    print(f'G1 and G2 are isomorphic: {GM.is_isomorphic()}')
    print(f'G2 is embeddable in G1: {GM.is_embeddable()}')

    G1 = kt.KappaComplex('A(x[1],y[2]), A(x[1],y[3]), C(x[2]), C(x[3]{p})', normalize=False)
    G2 = kt.KappaComplex('A(x[1],y[2]), A(x[1],y[3]), C(x[2]), C(x[3])', normalize=False)
    GM = GraphMatcher(G1, G2)
    print(f'G2 is embeddable in G1: {GM.is_embeddable()}')
    G2.name_list = kt.shift(G2.name_list)
    GM = GraphMatcher(G1, G2)
    print(f'G1 and G2 are isomorphic: {GM.is_isomorphic()}')

    G1 = kt.KappaComplex('A(x[.] y[2]), A(x[2] y[3]), A(x[3] y[4]), A(x[4] y[1]), B(x[1])')
    G2 = kt.KappaComplex('A(x y[2]), A(x[2] y)')
    print('G1:')
    print(G1.show())
    print('G2:')
    print(G2.show())
    maps = all_embeddings(G1, G2)
    print(f'number of embeddings from G2 into G1: {len(maps)} ')
    print_map(maps)

    G1 = kt.KappaComplex('A(b[1] a[2]), A(b[3] a[2]), B(a[1] x{p}), B(a[3] x{u})')
    G2 = kt.KappaComplex('A(b[2]), B(a[2] x)')
    maps = all_embeddings(G1, G2)
    print(f'number of embeddings from G2 into G1: {len(maps)} ')
    print_map(maps)
