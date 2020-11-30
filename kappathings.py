# Walter Fontana, 2020

import re
import copy
import json
import random
import networkx as nx


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def shift(in_list, n=1):
    new_list = []
    for i in range(n, len(in_list)):
        new_list += [in_list[i]]
    for i in range(0, n):
        new_list += [in_list[i]]
    return new_list


def cut_and_cup(c1, c2):
    """
    returns intersection and union of agent and bond set
    """
    agent_set1 = set()
    agent_set2 = set()
    bond_set1 = set()
    bond_set2 = set()
    for (agent1, site1), (agent2, site2) in c1.bonds:
        a1 = c1.info[agent1]['sID'] + c1.info[agent1]['type']
        a2 = c1.info[agent2]['sID'] + c1.info[agent2]['type']
        agent_set1.add(a1)
        agent_set1.add(a2)
        b = sorted([(a1, site1), (a2, site2)], key=lambda i: i[0])
        bond_set1.add(tuple(b))
    for (agent1, site1), (agent2, site2) in c2.bonds:
        a1 = c2.info[agent1]['sID'] + c2.info[agent1]['type']
        a2 = c2.info[agent2]['sID'] + c2.info[agent2]['type']
        agent_set2.add(a1)
        agent_set2.add(a2)
        b = sorted([(a1, site1), (a2, site2)], key=lambda i: i[0])
        bond_set2.add(tuple(b))
    return agent_set1 & agent_set2, agent_set1 | agent_set2, bond_set1 & bond_set2, bond_set1 | bond_set2


class KappaComplex:
    """
    :param data: a JSON formatted input -- string or list -- with schema [  list of agents  ] as defined in
                 kappasnap.KappaSnapShot OR a Kappa string
    :param count: (optional) the count associated with the complex

        constructs representations of the form

               self.agents[name] = { site_name: {'state': state, 'bond': bond_stub} }
            self.adjacency[name] = [ agent_names ]
                 self.info[name] = { 'id': identifier
                                     'type': agent_type
                                     'sID': optional KaSim-assigned identifier
                                     'degree': int n }
                    self.bonds   = { ( (agent1, site1), (agent2, site2) ) }  # a set!

            where everything, except e, is of type string, and
            * the interface dictionary of an agent is sorted by site name (needs Python 3.7+)
            * agent names are unique, consisting of type + identifier, eg Axin.42. (including the last dot), where
            the right and left separators (dots by default) can be chosen.
            * self.bonds is a list of unique tuples: (agent1, site1), (agent2, site2), lexicographically sorted
            on agent.
            * bonds are stubs of the form name@site indicating the name of the agent and site
            that anchors the other end of the bond
            * self.name_list is a list of agent names for the purpose of "sorting" the agent dictionary self.agents.
            A dictionary has no order by construction, but we can fake an order by iterating through it using an
            ordered list of its keys, whenever order is desired (such as in re-assigning identifiers or
            pretty printing)
            * All types are string, except when otherwise noted.

            Other properties:
            (2)
            self.size
            (3)
            self.composition (that's the sum formula, sorted by agent name)
    """

    def __init__(self, data, count=0, normalize=False, randomize=False, minimal=False):
        """
        :param data: the Kappa expression (JSON or kappa)
        :param count:
        :param normalize:
        :param randomize:
        :param minimal: set to True, if you don't use the representation for graph matching
        """
        # change these defs only if you know what you are doing
        self.symbols = r'[_~][a-zA-Z0-9_~+-]+|[a-zA-Z][a-zA-Z0-9_~+-]*'
        self.sID = r'x[0-9]+:'
        self.bondsep = '@'
        self.idsep = ('.', '.')  # not any of '(?:[_~][a-zA-Z0-9_~+-]+|[a-zA-Z][a-zA-Z0-9_~+-]*)'

        # properties of the complex
        self.count = count
        self.size = 0
        self.is_pattern = False
        self.composition = {}
        self.sum_formula = ''
        self.rarest_type = ''

        # data structures representing the complex; some redundancy here for convenience
        self.agents = {}
        self.adjacency = {}
        self.info = {}
        self.bonds = set()
        self.nbonds = {}  # nbonds[(a1, a2)] = n of bonds between a1 and a2 for a1 <= a2
        self.name_list = []
        self.type_slice = {}  # dict of pairs
        self.navigation = {}

        # auxiliary variables
        self.counter = 0
        self.next = 0

        # build regex'es
        site_name_re = r'(' + self.symbols + r')'
        internal_state_re = r'({(?:' + self.symbols + r'|[#]' + r')})'
        # binding_id_re = r'(\[(?:.|\d+)\])'
        binding_re = r'(\[(?:.*)\])'  # we still need to parse the bond expression
        # this can be . | # | number | site.agent (a stub)
        self.binding_state_re = \
            re.compile(r'^' + r'(?:\.|_|#|\d+)' + r'|(?:' + self.symbols + r')\.(?:' + self.symbols + r')')
        # using optional lookahead, since the internal state is optional and there is no prescribed order.
        # (gobble up the string with .*)
        self.site_re = \
            re.compile(r'^' + site_name_re + r'(?=.*' + internal_state_re + r')?' + r'(?=.*' + binding_re + r')?.*')
        agent_name_re = r'(' + self.symbols + r')'
        sID_re = r'(' + self.sID + r')?'
        agent_interface_re = r'\(([^()]*)\)'
        # to dissect agents
        self.agent_re = \
            re.compile(r'^' + sID_re + agent_name_re + agent_interface_re + r'$')
        # to find all agents (hence groups are non-capturing), before dissecting them
        self.agents_re = \
            re.compile(r'(?:' + self.sID + r')?(?:' + self.symbols + r')' + r'\([^()]*\)')

        # produce an internal representation of the complex from json or kappa
        if data:
            # determine input type
            if type(data) is list or 'node_type' in data:
                self.get_structure_from_json(data)
            else:
                self.get_structure_from_kappa(data)

        # size
        self.size = len(self.agents)
        # get the names list (for accessing the agent dictionary in a desired order)
        self.name_list = [k for k in self.agents]
        # get the composition
        self.get_composition()
        self.rarest_type = next(iter(self.composition))
        # sort name list by rarity (primary) and type (secondary)
        self.name_list = sorted(self.name_list,
                                key=lambda i: (self.composition[self.info[i]['type']], self.info[i]['type']))

        if randomize:
            # this is to re-index the same complex in different ways;
            # chiefly for the purpose of testing graph isomorphism
            self.permute_ids(self.randomize_ids())
        elif normalize:
            # agent list will be lexicographically sorted /and has increasing identifier/
            self.permute_ids(self.normalize_ids())  # this also makes the bond list

        if not minimal:
            # get the slice indices for each type in the sorted name_list
            start = 0
            for a_type in self.composition.keys():
                abundance = self.composition[a_type]
                self.type_slice[a_type] = (start, start + abundance)
                start += abundance

            # construct adjacency lists
            self.get_adjacency_lists()
            # assemble the navigation list (for rigid embeddings with kappamorph.py)
            self.make_navigation_list()

    def get_structure_from_json(self, data):
        """
        Given a complex in JSON snapshot format, construct the internal representation as in the class documentation.
        :param data: a JSON format (string or list) with schema [  list of agents  ] as defined in
                     kappasnap.KappaSnapShot
        This will not parse KaSim-supplied IDs. Use _from_kappa() for that!
        """

        def parse_link(link):
            """
            Parse the 'port_links' field in a way that deals with both old and
            and new JSON formats (gah).

            In new format we have extraneous viz information in links, e.g.:
              [[[0, 1], 0]]
            in old format we don't, i.e., we'd instead have:
              [[1, 0]]
            the problem in the current code is that the JSON parsing (and hence
            he actual spec being followed) is split between this file (kappathings.py)
            and kappasnap.py. All of this should be rewritten to use the snapshot
            parser in kappy, which correctly handles all the supported format.
            everything in KappaUtilities can easily use that representation. -YK
            (In new format, we have this (more redundant) representation
            link -> [ [ [ 0, list index of agent to which this agent is connected ],
                  index of connecting site of connected agent ] ] )
            """
            if type(link[0][0]) is list:
                # new links format
                partner_ind = link[0][0][1]
                site_ind = link[0][1]
            else:
                # old links format
                partner_ind = link[0][0]
                site_ind = link[0][1]
            return partner_ind, site_ind

        if type(data) is str:  # otherwise data is already parsed
            data = json.loads(data)
        self.bonds = set()  # set of tuples ((agent1, site1), (agent2, site2))
        this_complex = self.assign_identifiers(data)  # a list of 'agent' as per json schema
        for a in this_complex:
            agent_name = a['node_type']
            agent_type, identifier = self.extract_identifier(a['node_type'])
            degree = 0
            interface = {}
            for s in a['node_sites']:  # the sites
                site_name = s['site_name']
                link = s['site_type'][1]['port_links']  # site link (type list)
                if not link:
                    # agent isn't bound to anything
                    bond = '.'
                else:
                    degree += 1
                    # indices of the agent that's the binding partner, and of the
                    # site on that agent through which the bond is made -YK
                    partner_ind, site_ind = parse_link(link)
                    partner = this_complex[partner_ind]
                    agent2 = partner['node_type']
                    site2 = partner['node_sites'][site_ind]['site_name']
                    bond = agent2 + self.bondsep + site2
                    b = sorted([(a['node_type'], site_name), (agent2, site2)], key=lambda i: i[0])
                    self.bonds.add(tuple(b))  # unique bonds
                state_ = s['site_type'][1]['port_states']  # site state
                if not state_:
                    state = ''
                else:
                    state = state_[0]
                interface[site_name] = {'state': state, 'bond': bond}
            self.agents[agent_name] = interface
            self.info[agent_name] = {'id': identifier, 'type': agent_type, 'degree': degree}

    def get_structure_from_kappa(self, data):
        """
        Given a complex or pattern in kappa format, construct the internal representation
        as per class documentation.

        :param data: a Kappa expression of the complex
        """
        self.counter = 0
        expression = re.sub(r'\s+|\t+|\n+', ' ', data)  # remove line breaks and white matter
        # capture all agents
        match = self.agents_re.findall(expression)

        for agent in match:
            agent_type, identifier, agent_sID, interface = self.parse_agent(agent)
            if agent_sID is None: agent_sID = ''
            agent_name = self.attach_identifier(agent_type, identifier)
            self.agents[agent_name] = interface
            self.info[agent_name] = {'id': identifier, 'type': agent_type, 'sID': agent_sID, 'degree': -1}

        # replace bond ids with stub notation and get the set of bonds
        self.stubbify_bonds()

    def parse_agent(self, agent_expression):
        match = self.agent_re.match(agent_expression)
        if not match:
            exit('Invalid agent declaration <' + agent_expression + '>')
        agent_sID = match.group(1)
        agent_type = match.group(2)
        self.counter += 1
        identifier = str(self.counter)
        interface = {}

        # parse the agent interface
        iface = match.group(3)
        # since Kappa allows commas or whitespace as separators,
        # swap all commas for spaces and split by whitespace
        sites = iface.replace(',', ' ').split()
        for item in sites:
            try:
                site_name, state, bond = self.parse_site(item)
                interface[site_name] = {'state': state, 'bond': bond}
                # sort interface by key
                # interface = dict(sorted(interface.items()))
            except:
                exit('Could not parse site ' + item + ' in ' + agent_expression)

        return agent_type, identifier, agent_sID, interface

    def parse_site(self, site_expression):
        match = self.site_re.match(site_expression)
        if not match:
            exit('Could not parse site ' + site_expression)
        # return site name, internal state and binding state (without parentheses)
        site_name = match.group(1)
        if match.group(2):  # the modification state; it may be absent, so we need to check
            internal_state = match.group(2)[1:-1]  # remove parens
        else:
            internal_state = '#'  # don't care

        binding_state = '#'  # don't care (absent) by default
        if match.group(3):  # there is an explicit binding state
            binding_expression = match.group(3)[1:-1]  # remove parens
            match = self.binding_state_re.match(binding_expression)  # continue parsing
            if match:
                binding_state = match.group(0)  # either '.' or '#' or number or stub
                # warning: if the site name starts with '_' we have a problem; fix later...
            else:
                exit('Could not parse binding state ' + binding_expression)

        return site_name, internal_state, binding_state

    def stubbify_bonds(self):
        """
        replace numeric bond labels with unique bond stubs
        generate the bond set (self.bonds)
        """
        # If we are dealing with an object that contains a bond pattern, the degree of a node has no meaning.
        # The degree is used only for VF2 isomorphism checking, but not for pattern embeddings.
        self.bonds = set()
        bonds = {}
        for name in self.agents:
            degree = 0
            for site in self.agents[name]:
                link = self.agents[name][site]['bond']
                if link != '.':
                    if is_number(link):
                        degree += 1
                        if link in bonds:
                            [(name1, site1)] = bonds[link]
                            # self.agents[name1][site1]['bond'] = name + self.bondsep + site
                            # self.agents[name][site]['bond'] = name1 + self.bondsep + site1
                            self.agents[name1][site1]['bond'] = ''.join([name, self.bondsep, site])
                            self.agents[name][site]['bond'] = ''.join([name1, self.bondsep, site1])
                            b = sorted([(name1, site1), (name, site)], key=lambda i: i[0])
                            self.bonds.add(tuple(b))  # collect unique bonds
                        else:
                            bonds[link] = [(name, site)]
                    elif self.bondsep in self.agents[name][site]['bond']:
                        degree += 1
                    else:
                        # bond state is a ghost, or '_', or '#'
                        degree = -1  # reset and flag, just in case
                        self.is_pattern = True

            self.info[name]['degree'] = degree

    def kappa_expression(self):
        """
        convert internal rep into a kappa string
        """
        # If we are dealing with an object that contains a bond pattern, the degree of a node has no meaning.
        # The degree is used only for VF2 isomorphism checking, but not for pattern embeddings.

        i = 1
        num = {}
        for (p, q) in self.bonds:
            num[p] = i
            num[q] = i
            i += 1

        s = ''
        for name in self.agents:
            s += f'{self.info[name]["type"]}('
            for site in self.agents[name]:
                s += f'{site}'
                state = self.agents[name][site]['state']
                if state != '#':
                    s += '{' + f'{state}' + '}'
                link = self.agents[name][site]['bond']
                if link == '.':
                    s += '[.] '
                elif link == '#':
                    # s += '[#] '
                    s += ''
                elif self.bondsep in link:
                    ag, ste = link.split(self.bondsep)
                    s += f'[{num[(ag, ste)]}] '
                else:
                    s += f'[{link}] '
            s = s[:-1] + '), '
        return s[:-2]

    def get_adjacency_lists(self):
        """
        construct the adjacency views
        """
        for name1 in self.agents:
            adjacency = []
            for s1 in self.agents[name1]:
                if self.bondsep in self.agents[name1][s1]['bond']:
                    name2, x, x = self.agents[name1][s1]['bond'].partition(self.bondsep)
                    adjacency.append(name2)
            # sort on rarity and type
            # adj = sorted(adjacency, key=lambda i: (self.composition[self.info[i]['type']], self.info[i]['type']))
            self.adjacency[name1] = adjacency

    def get_bond_numbers(self):  # fuse with adjacency list function ?
        """
        fills in the nbonds dict, which reports the number of bonds between n1 and any other node it is bound to
        """
        s = set()
        for (a1, s1), (a2, s2) in self.bonds:
            if (a1, a2) in s:
                self.nbonds[(a1, a2)] += 1
            else:
                s.add((a1, a2))
                self.nbonds[(a1, a2)] = 1

    def extract_identifier(self, name):
        agent_type, identifier = name.split(self.idsep[0])[:2]
        if self.idsep[1] and self.idsep[0] != self.idsep[1]:
            identifier = identifier[:-1]
        return agent_type, identifier

    def attach_identifier(self, agent_type, id):
        return agent_type + self.idsep[0] + id + self.idsep[1]

    def make_navigation_list(self):
        # self.navigation[a1][a2] contains a site of a1 that anchors a bond to a2
        # (For the purpose of this array, we don't care about multiple bonds between the same agents.)
        # This is similar to self.bonds, but organized as a dictionary for convenience.
        self.navigation = {}
        for (a1, s1), (a2, s2) in self.bonds:  # names a1 and a2 in bonds have id attached
            self.navigation[(a1, a2)] = s1
            self.navigation[(a2, a1)] = s2

    def assign_identifiers(self, complex):  # complex is in JSON format
        """
        Generate unique names for the agents of a complex. The names are of the form 'type LDEL identifier RDEL',
        with {L,R}DEL delimiters, by default '.' (dot). For example: APC.2., which means the agent is of type APC
        with identifier 2 (because the complex contains multiple APCs).

        :param complex: a list of the JSON form [  list of agents  ]  (see schema)
        :return: a new dictionary of agents of the format { "node_type" : string, "node_sites" : [ list of sites ] }
        with altered string values. The original is not modified.
        """
        new_complex = copy.deepcopy(complex)  # make a "real" copy, since we're going to change it!
        self.counter = 0
        for a in new_complex:
            self.counter += 1
            a['node_type'] = self.attach_identifier(a['node_type'], str(self.counter))
        return new_complex

    def normalize_ids(self):
        permute = {}
        for i in range(0, self.size):
            name = self.name_list[i]
            permute[self.info[name]['id']] = str(i + 1)
        return permute

    def randomize_ids(self):
        l = [i for i in range(1, len(self.agents) + 1)]
        random.shuffle(l)
        permute = {}
        for i in range(0, self.size):
            name = self.name_list[i]
            permute[self.info[name]['id']] = str(l[i])
        return permute

    def permute_ids(self, permute):
        """
        We (re)assign identifiers according to
        :param permute
        Note: nodes must be lexicographically sorted by name.
        """
        self.bonds = set()  # reset

        # apply permutation
        info = {}
        renamed = {}
        for i in range(0, self.size):
            name1 = self.name_list[i]
            id1 = self.info[name1]['id']
            type1 = self.info[name1]['type']
            sID = self.info[name1]['sID']
            new_id1 = permute[id1]
            new_name1 = self.attach_identifier(type1, new_id1)
            self.name_list[i] = new_name1
            renamed[new_name1] = {}
            info[new_name1] = {'id': new_id1, 'type': type1, 'sID': sID, 'degree': self.info[name1]['degree']}
            for site1 in self.agents[name1]:
                renamed[new_name1][site1] = {}
                renamed[new_name1][site1]['state'] = self.agents[name1][site1]['state']
                if self.bondsep in self.agents[name1][site1]['bond']:
                    agent2, site2 = self.agents[name1][site1]['bond'].split(self.bondsep)  # name.oldid
                    type2, id2 = self.extract_identifier(agent2)
                    new_name2 = self.attach_identifier(type2, permute[id2])
                    renamed[new_name1][site1]['bond'] = new_name2 + self.bondsep + site2
                    b = sorted([(new_name1, site1), (new_name2, site2)], key=lambda i: i[0])
                    self.bonds.add(tuple(b))
                else:
                    renamed[new_name1][site1]['bond'] = self.agents[name1][site1]['bond']
        self.agents = renamed
        self.info = info

    def get_composition(self):
        """
        Get the "sum formula" of a complex. Agents are ordered by increasing abundance within the complex.
        """
        comp = {}
        for a in self.agents:
            a_type = self.info[a]['type']
            if a_type in comp.keys():
                comp[a_type] += 1
            else:
                comp[a_type] = 1

        # sort the dict by value and then key:
        self.composition = {k: v for k, v in sorted(comp.items(), key=lambda item: (item[1], item[0]))}

        self.sum_formula = ''
        for a_type in self.composition.keys():
            self.sum_formula += (a_type + '{' + str(self.composition[a_type]) + '}')

    def is_multigraph(self):
        """
        Test if the set of bonds implies a multigraph.
        :return: True / False
        """
        s = set()
        for (a1, s1), (a2, s2) in self.bonds:
            if (a1, a2) in s:
                return True
            else:
                s.add((a1, a2))
        return False

    def nodes(self):
        """
        This emulates the networkx' G.nodes() method returning a list of node names.
        """
        # return [k for k in self.agents]
        return self.name_list

    def order(self):
        """
        This works like __len__. For compatibility with networkx representation.
        """
        return self.size

    def degree(self):
        """
        This emulates networkx G.degree(), returning a list of (node, degree) pairs
        """
        l = []
        for name in self.agents:
            if self.info[name]['degree'] == -1:  # an unspecified degree suggests a pattern
                return []
            l += [(name, self.info[name]['degree'])]
        return l

    def number_of_edges(self, name1, name2):
        """
        Returns the number of bonds between two nodes; used in VF2 isomorphism.
        """
        if name1 > name2:
            pair = (name2, name1)
        else:
            pair = (name1, name2)

        if pair in self.nbonds:
            return self.nbonds[pair]
        else:
            return 0

    def __str__(self):
        """
        exhaustive pretty print
        """
        info = f'expression has {self.size} agents:\n'
        if self.is_pattern:
            info += f'expression is a pattern\n'
        for i in range(0, self.size):
            name = self.name_list[i]
            interface = ''
            iface = self.agents[name]
            for s in iface:
                interface += s + '{' + iface[s]['state'] + '}' + '[' + iface[s]['bond'] + '] '
            if not self.is_pattern:
                info += f"[d: {self.info[name]['degree']}] "
            info += self.info[name]['sID'] + name + '(' + interface[:-1] + ')\n'
        info += f'{len(self.bonds)} specific bonds:\n'
        for (a1, s1), (a2, s2) in self.bonds:
            info += a1 + self.bondsep + s1 + ' <-> ' + a2 + self.bondsep + s2 + '\n'
        info += 'composition:\n'
        info += self.sum_formula + '\n'
        return info

    def show(self, internal=False):
        """
        print representation
        if internal=False, print standard kappa expression
        """
        if internal:
            info = ''
            for i in range(0, self.size):
                name = self.name_list[i]
                interface = ''
                iface = self.agents[name]
                for s in iface:
                    interface += s + '{' + iface[s]['state'] + '}' + '[' + iface[s]['bond'] + '] '
                info += self.info[name]['sID'] + name + '(' + interface[:-1] + '), '
            print(info[:-2])  # remove last comma+blank
        else:
            print(self.kappa_expression())

    def __iter__(self):
        return iter(self.name_list)

    # Implementing __next__ as below does not work... The reason is that breaking out of a loop does not
    # reset self.next! Use Python's iter() function as above.
    #
    # def __next__(self):
    #     """
    #     Make KappaComplex an iterable, using the name_list. Ie, iterate through self.name_list and use the returned
    #     list element to access self.agents (or any other name-indexed data structure).
    #     """
    #     this = self.next
    #     if this == self.size:
    #         self.next = 0  # reset
    #         raise StopIteration
    #     else:
    #         self.next += 1
    #         return self.name_list[this]

    def __len__(self):
        return self.size

    def __getitem__(self, name):
        """
        Makes C[name] return the list of neighbors of node name; emulates the adjacency view of networkx
        """
        return self.adjacency[name]


# -----------------------------------------------------------------

if __name__ == '__main__':
    import time

    # e1 = " x222667:P(a1[.] a2[.] a3[.] d[1]), x251065:P(a1[.] a2[.] a3[.] d[1])"
    # e2 = " x22266:P(a1[.] a2[.] a3[.] d[1]), x251065:P(a1[.] a2[.] a3[.] d[1])"
    # agents, bonds = intersect_complexes(KappaComplex(e1), KappaComplex(e2))
    # print(agents)
    # print(bonds)

    # usage scenarios

    # print("from Kappa:")
    # # input a simple Kappa string to create the KappaComplex representation
    # data = ' A(o[1], p[2] t{p}[3]), B(x[1] y[2] z[.]), C(w[3], y[z.B])'
    # c = KappaComplex(data, count=175, normalize=True)
    # # print it
    # print(c)
    # # short version
    # c.show()
    #
    # print(f'count: {c.count}')
    # print(f'size: {c.size}')
    # print(f'sum formula: {c.sum_formula}')
    # print('')
    #
    # print("from JSON:")
    # # input a JSON string of the same complex as above
    # data = '[{"node_type":"A","node_sites":[{"site_name":"o","site_type":["port",{"port_links":[[[0,1],0]],' \
    #        '"port_states":[]}]},{"site_name":"p","site_type":["port",{"port_links":[[[0,1],1]],"port_states":[]}]},' \
    #        '{"site_name":"t","site_type":["port",{"port_links":[[[0,2],0]],"port_states":["p"]}]}]},' \
    #        '{"node_type":"B","node_sites":[{"site_name":"x","site_type":["port",{"port_links":[[[0,0],0]],' \
    #        '"port_states":[]}]},{"site_name":"y","site_type":["port",{"port_links":[[[0,0],1]],' \
    #        '"port_states":[]}]},{"site_name":"z","site_type":["port",{"port_links":[],"port_states":[]}]}]},' \
    #        '{"node_type":"C","node_sites":[{"site_name":"w","site_type":["port",{"port_links":[[[0,0],2]],' \
    #        '"port_states":[]}]}]}]'
    # c = KappaComplex(data, normalize=True)
    # print(c)
    # #
    # print("from file:")
    # # input a file containing one (large) Kappa string
    # line = open('TestData/bigly.ka', 'r').read()
    # # remove newlines that might occur in the file
    # line = re.sub(r'\n+', ' ', line)
    # # create a KappaComplex with whatever assignment of node identifiers arises
    # # (that's the normalize=False flag).
    # c1 = KappaComplex(line)
    # print(c1)
    # print("list of nodes:")
    # # a list of node names of the complex
    # print(c1.nodes())
    # # the complex appears as an iterable
    # print('iterating through the nodes of the complex and showing adjacency views')
    # # testing __iter__ and __getitem__ methods
    # for n in c1:
    #     print(f'adjacency of {n}: {c1[n]}')
    # print('')
    # #
    # print("normalize identifiers:")
    # # "normalize" the identifiers, which means: when node types are sorted lexicographically,
    # # nodes are assigned successively increasing identifiers.
    # c2 = KappaComplex(line, normalize=True)
    # print(c2)
    #
    # print("randomly permute identifiers:")
    # # assign identifiers randomly; this is useful for testing the isomorphism implementation
    # c3 = KappaComplex(line, randomize=True)
    # print(c3)
    #
    # print("list of nodes:")
    # # a list of node names of the complex
    # print(c3.nodes())
    # # the complex appears as an iterable
    # print('iterating through the nodes of the complex and showing adjacency views')
    # # testing __iter__ and __getitem__ methods
    # for n in c3:
    #     print(f'adjacency of {n}: {c3[n]}')
    # print('')
    # print(c3.sum_formula)
    # print(c3.type_slice)
    # for a_type in c3.type_slice.keys():
    #     (start, stop) = c3.type_slice[a_type]
    #     print(c3.name_list[start:stop])

    # print("big stuff from file:")
    # line = open('TestData/monster.ka', 'r').read()
    # line = re.sub(r'\n+', ' ', line)
    # start = time.process_time()
    # c1 = KappaComplex(line)
    # end = time.process_time()
    # print(f'seconds: {end - start}')
