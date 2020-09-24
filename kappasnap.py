import kappathings as kt
import kappy
import re
import gzip
import json
import os
import sys


class SnapShot:
    """
    Provides class Snapshot.
    Read a snapshot from a JSON or .ka file and generate internal representations of complexes.
    The complexes are stored as a list self.complex[] of objects of the kappathings.KappaComplex class.

    Basic usage:
        snapshot = SnapShot('TestData/snap.ka')
        print(f'{snapshot.number_of_distinct_complexes} complexes')
        for c in snapshot.complex:
            print(c)
    """

    def __init__(self, file, use_kappy=False):
        self.file = file
        self.time = 0.
        self.event = 0
        self.number_of_distinct_complexes = 0
        self.complexes = []  # list of 'KappaComplex'es

        if self.file.endswith(".json") or self.file.endswith(".json.gz"):
            S = JsonSnapShot(file, use_kappy=use_kappy)
        elif self.file.endswith(".ka"):
            S = KappaSnapShot(file, use_kappy=use_kappy)
        else:
            raise Exception("Unknown file extension %s" % self.file)
        self.time = S.time
        self.event = S.event
        self.number_of_distinct_complexes = S.number_of_distinct_complexes
        self.complexes = S.complexes

    def get_size_distribution(self, dictionary=False):
        """
        Generates the size distribution present in the snapshot

        :return: sorted list of pairs (length, count)
                 eg. [(1, 8583), (2, 2642), (7, 836), (4, 253), (103, 82)]
        """
        length = {}
        for c in self.complexes:
            if c.size in length.keys():
                length[c.size] += c.count
            else:
                length[c.size] = c.count

        if dictionary:
            d = {'length': [], 'count': []}
            for l, c in sorted(length.items(), key=lambda i: i[0], reverse=False):
                d['length'].append(l)
                d['count'].append(c)
            return d
        else:
            return [(l, c) for l, c in sorted(length.items(), key=lambda i: i[0], reverse=False)]


class KappaSnapShot:
    def __init__(self, kappa_file, use_kappy):
        # aux vars
        self.data = None
        self.use_kappy = use_kappy
        self.current_line = ''

        self.kappa_file = kappa_file
        self.time = 0.
        self.event = 0
        self.number_of_distinct_complexes = 0
        self.complexes = []  # list of 'KappaComplex'es

        self.load_and_unpack()
        self.data = None  # clear

    def load_and_unpack(self):
        """
        Load a Kappa snapshot file.
        """
        if not os.path.isfile(self.kappa_file):
            raise Exception("Cannot find snapshot file %s" % self.kappa_file)
        else:
            with open(self.kappa_file, "r") as data:
                self.data = data
                self.event = float(data.readline().split('Event:')[1][:-2].strip())
                data.readline()
                t = data.readline().split('T0')[1][:-2]
                self.time = float(re.sub(r'"', ' ', t).strip())
                data.readline()
                self.current_line = data.readline()[:-1]  # this should be the first line of the first complex

                while True:
                    entry = self.next_complex_from_file()
                    if not entry:
                        break
                    if self.use_kappy:
                        komplex = None # until we can repair this
                        sys.exit('Avoid kappy for now')
                        # komplex = kappy.KappaComplex.from_string(entry)
                    else:
                        # DIY...
                        # parse the entry
                        match = re.findall(r'%init: (.*?) \/\*(.*?) agents\*\/ (.*?)$', entry)[0]
                        # build the internal representation
                        komplex = kt.KappaComplex(match[2].strip(), count=int(match[0]))
                    self.complexes.append(komplex)
            self.number_of_distinct_complexes = len(self.complexes)

    def next_complex_from_file(self):
        entry = self.current_line
        line = ''
        while True:
            line = self.data.readline()
            if not line:
                self.current_line = ''
                break
            elif '%init' in line:
                if 'agents*/' in line:
                    self.current_line = line[:-1]
                else:
                    self.current_line = ''
                break
            entry = entry + line[:-1]  # remove \n
        return entry


class JsonSnapShot:
    def __init__(self, json_file, use_kappy):
        """
        :param json_file (type string): snapshot file in JSON format

        The informal schema of the json snapshot file (as of May 2020):

        {
        "snapshot_tokens": [ list ]
        "snapshot_file": string,
        "snapshot_event": int,
        "snapshot_time": float,
        "snapshot_agents": [ list of complexes ]     // this is really a misnomer: should be snapshot_complexes...

        ===> will be: "snapshot_complexes": [ list of complexes ]

        complex -> [ count, [  [  list of agents  ]  ]  ]
        agent ->  { "node_type" : string, "node_sites" : [ list of sites ] }

        ===> will be: agent ->  { "type" : string, "sites" : [ list of sites ] }

        site ->  { "site_name" : string, "site_type" : [ list of site states ] }

        ===> will be: site ->  { "name" : string, "type" : [ list of site states ] }

        site state -> [ "port", { "port_links" : [ list of links ], "port_states" : [ list of internal states ] } ]

        ===> will be: site state -> [ "port", { "links" : [ list of links ], "states" : [ list of internal states ] } ]

        link -> [ [ [ 0, list index of agent to which this agent is connected ],
                  index of connecting site of connected agent ] ]
        internal state -> [ state ]
        }

        establishes:
            self.complex = [ a list of 'kappathings.KappaComplex' ]

        methods:
            self.get_size_distribution
        """
        self.json_file = json_file
        self.use_kappy = use_kappy
        self.snap_name = None  # this is the internal name of the snapshot
        self.data = None
        self.time = 0.
        self.event = 0
        self.number_of_distinct_complexes = 0
        self.complexes = []  # list of 'KappaComplex'es

        self.load()
        self.unpack()
        self.data = None  # garbage collect the raw JSON data

    def load(self):
        """
        Load a JSON snapshot file. Can be a plain JSON file (.json) or
        a compressed file (.json.gz).
        """
        if not os.path.isfile(self.json_file):
            raise Exception("Cannot find snapshot file %s" % self.json_file)
        if self.json_file.endswith(".json.gz"):  # uncompress first
            with gzip.open(self.json_file, "rb") as data:
                try:
                    self.data = json.load(data)
                except:
                    raise Exception('Data format issue?')
        else:
            # read uncompressed file
            with open(self.json_file, "r") as data:
                try:
                    self.data = json.load(data)
                except:
                    raise Exception('Data format issue?')

    def unpack(self):
        #
        # Reminder:
        # data['snapshot_agents'] [i] [1] [0] [j]
        #                            |
        #                            ith complex
        #                                        |
        #                                        jth agent
        #

        self.snap_name = self.data['snapshot_file']
        self.time = float(self.data['snapshot_time'])
        self.event = int(self.data['snapshot_event'])
        self.number_of_distinct_complexes = len(self.data['snapshot_agents'])
        # this function encodes part of the JSON spec that's being parsed
        # there should instead be just one spot where JSON (or KA) is parsed
        # and converted into an internal representation. ideally this should be
        # done by kappy and this package should be agnostic to the
        # file format spec -YK
        for c in self.data['snapshot_agents']:
            abundance, comp_info = c[0], c[1]
            # deal with old JSON format which didn't have
            # here a list of lists (ugh..) -YK
            if type(comp_info[0]) != list:
                comp_info = [comp_info]
            if self.use_kappy:
                komplex = None # until we can repair this
                sys.exit('Avoid kappy for now')
                # komplex = kappy.KappaComplex.from_string(comp_info[0])
            else:
                # DYI...
                komplex = kt.KappaComplex(comp_info[0], count=c[0])
            self.complexes.append(komplex)


# -----------------------------------------------------------------


if __name__ == '__main__':

    import kappamorph as km
    import kappaviz as viz

    # Yarden tests (these take 10 minutes for rigid and 43 minutes for VF2)

    snap1_obj = SnapShot('TestData/snap_large.ka')
    snap1_size = len(snap1_obj.complexes)

    intersection = []
    SGM = km.SiteGraphMatcher()
    for complex1 in snap1_obj.complexes:
        for complex2 in snap1_obj.complexes:
            if SGM.isomorphic(complex1, complex2):
                intersection.append(complex1)
                break
    print(f'size of intersection: {len(intersection)}')

    # intersection_vf = []
    # for complex1 in snap1_obj.complex:
    #     for complex2 in snap2_obj.complex:
    #         if km.isomorphic_vf2(complex1, complex2):
    #             intersection_vf.append(complex1)
    #             break
    # print(f'size of intersection: {len(intersection_vf)}')
    # print(f'intersections are equal: {intersection == intersection_vf}')

    # usage scenarios

    # print('from Kappa')
    # snapshot1 = SnapShot('TestData/snap.ka')
    # print(f'{snapshot1.number_of_distinct_complexes} complexes')
    # for c in snapshot1.complex:
    #     print(c)
    #
    # print('from JSON')
    # snapshot2 = SnapShot('TestData/snap.json')
    # print(f'{snapshot2.number_of_distinct_complexes} complexes')
    # for c in snapshot2.complex:
    #     print(c)
    #
    # print('large file')
    # snapshot = SnapShot('TestData/snap_large.ka')
    # print(f'{snapshot.number_of_distinct_complexes} complexes')
    # max = 10
    # i = 0
    # while True:
    #     print(snapshot.complex[i])
    #     i += 1
    #     if i == max:
    #         break
    #
    # # size distribution
    # dist = snapshot.get_size_distribution()
    # for item in dist:
    #     print(f'size: {item[0]}   count: {item[1]}')
    # # viz.show_ranked_complexes(snapshot, prog='sfdp', sort='size', cutoff=4, cols=2, rows=2, node_size=20, font_size=4)
    # # viz.show_ranked_complexes(snapshot, sort='count', cutoff=4, cols=2, rows=2, node_size=20)
