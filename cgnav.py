from sets import Set
from time import time
from collections import OrderedDict
import MySQLdb
import sys
import curses


class PatternConcept(object):
    def __init__(self,extent,intent):
        print("extent={0}".format(extent))
        self.extent = extent
        self.intent = intent
        self.upper = []
        self.lower = []

    def __str__(self):
        return self.extent.__str__()

    def __repr__(self):
        return self.extent.__repr__()

    def lessequal(self,concept):
        return set(self.extent) <= set(concept.extent)


class PatternStructure(object):

    def __init__(self,delta,M):
        self.delta = delta
        self.M = M
        self.__generators = None

    def generators(self):

        if self.__generators:
            return self.__generators

        generators = []
        components = self.decompose(self.delta)
        for comp in components:
            for n in comp.E[0]:
                generators.append( (n.name,self.copy((n,comp))) )

        self.__generators = generators
        return generators


    def windowed_product(self,pattern1,pattern2):
        return self.__windowed_product(pattern1[0],pattern2[0])


    def product(self,graph1,graph2):

        matrix = {}

        components = []
        for n1 in graph1.E[0]:
            for n2 in graph2.E[0]:
                if not n1.node_id + n2.node_id in matrix:
                    components.append(self.__windowed_product(n1,n2,matrix)[1])
        return components


    def __windowed_product(self,node1,node2,matrix=None):

        if matrix is None:
            matrix = {}

        node0 = self.sup(node1,node2)
        matrix[node1.node_id + node2.node_id] = node0
        queue = [(node1,node2)]

        E = [ [] for i in range(len(self.M)) ]
        E[0].append(node0)

        while queue:
            node1,node2 = queue.pop(0)

            node = matrix[node1.node_id + node2.node_id]
            for n in range(1,len(self.M)):
                for k in range(n):
                    for r1 in node1.incidences.get((k,n),[]):
                        for r2 in node2.incidences.get((k,n),[]):
                            star = []
                            for j in range(n):
                                if j != k:
                                    n1 = r1.nodes[j]
                                    n2 = r2.nodes[j]
                                    pos = n1.node_id + n2.node_id
                                    if pos in matrix:
                                        neighbor = matrix[pos]
                                    else:
                                        neighbor = self.sup(n1,n2)
                                        E[0].append(neighbor)
                                        matrix[pos] = neighbor
                                        queue.append((n1,n2))
                                else:
                                    neighbor = node
                                star += [neighbor]
                            if k == 0:
                                r = RNode(star)
                                r.atts = r1.atts & r2.atts
                                E[n].append(r)

        return (node0,AbstractConceptGraph(E))

    def terminal(self):
        n0 = Node((),"UNIT")
        n0.atts = set(self.M[0])
        E = [[n0]]
        for i in range(1,len(self.M)):
            r = RNode([n0]*i)
            r.atts = set(self.M[i])
            E += [[r]]
        return (n0,AbstractConceptGraph(E))

    def copy(self,pattern):
        return self.windowed_product(pattern,self.terminal())

    def decompose(self,graph):
        return self.product(graph,self.terminal()[1])

    def sup(self,node1,node2):
        node_id = node1.node_id + node2.node_id
        name = "{0}:{1}".format(node1.name,node2.name).strip(":")
        node = Node(node_id,name)
        node.atts = node1.atts & node2.atts
        return node

    def find_morphisms(self,pattern1,pattern2):

        n1 = self.copy(pattern1)[0]
        n2 = self.copy(pattern2)[0]
        if not n1.atts <= n2.atts:
            return False
        r1 = self.create_neighborhoods(n1) 
        if r1 is None:
            return True
        n1.target = n2
        # print "matching " + n1.__str__() + " with " + n2.__str__()
        if self.find_match(r1):
            return True
        else:
            return False

    def create_neighborhoods(self,n0):

        n0.distance = 0
        sortlist = [n0]
        rnodes = []
        i = 0

        while i < len(sortlist):
            node = sortlist[i]
            d = node.distance
            for n in range(len(self.M)):
                for k in range(n):
                    for rnode in node.incidences.get((k,n),[]):
                        if rnode.origin is None:
                            rnode.origin = (k,n)
                            rnodes.append(rnode)
                            for j,neighbor in enumerate(rnode.nodes):
                                if neighbor.distance is None:
                                    neighbor.distance = d+1
                                    sortlist.append(neighbor)
            i = i+1

        for i in range(len(rnodes)-1):
            rnodes[i].next = rnodes[i+1]
            rnodes[i+1].prev = rnodes[i]

        if not rnodes:
            return None

        return rnodes[0]

    def find_match(self,r1):
        if r1.origin is None:
            print "WARNING: TRIPLE NOT CLASSIFIED"
            return False

        for r2 in r1.nodes[r1.origin[0]].incidences.get(r1.origin,[]):
            if self.matches(r1,r2):

                undo = [ r1.nodes[i].target for i in range(len(r1.nodes)) ] 
                for k in range(len(r1.nodes)):
                    r1.nodes[k].target = r2.nodes[k]

                if r1.next is None:

                    return True

                elif self.find_match(r1.next):
                    return True

                for k in range(len(r1.nodes)):
                    r1.nodes[k].target = undo[k]

        return False

    def matches(self,r1,r2):
        for k in range(len(r1.nodes)):
            if r1.nodes[k].target not in [None,r2.nodes[k]]:
                return False
        for k in range(len(r1.nodes)):
            if not r1.nodes[k].atts <= r2.nodes[k].atts:
                return False
        return True


class LatticeRunner(object):
    def __init__(self,lattice):
        self.lattice = lattice
        self.location = lattice[-1]

    def where(self):
        print self.location

    def lower(self):
        for i,c in enumerate(self.location.lower):
            print i, ") ", c

    def upper(self):
        for i,c in enumerate(self.location.upper):
            print i, ") ", c

    def up(self,i):
        self.location = self.location.upper[i]

    def down(self,i):
        self.location = self.location.lower[i]

# depends on: PatternStructure.{terminal,generators,windowed_product,find_morphisms}
# depends on: PatternConcept.{__init__,intent,extent,upper,lower,lessequal}
def build(pst):

    unit = pst.terminal()

    stack,lattice = [],[]
    unit_extent,unit_additions = [],[]

    for g,p in pst.generators():
        if pst.find_morphisms(unit,p):
            unit_extent.append(g)
        else:
            unit_additions.append((g,p))

    c_zero = PatternConcept(unit_extent,unit)
    stack.append((c_zero,unit_additions))
    lattice.append(c_zero)


    while stack:

        concept,additions = stack[-1]

        if not additions:
            stack.pop()
            continue

        g0,p0 = additions.pop()
        next_intent = pst.windowed_product(concept.intent,p0)

        next_extent,next_additions = [],[]
        j,approved = 0,False

        for g,p in pst.generators():
            if j<len(concept.extent) and concept.extent[j] == g:
                next_extent.append(g)
                j += 1
            elif g == g0:
                next_extent.append(g)
                approved = True
            else:
                if approved:
                    if pst.find_morphisms(next_intent,p):
                        next_extent.append(g)
                    else:
                        next_additions.append((g,p))
                else:
                    if pst.find_morphisms(next_intent,p):
                        break
                    else:
                        continue

        if approved:
            next_concept = PatternConcept(next_extent,next_intent)
            stack.append((next_concept,next_additions))

            for c in lattice:
                if c.lessequal(next_concept):
                    if not any(cu.lessequal(next_concept) for cu in c.upper):
                        c.upper.append(next_concept)
                        next_concept.lower.append(c)

            lattice.append(next_concept)

    return lattice

class ParseError(Exception):
    pass

class Node(object):

    def __init__(self,node_id,name):
        self.node_id = node_id
        self.name = name
        self.atts = set()
        self.incidences = {}
        self.distance = None
        self.target = None

class RNode(object):

    def __init__(self,nodes):
        n = len(nodes)
        for i,node in enumerate(nodes):
            node.incidences.setdefault((i,n),[]).append(self)
        self.atts = set()
        self.nodes = nodes
        self.origin = None
        self.next = None
        self.prev = None

class AbstractConceptGraph(object):

    def __init__(self,E):
        self.E = E

class PCFReader(object):

    def read(self,filename):

        with open(filename) as pcf_file:
            lines = [ line.strip() for line in pcf_file.readlines() ]

        it = iter(lines)
        if not it.next() == "B-PCF":
            raise ParseError("Wrong File Format")

        it.next()
        n_objs = int(it.next())
        n_atts = int(it.next())

        it.next()
        nodeseq = [ Node((i,),it.next()) for i in range(n_objs) ]
        nodes = { node.name:node for node in nodeseq }
        atts = [ it.next() for k in range(n_atts) ]

        for k in range(n_objs):
            nodeseq[k].atts = set([ atts[i] for i,char in enumerate(it.next()) if char in ["x","X"] ])

        E = [ nodeseq ]
        M = [ atts ]

        while True:
            try:
                it.next()
                n_objs = int(it.next())
                n_atts = int(it.next())

                it.next()
                rnodes = [ RNode([ nodes[name] for name in it.next().split(",") ]) for i in range(n_objs) ]
                atts = [ it.next() for k in range(n_atts) ]

                for k in range(n_objs):
                    rnodes[k].atts = set([ atts[i] for i,char in enumerate(it.next()) if char in ["x","X"] ])

                E += [ rnodes ]
                M += [ atts ]

            except StopIteration:
                break

        return PatternStructure( AbstractConceptGraph(E),M )

if __name__ == '__main__':

    pst = PCFReader().read("family_relations.pcf")

    lattice = build(pst)
    runner = LatticeRunner(lattice)
    state = 0

    while True:
        print('\n\n')
        runner.where()
        if(state == 0):
            print "Lower Neighbors:"
            runner.lower()
        elif(state == 1):
            print "Upper Neighbors:"
            runner.upper()
        choice = raw_input("What do you want to do now? ")

        if choice == 'u':
            state = 1

        elif choice == 'l':
            state = 0

        else:
            if(state == 0):
                runner.down(int(choice))
            elif(state == 1):
                runner.up(int(choice))

