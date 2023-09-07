import pygraphviz


class FormalContext(object):

    def __init__(self,M):
        self.__objects = []
        self.__intent = {}

        self.__attributes = M
        self.__extent = { m:set() for m in M }

    def _repr_html_(self):
        header_cells = [ "<td>{0}</td>".format(m) for m in self.__attributes ]
        rows = [ "<tr><td></td>{0}</tr>".format("".join(header_cells)) ]

        for g in self.__objects:
            row_cells = [ "<td>{0}</td>".format("x" if m in self.__intent[g] else "") for m in self.__attributes ]
            rows += [ "<tr><td>{0}</td>{1}</tr>".format(g,"".join(row_cells)) ]

        table = "<table style='display:inline-block'>{0}</table>".format("".join(rows))
        return table

    def objects(self):
        return list(self.__objects)

    def attributes(self):
        return list(self.__attributes)

    def add_object(self,objname,intent):

        self.__objects.append(objname)
        self.__intent[objname] = set(intent)

        for m in self.__attributes:
            if m in intent:
                self.__extent[m].add(objname)
            else:
                self.__extent[m].discard(objname)

    def extent(self,atts):

        if not atts:
            return set(self.__objects)

        return set.intersection(*[self.__extent[m] for m in atts])

    def intent(self,objs):

        if not objs:
            return set(self.__attributes)

        return set.intersection(*[self.__intent[g] for g in objs])


class PowerContextFamily(object):

    def __init__(self):
        self.K = []
        self.M = []

    def __getitem__(self,key):
        return self.K[key]

    # def __setitem__(self,key,value):
    #     self.K[key] = value

    def add_context(self,ctx):
        self.K += [ ctx ]
        self.M += [ ctx.attributes() ]

    def schema(self):
        return AttributeSystem(self.M)

    def _repr_html_(self):
        return "&nbsp;&nbsp;".join(ctx._repr_html_() for ctx in self.K )

class PCFFile(object):

    class ParseError(Exception):
        pass

    def __init__(self,filename=None):

        if filename:
            with open(filename) as pcf_file:
                self.__lines = filter(None, [ line.strip() for line in pcf_file.readlines() ])
        else:
            self.__lines = ["B-PCF","0","0"]

        self.__it = iter(self.__lines)


    def __parse_context(self,relational=False):

        n_objs = int(self.__it.next())
        n_atts = int(self.__it.next())

        if relational:
            objs = [ tuple(self.__it.next().split(",")) for g in range(n_objs) ]
        else:
            objs = [ self.__it.next() for g in range(n_objs) ]
        atts = [ self.__it.next() for m in range(n_atts) ]

        ctx = FormalContext(atts)
        for k in range(n_objs):
            obj_intent = set([ atts[i] for i,char in enumerate(self.__it.next()) if char in ["x","X"] ])
            ctx.add_object(objs[k],obj_intent)

        return ctx


    def parse(self):

        if not self.__it.next() == "B-PCF":
            raise self.ParseError("Wrong File Format")

        pcf = PowerContextFamily()
        pcf.add_context(self.__parse_context())

        while True:
            try:
                pcf.add_context(self.__parse_context(True))
            except StopIteration:
                break

        return pcf


class AttributeSystem(object):

    def __init__(self,M):
        self.__M = [ set(S) for S in M ]
        while len(self.__M) <= 2:
            self.__M += [ set() ]

    def __len__(self):
        return len(self.__M)

    def __getitem__(self,key):
        return self.__M[key]

    def __setitem__(self,key,value):
        self.__M[key] = value

    def clone(self):
        return AttributeSystem([ set(S) for S in self.__M ])

    def union(self,other):
        n = max(len(self),len(other))
        M = []
        for i in range(n):
            if i >= len(self):
                M += [ set(other[i]) ]
            elif i >= len(other):
                M += [ set(self[i]) ]
            else:
                M += [ (self[i] | other[i]) ]
        return AttributeSystem(M)


class Node(object):

    def __init__(self,node_id,name="",atts=set()):
        self.id = node_id
        self.name = name
        self.atts = atts
        self.incidences = {}
        self.distance = None
        self.target = None

    def __str__(self):
        return str(self.id)

    def __repr__(self):
        return str(self.id)

    def inf(self,other):
        node_id = "{0}:{1}".format(self.id,other.id).strip(":") or ":"
        node = Node(node_id)
        node.atts = self.atts & other.atts
        return node

    def sup(self,other):
        node_id = "{0}:{1}".format(self.id,other.id).strip(":") or ":"
        node = Node(node_id)
        node.atts = self.atts | other.atts
        return node

class RNode(object):

    def __init__(self,nodes,atts):
        n = len(nodes)
        for i,node in enumerate(nodes):
            node.incidences.setdefault((i,n),[]).append(self)
        self.atts = set(atts)
        self.nodes = nodes
        self.origin = None
        self.next = None
        self.prev = None

    def __str__(self):
        return "({0})".format(",".join(self.nodes))

    def __repr__(self):
        return "({0})".format(",".join([ node.__repr__() for node in self.nodes]))

    def id(self):
        return tuple(v.id for v in self.nodes)


class Table(object):

    def __init__(self,header):
        self.header = header
        self.data = []

    def append(self,el):
        if isinstance(el,Table):
            self.data += el.data
        else:
            self.data += [el]

    def _repr_html_(self):
        header_cells = [ "<td bgcolor=gray>{0}</td>".format(c) for c in self.header ]
        rows = [ "<tr>{0}</tr>".format("".join(header_cells)) ]

        for row in self.data:
            row_cells = [ "<td>{0}</td>".format(row[c]) for c in self.header ]
            rows += [ "<tr>{0}</tr>".format("".join(row_cells)) ]

        table = "<table style='display:inline-block'>{0}</table>".format("".join(rows))
        return table


class IntensionGraph(object):

    def __init__(self,schema):
        self.E = [ [] for i in range(len(schema)) ]
        self.schema = schema.clone()
        self.window = {}


    def __getitem__(self,key):
        return self.E[key]


    def __setitem__(self,key,value):
        self.E[key] = value


    def __len__(self):
        return len(self.schema)


    def __str__(self):
        return str(self.E)


    def _repr_svg_(self):
        G = pygraphviz.AGraph(directed=True)

        if len(self.schema) <= 3:
            for v in self.E[0]:
                G.add_node(v.id,xlabel=",".join(v.atts))
            for e in self.E[1]:
                G.add_node(e.id(),label=",".join(e.atts),shape="rect")
                G.add_edge(e.id(),e.nodes[0].id)
            for e in self.E[2]:
                G.add_edge(e.id(),label=",".join(e.atts))

        else:
            for v in self.E[0]:
                G.add_node(v.id,xlabel=",".join(v.atts),fillcolor="white",style="filled")
            for k in range(1,len(self.schema)):
                for e in self.E[k]:
                    G.add_node(e.id(),label=",".join(e.atts),fillcolor="white",style="filled",shape="rect")
                    for i,v in enumerate(e.nodes):
                        G.add_edge(e.id(),v.id,xlabel=str(i+1))

        for v in self.window.values():
            n = G.get_node(v.id)
            n.attr["style"]="filled"
            n.attr["fillcolor"]="yellow"

        return G.draw(path=None,format='svg',prog='dot')


    def __add__(self,other):
        return self.sum(other)


    def __mul__(self,other):
        return self.product(other)


    def clone(self):
        C = IntensionGraph(self.schema)
        C[0] = [ Node(n.id,n.name,n.atts) for n in self[0] ]

        for k in range(1,len(self.schema)):
            C[k] = [ RNode([C.nodesByID(n.id) for n in rnode.nodes],rnode.atts) for rnode in self[k] ]

        C.window = { name:C.nodesByID(n.id) for name,n in self.window.iteritems() }

        return C

    def components(self):
        return self.product(IntensionGraph.terminal(self.schema),components=True)


    def nodesByID(self,node_id):
        for v in self[0]:
            if v.id == node_id:
                return v

        return None


    @staticmethod
    def initial(window_keys=[]):
        I = IntensionGraph(AttributeSystem([set()]))
        I.window = {}

        for k in window_keys:
            n = Node(k)
            I[0] += [n]
            I.window[k] = n

        return I


    @staticmethod
    def terminal(schema,window_keys=[]):
        T = IntensionGraph(schema)

        n0 = Node(":","")
        n0.atts = set(schema[0])
        T[0] += [n0]

        for i in range(1,len(schema)):
            if schema[i]:
                T[i] += [ RNode([n0]*i,set(schema[i])) ]

        T.window = { k:n0 for k in window_keys }
        return T


    @staticmethod
    def star(atts,*names):
        M = [ set() ] + [ set() for name in names ]
        M[len(names)] = set(atts)

        S = IntensionGraph(AttributeSystem(M))
        S[0] = [ Node(name,name) for name in names ]
        S[len(names)] = [ RNode(S[0],atts) ]
        S.window = { name:S[0][i] for i,name in enumerate(names) }

        return S


    @staticmethod
    def node(atts,name):
        N = IntensionGraph(AttributeSystem([ set(atts) ]))
        n0 = Node(name,name)
        n0.atts = set(atts)

        N[0] = [n0]
        N.window = {name:n0}

        return N


    @staticmethod
    def ig(pcf):

        ig = IntensionGraph(pcf.schema())

        nodes = []
        for g in pcf.K[0].objects():
            v = Node(g,g)
            v.atts = pcf.K[0].intent([g])
            nodes.append(v)
        ig[0] += nodes

        nodeByID = { v.id:v for v in nodes }

        for k in range(1,len(pcf.K)):
            rnodes = []
            for t in pcf.K[k].objects():
                rnodes += [ RNode([ nodeByID[g] for g in t ], pcf.K[k].intent([t])) ]
            ig[k] += rnodes

        return ig


    @staticmethod
    def objects(pcf):

        D = IntensionGraph.ig(pcf)
        complist = D.components()
        objects = []
        for C in complist:
            objects += [ C.windowed(n.id) for n in C[0] ]

        return objects


    def pcf(self):

        ctx = FormalContext(self.schema[0])
        for v in self.E[0]:
            ctx.add_object(v.id,v.atts)

        pcf = PowerContextFamily()
        pcf.add_context(ctx)

        for k in range(1,len(self.E)):
            rctx = FormalContext(self.schema[k])
            for e in self.E[k]:
                rctx.add_object(e.id(),e.atts)
            pcf.add_context(rctx)

        return pcf


    def rename(self,f):
        R = self.clone()
        for v in R[0]:
            v.id = f(v.id)

        return R


    def windowed(self,*args,**kwargs):
        W = self.clone()
        W.window = {}
        for i,nid in enumerate(args):
            W.window[str(i)] = self.nodesByID(str(nid))
        for name,nid in kwargs.iteritems():
            W.window[name] = self.nodesByID(str(nid))

        return W


    def quotient(self,relation):

        cls = { node.id:set([node.id]) for node in self[0] }

        for i,j in relation:
            union = cls[i] | cls[j]
            for k in union:
                cls[k] = union

        Q = IntensionGraph(self.schema)
        ids = cls.keys()

        for nid in ids:
            nclass = cls[nid]
            if not isinstance(nclass,Node):
                node = Node(nid)
                node.atts = set.union(*[ self.nodesByID(i).atts for i in nclass ])
                for i in nclass:
                    cls[i] = node
                Q[0] += [node]

        for k in range(1,len(self.schema)):
            edgesByID = {}
            for e in self[k]:
                t = tuple([ cls[nid].id for nid in e.id() ])
                if t in edgesByID:
                    edgesByID[t].atts |= e.atts
                else:
                    rnode = RNode([ cls[nid] for nid in e.id() ],e.atts)
                    edgesByID[t] = rnode
                    Q[k] += [rnode]

        Q.window = { name:cls[node.id] for name,node in self.window.iteritems() }
        return Q


    def sum(self,other):

        S1 = self.rename(lambda s:"0:{0}".format(s))
        S2 = other.rename(lambda s:"1:{0}".format(s))

        S = IntensionGraph(S1.schema.union(S2.schema))
        n = max(len(S1.schema),len(S2.schema))
        for k in range(n):
            if k >= len(S1.schema):
                S[k] = list(S2[k])
            elif k >= len(S2.schema):
                S[k] = list(S1[k])
            else:
                S[k] = S1[k] + S2[k]
        S.window = S1.window.copy()
        S.window.update(S2.window)

        common_keys = set(S1.window.keys()) & set(S2.window.keys())
        relation = [ (S1.window[k].id,S2.window[k].id) for k in common_keys ]

        S = S.quotient(relation)

        proj = lambda s:s.split(":",1)[1]
        proj_ids = set( proj(v.id) for v in S[0] )

        if len(proj_ids) == len(S[0]):
            S = S.rename(proj)

        return S


    def product(self,other,components=False,matrix=None):
        common_keys = set(self.window.keys()) & set(other.window.keys())

        if not common_keys:
            matrix = {}

            # hack to ensure self != other (because windows must be set independently on both operands).
            # this is error-prone because cloning relies on multiplication.
            save1 = self.window
            save2 = other.window
            if self == other:
                other = other.clone()

            complist = []
            for nd1 in self[0]:
                for nd2 in other[0]:
                    if not (nd1.id,nd2.id) in matrix:
                        self.window = {"0":nd1}
                        other.window = {"0":nd2}
                        prod = self.product(other,matrix=matrix)
                        prod.window = {}
                        complist.append(prod)

            self.window = save1
            other.window = save2

            if components:
                return complist
            else:
                G = IntensionGraph(self.schema.union(other.schema))
                n = max(len(self.schema),len(other.schema))
                for C in complist:
                    for k in range(n):
                        G[k] += list(C[k])

                return G


        else:
            if matrix is None:
                matrix = {}

            key = iter(common_keys).next()

            node1 = self.window[key]
            node2 = other.window[key]

            node0 = node1.inf(node2)
            matrix[(node1.id,node2.id)] = node0
            queue = [(node1,node2)]

            P = IntensionGraph(self.schema.union(other.schema))
            P[0] += [node0]

            while queue:
                node1,node2 = queue.pop(0)

                node = matrix[(node1.id,node2.id)]
                for n in range(1,min(len(self.schema),len(other.schema))):
                    for k in range(n):
                        for r1 in node1.incidences.get((k,n),[]):
                            for r2 in node2.incidences.get((k,n),[]):
                                if r1.atts & r2.atts:
                                    star = []
                                    for j in range(n):
                                        if j != k:
                                            n1 = r1.nodes[j]
                                            n2 = r2.nodes[j]
                                            pos = (n1.id,n2.id)
                                            if pos in matrix:
                                                neighbor = matrix[pos]
                                            else:
                                                neighbor = n1.inf(n2)
                                                P[0].append(neighbor)
                                                matrix[pos] = neighbor
                                                queue.append((n1,n2))
                                        else:
                                            neighbor = node
                                        star += [neighbor]
                                    if k == 0:
                                        P[n] += [ RNode(star,r1.atts & r2.atts) ]

            P.window = { i:matrix.get((self.window[i].id,other.window[i].id)) for i in common_keys }
            if None not in P.window.values():
                return P
            else:
                return IntensionGraph.initial(common_keys)


    def leq(self,other):
        return self.morphisms(other,False)


    def morphisms(self,other,all=True):

        G1 = self.clone()
        G2 = other.clone()

        header = [ node.id for node in G1[0] ]

        if not G1.window:
            if all:
                table = Table(header)
                if len(G1[0])==0:
                    table.append({})
                    return table
                else:
                    n1 = G1[0][0]
                    if len(G2[0])==0:
                        return table
                    else:
                        G1.window = { "0": n1 }
                        for n2 in G2[0]:
                            G2.window = { "0": n2 }
                            table.append(G1.morphisms(G2,True))
                        return table

            else:
                if len(G1[0])==0:
                    return True
                else:
                    n1 = G1[0][0]
                    if len(G2[0])==0:
                        return False
                    else:
                        G1.window = { "0": n1 }
                        for n2 in G2[0]:
                            G2.window = { "0": n2 }
                            if G1.morphisms(G2,False):
                                return True
                        return False

        else:
            name,n1 = G1.window.iteritems().next()
            if name not in G2.window:
                return Table(header) if all else False
            n2 = G2.window[name]

            if not n1.atts <= n2.atts:
                return Table(header) if all else False
            r1 = G1.create_neighborhoods(n1)
            if r1 is None:
                if all:
                    table = Table(header)
                    table.append({n1.id:n2.id})
                    return table
                else:
                    return True
            n1.target = n2
            # print "matching " + n1.__str__() + " with " + n2.__str__()
            if all:
                table = Table(header)
                G1.find_match(r1,table)
                return table
            else:
                if G1.find_match(r1):
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
            for n in range(len(self.schema)):
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


    def find_match(self,r1,table=None):
        if r1.origin is None:
            print "WARNING: TRIPLE NOT CLASSIFIED"
            return False

        for r2 in r1.nodes[r1.origin[0]].target.incidences.get(r1.origin,[]):
            if self.matches(r1,r2):

                undo = [ r1.nodes[i].target for i in range(len(r1.nodes)) ]
                for k in range(len(r1.nodes)):
                    r1.nodes[k].target = r2.nodes[k]

                if table:
                    if r1.next is None:
                        table.append({node.id:node.target.id for node in self[0]})
                    else:
                        self.find_match(r1.next,table)

                else:
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


    def minimize(self):
        table = self.morphisms(self)
        min_size = len(table.header)
        fmin = table.data[0]

        for endo in table.data:
            m = len(set(endo.itervalues()))
            if m < min_size:
                min_size = m
                fmin = endo

        core_ids = fmin.values()

        C = self.clone()
        C[0] = [ n for n in C[0] if n.id in core_ids ]
        for k in range(1,len(C.schema)):
            C[k] = [ e for e in C[k] if all(n in C[0] for n in e.nodes) ]

        return C


star = IntensionGraph.star
node = IntensionGraph.node
terminal = IntensionGraph.terminal
initial = IntensionGraph.initial
ig = IntensionGraph.ig


class PatternStructure(object):

    def __init__(self,pattern_class,pcf):
        self.cls = pattern_class
        self.pcf = pcf


    def build(self):

        unit = self.cls.terminal(self.pcf.schema(),["0"])

        stack,lattice = [],[]
        unit_extent,unit_additions = [],[]

        generators = [ (pattern.window["0"].id,pattern.minimize()) for pattern in self.cls.objects(self.pcf) ]

        for g,pg in generators:
            if unit.leq(pg):
                unit_extent.append(g)
            else:
                unit_additions.append((g,pg))

        c_zero = PatternConcept(unit_extent,unit)
        stack.append((c_zero,unit_additions))
        lattice.append(c_zero)


        while stack:
            #print("NEXT STEP")
            concept,additions = stack[-1]
            #print("current extent: {0}".format(",".join(concept.extent)))
            #print("     additions: {0}".format(",".join([ p[0] for p in additions ])))
            if not additions:
                stack.pop()
                continue

            g0,pg0 = additions.pop()
            #print("extent + {0}".format(g0))
            next_intent = (concept.intent * pg0).minimize()

            next_extent,next_additions = [],[]
            j,approved = 0,False

            for g,pg in generators:

                if j<len(concept.extent) and concept.extent[j] == g:
                    next_extent.append(g)
                    j += 1
                elif g == g0:
                    next_extent.append(g)
                    approved = True
                else:
                    if approved:
                        if next_intent.leq(pg):
                            next_extent.append(g)
                        else:
                            next_additions.append((g,pg))
                    else:
                        if next_intent.leq(pg):
                            #print("rejected! {0} -> {1}".format(g0,g))
                            break
                        else:
                            continue

            if approved:
                #print("next concept: extent + {0}".format(g0))
                next_concept = PatternConcept(next_extent,next_intent)
                stack.append((next_concept,next_additions))

                for c in lattice:
                    if c.lessequal(next_concept):
                        if not any(cu.lessequal(next_concept) for cu in c.upper):
                            c.upper.append(next_concept)
                            next_concept.lower.append(c)

                lattice.append(next_concept)

        return lattice


# class PatternStructure(object):
#
#    def __init__(self,delta):
#        self.delta = delta
#        self.__generators = None
#
#    def generators(self):
#
#        if self.__generators:
#            return self.__generators
#
#        generators = []
#        components = self.components(self.delta)
#        for comp in components:
#            for n in comp.E[0]:
#                generators.append( (n.id,self.clone((n,comp))) )
#
#        self.__generators = generators
#        return generators
#
#
## depends on: PatternStructure.{terminal,generators,product,morphisms}
## depends on: PatternConcept.{__init__,intent,extent,upper,lower,lessequal}
#    def build(pst):
#
#        unit = pst.terminal()
#
#        stack,lattice = [],[]
#        unit_extent,unit_additions = [],[]
#
#        for g,p in pst.generators():
#            if pst.morphisms(unit,p):
#                unit_extent.append(g)
#            else:
#                unit_additions.append((g,p))
#
#        c_zero = PatternConcept(unit_extent,unit)
#        stack.append((c_zero,unit_additions))
#        lattice.append(c_zero)
#
#
#        while stack:
#            # print("NEXT STEP")
#            concept,additions = stack[-1]
#            # print("current extent: {0}".format(",".join(concept.extent)))
#            # print("     additions: {0}".format(",".join([ p[0] for p in additions ])))
#            if not additions:
#                stack.pop()
#                continue
#
#            g0,p0 = additions.pop()
#            # print("extent + {0}".format(g0))
#            next_intent = pst.product(concept.intent,p0)
#
#            next_extent,next_additions = [],[]
#            j,approved = 0,False
#
#            for g,p in pst.generators():
#
#                if j<len(concept.extent) and concept.extent[j] == g:
#                    next_extent.append(g)
#                    j += 1
#                elif g == g0:
#                    next_extent.append(g)
#                    approved = True
#                else:
#                    if approved:
#                        if pst.morphisms(next_intent,p):
#                            next_extent.append(g)
#                        else:
#                            next_additions.append((g,p))
#                    else:
#                        if pst.morphisms(next_intent,p):
#                            # print("rejected! {0} -> {1}".format(g0,g))
#                            break
#                        else:
#                            continue
#
#            if approved:
#                # print("next concept: extent + {0}".format(g0))
#                next_concept = PatternConcept(next_extent,next_intent)
#                stack.append((next_concept,next_additions))
#
#                for c in lattice:
#                    if c.lessequal(next_concept):
#                        if not any(cu.lessequal(next_concept) for cu in c.upper):
#                            c.upper.append(next_concept)
#                            next_concept.lower.append(c)
#
#                lattice.append(next_concept)
#
#        return lattice


class PatternConcept(object):
    def __init__(self,extent,intent):
        # print("extent={0}".format(extent))
        self.extent = extent
        self.intent = intent
        self.upper = []
        self.lower = []

    def extent(self):
        return self.extent

    def intent(self):
        return self.intent

    def lessequal(self,concept):
        return set(self.extent) <= set(concept.extent)


class LatticeWalker(object):

    def __init__(self,pcffile):
        pst = PCFReader().read(pcffile)
        self.lattice = build(pst)
        self.concept = self.lattice[-1]
        self.direction = 0

    def turn(self):
        self.direction = 1 - self.direction

    def go(self,edge):
        if(self.direction == 0):
            self.concept = self.concept.lower[edge]
            self.where()
            self.lower()
        elif(self.direction == 1):
            self.concept = self.concept.upper[edge]
            self.where()
            self.upper()

    def where(self):
        return self.concept

    def lower(self):
        print("Lower Neighbors:")
        for i,c in enumerate(self.concept.lower):
            print i, ") ",c

    def upper(self):
        print("Upper Neighbors:")
        for i,c in enumerate(self.location.upper):
            print i, ") ",c

