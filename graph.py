import copy
import numpy as np
class Graph():

    def __init__(self, param, vertices):
        self.V = vertices
        self.graph = param.AllLine2Bus
        self.param = param
    def get_list_lines_in_loop(self,lineOn,lineOff):
        self.lineOff = lineOff
        self.lineC = copy.deepcopy(self.param.AllLine2Bus)
        #
        self.busC = copy.deepcopy(self.param.BUSC)

        for li in self.lineOff:
            for bi in self.lineC[li]:
                self.busC[bi].remove(li)
        #
        self.bus2bus = copy.deepcopy(self.param.bus2bus)
        for li in self.lineOff:
            bi = self.lineC[li]
            self.bus2bus[bi[0]].remove(bi[1])
            self.bus2bus[bi[1]].remove(bi[0])
        #
        for li in self.lineOff:
            self.lineC.pop(li)
        node = self.lineC[lineOn][0]
        loop = self.getloop(node)
        self.line_in_loop = []
        for i in range(len(loop[1])):
            line = {loop[1][i-1],loop[1][i]}
            position = self.param.busesin1line.index(line)
            self.line_in_loop.append(self.param.listline[position])
        
        return self.line_in_loop
         
    def dfs(self, node, parent, path):
        if node in path:
            start_idx = path.index(node)
            loop_found = path[start_idx:]
            set_loop_found = set (loop_found)
            if set_loop_found not in self.setloops.values():
                l = len(self.loops)
                self.loops[l+1] = loop_found
                self.setloops[l+1] = set_loop_found
            return 
        path.append(node)
        for neighbor in self.bus2bus[node]:
            if neighbor != parent:
                self.dfs(neighbor, node, path.copy())
        path.pop()
    def getloop(self,node):
        self.loops = dict()
        self.setloops = dict()
        
        self.dfs(node, None, [])
        return self.loops
    
    
    


    # A utility function to find set of an element i
    # (truly uses path compression technique)
    def find(self, parent, i):
        if parent[i-1] != i:
            # Reassignment of node's parent
            # to root node as
            # path compression requires
            parent[i-1] = self.find(parent, parent[i-1])
        return parent[i-1]

    # A function that does union of two sets of x and y
    # (uses union by rank)
    def union(self, parent, rank, x, y):

        # Attach smaller rank tree under the root of
        # the high rank tree (Union by Rank)
        if rank[x-1] < rank[y-1]:
            parent[x-1] = y
        elif rank[x-1] > rank[y-1]:
            parent[y-1] = x

        # If ranks are the same, then make one as root
        # and increment its rank by one
        else:
            parent[y-1] = x
            rank[x-1] += 1

    # The main function to construct MST
    # using Kruskal's algorithm
    def KruskalMST(self):

        # This will store the resultant MST
        result = []

        # An index variable, used for sorted edges
        i = 0

        # An index variable, used for result[]
        e = 0

        parent = []
        rank = []

        # Create V subsets with single elements
        for node in self.param.BUS.keys():
            parent.append(node)
            rank.append(0)

        # Number of edges to be taken is less than V-1
        while e < self.V - 1:

            # Pick the edge and increment
            # the index for the next iteration
            u, v = self.graph[self.param.listline[i]]
            i = i + 1
            x = self.find(parent, u)
            y = self.find(parent, v)

            # If including this edge doesn't
            # cause a cycle, then include it in result
            # and increment the index of result
            # for the next edge
            if x != y:
                e = e + 1
                result.append({u, v})
                self.union(parent, rank, x, y)
            # Else discard the edge

        
        lineOn = set()
        for li in result:
            position = self.param.busesin1line.index(li)
            lineOn.add(self.param.listline[position]) 
        AllLine = set(self.param.AllLine2Bus.keys())
        lineOff = AllLine.difference(lineOn)
        
        return list(lineOff)
    
    def get_all_loops(self):
        fundamental_loops = dict()
        an_init_tree = self.KruskalMST()
        count = 1
        #iterate lineOFFs in initial tree
        for li in an_init_tree:
            #create a config with one fundamental loop
            new_config = an_init_tree.copy()
            new_config.remove(li)
            loop_found = self.get_list_lines_in_loop(li,new_config)
            fundamental_loops[count] = loop_found
            count += 1
        set_fundamental_loops = dict()
        for k,v in fundamental_loops.items():
            set_fundamental_loops[k] = set(v)
        return fundamental_loops, set_fundamental_loops
    
    def init_pop(self,no_chromo):
        an_init_tree = self.KruskalMST()
        leng_chromo = len(an_init_tree)
        first_pop = []
        first_pop.append(an_init_tree)
        for chromo_idx in range(no_chromo):
            random_gene_idx = np.random.choice(range(leng_chromo))
            exchg_branch = first_pop[chromo_idx][random_gene_idx]
            # create a copy version of chromo
            new_config = first_pop[chromo_idx].copy()
            new_config.remove(exchg_branch)
            # use dfs to find loop 
            loop_found = self.get_list_lines_in_loop(exchg_branch,new_config)
            random_line_in_loop_idx = np.random.choice(range(len(loop_found)))
            exchg_branch2 = loop_found[random_line_in_loop_idx]
            new_config.append(exchg_branch2)
            first_pop.append(new_config)
        return np.array(first_pop)
    
    def detect_cocycle(self,lineOff,newlineOff,mutual):
        lineC = copy.deepcopy(self.param.AllLine2Bus)
        busC = self.param.BUSC
        busC_copy = copy.deepcopy(busC)
        for li in lineOff:
            for bi in lineC[li]:
                busC_copy[bi].remove(li)
        for bi in lineC[newlineOff]:
            busC_copy[bi].discard(newlineOff)
        for k in busC_copy.keys():
            if len(busC_copy[k]) == 0:
                island_bus = k
                break
        p1_lineoff = lineC[newlineOff][0]
        p2_lineoff = lineC[newlineOff][1]
        if "island_bus" in locals():
            setfound = busC[island_bus].intersection(mutual)
            print(setfound)
            linefound = setfound.pop()
        else:
            # doesnt find any cocycle
            if newlineOff in lineOff:
                linefound = newlineOff
            else:
                for li in mutual:
                    for bi in lineC[li]:
                        if bi == p1_lineoff or bi == p2_lineoff:
                            linefound = li
                            break
        if "linefound" not in locals():
            linefound = mutual.pop()

        return linefound
    
    def is_connected(self,graph,start,end):
        self.graph = graph
        self.start = start
        self.end = end
        self.visited = set()
        return self.newdfs(self.start)
    
    def newdfs(self,node):
        self.visited.add(node)
        for neighbor in self.graph.get(node, []):
            if neighbor == self.end:
                return True
            if neighbor not in self.visited:
                if self.newdfs(neighbor):
                    return True
        return False