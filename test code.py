import numpy as np
def find_fundamental_loops_undirected(graph):
    loops = dict()
    setLoops = dict() 

    def dfs(node, parent, path):
        if node in path:
            
            start_idx = path.index(node)
            loop_found = path[start_idx:]
            set_loop_found = set(loop_found)
            if set_loop_found not in setLoops.values():
                l = len(loops)
                loops[l+1] = loop_found
                setLoops[l+1] = set_loop_found
            return

        path.append(node)
        for neighbor in graph[node]:
            if neighbor != parent:
                dfs(neighbor, node, path.copy())
        path.pop()

    for node in graph:
        dfs(node, None, [])

    return loops,setLoops

# Example undirected graph represented as an adjacency list
graph = {
    1: (2,),
    2: (3, 1, 7),
    3: (4, 2),
    4: (5,3,9),
    5: (4,6,11),
    6: (5,),
    7: (8,9,2),
    8: (7,10),
    9: (7,4,11,10),
    10:(12,8,11,9),
    11: (9,10,12,5),
    12: (11,10)
}
lineC = {
    1:[2,1],
    2:(2,3),
    3:(3,4),
    4:(4,5),
    5:(5,6),
    6:(2,7),
    7:(7,8),
    8:(9,4),
    9:(9,10),
    10:(10,12),
    11:(5,11),
    12:(7,9),
    13:(9,11),
    14:(8,10),
    15:(10,11),
    16:(11,12)
}
loop_count = 0
fundamental_loops,fundamental_loops_check = find_fundamental_loops_undirected(graph)
# print(fundamental_loops)
for k,v in fundamental_loops.items():
    print(v)
    
setpop = set()
for v1 in fundamental_loops_check.values():
    for k,v2 in fundamental_loops_check.items():
        if v2 != v1:
            if v1.issubset(v2):
                setpop.add(k)
for k in setpop:
    fundamental_loops_check.pop(k)
print(fundamental_loops_check)
busC = {k:[] for k in graph.keys()}
for k,v in lineC.items():
    for bi in v:
        busC[bi].append(k)
a = {1,2,3,4}
a.discard(0)
print(a)
