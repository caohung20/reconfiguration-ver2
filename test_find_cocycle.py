import calendar
import numpy as np
import time
import pygad
import copy
from config import Configuration
from parameter import Parameter
from powerflow import PowerFlow
from graph import Graph   
    



class GA():
    def __init__(self):
        global param
        self.g = Graph(param,param.nBus)
        


    def fitness_func(self, reconfiguration, solution, solution_idx):
        global param
        config = Configuration(param=param, lineOff=solution)
        pf = PowerFlow(config=config)
        r1 = pf.run1Config_WithObjective()
        fitness = -r1['Objective']
        return fitness


    def crossover_func(self, parents, offspring_size, reconfiguration):
        global param
        
        
        parent1 = [6, 12, 24, 35, 34]
        parent2 =  [8, 6, 20, 36, 5]

        set_par2 = set(parent2)
        random_chosen_point = 4
        #this is exchange line
        exchg_branch = parent1[random_chosen_point]
        parent1.remove(exchg_branch)
        # use dfs to find loop
        print(exchg_branch,parent1) 
        loop_found = self.g.get_list_lines_in_loop(exchg_branch,parent1)
        
        set_loop_found = set(loop_found)
        print(set_loop_found)
        mutual = set_loop_found.intersection(set_par2)
        print(mutual)
        if len(mutual) == 1:
            exchg_branch2 = mutual.pop()
        else:
            exchg_branch2 = self.g.detect_cocycle(set_par2,exchg_branch,mutual)
        #crossover
        parent1.append(exchg_branch2)
        position = parent2.index(exchg_branch2)
        parent2[position] = exchg_branch
        return
        
        

        
if __name__ == '__main__':
    fi = 'tromvia.xlsx'
    param = Parameter(fi)
    start = time.time()
    ga = GA()
    ga.crossover_func("a","b","calendar")
    """
    graph = Graph(param, param.nBus)

    initial_population = graph.init_pop(20) 
    lineOff = {36, 15, 33, 35, 34}  
    lineC = copy.deepcopy(param.AllLine2Bus)
    bus2bus = copy.deepcopy(param.bus2bus)
    newlineOff = 15
    #
    for li in lineOff:
        v = lineC[li]
        bus2bus[v[0]].discard(v[1])
        bus2bus[v[1]].discard(v[0])
    while True:
    
        count = 0
        for pi in lineC[newlineOff]:
            check = graph.isReachable(1, pi, bus2bus)
            if not check:
                break
            count += 1
            if count == 2:
                print("da print check:",check)
        break
    print("khong ket noi")
    """

    
