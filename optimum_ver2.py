import numpy as np
import time
import pygad
from config import Configuration
from parameter import Parameter
from powerflow import PowerFlow
from graph import Graph
from new_kernel import Configuration as CF

    



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
        global no_ex_gene
        list_parents = parents.tolist()
        offspring = []
        idx = 0
        while idx != offspring_size[0]:
            parent1 = list_parents[idx % parents.shape[0]]
            parent2 = list_parents[(idx + 1) % parents.shape[0]]
            set_par2 = set(parent2)
            random_chosen_points = np.random.randint(0, offspring_size[1]-1, no_ex_gene)
            for random_chosen_point in random_chosen_points:
                #this is exchange line
                exchg_branch = parent1[random_chosen_point]
                parent1.remove(exchg_branch)
                # use dfs to find loop 
                loop_found = self.g.get_list_lines_in_loop(exchg_branch,parent1)
                set_loop_found = set(loop_found)
                mutual = set_loop_found.intersection(set_par2)
                if len(mutual) == 1:
                    exchg_branch2 = mutual.pop()
                else:
                    exchg_branch2 = self.g.detect_cocycle(set_par2,exchg_branch,mutual)
                #crossover
                parent1.append(exchg_branch2)
                position = parent2.index(exchg_branch2)
                parent2[position] = exchg_branch
            idx += 1
            offspring.append(parent1)

        return np.array(offspring)

    def mutation_func(self,offspring,reconfiguration):
        global param
        for chromosome_idx in range(offspring.shape[0]):
            random_gene_idx = np.random.choice(range(offspring.shape[1]))
            exchg_branch = offspring[chromosome_idx,random_gene_idx]
            # create a copy version of chromo
            new_config = list(offspring[chromosome_idx,:])
            new_config.remove(exchg_branch)
            # use dfs to find loop 
            loop_found = self.g.get_list_lines_in_loop(exchg_branch,new_config)
            random_line_in_loop_idx = np.random.choice(range(len(loop_found)))
            exchg_branch2 = loop_found[random_line_in_loop_idx]
            offspring[chromosome_idx,random_gene_idx] = exchg_branch2
        return offspring


if __name__ == '__main__':
    fi = 'tromvia.xlsx'
    param = Parameter(fi)
    start = time.time()
    ga = GA()
    
    graph = Graph(param, param.nBus)

    initial_population = graph.init_pop(20)

    no_ex_gene = 2

    reconfiguration = pygad.GA(num_generations=500, num_parents_mating=10,
                        fitness_func=ga.fitness_func,
                        initial_population=initial_population,gene_type=int,
                        crossover_type=ga.crossover_func,
                        mutation_type=ga.mutation_func,mutation_num_genes=1,
                        crossover_probability=0.8,
                        mutation_probability=0.005,allow_duplicate_genes=False)
    
    reconfiguration.run()

    # Retrieve the best solution and its fitness value
    best_solution = reconfiguration.best_solution()

    # Print the best solution and its fitness value
    print("Best Solution: ", best_solution)
    end = time.time()
    print(end-start)
    reconfiguration.summary()
    #
    reconfiguration.plot_fitness()
    #

    
    