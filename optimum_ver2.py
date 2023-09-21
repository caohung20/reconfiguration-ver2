import numpy as np
import time
import pygad
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
        offspring = []
        
        idx = 0
        while len(offspring) != offspring_size[0]:
            
            parent1 = list(parents[idx % parents.shape[0], :])
            parent2 = parents[(idx + 1) % parents.shape[0], :]
            set_par2 = set(parent2)
            random_chosen_point = np.random.choice(range(offspring_size[1]))
            #this is exchange line
            exchg_branch = parent1[random_chosen_point]
            parent1.remove(exchg_branch)
            # use dfs to find loop 
            loop_found = self.g.get_list_lines_in_loop(exchg_branch,parent1)
            
            set_loop_found = set(loop_found)
            mutual = set_loop_found.intersection(set_par2)
            exchg_branch2 = mutual.pop()
            parent1.append(exchg_branch2) 

            offspring.append(parent1)

            idx += 1

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
    print(initial_population)
    reconfiguration = pygad.GA(num_generations=200, num_parents_mating=10,
                        fitness_func=ga.fitness_func,
                        initial_population=initial_population,gene_type=int,
                        crossover_type=ga.crossover_func,
                        mutation_type=ga.mutation_func,mutation_num_genes=2,
                        crossover_probability=0.9,
                        mutation_probability=0.005)
                        #parallel_processing=4 ,mutation_num_genes=10)
    reconfiguration.run()

    # Retrieve the best solution and its fitness value
    best_solution = reconfiguration.best_solution()

    # Print the best solution and its fitness value
    print("Best Solution: ", best_solution)
    end = time.time()
    print(end-start)
    reconfiguration.summary()
    
    
    
