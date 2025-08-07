# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 17:29:51 2024

@author: shun
"""

import trip_chain_simulator

simulator = trip_chain_simulator.Trip_chain_simulator()

simulator.read_input_data("input settings.csv", "input od cost.csv", "input trip chain data.csv", "trip chain", "input initial parameter values.csv")

#simulator.optimize_parameters(1, 50, 500, 10, "output estimation result.csv")
#simulator.print_fitness_of_multiple_cases(1, 100, "output fitness.csv")

#simulator.print_statistics(100, "output summary trip chain.csv", "output summary od freq.csv", "output summary visits.csv", order_insensitive = False, count_unobserved = False)

#simulator.read_input_data("input settings.csv", "input od cost 2.csv", "input trip chain data.csv", "trip chain", "input initial parameter values.csv")
simulator.set_random_seed(0)
simulator.print_statistics(100, "output summary trip chain 2a.csv", "output summary od freq 2a.csv", "output summary visits 2a.csv", order_insensitive = True, count_unobserved = False)
#simulator.print_statistics(100, "output summary trip chain 2b.csv", "output summary od freq 2b.csv", "output summary visits 2b.csv", order_insensitive = True, count_unobserved = False)