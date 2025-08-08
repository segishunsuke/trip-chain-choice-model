# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 17:29:51 2024

@author: shun
"""

import trip_chain_simulator_nested_logit

simulator = trip_chain_simulator_nested_logit.Trip_chain_simulator()

simulator.read_input_data("input settings.csv", "input od cost.csv", "input trip chain data.csv", "trip chain", "input initial parameter values nested logit.csv")

# Estimate the parameters by MLE.
# The argument of this method is the name of the CSV file in which to output the estimation results.
# After execution, the parameters in the object are updated to the estimated values.
simulator.optimize_parameters_by_mle("output estimation result nested logit.csv")

# You can also use the Poisson pseudo-likelihood to estimate the parameters.
#simulator.optimize_parameters(1, 13, 130, 10, "output estimation result nesteed logit f.csv")

simulator.print_fitness_of_multiple_cases(1, 100, "output fitness nested logit.csv")

#simulator.set_random_seed(0)

simulator.print_statistics(100, "output summary trip chain nested logit.csv", "output summary od freq nested logit.csv", "output summary visits nested logit.csv", order_insensitive = False, count_unobserved = False)
