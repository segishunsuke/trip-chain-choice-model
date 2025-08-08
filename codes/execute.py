# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 17:29:51 2024

@author: shun
"""

import trip_chain_simulator

# Create an object for parameter estimation and prediction.
simulator = trip_chain_simulator.Trip_chain_simulator()

# Read the configuration and data files.
# The arguments of this method, in order, are: the configuration file, the OD travel cost data file,
# the observed trip chain data file, the column name that contains the trip chain in the data file,
# and the parameter configuration file.
simulator.read_input_data("input settings.csv", "input od cost.csv", "input trip chain data.csv", "trip chain", "input initial parameter values.csv")

# Estimate the parameters.
# The arguments of this method, in order, are: the number of simulations used to evaluate the objective function,
# the number of individuals per GA generation, the total number of function evaluations in the GA,
# the number of random number streams used, and the CSV file to output the estimation results.
# After execution, the parameters in the object are updated to the estimated values.
simulator.optimize_parameters(1, 13, 130, 10, "output estimation result.csv")

# Outputs the results of evaluating the objective function (goodness of fit).
# The arguments of this method, in order, are: the number of simulations to evaluate the objective function,
# the number of random number streams, and the output CSV file.
simulator.print_fitness_of_multiple_cases(1, 100, "output fitness.csv")

# Use this method to specify the seed for the random number generator.
# The argument should be a long integer.
#simulator.set_random_seed(0)

# Performs simulations to compute various predicted values under the current parameter values and outputs the results.
# The arguments of this method, in order, are: the number of simulations to compute predicted values,
# the CSV file for predicted trip chain frequencies, the CSV file for predicted OD flows,
# the CSV file for predicted number of visitors per place, whether to distinguish trip chains by place order,
# and whether to include unobserved trip chains in the output.
simulator.print_statistics(100, "output summary trip chain.csv", "output summary od freq.csv", "output summary visits.csv", order_insensitive = False, count_unobserved = False)

