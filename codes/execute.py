# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 17:29:51 2024

@author: shun
"""

import trip_chain_simulator

# Create an object for estimation and prediction.
simulator = trip_chain_simulator.Trip_chain_simulator()

# Read the configuration and data files.
# The method's arguments are, in order, the configuration file name, the data file name for OD travel costs, the data file name for the observed trip chain, the name of the column in which the trip chain is described in the data file, and the parameter configuration file name.
simulator.read_input_data("input settings.csv", "input od cost.csv", "input trip chain data.csv", "trip chain", "input initial parameter values.csv")

# Estimate the parameters.
# The method's arguments are, in order, the number of simulations used to evaluate the objective function, the number of individuals per generation of the GA, the number of times the objective function is evaluated in the GA, the number of random number series used, and the CSV file name to output the optimization results.
# After completion, the parameters in the object are updated to the estimated values.
simulator.optimize_parameters(1, 12, 120, 10, "output estimation result.csv")

# Outputs the results of evaluating the objective function (goodness of fit).
# The method's arguments are, in order, the number of simulations used to evaluate the objective function, the number of random number series used, and the CSV file to output the results.
simulator.print_fitness_of_multiple_cases(1, 100, "output fitness.csv")

# Use this method to specify the random number sequence to use.
# The argument is a long int type number.
#simulator.set_random_seed(0)

# Performs simulations to find various predicted values under the current parameter values and outputs the results.
# The method's arguments are, in order, the number of simulations to use in calculating the predicted values, the CSV file name to output the predicted values of the trip chain frequencies, the CSV file name to output the predicted values of the OD traffic volume, the CSV file name to output the predicted values of the number of visitors for individual places, whether to distinguish trip chains by the order in which places are visited, and whether to record the predicted values of unobserved trip chains.
simulator.print_statistics(100, "output summary trip chain.csv", "output summary od freq.csv", "output summary visits.csv", order_insensitive = False, count_unobserved = False)
