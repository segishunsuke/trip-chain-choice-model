# -*- coding: utf-8 -*-

ctypedef struct Trip_chain_data:
    int n_ad_pairs  # In this code, n_ad_pairs = n_ports * n_ports, where ports represent entry and exit points.
    int *n_trip_chains  # The number of trip chain patterns for each arrival-departure pair
    int **len_trip_chains  # Length of each stored trip chain pattern
    int ***trip_chains     # The trip chain patterns themselves
    int **freq_observation  # Observed frequencies of each trip chain pattern
    int **freq_simulation   # Simulated frequencies of each trip chain pattern
    int freq_simulation_not_observed  # Frequency count of simulated chains not found in the observed set

cdef Trip_chain_data* init_trip_chain_data(int n_ad_pairs) noexcept nogil
cdef void append_trip_chain_data(Trip_chain_data *this, int id_ad_pair, int *trip_chain, int len_trip_chain, int freq_observation) noexcept nogil
cdef int search_trip_chain_data(Trip_chain_data *this, int id_ad_pair, int *trip_chain, int len_trip_chain) noexcept nogil
cdef void sort_trip_chain_data(Trip_chain_data *this) noexcept nogil
cdef void sort_trip_chain_data_observed_frequency(Trip_chain_data *this) noexcept nogil
cdef void delete_trip_chain_data(Trip_chain_data *this) noexcept nogil

cdef class Trip_chain_simulator:
    cdef int n_places
    cdef int n_ports
    cdef int n_all_places
    cdef int n_sample
    cdef int n_iteration
    cdef unsigned long random_seed
    cdef double poisson_shift

    cdef Trip_chain_data *trip_chain_data

    cdef int[::1] arrival_ports
    cdef int[::1] depart_ports
    
    cdef int[:,::1] od_freq_observation
    cdef int[:,::1] od_freq_simulation
    cdef int[:,::1] visits_observation
    cdef int[:,::1] visits_simulation
    
    cdef Trip_chain_data *trip_chain_data_for_statistics
    
    cdef double[:,::1] od_cost
    cdef double[::1] alpha
    cdef double beta
    cdef double sigma_t

    cdef bint input_data_have_been_read

    cdef void set_parameters(self, double* theta)
    cdef double total_utility_of_trip_chain(self, int *trip_chain, int *len_trip_chain, double *epsilon, double *epsilon_t)
    cdef void heuristic_search(self, int *trip_chain, int *len_trip_chain, double *epsilon, double *epsilon_t, int arrival_port, int depart_port)
    cdef void simulate(self)
    cdef double compute_v_fitness(self)
    cdef void simulate_for_statistics(self, bint order_insensitive, bint count_unobserved)

cdef double fitness(double* theta)