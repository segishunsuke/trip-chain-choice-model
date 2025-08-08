# -*- coding: utf-8 -*-
# cython: cdivision=True
# cython: boundscheck=False

import numpy as np
import csv
from libc.stdio cimport printf
from libc.stdlib cimport malloc, realloc, free, exit
from libc.math cimport exp, log, sqrt
from libc.time cimport time
cimport geneticr
import myoptimize

cdef extern from "mt19937ar.c":
    MT_state *init_mt_state() nogil
    void dealloc_mt_state(MT_state *mt_state) nogil
    void init_genrand(MT_state *mt_state, unsigned long s) nogil
    double genrand_real2(MT_state *mt_state) nogil
    double snormal_rand(MT_state *mt_state) nogil

cdef extern from "stdlib.h":
    void qsort(void *array, size_t count, size_t size, int (*compare)(const void*, const void*))

cdef int double_compare(const void *a, const void *b) noexcept:
    cdef double ia, ib
    ia = (<double*>a)[0]
    ib = (<double*>b)[0]
    if ia < ib:
        return -1
    elif ia == ib:
        return 0
    else:
        return 1

cdef int* imalloc(int n) noexcept nogil:
    cdef:
        int *buf
    buf = <int*>malloc(n * sizeof(int))
    if buf == NULL:
        printf("Memory error\n")
        exit(1)
    for i in range(n):
        buf[i] = 0
    return buf

cdef int** i2malloc(int n) noexcept nogil:
    cdef:
        int **buf
        int i
    buf = <int**>malloc(n * sizeof(int*))
    if buf == NULL:
        printf("Memory error\n")
        exit(1)
    for i in range(n):
        buf[i] = NULL
    return buf

cdef int*** i3malloc(int n) noexcept nogil:
    cdef:
        int ***buf
        int i
    buf = <int***>malloc(n * sizeof(int**))
    if buf == NULL:
        printf("Memory error\n")
        exit(1)
    for i in range(n):
        buf[i] = NULL
    return buf

cdef int* irealloc(int* buf, int n, int value) noexcept nogil:
    buf = <int*>realloc(buf, n * sizeof(int))
    if buf == NULL:
        printf("Memory error\n")
        exit(1)
    buf[n-1] = value
    return buf

cdef int** i2realloc(int** buf, int n1, int n2) noexcept nogil:
    buf = <int**>realloc(buf, n1 * sizeof(int*))
    if buf == NULL:
        printf("Memory error\n")
        exit(1)
    buf[n1-1] = <int*>malloc(n2 * sizeof(int))
    if buf[n1-1] == NULL:
        printf("Memory error\n")
        exit(1)
    return buf

# The Trip_chain_data object stores trip chain data for all arrival-departure port pairs.
cdef Trip_chain_data* init_trip_chain_data(int n_ad_pairs) noexcept nogil:
    cdef:
        Trip_chain_data *this
    this = <Trip_chain_data*>malloc(sizeof(Trip_chain_data))
    this.n_ad_pairs = n_ad_pairs
    this.n_trip_chains = imalloc(n_ad_pairs)
    this.len_trip_chains = i2malloc(n_ad_pairs)
    this.trip_chains = i3malloc(n_ad_pairs)
    this.freq_observation = i2malloc(n_ad_pairs)
    this.freq_simulation = i2malloc(n_ad_pairs)
    this.freq_simulation_not_observed = 0
    return this

# Appends a new trip chain pattern to the specified arrival-departure pair (id_ad_pair).
# After appending, sorts the trip chains for that pair in ascending order of length.
cdef void append_trip_chain_data(Trip_chain_data *this, int id_ad_pair, int *trip_chain, int len_trip_chain, int freq_observation) noexcept nogil:
    cdef:
        int id_trip_chain, j
        bint swap
    this.n_trip_chains[id_ad_pair] += 1
    this.len_trip_chains[id_ad_pair] = irealloc(this.len_trip_chains[id_ad_pair], this.n_trip_chains[id_ad_pair], len_trip_chain)
    this.trip_chains[id_ad_pair] = i2realloc(this.trip_chains[id_ad_pair], this.n_trip_chains[id_ad_pair], len_trip_chain)
    this.freq_observation[id_ad_pair] = irealloc(this.freq_observation[id_ad_pair], this.n_trip_chains[id_ad_pair], freq_observation)
    this.freq_simulation[id_ad_pair] = irealloc(this.freq_simulation[id_ad_pair], this.n_trip_chains[id_ad_pair], 0)
    
    id_trip_chain = this.n_trip_chains[id_ad_pair] - 1
    for j in range(len_trip_chain):
        this.trip_chains[id_ad_pair][id_trip_chain][j] = trip_chain[j]
    
    for id_trip_chain in range(this.n_trip_chains[id_ad_pair] - 1, 0, -1):
        swap = False
        if this.len_trip_chains[id_ad_pair][id_trip_chain] < this.len_trip_chains[id_ad_pair][id_trip_chain-1]:
            swap = True
        
        if swap:
            this.len_trip_chains[id_ad_pair][id_trip_chain], this.len_trip_chains[id_ad_pair][id_trip_chain-1] = this.len_trip_chains[id_ad_pair][id_trip_chain-1], this.len_trip_chains[id_ad_pair][id_trip_chain]
            this.trip_chains[id_ad_pair][id_trip_chain], this.trip_chains[id_ad_pair][id_trip_chain-1] = this.trip_chains[id_ad_pair][id_trip_chain-1], this.trip_chains[id_ad_pair][id_trip_chain]
            this.freq_observation[id_ad_pair][id_trip_chain], this.freq_observation[id_ad_pair][id_trip_chain-1] = this.freq_observation[id_ad_pair][id_trip_chain-1], this.freq_observation[id_ad_pair][id_trip_chain]
            this.freq_simulation[id_ad_pair][id_trip_chain], this.freq_simulation[id_ad_pair][id_trip_chain-1] = this.freq_simulation[id_ad_pair][id_trip_chain-1], this.freq_simulation[id_ad_pair][id_trip_chain]
        else:
            break

# Searches for a matching trip chain within the specified arrival-departure pair.
# Returns the index of the matched pattern, or -1 if not found.
# Comparison is done by internal places only (excluding arrival and departure ports).
cdef int search_trip_chain_data(Trip_chain_data *this, int id_ad_pair, int *trip_chain, int len_trip_chain) noexcept nogil:
    cdef:
        bint equal
        int id_trip_chain, j
    for id_trip_chain in range(this.n_trip_chains[id_ad_pair]):
        if this.len_trip_chains[id_ad_pair][id_trip_chain] == len_trip_chain:
            equal = True
            for j in range(1, len_trip_chain - 1):
                if this.trip_chains[id_ad_pair][id_trip_chain][j] != trip_chain[j]:
                    equal = False
                    break
            if equal:
                return id_trip_chain
        if this.len_trip_chains[id_ad_pair][id_trip_chain] > len_trip_chain:
            break  # Sorted by length: no need to continue
    return -1

# Sort the trip chains for each Arrival-Departure pair in ascending order of length.
# When lengths are equal, chains are sorted by descending observed frequency.
cdef void sort_trip_chain_data(Trip_chain_data *this) noexcept nogil:
    cdef:
        int id_ad_pair, id_trip_chain, j
        bint swap, swapped
    for id_ad_pair in range(this.n_ad_pairs):
        while True:
            swapped = False
            for id_trip_chain in range(this.n_trip_chains[id_ad_pair] - 1, 0, -1):
                swap = False
                if this.len_trip_chains[id_ad_pair][id_trip_chain] < this.len_trip_chains[id_ad_pair][id_trip_chain-1]:
                    swap = True
                elif this.len_trip_chains[id_ad_pair][id_trip_chain] == this.len_trip_chains[id_ad_pair][id_trip_chain-1]:
                    if this.freq_observation[id_ad_pair][id_trip_chain] > this.freq_observation[id_ad_pair][id_trip_chain-1]:
                        swap = True
                
                if swap:
                    this.len_trip_chains[id_ad_pair][id_trip_chain], this.len_trip_chains[id_ad_pair][id_trip_chain-1] = this.len_trip_chains[id_ad_pair][id_trip_chain-1], this.len_trip_chains[id_ad_pair][id_trip_chain]
                    this.trip_chains[id_ad_pair][id_trip_chain], this.trip_chains[id_ad_pair][id_trip_chain-1] = this.trip_chains[id_ad_pair][id_trip_chain-1], this.trip_chains[id_ad_pair][id_trip_chain]
                    this.freq_observation[id_ad_pair][id_trip_chain], this.freq_observation[id_ad_pair][id_trip_chain-1] = this.freq_observation[id_ad_pair][id_trip_chain-1], this.freq_observation[id_ad_pair][id_trip_chain]
                    this.freq_simulation[id_ad_pair][id_trip_chain], this.freq_simulation[id_ad_pair][id_trip_chain-1] = this.freq_simulation[id_ad_pair][id_trip_chain-1], this.freq_simulation[id_ad_pair][id_trip_chain]
                    swapped = True
            
            if not swapped:
                break

# Sort the trip chains for each Arrival-Departure pair in descending observed frequency.
cdef void sort_trip_chain_data_observed_frequency(Trip_chain_data *this) noexcept nogil:
    cdef:
        int id_ad_pair, id_trip_chain, j
        bint swap, swapped
    for id_ad_pair in range(this.n_ad_pairs):
        while True:
            swapped = False
            for id_trip_chain in range(this.n_trip_chains[id_ad_pair] - 1, 0, -1):
                swap = False
                if this.freq_observation[id_ad_pair][id_trip_chain] > this.freq_observation[id_ad_pair][id_trip_chain-1]:
                    swap = True
                
                if swap:
                    this.len_trip_chains[id_ad_pair][id_trip_chain], this.len_trip_chains[id_ad_pair][id_trip_chain-1] = this.len_trip_chains[id_ad_pair][id_trip_chain-1], this.len_trip_chains[id_ad_pair][id_trip_chain]
                    this.trip_chains[id_ad_pair][id_trip_chain], this.trip_chains[id_ad_pair][id_trip_chain-1] = this.trip_chains[id_ad_pair][id_trip_chain-1], this.trip_chains[id_ad_pair][id_trip_chain]
                    this.freq_observation[id_ad_pair][id_trip_chain], this.freq_observation[id_ad_pair][id_trip_chain-1] = this.freq_observation[id_ad_pair][id_trip_chain-1], this.freq_observation[id_ad_pair][id_trip_chain]
                    this.freq_simulation[id_ad_pair][id_trip_chain], this.freq_simulation[id_ad_pair][id_trip_chain-1] = this.freq_simulation[id_ad_pair][id_trip_chain-1], this.freq_simulation[id_ad_pair][id_trip_chain]
                    swapped = True
            
            if not swapped:
                break

cdef void delete_trip_chain_data(Trip_chain_data *this) noexcept nogil:
    cdef:
        int id_ad_pair, id_trip_chain, j
    if this != NULL:
        for id_ad_pair in range(this.n_ad_pairs):
            for id_trip_chain in range(this.n_trip_chains[id_ad_pair]):
                free(this.trip_chains[id_ad_pair][id_trip_chain])
            free(this.len_trip_chains[id_ad_pair])
            free(this.trip_chains[id_ad_pair])
            free(this.freq_observation[id_ad_pair])
            free(this.freq_simulation[id_ad_pair])
        free(this.n_trip_chains)
        free(this.len_trip_chains)
        free(this.trip_chains)
        free(this.freq_observation)
        free(this.freq_simulation)
        free(this)

cdef class Trip_chain_simulator:
    def __init__(self):
        self.trip_chain_data = NULL
        self.trip_chain_data_for_statistics = NULL
        self.n_iteration = 1
        self.random_seed = <unsigned long>time(NULL)
        self.beta = 10.0
        self.gam = 0.75
        self.poisson_shift = 1.0e-4
        self.input_data_have_been_read = False

    def __dealloc__(self):
        if self.trip_chain_data != NULL:
            delete_trip_chain_data(self.trip_chain_data)
            self.trip_chain_data = NULL
        global simulator_context
        if simulator_context is self:
            simulator_context = None

    def revise_random_seed(self):
        cdef unsigned long random_seed_prev
        random_seed_prev = self.random_seed
        while random_seed_prev == self.random_seed:
            self.random_seed = <unsigned long>time(NULL)

    def set_random_seed(self, unsigned long random_seed):
        self.random_seed = random_seed

    def set_poisson_shift(self, double poisson_shift):
        if poisson_shift <= 1.0e-15:
            raise ValueError("poisson_shift must be greater than 1.0e-15.")
        self.poisson_shift = poisson_shift

    cdef void set_parameters(self, double* theta):
        cdef int k
        for k in range(self.n_places):
            self.alpha[k] = theta[k]
        self.beta = theta[self.n_places]
        self.gam = theta[self.n_places + 1]

    def set_as_context(self):
        global simulator_context
        simulator_context = self

    def read_input_data(self, file_settings, file_od_cost, file_trip_chain, property_trip_chain, file_initial_parameter_values):
        cdef:
            int i, j, k, it, arrival_port, depart_port, id_ad_pair, id_trip_chain, len_trip_chain
            double od_cost_norm
            int *trip_chain = NULL

        self.input_data_have_been_read = False

        if self.trip_chain_data != NULL:
            delete_trip_chain_data(self.trip_chain_data)
            self.trip_chain_data = NULL

        with open(file_settings, "r") as fin:
            reader = csv.reader(fin)
            data_settings = [row for row in reader]

        self.n_places = int(data_settings[0][1])
        if self.n_places <= 0:
            raise ValueError(f"Number of places specified in {file_settings} must be positive.")
        self.n_ports = int(data_settings[1][1])
        if self.n_ports <= 0:
            raise ValueError(f"Number of ports specified in {file_settings} must be positive.")
        self.n_all_places = self.n_places + self.n_ports
        self.poisson_shift = float(data_settings[2][1])
        if self.poisson_shift <= 1.0e-15:
            raise ValueError(f"Shift parameter of Poisson likelihood specified in {file_settings} must be greater than 1.0e-15.")
        od_cost_norm = float(data_settings[3][1])

        with open(file_od_cost, "r") as fin:
            reader = csv.reader(fin)
            data_od_cost = [row for row in reader]

        self.od_cost = np.zeros((self.n_all_places, self.n_all_places))
        for it in range(1, len(data_od_cost)):
            i = int(data_od_cost[it][0])
            j = int(data_od_cost[it][1])
            self.od_cost[i,j] = float(data_od_cost[it][2])
        if od_cost_norm <= 1.0e-15:
            od_cost_norm = np.percentile(self.od_cost.base, 95)
        for i in range(self.n_all_places):
            for j in range(self.n_all_places):
                self.od_cost[i,j] /= od_cost_norm
        print(f"OD costs are divided by {od_cost_norm}.")

        with open(file_trip_chain, "r") as fin:
            reader = csv.DictReader(fin)
            dict_trip_chain = [row for row in reader]

        self.n_sample = len(dict_trip_chain)
        self.arrival_ports = np.zeros(self.n_sample, dtype=int)
        self.depart_ports = np.zeros(self.n_sample, dtype=int)
        
        self.od_freq_observation = np.zeros((self.n_all_places, self.n_all_places), dtype=int)
        self.od_freq_simulation = np.zeros((self.n_all_places, self.n_all_places), dtype=int)
        self.visits_observation = np.zeros((self.n_ports * self.n_ports, self.n_places), dtype=int)
        self.visits_simulation = np.zeros((self.n_ports * self.n_ports, self.n_places), dtype=int)

        trip_chain = imalloc(self.n_places + 2)
        self.trip_chain_data = init_trip_chain_data(self.n_ports * self.n_ports)

        for i in range(self.n_sample):
            trip_chain_py = np.fromstring(dict_trip_chain[i][property_trip_chain].strip("[]"), sep=" ", dtype=int)
            len_trip_chain = len(trip_chain_py)
            if len_trip_chain > self.n_places + 2:
                free(trip_chain)
                raise ValueError(f"Trip chain {i} is too long.")
            for j in range(len_trip_chain):
                trip_chain[j] = trip_chain_py[j]
                if not (0 <= trip_chain[j] < self.n_all_places):
                    free(trip_chain)
                    raise ValueError(f"Trip chain {i} includes an undefined place.")

            arrival_port = trip_chain[0]
            depart_port = trip_chain[len_trip_chain - 1]
            id_ad_pair = (arrival_port - self.n_places) * self.n_ports + (depart_port - self.n_places)
            if (not (self.n_places <= arrival_port < self.n_all_places)) or (not (self.n_places <= depart_port < self.n_all_places)):
                free(trip_chain)
                raise ValueError(f"Trip chain {i} includes an undefined port.")

            self.arrival_ports[i] = arrival_port
            self.depart_ports[i] = depart_port
            
            for j in range(len_trip_chain - 1):
                self.od_freq_observation[trip_chain[j], trip_chain[j + 1]] += 1
            for j in range(1, len_trip_chain-1):
                self.visits_observation[id_ad_pair, trip_chain[j]] += 1

            id_trip_chain = search_trip_chain_data(self.trip_chain_data, id_ad_pair, trip_chain, len_trip_chain)
            if id_trip_chain == -1:
                append_trip_chain_data(self.trip_chain_data, id_ad_pair, trip_chain, len_trip_chain, 1)
            else:
                self.trip_chain_data.freq_observation[id_ad_pair][id_trip_chain] += 1
        free(trip_chain)
        
        sort_trip_chain_data(self.trip_chain_data)
        
        self.alpha = np.zeros(self.n_places)
        with open(file_initial_parameter_values, "r") as fin:
            reader = csv.reader(fin)
            data_parameters = [row for row in reader]
        for k in range(self.n_places):
            self.alpha[k] = float(data_parameters[k][1])
        self.beta = float(data_parameters[self.n_places][1])
        if self.beta <= 0.0:
            raise ValueError(f"beta specified in {file_initial_parameter_values} must be positive.")
        self.gam = float(data_parameters[self.n_places + 1][1])
        if self.gam <= 0.0:
            raise ValueError(f"gam specified in {file_initial_parameter_values} must be positive.")

        self.input_data_have_been_read = True
    
    cdef void search(self, int *trip_chain, int *len_trip_chain, MT_state *mt_state, int arrival_port, int depart_port):
        cdef:
            int i, j, k, current_place, next_place
            double max_v, sum_p, v_continue, p_continue
            double *v
            double *p
            int *visited
            double z
        
        v = <double*>malloc(sizeof(double) * (self.n_places+1))
        p = <double*>malloc(sizeof(double) * (self.n_places+1))
        visited = <int*>malloc(sizeof(int) * self.n_places)
        
        for k in range(self.n_places):
            visited[k] = False
        
        j = 0
        trip_chain[0] = arrival_port
        
        while True:
            current_place = trip_chain[j]
            
            for k in range(self.n_places):
                if not visited[k]:
                    v[k] = self.alpha[k] - self.beta * self.od_cost[current_place, k]
            if j > 0:
                v[self.n_places] = - self.beta * self.od_cost[current_place, depart_port]
            
            max_v = -1.0e100
            for k in range(self.n_places):
                if not visited[k]:
                    if max_v < v[k]:
                        max_v = v[k]
            
            for k in range(self.n_places):
                if not visited[k]:
                    p[k] = exp(v[k] - max_v)
            
            sum_p = 0.0
            for k in range(self.n_places):
                if not visited[k]:
                    sum_p += p[k]
            
            for k in range(self.n_places):
                if not visited[k]:
                    p[k] /= sum_p
            
            if j > 0:
                if sum_p > 0.0:
                    v_continue = max_v + log(sum_p)
                    
                    if v_continue > v[self.n_places]:
                        max_v = v_continue
                    else:
                        max_v = v[self.n_places]
                    
                    p_continue = exp( self.gam * (v_continue - max_v) )
                    p[self.n_places] = exp( self.gam * (v[self.n_places] - max_v) )
                    
                    sum_p = p_continue + p[self.n_places]
                    
                    p_continue /= sum_p
                    p[self.n_places] /= sum_p
                else:
                    p_continue = 0.0
                    p[self.n_places] = 1.0
            else:
                p_continue = 1.0
                p[self.n_places] = 0.0
            
            z = genrand_real2(mt_state)
            sum_p = 0.0
            next_place = -1
            for k in range(self.n_places):
                if not visited[k]:
                    sum_p += p_continue * p[k]
                    if sum_p >= z:
                        next_place = k
                        visited[next_place] = True
                        break
            if next_place == -1:
                next_place = depart_port
            
            j += 1
            trip_chain[j] = next_place
            
            if next_place == depart_port:
                len_trip_chain[0] = j + 1
                break
        
        free(v)
        free(p)
        free(visited)

    cdef void simulate(self):
        cdef:
            int i, j, k, k2, id_trip_chain, i_iter
            int arrival_port, depart_port, id_ad_pair
            int *trip_chain
            int len_trip_chain[1]
            MT_state *mt_state

        mt_state = init_mt_state()
        init_genrand(mt_state, self.random_seed)
        trip_chain = <int*>malloc(sizeof(int) * (self.n_places + 2))

        for arrival_port in range(self.n_places, self.n_all_places):
            for depart_port in range(self.n_places, self.n_all_places):
                id_ad_pair = (arrival_port - self.n_places) * self.n_ports + (depart_port - self.n_places)
                for id_trip_chain in range(self.trip_chain_data.n_trip_chains[id_ad_pair]):
                    self.trip_chain_data.freq_simulation[id_ad_pair][id_trip_chain] = 0
        self.trip_chain_data.freq_simulation_not_observed = 0

        for i_iter in range(self.n_iteration):
            for i in range(self.n_sample):
                self.search(trip_chain, len_trip_chain, mt_state, self.arrival_ports[i], self.depart_ports[i])
                id_ad_pair = (self.arrival_ports[i] - self.n_places) * self.n_ports + (self.depart_ports[i] - self.n_places)
                id_trip_chain = search_trip_chain_data(self.trip_chain_data, id_ad_pair, trip_chain, len_trip_chain[0])
                if id_trip_chain == -1:
                    self.trip_chain_data.freq_simulation_not_observed += 1
                else:
                    self.trip_chain_data.freq_simulation[id_ad_pair][id_trip_chain] += 1

        dealloc_mt_state(mt_state)
        free(trip_chain)

    cdef double compute_v_fitness(self):
        cdef:
            int arrival_port, depart_port, id_ad_pair, id_trip_chain
            double v_fitness
        v_fitness = 0.0
        for arrival_port in range(self.n_places, self.n_all_places):
            for depart_port in range(self.n_places, self.n_all_places):
                id_ad_pair = (arrival_port - self.n_places) * self.n_ports + (depart_port - self.n_places)
                for id_trip_chain in range(self.trip_chain_data.n_trip_chains[id_ad_pair]):
                    v_fitness += self.trip_chain_data.freq_observation[id_ad_pair][id_trip_chain] * log(<double>(self.trip_chain_data.freq_simulation[id_ad_pair][id_trip_chain]) / <double>(self.n_iteration) + self.poisson_shift) - (<double>(self.trip_chain_data.freq_simulation[id_ad_pair][id_trip_chain]) / <double>(self.n_iteration) + self.poisson_shift)
        v_fitness += - (<double>(self.trip_chain_data.freq_simulation_not_observed) / <double>(self.n_iteration) + self.poisson_shift)
        return v_fitness

    cdef void simulate_for_statistics(self, bint order_insensitive, bint count_unobserved):
        cdef:
            int i, j, k, k2, id_trip_chain, i_iter
            int arrival_port, depart_port, id_ad_pair
            int *trip_chain
            int len_trip_chain[1]
            bint swapped

        mt_state = init_mt_state()
        init_genrand(mt_state, self.random_seed)
        trip_chain = <int*>malloc(sizeof(int) * (self.n_places + 2))

        for arrival_port in range(self.n_places, self.n_all_places):
            for depart_port in range(self.n_places, self.n_all_places):
                id_ad_pair = (arrival_port - self.n_places) * self.n_ports + (depart_port - self.n_places)
                for id_trip_chain in range(self.trip_chain_data_for_statistics.n_trip_chains[id_ad_pair]):
                    self.trip_chain_data_for_statistics.freq_simulation[id_ad_pair][id_trip_chain] = 0
        self.trip_chain_data_for_statistics.freq_simulation_not_observed = 0

        for k in range(self.n_all_places):
            for k2 in range(self.n_all_places):
                self.od_freq_simulation[k, k2] = 0
        for id_ad_pair in range(self.n_ports * self.n_ports):
            for k in range(self.n_places):
                self.visits_simulation[id_ad_pair, k] = 0

        for i_iter in range(self.n_iteration):
            for i in range(self.n_sample):
                self.search(trip_chain, len_trip_chain, mt_state, self.arrival_ports[i], self.depart_ports[i])
                
                if order_insensitive:
                    while True:
                        swapped = False
                        for j in range(1, len_trip_chain[0]-2):
                            if trip_chain[j] > trip_chain[j+1]:
                                trip_chain[j], trip_chain[j+1] = trip_chain[j+1], trip_chain[j]
                                swapped = True
                        if not swapped:
                            break
                
                id_ad_pair = (self.arrival_ports[i] - self.n_places) * self.n_ports + (self.depart_ports[i] - self.n_places)
                id_trip_chain = search_trip_chain_data(self.trip_chain_data_for_statistics, id_ad_pair, trip_chain, len_trip_chain[0])
                if id_trip_chain == -1:
                    if count_unobserved:
                        append_trip_chain_data(self.trip_chain_data_for_statistics, id_ad_pair, trip_chain, len_trip_chain[0], 0)
                        id_trip_chain = search_trip_chain_data(self.trip_chain_data_for_statistics, id_ad_pair, trip_chain, len_trip_chain[0])
                        self.trip_chain_data_for_statistics.freq_simulation[id_ad_pair][id_trip_chain] = 1
                    else:
                        self.trip_chain_data_for_statistics.freq_simulation_not_observed += 1
                else:
                    self.trip_chain_data_for_statistics.freq_simulation[id_ad_pair][id_trip_chain] += 1

                for j in range(len_trip_chain[0]-1):
                    self.od_freq_simulation[trip_chain[j], trip_chain[j + 1]] += 1
                for j in range(1, len_trip_chain[0]-1):
                    self.visits_simulation[id_ad_pair, trip_chain[j]] += 1

        dealloc_mt_state(mt_state)
        free(trip_chain)

    def print_statistics(self, int n_iteration, file_summary_trip_chain, file_summary_od_freq, file_summary_visits, bint order_insensitive = False, bint count_unobserved = False):
        cdef:
            int k, k2, id_trip_chain, i, j
            int arrival_port, depart_port, id_ad_pair, sum_visits, len_trip_chain, freq_observation
            int *trip_chain
            bint swapped

        if not self.input_data_have_been_read:
            raise RuntimeError("Input data must be read first using read_input_data().")

        if n_iteration <= 0:
            raise ValueError("n_iteration must be positive.")
        self.n_iteration = n_iteration
        
        delete_trip_chain_data(self.trip_chain_data_for_statistics)
        self.trip_chain_data_for_statistics = NULL
        
        trip_chain = imalloc(self.n_places + 2)
        self.trip_chain_data_for_statistics = init_trip_chain_data(self.n_ports * self.n_ports)
        
        for arrival_port in range(self.n_places, self.n_all_places):
            for depart_port in range(self.n_places, self.n_all_places):
                id_ad_pair = (arrival_port - self.n_places) * self.n_ports + (depart_port - self.n_places)
                for id_trip_chain in range(self.trip_chain_data.n_trip_chains[id_ad_pair]):
                    len_trip_chain = self.trip_chain_data.len_trip_chains[id_ad_pair][id_trip_chain]
                    freq_observation = self.trip_chain_data.freq_observation[id_ad_pair][id_trip_chain]
                    for j in range(len_trip_chain):
                        trip_chain[j] = self.trip_chain_data.trip_chains[id_ad_pair][id_trip_chain][j]
                    
                    if order_insensitive:
                        while True:
                            swapped = False
                            for j in range(1, len_trip_chain-2):
                                if trip_chain[j] > trip_chain[j+1]:
                                    trip_chain[j], trip_chain[j+1] = trip_chain[j+1], trip_chain[j]
                                    swapped = True
                            if not swapped:
                                break
                    
                    id_trip_chain = search_trip_chain_data(self.trip_chain_data_for_statistics, id_ad_pair, trip_chain, len_trip_chain)
                    if id_trip_chain == -1:
                        append_trip_chain_data(self.trip_chain_data_for_statistics, id_ad_pair, trip_chain, len_trip_chain, freq_observation)
                    else:
                        self.trip_chain_data_for_statistics.freq_observation[id_ad_pair][id_trip_chain] += freq_observation
        free(trip_chain)
        
        self.simulate_for_statistics(order_insensitive, count_unobserved)
        
        sort_trip_chain_data_observed_frequency(self.trip_chain_data_for_statistics)

        with open(file_summary_trip_chain, "w") as fout:
            fout.write("trip chain,fq observed,fq predicted\n")
            for arrival_port in range(self.n_places, self.n_all_places):
                for depart_port in range(self.n_places, self.n_all_places):
                    id_ad_pair = (arrival_port - self.n_places) * self.n_ports + (depart_port - self.n_places)
                    for id_trip_chain in range(self.trip_chain_data_for_statistics.n_trip_chains[id_ad_pair]):
                        fout.write("[")
                        for j in range(self.trip_chain_data_for_statistics.len_trip_chains[id_ad_pair][id_trip_chain]):
                            fout.write(str(self.trip_chain_data_for_statistics.trip_chains[id_ad_pair][id_trip_chain][j]))
                            if j < self.trip_chain_data_for_statistics.len_trip_chains[id_ad_pair][id_trip_chain] - 1:
                                fout.write(" ")
                            else:
                                fout.write("],")
                        fout.write(str(self.trip_chain_data_for_statistics.freq_observation[id_ad_pair][id_trip_chain]) + ",")
                        fout.write(str(<double>(self.trip_chain_data_for_statistics.freq_simulation[id_ad_pair][id_trip_chain]) / <double>(self.n_iteration)) + "\n")

        with open(file_summary_od_freq, "w") as fout:
            fout.write("origin,destination,fq observed,fq predicted\n")
            for i in range(self.n_all_places):
                for j in range(self.n_all_places):
                    fout.write(str(i) + ",")
                    fout.write(str(j) + ",")
                    fout.write(str(self.od_freq_observation[i, j]) + ",")
                    fout.write(str(<double>(self.od_freq_simulation[i, j]) / <double>(self.n_iteration)) + "\n")
        
        with open(file_summary_visits, "w") as fout:
            fout.write("arrival,departure,place,visits observed,visits predicted\n")
            for arrival_port in range(self.n_places, self.n_all_places):
                for depart_port in range(self.n_places, self.n_all_places):
                    id_ad_pair = (arrival_port - self.n_places) * self.n_ports + (depart_port - self.n_places)
                    
                    sum_visits = 0
                    for k in range(self.n_places):
                        sum_visits += self.visits_observation[id_ad_pair, k]
                    
                    if sum_visits > 0:
                        for k in range(self.n_places):
                            fout.write(str(arrival_port) + ",")
                            fout.write(str(depart_port) + ",")
                            fout.write(str(k) + ",")
                            fout.write(str(self.visits_observation[id_ad_pair, k]) + ",")
                            fout.write(str(<double>(self.visits_simulation[id_ad_pair, k]) / <double>(self.n_iteration)) + "\n")

    def print_fitness_of_multiple_cases(self, int n_iteration, int n_sample_cases, file_fitness_statistics):
        cdef:
            int id_sample_case
            double v_fitness
            double average_v_fitness, average_v_fitness_2

        if not self.input_data_have_been_read:
            raise RuntimeError("Input data must be read first using read_input_data().")

        if n_iteration <= 0:
            raise ValueError("n_iteration must be positive.")
        self.n_iteration = n_iteration

        if n_sample_cases <= 0:
            raise ValueError("n_sample_cases must be positive.")

        fout = open(file_fitness_statistics, "w")
        fout.write("Case,fitness\n")

        average_v_fitness = 0.0
        average_v_fitness_2 = 0.0

        for id_sample_case in range(n_sample_cases):
            self.revise_random_seed()
            self.simulate()
            v_fitness = self.compute_v_fitness()
            average_v_fitness += v_fitness
            average_v_fitness_2 += v_fitness ** 2
            fout.write(str(id_sample_case) + ",")
            fout.write(str(v_fitness) + "\n")

        average_v_fitness /= n_sample_cases
        average_v_fitness_2 = (average_v_fitness_2 - average_v_fitness ** 2 * n_sample_cases) / (n_sample_cases - 1)

        fout.write("average," + str(average_v_fitness) + "\n")
        fout.write("stddev," + str(sqrt(average_v_fitness_2)) + "\n")
        fout.close()

    def optimize_parameters(self, int n_iteration, int n_gene, int n_func_evaluation, int n_solutions, file_estimation_result):
        cdef:
            int k, i_solution
            double *minr
            double *maxr
            double **solutions
            double *buf
            double best_fitness, median
            geneticr.Population *population

        if not self.input_data_have_been_read:
            raise RuntimeError("Input data must be read first using read_input_data().")

        if n_iteration <= 0:
            raise ValueError("n_iteration must be positive.")
        self.n_iteration = n_iteration

        if n_gene <= 0:
            raise ValueError("n_gene must be positive.")

        if n_func_evaluation <= 0:
            raise ValueError("n_func_evaluation must be positive.")

        if n_solutions <= 0:
            raise ValueError("n_solutions must be positive.")

        minr = <double*>malloc(sizeof(double) * (self.n_places + 2))
        maxr = <double*>malloc(sizeof(double) * (self.n_places + 2))
        solutions = <double**>malloc(sizeof(double*) * n_solutions)
        for i_solution in range(n_solutions):
            solutions[i_solution] = <double*>malloc(sizeof(double) * (self.n_places + 3))
        buf = <double*>malloc(sizeof(double) * n_solutions)

        self.set_as_context()
        population = geneticr.init_population(n_gene, self.n_places + 2, minr, maxr, fitness)

        best_fitness = -1.0e100
        i_solution = 0
        while i_solution < n_solutions:
            self.revise_random_seed()
            for k in range(self.n_places):
                minr[k] = self.alpha[k] - 0.1
                maxr[k] = self.alpha[k] + 0.1
            minr[self.n_places] = self.beta - 1.0
            maxr[self.n_places] = self.beta + 1.0
            minr[self.n_places + 1] = self.gam - 0.01
            maxr[self.n_places + 1] = self.gam + 0.01

            for k in range(self.n_places):
                population.genes[0].arrayr[k] = self.alpha[k]
            population.genes[0].arrayr[self.n_places] = self.beta
            population.genes[0].arrayr[self.n_places + 1] = self.gam

            geneticr.simulate_evolution(population, n_func_evaluation, True)

            for k in range(self.n_places):
                self.alpha[k] = population.genes[0].arrayr[k]
            self.beta = population.genes[0].arrayr[self.n_places]
            self.gam = population.genes[0].arrayr[self.n_places + 1]

            for k in range(self.n_places):
                solutions[i_solution][k] = self.alpha[k]
            solutions[i_solution][self.n_places] = self.beta
            solutions[i_solution][self.n_places + 1] = self.gam
            solutions[i_solution][self.n_places + 2] = population.genes[0].fit

            if population.genes[0].fit > best_fitness:
                best_fitness = population.genes[0].fit
                i_solution = 0
            else:
                i_solution += 1

        with open(file_estimation_result, "w") as fout:
            fout.write(",median")
            for i_solution in range(n_solutions):
                fout.write(",Solution_" + str(i_solution))
            fout.write("\n")

            for k in range(self.n_places + 2):
                for i_solution in range(n_solutions):
                    buf[i_solution] = solutions[i_solution][k]
                qsort(buf, n_solutions, sizeof(double), double_compare)
                if n_solutions % 2 == 0:
                    median = 0.5 * (buf[n_solutions // 2 - 1] + buf[n_solutions // 2])
                else:
                    median = buf[n_solutions // 2]

                if k < self.n_places:
                    fout.write("alpha[" + str(k) + "]," + str(median))
                    self.alpha[k] = median
                elif k == self.n_places:
                    fout.write("beta," + str(median))
                    self.beta = median
                else:
                    fout.write("gam," + str(median))
                    self.gam = median

                for i_solution in range(n_solutions):
                    fout.write("," + str(solutions[i_solution][k]))
                fout.write("\n")

            fout.write("fitness,")
            for i_solution in range(n_solutions):
                fout.write("," + str(solutions[i_solution][self.n_places + 2]))
            fout.write("\n")

        free(minr)
        free(maxr)
        for i_solution in range(n_solutions):
            free(solutions[i_solution])
        free(solutions)
        free(buf)
        geneticr.dealloc_population(population)
    
    cdef double compute_log_likelihood(self):
        cdef:
            int i, j, k, k2, id_trip_chain, current_place, next_place
            int arrival_port, depart_port, id_ad_pair
            double value, max_v, sum_p, v_continue, p_continue
            double *v
            double *p
            bint *visited
        
        v = <double*>malloc(sizeof(double) * (self.n_places+1))
        p = <double*>malloc(sizeof(double) * (self.n_places+1))
        visited = <bint*>malloc(sizeof(bint) * self.n_places)
        
        value = 0.0
        for arrival_port in range(self.n_places, self.n_all_places):
            for depart_port in range(self.n_places, self.n_all_places):
                id_ad_pair = (arrival_port - self.n_places) * self.n_ports + (depart_port - self.n_places)
                for id_trip_chain in range(self.trip_chain_data.n_trip_chains[id_ad_pair]):
                    depart_port = self.trip_chain_data.trip_chains[id_ad_pair][id_trip_chain][self.trip_chain_data.len_trip_chains[id_ad_pair][id_trip_chain]-1]
                    
                    for k in range(self.n_places):
                        visited[k] = False
                    
                    for j in range(self.trip_chain_data.len_trip_chains[id_ad_pair][id_trip_chain]-1):
                        current_place = self.trip_chain_data.trip_chains[id_ad_pair][id_trip_chain][j]
                        next_place = self.trip_chain_data.trip_chains[id_ad_pair][id_trip_chain][j+1]
                        
                        if next_place < self.n_places:
                            if visited[next_place]:
                                import sys
                                sys.exit("The same place is visited twice.")
                        
                        for k in range(self.n_places):
                            if not visited[k]:
                                v[k] = self.alpha[k] - self.beta * self.od_cost[current_place, k]
                        if j > 0:
                            v[self.n_places] = - self.beta * self.od_cost[current_place, depart_port]
                        
                        max_v = -1.0e100
                        for k in range(self.n_places):
                            if not visited[k]:
                                if max_v < v[k]:
                                    max_v = v[k]
                        
                        for k in range(self.n_places):
                            if not visited[k]:
                                p[k] = exp(v[k] - max_v)
                        
                        sum_p = 0.0
                        for k in range(self.n_places):
                            if not visited[k]:
                                sum_p += p[k]
                        
                        for k in range(self.n_places):
                            if not visited[k]:
                                p[k] /= sum_p
                        
                        if j > 0:
                            if sum_p > 0.0:
                                v_continue = max_v + log(sum_p)
                                
                                if v_continue > v[self.n_places]:
                                    max_v = v_continue
                                else:
                                    max_v = v[self.n_places]
                                
                                p_continue = exp( self.gam * (v_continue - max_v) )
                                p[self.n_places] = exp( self.gam * (v[self.n_places] - max_v) )
                                
                                sum_p = p_continue + p[self.n_places]
                                
                                p_continue /= sum_p
                                p[self.n_places] /= sum_p
                            else:
                                p_continue = 0.0
                                p[self.n_places] = 1.0
                        else:
                            p_continue = 1.0
                            p[self.n_places] = 0.0
                        
                        if next_place >= self.n_places:
                            value += self.trip_chain_data.freq_observation[id_ad_pair][id_trip_chain] * log(p[self.n_places])
                        else:
                            value += self.trip_chain_data.freq_observation[id_ad_pair][id_trip_chain] * log(p_continue * p[next_place])
                        
                        if next_place < self.n_places:
                            visited[next_place] = True
        
        free(v)
        free(p)
        free(visited)
        return value
    
    def optimize_parameters_by_mle(self, file_estimation_result):
        cdef:
            int k
        
        if not self.input_data_have_been_read:
            raise RuntimeError("Input data must be read first using read_input_data().")
        
        self.set_as_context()
        
        theta = np.ones(self.n_places+2)
        
        for k in range(self.n_places):
            theta[k] = self.alpha[k]
        theta[self.n_places] = self.beta
        theta[self.n_places+1] = self.gam
        
        typx = np.ones(self.n_places+2)
        typx[self.n_places] = 0.0
        typx[self.n_places+1] = 0.0
        res = myoptimize.optimize(log_likelihood, theta, grad = grad_log_likelihood, typx = typx, tol = 1.0e-5)
        theta[:] = res.x
        log_likelihood(theta)
        
        with open(file_estimation_result, "w") as fout:
            for k in range(self.n_places):
                fout.write("alpha["+str(k)+"],"+str(self.alpha[k])+"\n")
            fout.write("beta,"+str(self.beta)+"\n")
            fout.write("gam,"+str(self.gam)+"\n")
            fout.write("log-likelihood,"+str(ll_current))

# For function dispatching via a global variable reference to the current simulator instance
cdef Trip_chain_simulator simulator_context = None

cdef double fitness(double* theta):
    if simulator_context is None:
        raise RuntimeError("Simulator context is not set.")
    simulator_context.set_parameters(theta)
    simulator_context.simulate()
    return simulator_context.compute_v_fitness()

cdef double ll_current = 0.0

def log_likelihood(theta):
    global ll_current
    cdef:
        int i
        double *theta_c
    if simulator_context is None:
        raise RuntimeError("Simulator context is not set.")
    theta_c = <double*>malloc(sizeof(double) * theta.shape[0])
    for i in range(theta.shape[0]):
        theta_c[i] = theta[i]
    simulator_context.set_parameters(theta_c)
    free(theta_c)
    ll_current = simulator_context.compute_log_likelihood()
    return - ll_current

def grad_log_likelihood(theta):
    print(f"thata: {theta}")
    print(f"log likelihood: {ll_current}\n")
    
    grad = theta.copy()
    f_current = - ll_current
    for i in range(theta.size):
        theta_p = theta.copy()
        if i in [simulator_context.n_places, simulator_context.n_places + 1]:
            dtheta = theta[i] * 1.0e-7
        else:
            dtheta = 1.0e-7
        theta_p[i] = theta[i] + dtheta
        f_p = log_likelihood(theta_p)
        grad[i] = (f_p - f_current) / dtheta
    return grad


