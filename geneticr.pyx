# cython: cdivision=True
# cython: boundscheck=False

from libc.stdlib cimport malloc, realloc, free, exit
from libc.stdio cimport printf
import sys
import time

cdef extern from "mt19937ar.c":
    MT_state *init_mt_state()
    void dealloc_mt_state(MT_state *mt_state)
    void init_genrand(MT_state *mt_state, unsigned long s)
    double genrand_real1(MT_state *mt_state)
    double genrand_real2(MT_state *mt_state)
    double snormal_rand(MT_state *mt_state)

cdef extern from "stdlib.h":
    void qsort(void *array, size_t count, size_t size, int (*compare)(const void*, const void*))

cdef int gene_compare(const void *a, const void *b) noexcept:
    cdef Gene ia, ib
    ia = (<Gene*>a)[0]
    ib = (<Gene*>b)[0]
    if ia.fit > ib.fit:
        return -1
    elif ia.fit == ib.fit:
        return 0
    else:
        return 1

cdef double dmin(double a, double b):
    if a < b:
        return a
    else:
        return b

cdef double dmax(double a, double b):
    if a > b:
        return a
    else:
        return b

cdef void init_gene(Gene *gene, int size_arrayr):
    gene.arrayr = <double*>malloc(size_arrayr * sizeof(double))
    if gene.arrayr == NULL:
        printf('Memory error\n')
        exit(1)
    gene.fit = 0.0

cdef void dealloc_gene(Gene *gene):
    free(gene.arrayr)

cdef Population *init_population(int n_gene, int size_arrayr, double *minr, double *maxr, \
    double (*func)(double*), double mutation_rater = 0.25):
    
    cdef Population *sself
    cdef int g
    cdef MT_state *mt_state
    
    if n_gene%2 != 0:
        printf('n_gene must be an even number.\n')
        exit(1)
    
    sself = <Population*>malloc(sizeof(Population))
    if sself == NULL:
        printf('Memory error\n')
        exit(1)
    
    sself.n_gene = n_gene
    sself.size_arrayr = size_arrayr
    sself.minr = minr
    sself.maxr = maxr
    sself.func = func
    sself.mutation_rater = mutation_rater
    
    sself.genes = <Gene*>malloc((n_gene+n_gene/2) * sizeof(Gene))
    if sself.genes == NULL:
        printf('Memory error\n')
        exit(1)
    
    for g in range(n_gene+n_gene/2):
        init_gene(&sself.genes[g], size_arrayr)
    
    sself.mt_state = init_mt_state()
    init_genrand(sself.mt_state, int(time.time()))
    
    return sself

cdef void dealloc_population(Population *sself):
    cdef int g
    
    for g in range(sself.n_gene+sself.n_gene/2):
        dealloc_gene(&sself.genes[g])
    
    dealloc_mt_state(sself.mt_state)
    
    free(sself)

cdef void initialize_genes(Population *sself, int g_start):
    cdef int g, i
    
    for g in range(g_start, sself.n_gene):
        for i in range(sself.size_arrayr):
            sself.genes[g].arrayr[i] = sself.minr[i] \
                + (sself.maxr[i] - sself.minr[i])*genrand_real1(sself.mt_state)

cdef void compute_fitness(Population *sself, int g_start, int g_end):
    cdef int g
    
    for g in range(g_start, g_end):
        sself.genes[g].fit = sself.func(sself.genes[g].arrayr)

cdef void sort_genes(Population *sself, int g_end):
    qsort(sself.genes, g_end, sizeof(Gene), gene_compare)

cdef void selection(Population *sself, int *iparent1, int *iparent2):
    cdef double z, summ
    cdef int g, parent1, parent2
    
    z = genrand_real2(sself.mt_state) * 0.5*sself.n_gene*(sself.n_gene+1)
    summ = 0.0
    for g in range(sself.n_gene):
        summ += sself.n_gene - g
        if summ > z:
            break
    parent1 = g
    
    z = genrand_real2(sself.mt_state) * 0.5*sself.n_gene*(sself.n_gene-1)
    summ = 0.0
    for g in range(sself.n_gene):
        if g < parent1:
            summ += (sself.n_gene-1) - g
        elif g > parent1:
            summ += (sself.n_gene-1) - (g-1)
        if summ > z:
            break
    parent2 = g
    
    iparent1[0], iparent2[0] = parent1, parent2

cdef void crossover_and_mutation(Population *sself, int child, int parent1, int parent2):
    cdef int i
    cdef double m, M, z
    
    for i in range(sself.size_arrayr):
        m = dmin(sself.genes[parent1].arrayr[i], sself.genes[parent2].arrayr[i])
        M = dmax(sself.genes[parent1].arrayr[i], sself.genes[parent2].arrayr[i])
        z = genrand_real1(sself.mt_state)
        sself.genes[child].arrayr[i] = m + (M - m)*z
        
        z = genrand_real2(sself.mt_state)
        if z < sself.mutation_rater:
            sself.genes[child].arrayr[i] += (M - m)*snormal_rand(sself.mt_state)

cdef int replacement(Population *sself):
    cdef int pair, parent1, parent2, replacement_occured
    cdef double min_fit_parent, max_fit_child
    
    for pair in range(sself.n_gene/2):
        selection(sself, &parent1, &parent2)
        crossover_and_mutation(sself, sself.n_gene+pair, parent1, parent2)
    
    compute_fitness(sself, sself.n_gene, sself.n_gene+sself.n_gene/2)
    
    min_fit_parent = 1.0e100
    for g in range(0, sself.n_gene):
        if min_fit_parent > sself.genes[g].fit:
            min_fit_parent = sself.genes[g].fit
    
    max_fit_child = -1.0e100
    for g in range(sself.n_gene, sself.n_gene+sself.n_gene/2):
        if max_fit_child < sself.genes[g].fit:
            max_fit_child = sself.genes[g].fit
        
    if max_fit_child > min_fit_parent:
        replacement_occured = True
        sort_genes(sself, sself.n_gene+sself.n_gene/2)
    else:
        replacement_occured = False
    
    return replacement_occured

cdef void simulate_evolution(Population *sself, int n_func_evaluation, int initial_solution):
    cdef int i_func_evaluation, i, n_no_replacement, replacement_occured
    
    i_func_evaluation = 0
    
    if not initial_solution:
        initialize_genes(sself, 0)
    else:
        initialize_genes(sself, 1)
    
    compute_fitness(sself, 0, sself.n_gene)
    i_func_evaluation += sself.n_gene
    sort_genes(sself, sself.n_gene)
    
    n_no_replacement = 0
    while i_func_evaluation < n_func_evaluation:
        replacement_occured = replacement(sself)
        i_func_evaluation += sself.n_gene/2
        if replacement_occured:
            n_no_replacement = 0
        else:
            n_no_replacement += 1
        
        if n_no_replacement >= 0.1875*sself.n_gene:
            printf('Trapped into a local optimum. Reinitialzie the population.\n')
            initialize_genes(sself, 1)
            compute_fitness(sself, 1, sself.n_gene)
            i_func_evaluation += sself.n_gene - 1
            sort_genes(sself, sself.n_gene)
        
        printf('# Eval: %d/%d, Best fit: %lf\n', i_func_evaluation, n_func_evaluation, sself.genes[0].fit)
