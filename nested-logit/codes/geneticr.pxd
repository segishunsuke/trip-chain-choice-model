cdef extern from "mt19937ar.c":
    ctypedef struct MT_state

cdef struct Gene:
    double *arrayr
    double fit

cdef struct Population:
    int n_gene
    int size_arrayr
    double *minr
    double *maxr
    double (*func)(double*)
    double mutation_rater
    Gene *genes
    double *initial_arrayr
    MT_state *mt_state

cdef Population *init_population(int n_gene, int size_arrayr, double *minr, double *maxr, \
    double (*func)(double*), double mutation_rater = ?)
cdef void dealloc_population(Population *sself)
cdef void simulate_evolution(Population *sself, int n_func_evaluation, int initial_solution)

cdef int gene_compare(const void *a, const void *b)
cdef double dmin(double a, double b)
cdef double dmax(double a, double b)
cdef void init_gene(Gene *gene, int size_arrayr)
cdef void dealloc_gene(Gene *gene)
cdef void initialize_genes(Population *sself, int g_start)
cdef void compute_fitness(Population *sself, int g_start, int g_end)
cdef void sort_genes(Population *sself, int g_end)
cdef void selection(Population *sself, int *iparent1, int *iparent2)
cdef void crossover_and_mutation(Population *sself, int child1, int parent1, int parent2)
cdef int replacement(Population *sself)
