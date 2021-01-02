# Recurrence network analysis measures

# Already defined: degree, density, local_clustering coefficient,
# closeness_centrality, laplacian_matrix, laplacian_spectrum
# Lacks: local efficiency

"""
    averageclustering(graph)

Returns the global clustering coefficient of a graph that
represents a recurrence network, defined as the average of
local clustering coefficients across all the graph's nodes.

# Description

For a given vertex `i` of the graph, the local clustering coefficient
measures the fraction of pairs neighboring vertices of `i` that are
mutually close, i.e.:[^1] 

```math
\\mathcal{C}_i = \\frac{\\textrm{number of triangles including vertex i}}{\\textrm{number of triples centred on vertex i}}
```
or given the adjacency matrix `A`:

```math
\\mathcal{C}_i = \\frac{\\sum_{j,k}^N {A_{jk}A_{ij}A_{ik}}}{\\sum_{j,k}^N {A_{ij}A_{ik}(1-\\delta_{jk})}}
```

(with ``\$\\mathcal{C}_i\$ = 0`` for vertices with less than two neighbors)[^2].

The global clustering coefficient introduced by Watts and Strogatz [^3]
is defined as the arithmetic mean of the local clustering coefficients
taken over all vertices of the network:

```math
\\hat{\\mathcal{C}} = \\frac{1}{N} \\sum_{i=1}^N{\\mathcal{C_i}}
```

This measure gives equal weights to all vertices, regardless of how much
connected they are to the rest of the network; if the connectivity is
heterogeneous across the network, the result will be dominated by contributions
from the most abundant type of vertices. In the case of very sparsely
connected networks, with a significant number of nodes with one or no connection,
this can lead to an underestimation of the actual fraction of triangles in the network.

For an alternative calculation, see [`transitivity`](@ref) -- which is
equivalent to [`global_clustering_coefficient`] in `LightGraphs.jl`

# References

[^1]: R.V. Donner *et al.*, [The geometry of chaotic dynamics — a complex network perspective, *Eur. Phys. J.* B 84, 653–672 (2011)](https://doi.org/10.1140/epjb/e2011-10899-1)

[^2] R.V. Donner *et al.* "Complex Network Analysis of Recurrences", in:
Webber, C.L. & Marwan N. (eds.) *Recurrence Quantification Analysis. Theory and Best Practices*, Springer, pp. 101-165 (2015).

[^3]: D.J. Watts & S.H. Strogatz, ["Collective dynamics of 'small-world' networks", *Nature 393*(6684), 440–442 (1998)](https://doi.org/10.1038%2F30918)
"""
averageclustering(graph) = mean(local_clustering_coefficient(graph))

"""
    transitivity(graph)

Returns the global clustering coefficient of a graph that
represents a recurrence network, following Barrat's and Weigt's
formulation, frequently referred to as *network transitivity*.

# Description

In graph theory the term *transitivity* is related to the relationship
of mutual adjacency in triples of vertices, such that if the
vertex `i` is connected with both `j` and `k`, then `j` and `k`
are also connected between them. The degree of transitivity for
a given vertex is defined as the relative frequency of such closed 3-loops
or "triangles", which is quantified by the *local clustering coefficient*:[^1] 

```math
\\mathcal{C}_i = \\frac{\\textrm{number of triangles including vertex i}}{\\textrm{number of triples centred on vertex i}}
```

The network transitivity measures the fraction of closed triangles with respect to
the number of linked triples of vertices in the whole network:

```math
\\mathcal{T} = \\frac{3 \\times \\textrm{number of triangles in the network}}{\\textrm{number of linked triples of vertices}}
```

This coincides with Barrat's and Weigt's definition of
*global clustering coefficient*[^2], and actually it is equivalent to
the function `global_clustering_coefficient` from the package `LightGraphs`.
It should not be confused with measure proposed by Watts and Strogatz based on
the average of the local clustering coefficient (see [`averageclustering`](@ref)).
The transitivity coefficient characterizes the effective global dimensionality
of the system, giving equal weight to all triangles in the network[^1],[^3].


# References

[^1]: R.V. Donner *et al.*, [The geometry of chaotic dynamics — a complex network perspective, *Eur. Phys. J.* B 84, 653–672 (2011)](https://doi.org/10.1140/epjb/e2011-10899-1)

[^2]: A. Barrat & M. Weight, ["On the properties of small-world network models", *The European Physical Journal B* 13, 547–560 (2000)](https://doi.org/10.1007/s100510050067)

[^3] R.V. Donner *et al.* "Complex Network Analysis of Recurrences", in:
Webber, C.L. & Marwan N. (eds.) *Recurrence Quantification Analysis. Theory and Best Practices*, Springer, pp. 101-165 (2015).
"""
transitivity(graph) = global_clustering_coefficient(graph)

"""
    averagepath(graph)

Returns the average path length of a graph that
represents a recurrence network, defined as the mean value of the
shortest path lengths taken over all pairs of vertices in the network.

For disconnected pairs of vertices, the shortest path length is set to zero by definition,
but in most practical applications, this has no major impact on the corresponding statistics.[^1]

This measure is calculated as the arithmetic mean of of the inverse
closeness centrality.[^2]

# References

[^1]: R.V. Donner *et al.*, ["Recurrence networks - a novel paradigm for nonlinear time series analysis"
*New Journal of Physics* 12,  033025 (2010)](https://doi.org/10.1088/1367-2630/12/3/033025)

[^2]: R.V. Donner *et al.* "Complex Network Analysis of Recurrences", in:
Webber, C.L. & Marwan N. (eds.) *Recurrence Quantification Analysis. Theory and Best Practices*, Springer, pp. 101-165 (2015).

"""
averagepath(graph) = mean(1 ./ closeness_centrality(graph))

# global efficiency pending on harmonic centrality
