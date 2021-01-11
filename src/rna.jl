"""
    rna(R)
    rna(args...; kwargs...)

Calculate a set of Recurrence Network parameters.
The input `R` can be a symmetric recurrence matrix that is
interpreted as the adjacency matrix of an undirected complex network,
such that linked vertices are neighboring points in the phase space.

Alternatively, the inputs can be a graph object or any 
valid inputs to the `SimpleGraph` constructor of the
[LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl) package.

## Return

The returned value is a dictionary that contains the following entries,
with the corresponding global network properties[1, 2]:
* `:density`: edge density, approximately equivalent to the global recurrence rate in the phase space.
* `:transitivity`: network transitivity, which describes the
global clustering of points following Barrat's and Weigt's formulation [3].
* `:averagepath`: mean value of the shortest path lengths taken over
all pairs of connected vertices, related to the average separation
between points in the phase.
* `:diameter`: maximum value of the shortest path lengths between
pairs of connected vertices, related to the phase space diameter.

## References

[1]: R.V. Donner *et al.* "Recurrence networks — a novel paradigm for nonlinear time series analysis",
*New Journal of Physics* 12, 033025 (2010)
[DOI:10.1088/1367-2630/12/3/033025](https://doi.org/10.1088/1367-2630/12/3/033025)

[2]: R.V. Donner *et al.*, The geometry of chaotic dynamics — a complex network perspective,
*Eur. Phys. J.* B 84, 653–672 (2011)
[DOI:10.1140/epjb/e2011-10899-1](https://doi.org/10.1140/epjb/e2011-10899-1)

[3]: A. Barrat & M. Weight, "On the properties of small-world network models",
*The European Physical Journal B* 13, 547–560 (2000)
[DOI:10.1007/s100510050067](https://doi.org/10.1007/s100510050067)
"""
function rna(args...; kwargs...)
    graph = SimpleGraph(args...; kwargs...)
    return Dict{Symbol, Float64}(
        :density => density(graph),
        :transitivity => global_clustering_coefficient(graph),
        :averagepath => mean(1 ./ closeness_centrality(graph)),
        :diameter => diameter(graph)
    )
end
