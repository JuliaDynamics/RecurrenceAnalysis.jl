using Graphs

"""
    SimpleGraph(R::AbstractRecurrenceMatrix)

Create a recurrent network as a `SimpleGraph`, from the symmetric recurrence matrix `R`.
That matrix is typically the result of calculating a [`RecurrenceMatrix`](@ref)
or a [`JointRecurrenceMatrix`](@ref).

The recurrence structure of `R` is interpreted as the adjacency matrix of an
undirected complex network, where two different vertices are connected if they are
neighbors in the embedded phase space, i.e.

```math
A_{i,j} = R_{i,j} - \\delta_{i,j}
```

Following this definition, diagonal points of `R` are ommited, i.e.
the graph does not contain self-connected nodes.

See the package [Graphs.jl](https://juliagraphs.org/Graphs.jl/stable/)
for further options to work with `SimpleGraph` objects, besides the functions
for Recurrence Network Analysis provided in this package.

# References

[1] : R.V. Donner *et al.* "Recurrence networks — a novel paradigm for nonlinear time series analysis",
*New Journal of Physics* 12, 033025 (2010).
[DOI:10.1088/1367-2630/12/3/033025](https://doi.org/10.1088/1367-2630/12/3/033025)

[2] : R.V. Donner *et al.* "Complex Network Analysis of Recurrences", in:
Webber, C.L. & Marwan N. (eds.) *Recurrence Quantification Analysis.
Theory and Best Practices*, Springer, pp. 101-165 (2015).
"""
function Graphs.SimpleGraphs.SimpleGraph(R::AbstractRecurrenceMatrix)
    graph = SimpleGraph(R.data)
    delta = SimpleGraphFromIterator(Edge(v,v) for v in 1:oldsize(graph, 1))
    graph = difference(graph, delta)
    return graph
end
