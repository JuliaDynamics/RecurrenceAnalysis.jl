using LightGraphs

"""
    SimpleGraph(R::AbstractRecurrenceMatrix[, keepdiagonal=false])

Construct a `SimpleGraph` from the symmetric adjacency matrix defined by `R`.
That matrix is typically the result of calculating a [`RecurrenceMatrix`](@ref)
or a [`JointRecurrenceMatrix`](@ref).

Diagonal points, i.e. self-connected nodes of the graph, are ommited by default.
Set the optional argument `keepdiagonal` as `true` to keep them.
"""
function LightGraphs.SimpleGraphs.SimpleGraph(R::AbstractRecurrenceMatrix, keepdiagonal=false)
    graph = SimpleGraph(R.data)
    if !keepdiagonal
        delta = SimpleGraphFromIterator(Edge(v,v) for v=1:size(graph, 1))
        graph = difference(graph, delta)
    end
    return graph
end

"""
    SimpleDiGraph(R::AbstractRecurrenceMatrix)

Construct a `SimpleDiGraph` from the adjacency matrix defined by `R`. 
That matrix can be the result of calculating a [`RecurrenceMatrix`](@ref),
a [`JointRecurrenceMatrix`](@ref) or a [`CrossRecurrenceMatrix`](@ref).

The point `R[i,j]` of the matrix represents a connection from node `i` to `j`.
Points in the diagonal are taken as self-connections. 
"""
LightGraphs.SimpleGraphs.SimpleDiGraph(R::AbstractRecurrenceMatrix) = SimpleDiGraph(R.data)
