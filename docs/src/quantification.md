# Recurrence Quantification Analysis
A [`RecurrenceMatrix`](@ref) can be analyzed in several ways to yield information about the dynamics of the trajectory. All these various *measures* and functions are collectively called "Recurrence Quantification Analysis" (RQA).

To understand how each measure can be useful, we suggest to see the review articles listed in our documentation strings, namely:

1. N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems", *Phys. Reports 438*(5-6), 237-329 (2007).
2. N. Marwan & C.L. Webber, "Mathematical and computational foundations of recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).

You can also check the wikipedia page for [Recurrence quantification analysis](https://en.wikipedia.org/wiki/Recurrence_quantification_analysis).

The functions described in this page all accept a recurrence matrix (`x`), see [`RecurrenceMatrix`](@ref).

## RQA Measures

### All-in-one Bundle
In case you need all of the RQA-related functions (see below) and you don't want to write 10 lines of code to compute them all (since they are so many) we provide an all-in-one function that computes all of them and returns a dictionary with the results!
```@docs
rqa
```
---

!!! note "Return values for empty histograms"
    It may be the case that for a given recurrence matrix some structures do not exist at all. For example there are recurrence matrices that have no vertical lengths (or no vertical lengths with length less than `lmin`). In such cases the behavior of our RQA pipeline is the following:

    1. Quantities that represent maximum or average values are `0.0`.
    2. Quantities that represent entropies are `NaN`.

---

See also the [`@windowed`](@ref) macro for a windowed version of `rqa`.

## Example

For example, here are the RQA measures for the same example trajectory we used to make [Simple Recurrence Plots](@ref)
```@example MAIN
using RecurrenceAnalysis, DynamicalSystemsBase

# Create trajectory of Roessler system
@inbounds function roessler_rule(u, p, t)
    a, b, c = p
    du1 = -u[2]-u[3]
    du2 = u[1] + a*u[2]
    du3 = b + u[3]*(u[1] - c)
    return SVector(du1, du2, du3)
end
p0 = [0.15, 0.2, 10.0]
u0 = ones(3)
ro = CoupledODEs(roessler_rule, u0, p0)
N = 2000; Δt = 0.05
X, t = trajectory(ro, N*Δt; Δt, Ttr = 10.0)

# Make a recurrence matrix with fixed threshold
R = RecurrenceMatrix(X, 5.0)

# Compute RQA measures
rqa(R)
```

### Classical RQA Measures
```@docs
recurrencerate
determinism
dl_average
dl_max
dl_entropy
divergence
trend
```

### Extended RQA Measures
```@docs
laminarity
trappingtime
vl_average
vl_max
vl_entropy
```

### Recurrence Time Measures
```@docs
meanrecurrencetime
nmprt
rt_entropy
rt_average
```

## Keyword table
Since most of the above functions can be fined tuned with keyword arguments, here is a table summarizing them that could be of use:


| Argument  | Default   | Functions | Description |
| --------  | --------  | --------- | -----------
| `theiler` | 0 for `CrossRecurrenceMatrix`, 1 otherwise.  | `recurrencerate`, `determinism`, `*_average`, `*_max`, `*_entropy`, `divergence`, `trend`, `laminarity`, `trappingtime`,  `meanrecurrencetime`, `nmprt` | Theiler window: number of diagonals around the LOI **excluded** from the analysis. The value `0` means that the LOI is _included_ in the analysis. Use `1` to exclude the LOI. |
| `lmin`    | 2         | `determinism`, `*_average`, `*_max`, `*_entropy`, `divergence`, `laminarity`, `trappingtime`,  `meanrecurrencetime`, `nmprt` | Minimum length of the recurrent structures (diagonal or vertical) considered in the analysis. |
| `border`  | 10        | `trend`  | Number of diagonals excluded from the analysis near the border of the matrix. |

## Recurrence Structures Histograms
The functions that we list in this page internally compute histograms of some recurrence structures, like e.g. the vertical lengths.
You can access these values directly with the following function:
```@docs
recurrencestructures
```
