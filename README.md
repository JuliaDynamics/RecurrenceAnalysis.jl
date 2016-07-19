# RecurrenceAnalysis

**Julia tools for Recurrence Quantification Analysis (RQA)**

[![Build Status](https://travis-ci.org/heliosdrm/RecurrenceAnalysis.jl.svg?branch=master)](https://travis-ci.org/heliosdrm/RecurrenceAnalysis.jl)

[![RecurrenceAnalysis](http://pkg.julialang.org/badges/RecurrenceAnalysis_0.4.svg)](http://pkg.julialang.org/?pkg=RecurrenceAnalysis) [![RecurrenceAnalysis](http://pkg.julialang.org/badges/RecurrenceAnalysis_0.5.svg)](http://pkg.julialang.org/?pkg=RecurrenceAnalysis) 


## Basic functions

The following functions can be used to quantify distances and recurrences in time series:

| Function                | Description                                                              |
| --------                | -----------                                                              |
| `embed`                 | Embed a time series in *m* dimensions with a given delay                 |
| `distancematrix`        | Create a distance matrix from one or two (possibly embedded) time series |
| `recurrencematrix`      | Create a recurrence matrix from a (possibly embedded) time series        |
| `crossrecurrencematrix` | Create a cross recurrence matrix from two (possibly embeded) time series |
| `jointrecurrencematrix` | Create a joint recurrence matrix from two (possibly embeded) time series |

Functions for the calculation of [RQA parameters](http://www.recurrence-plot.tk/rqa.php) from recurrence matrices:

| Function         | Description                                                               |
| --------         | -----------                                                               |
| `recurrencerate` | Recurrence Rate (*RR*): ratio of recurrent points                         |
| `determinism`    | Determinism (*DET*): ratio of recurrent points forming diagonal lines     |
| `avgdiag`        | Average diagonal length (*L*)                                             |
| `maxdiag`        | Length of longest diagonal (*Lmax*)                                       |
| `divergence`     | Divergence (*DIV*): inverse of *Lmax*                                     |
| `entropy`        | Entropy (*ENT*) of the diagonal length histogram                          |
| `trend`          | Trend (*TND*) of paling recurrent points away from the LOI                |
| `laminarity`     | Laminarity (*LAM*): ratio of recurrent points forming vertical structures |
| `trappingtime`   | Trapping time (*TT*): average length of vertical structures               |
| `maxvert`        | Length of longest vertical structure (*Vmax*)                             |

A typical RQA workflow may be like this:

```julia
# 1. Transform time series x into an embedded vector
xe = embed(x, 5, 3)                 # Embedding dimension = 5, delay = 3 
# 2. Compute recurrence matrix
rmat = recurrencematrix(xe, 0.05)   # Threshold ε = 5% of maximum distance
# 3. Calculate RQA parameters
rr = recurrencerate(rmat)
detm = determinism(rmat)
ent = entropy(rmat)
...
```
### Plotting recurrence matrices

This package does not provide plotting tools. The matrices created by
`recurrencematrix`, `crossrecurrencematrix` and `jointrecurrencematrix` are
sparse matrices of boolean values, that can be converted to full matrices
of the appropriate type to display them on screen or save them as image files.

For instance, with the package [Winston](https://github.com/nolta/Winston.jl),
you can create a recurrence plot from the matrix `rmat`:

```julia
using Winston
colormap("grays")             # Set the colormap to have a B&W plot
rp = imagesc(0xff*full(rmat)) # Show the plot on the screen (in REPL mode)
savefig(rp, "rmat.png")       # Write the plot to disk (PNG)
```

(Note: the *y*-axis of such a plot is flipped with respect to the customary
orientation in recurrence plots.)

### RQA options:

Some functions can be tuned with options passed as keyword arguments:

| Argument  | Default   | Functions | Description |
| --------  | --------  | --------- | ----------- 
| `metric`  | `"max"`   | `distancematrix`<br/>`recurrencematrix`<br/>`crossrecurrencematrix`<br/>`jointrecurrencematrix` | Norm used to measure distances between points. Possible values: `max`, `inf` (infinity norm, same as `max`), and `euc` (Euclidean norm). |
| `scale`   | `maximum` | `recurrencematrix`<br/>`crossrecurrencematrix`<br/>`jointrecurrencematrix` | Function or fixed number to scale the distances between points before applying the given threshold. Use `maximum` (default) if the threshold is to be taken as a fraction of the maximum distance, `mean` if it is a fraction of the mean distance, etc., and `1` (identity scale) to use an absolute threshold. |
| `theiler` | `0` for `recurrencerate`,<br/>`1` for other functions  | `recurrencerate`<br/>`determinism`<br/>`avgdiag`<br/>`maxdiag`<br/>`divergence`<br/>`entropy`<br/>`trend` | 'Theiler' window: number of diagonals around the LOI excluded from the analysis. |
| `lmin`    | 2         | `determinism`<br/>`avgdiag`<br/>`maxdiag`<br/>`divergence`<br/>`entropy`<br/>`laminarity`<br/>`trappingtime`<br/>`maxvert` | Minimum length of the recurrent structures (diagonal or vertical) considered in the analysis. |
| `border`  | 10        | `trend`  | Number of diagonals excluded from the analysis near the border of the matrix. |

Note: In order to keep the functions simpler and avoid confusion with the scaling of the distances, there is no option to normalize the input, although that is a customary procedure in RQA. This can be done *prior* to using the functions of this package.


## Auxiliary functions

**Functions to estimate the optimal delay for a time series:**

| Function          | Description                                           |
| --------          | -----------                                           |
| `autocorrelation` | Autocorrelation for positive delays                   |
| `ami`             | Average mutual information                            |
| `gmi`             | Generalised mutual information based on Rényi entropy |

Examples:

```julia
autocorrelation(x) # Values from delay=0 to delay=length(x)-1
ami(x, 10)         # AMI for delay = 10
ami(x, 1:10)       # Same, for delays from 1 to 10
ami(x, (1,10))     # Same as above
gmi(x, 10, 0.1)    # GMI for delay = 10; recurrences for distances < 0.1
gmi(x, 1:10, 0.1)  # Same as with AMI ...
gmi(x, (1,10), 0.1)  
```
The AMI functions estimate the marginal and joint entropies of the original and delayed signals using histograms of their values. A third argument of the `ami` function can be used to tell how many bins those histograms must have. That argument must be either an integer, or one of the strings `"Sturges"` or `"FD"` (Freeman-Diaconis criterion). Sturges' criterion is used by default.

The GMI functions use recurrence matrices to estimate Rényi entropies (cf. [Marwan et al, *Phys Rep*:2007](http://www.recurrence-plot.tk/marwan_PhysRep2007.pdf)). In order to make the calculations more efficient, the options available for the basic functions are not available; the assumed settings are maximum norm and *no scaling* (because the scale would not be preserved for different delays); i.e. the time series should be normalised beforehand, or use a threshold referred to the absolute values of the time series.


**Functions to estimate the optimal embedding dimensions for a time series:**

| Function | Description                                                                                                       |
| -------- | -----------                                                                                                       |
| `fnn`    | False Nearest Neighbours (FNN, [Kennel et al., *Phys Rev A*:1992](http://dx.doi.org/10.1103/PhysRevA.45.3403))    |
| `afnn`   | Averaged FNN ([Cao, *Phys D Nonlin Phen*:1997](http://dx.doi.org/10.1016/S0167-2789(97)00118-8))                  |
| `ffnn`   | False First Nearesth Neighbours ([Krakovská et al., *J Complex Sys*:2015](http://dx.doi.org/10.1155/2015/932750)) |

These functions have the optional keyword argument `metric`, with the same meaning and default value as above.

# To-do list:

 * FAN method to define recurrence plots
 * Criteria to define the optimal radius
 * Windowed RQA
 * Recurrence Network analysis
