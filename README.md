[![Build Status](https://travis-ci.org/heliosdrm/RecurrenceAnalysis.jl.svg?branch=master)](https://travis-ci.org/heliosdrm/RecurrenceAnalysis.jl)

# RecurrenceAnalysis

Julia tools for Recurrence Quantification Analysis (RQA)

## Basic functions

The following functions can be used to quantify distances and recurrences in time series:

| Function                | Description                                                              |
| --------                | -----------                                                              |
| `embed`                 | Embed a time series in *m* dimensions with a given delay                 |
| `distancematrix`        | Create a distance matrix from one or two (possibly embedded) time series |
| `recurrencematrix`      | Create a recurrence matrix from a (possibly embedded) time series        |
| `crossrecurrencematrix` | Create a cross recurrence matrix from two (possibly embeded) time series |
| `jointrecurrencematrix` | Create a joint recurrence matrix from two (possibly embeded) time series |

Functions for the calculation of RQA parameters from recurrence matrices:

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

# Auxiliary functions

Functions to estimate the optimal delay for a time series:

| Function          | Description                                           |
| --------          | -----------                                           |
| `autocorrelation` | Autocorrelation for positive delays                   |
| `ami`             | Average mutual information                            |
| `gmi`             | Generalised mutual information based on Rényi entropy |


Functions to estimate the optimal embedding dimensions for a time series:

| Function | Description                                                            |
| -------- | -----------                                                            |
| `fnn`    | False Nearest Neighbours (FNN, Kennel et al., *Phys Rev A*:1992)         |
| `afnn`   | Averaged FNN (Cao, *Phys D Nonlin Phen*:1997)                            |
| `ffnn`   | False First Nearesth Neighbours (Krakovská et al., *J Complex Sys*:2015) |

# Plotting recurrence matrices

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

# To-do list:

 * FAN method to define recurrence plots
 * Criteria to define the optimal radius
 * Windowed RQA
 * Recurrence Network analysis
