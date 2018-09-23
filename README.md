# RecurrenceAnalysis

**Julia tools for Recurrence Plots and Recurrence Quantification Analysis (RQA)**

[![Build Status](https://travis-ci.org/heliosdrm/RecurrenceAnalysis.jl.svg?branch=master)](https://travis-ci.org/heliosdrm/RecurrenceAnalysis.jl)

[![RecurrenceAnalysis](http://pkg.julialang.org/badges/RecurrenceAnalysis_0.6.svg)](http://pkg.julialang.org/?pkg=RecurrenceAnalysis)

## Introduction

The typical work for the analysis of recurrences in time series consists of:

1. (Optionally) embed the time series in an *m*-dimensional phase space.
2. Create a recurrence plot of the (possibly embedded) time series, or a
cross-recurrence plot or a joint-recurrence plot if the analysis concerns a pair of time series.
3. Quantify the observed recurrence patterns with Recurrence Quantification Analysis (RQA) parameters.

This can be done with the package `RecurrenceAnalysis` as shown in the next example (see more details about each step in the following sections):

```julia
using RecurrenceAnalysis
# 1. Transform time series x into an embedded vector
xe = embed(x, 5, 3)              # Embedding dimension = 5, delay = 3
# 2. Compute recurrence matrix
rmat = recurrencematrix(xe, 1.5) # Threshold ε = 1.5
# 3. Calculate RQA parameters
rr = recurrencerate(rmat)
detm = determinism(rmat)
ent = entropy(rmat)
# 4. Alternative: all RQA parameters at once
rqa_par = rqa(rmat)              # rqa_par["RR"] == rr, etc.
...
```

## Embedding time series

To embed a time series `x` in `m` dimensions with a delay of `d` samples between dimensions, just do:

```julia
xe = embed(x, m, d)
```

The time series `x` can be a vector of numeric values, or a matrix with `n` columns representing a series of `n`-dimensional vectors. In such case, the resulting `xe` will be a matrix with `m*n` columns.

## Recurrence plots

The following functions can be used to quantify distances and recurrences in time series:

| Function                | Description                                                              |
| --------                | -----------                                                              |
| `distancematrix`        | Create a distance matrix from one or two (possibly embedded) time series |
| `recurrencematrix`      | Create a recurrence matrix from a (possibly embedded) time series        |
| `crossrecurrencematrix` | Create a cross recurrence matrix from two (possibly embeded) time series |
| `jointrecurrencematrix` | Create a joint recurrence matrix from two (possibly embeded) time series |

The arguments required by those functions are the time series and (except for `distancematrix`) the radius or threshold distance in the phase space to consider that two points are recurrent.

In addition they can take the following keyword arguments:

| Argument    | Default   | Functions | Description |
| --------    | --------  | --------- | ----------- 
| `metric`    | `"max"`   | All       | Norm used to measure distances between points. Possible values: `"max"`, (maximum or infinity norm, also identified as `"inf"`), `"euclidean"` (*L*<sub>2</sub> or Euclidean norm), or `"manhattan"` (*L*<sub>1</sub> or Manhattan norm, also `"cityblock"` or `"taxicab"`). |
| `scale`     | 1         | `recurrencematrix`<br/>`crossrecurrencematrix`<br/>`jointrecurrencematrix` | Function or fixed number to scale the threshold or radius that is used to identify recurrences. Use `maximum` if the threshold is to be taken as a fraction of the maximum distance, `mean` if it is a fraction of the mean distance, etc., and `1` (identity scale, applied by default) to keep the threshold without scaling. |
| `fixedrate` | `false`   | `recurrencematrix`<br/>`crossrecurrencematrix`<br/>`jointrecurrencematrix` | Flag that indicates if the radius should be taken as a fixed target rate (between 0 and 1) of recurrences. |

#### Visualization of recurrence plots

Recurrence plots can be visualized and saved as images with the tools of other packages for Julia, which allow to render numeric matrices as images or "heat maps" that represent numbers as grades of a gray scale (or other color scales). Those packages are independent of `RecurrenceAnalysis` and must be installed separately. Some packages that can be used for this purpose are:

* [Images](https://github.com/JuliaImages/Images.jl)
(with [ImageView](https://github.com/JuliaImages/ImageView.jl) for visualization
on screen).
* [PyPlot](https://github.com/JuliaPy/PyPlot.jl)
* [GR](https://github.com/jheinen/GR.jl)
* [Plotly](https://github.com/plotly/Plotly.jl) or [PlotlyJS](https://github.com/sglyon/PlotlyJS.jl)
* [Gadfly](https://github.com/GiovineItalia/Gadfly.jl)
* [Makie](https://github.com/JuliaPlots/Makie.jl)
* [Winston](https://github.com/JuliaGraphics/Winston.jl)

The function `recurrenceplot` can be used to transform the sparse matrices of boolean values created by `recurrencematrix`, `crossrecurrencematrix` or `jointrecurrencematrix` into "ready-to-print" (or "to-display") dense numeric matrices of adequate dimensions, as generally required by the mentioned plotting packages.

For instance, using the package ImageView, the matrix `rmat` returned by `recurrencematrix` might be displayed on the screen just with:

```julia
using ImageView
rplot = recurrenceplot(rmat)
imshow(rp)
```

In that basic example, the matrix `rplot` will differ from `rmat` in the orientation of the dimensions (i.e. the main diagonal of `rmat` will run from the bottom left to the upper right cells of `rplot`), and in the scale of the values, which will be coded as `0.0` to represent black points where `rmat` had values (in the pairs of recurrent time points), and `1.0` to represent white points in the empty cells.

The numeric endpoints of the color scale and the dimensions of the plot can be customized as in the next example:

```julia
rplot = recurrenceplot(rmat, (0xff, 0x00), width=300)
```

* The pair of values of the second optional argument are the numeric codes for black and white points, respectively &mdash; `(0.0, 1.0)` by default. In this example the scale contains integer values ranging from 0 for white to 255 for black, in hexadecimal notation (from 00 to FF),
as required by the function `imagesc` of the package Winston.
* The keyword argument `width=300` sets the number of pixels that the plot will contain from left to right (the other dimension is calculated automatically to keep the proportions of `rmat`). Alternatively it is possible to set the `height`, or both to define the maximum dimensions of the plot. (To make the plot of the exact size `width`×`height`, regardless of the proportions of the original matrix, add `exactsize=true` as extra keyword argument.)

If the image produced by `recurrenceplot` is smaller than the size of the original matrix, each pixel will contain a "grey" value (between the codes of white and black) proportional to the number of recurrent points in the block of cells contained by the pixel.

Each plotting package can use different functions to display the matrices, like `imshow`, `heatmap`, `spy`, `imagesc`, etc.
The mapping between the rows/columns of the matrix and the X-Y axes of the plot may also differ between functions or packages for plotting.

Some plotting packages can be used with a common syntax as backends of the "meta-package" [Plots](https://github.com/JuliaPlots/Plots.jl/).
In that case, the function `plot` can be used with the series type `:image` or `:heatmap` to visualize the recurrence plots (see the support of that feature by the different backends in <http://docs.juliaplots.org/latest/supported/>, and notice that the default orientations of the Y-axis are different in both series types.)

### Recurrence Quantification Analysis (RQA)

The following functions are used to calculate [RQA parameters](http://www.recurrence-plot.tk/rqa.php) from recurrence matrices:

| Function         | Abbreviation | Description                                                        |
| --------         | ------------ | -----------                                                        |
| `recurrencerate` | `"RR"`       | Recurrence Rate: ratio of recurrent points                         |
| `determinism`    | `"DET"`      | Determinism: ratio of recurrent points forming diagonal structures |
| `avgdiag`        | `"L"`        | Average diagonal length                                            |
| `maxdiag`        | `"Lmax"`     | Length of longest diagonal                                         |
| `divergence`     | `"DIV"`      | Divergence: inverse of *Lmax*                                      |
| `entropy`        | `"ENT"`      | Entropy of the diagonal length histogram                           |
| `trend`          | `"TND"`      | Trend of paling recurrent points away from the LOI                 |
| `laminarity`     | `"LAM"`      | Laminarity: ratio of recurrent points forming vertical structures  |
| `trappingtime`   | `"TT"`       | Trapping time: average length of vertical structures               |
| `maxvert`        | `"Vmax"`     | Length of longest vertical structure                               |

In addition, the function `rqa` may be used to calculate all the RQA parameters, which are returned in a dictionary with keys equal to the abbreviations of the previous table. This alternative is more efficient than calculating the parameters individually, since the diagonal and vertical structures are analysed only once.

Some functions can be tuned with options, passed as keyword arguments:

| Argument  | Default   | Functions | Description |
| --------  | --------  | --------- | ----------- 
| `metric`  | `"max"`   | `distancematrix`<br/>`recurrencematrix`<br/>`crossrecurrencematrix`<br/>`jointrecurrencematrix` | Norm used to measure distances between points. Possible values: `"max"`, (maximum or infinity norm, also identified as `"inf"`), `"euclidean"` (Euclidean norm), or `"manhattan"` (*L*<sub>1</sub> or Manhattan distance, also `"cityblock"` or `"taxicab"`). |
| `scale`   | 1         | `recurrencematrix`<br/>`crossrecurrencematrix`<br/>`jointrecurrencematrix` | Function or fixed number to scale the threshold or radius that is used to identify recurrences. Use `maximum` if the threshold is to be taken as a fraction of the maximum distance, `mean` if it is a fraction of the mean distance, etc., and `1` (identity scale, applied by default) to keep the threshold without scaling. |
| `theiler` | 0         | `recurrencerate`<br/>`determinism`<br/>`avgdiag`<br/>`maxdiag`<br/>`divergence`<br/>`entropy`<br/>`trend`<br/>`laminarity`<br/>`trappingtime`<br/>`maxvert` | 'Theiler' window: number of diagonals around the LOI excluded from the analysis. |
| `lmin`    | 2         | `determinism`<br/>`avgdiag`<br/>`maxdiag`<br/>`divergence`<br/>`entropy`<br/>`laminarity`<br/>`trappingtime`<br/>`maxvert` | Minimum length of the recurrent structures (diagonal or vertical) considered in the analysis. |
| `border`  | 10        | `trend`  | Number of diagonals excluded from the analysis near the border of the matrix. |

Note: In order to keep the functions simpler and avoid confusion with the scaling of the distances, there is no option to normalize the input, although that is a customary procedure in RQA. This can be done *prior* to using the functions of this package.

The function `rqa` also accepts all those keyword arguments, which are passed down to the corresponding elementary functions. That function also accepts other keyword arguments for the Theiler window and the minimum length of recurrent structures, which are applied only in the calculation of some parameters:

* `theilerdiag` overrides `theiler` in the calculation of parameters related to diagonal structures, i.e. DET, L, Lmax, DIV, ENT and TND.
* `theilervert` overrides `theiler` in the calculation of parameters related to vertical structures, i.e. LAM, TT and Vmax.
* `lmindiag` overrides `lmin` in the calculation of parameters related to diagonal structures.
* `lminvert` overrides `vmin` in the calculation of parameters related to vertical structures.

If only the RR and the parameters based on diagonal structures are required (DET, L, Lmax, DIV and TND), `rqa` can be called with the keyword argument `onlydiagonal=true`, to reduce the time of calculations.

### Comparison with other RQA software packages

After version 0.2.0 some RQA methods and the defaults of their options were changed to make their results comparable with those provided by other software packages. See the corresponding section of the [Wiki](https://github.com/heliosdrm/RecurrenceAnalysis.jl/wiki/Comparison-of-software-packages-for-RQA) for details.

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

The AMI functions estimate the marginal and joint entropies of the original and delayed signals using histograms of their values. A third argument of the `ami` function can be used to tell how many bins those histograms must have. That argument must be either an integer, or one of the strings `"Sturges"` or `"FD"` (Freeman-Diaconis criterion). Sturges&rsquo; criterion is used by default.

The GMI functions use recurrence matrices to estimate Rényi entropies (cf. [Marwan et al., *Phys Rep*:2007](http://www.recurrence-plot.tk/marwan_PhysRep2007.pdf)). In order to make the calculations more efficient, the options available for the basic functions are not available; the assumed settings are maximum norm and *no scaling* (because the scale would not be preserved for different delays); i.e. the time series should be normalised beforehand, or use a threshold referred to the absolute values of the time series.


**Functions to estimate the optimal embedding dimensions for a time series:**

| Function | Description                                                                                                       |
| -------- | -----------                                                                                                       |
| `fnn`    | False Nearest Neighbours (FNN, [Kennel et al., *Phys Rev A*:1992](http://dx.doi.org/10.1103/PhysRevA.45.3403))    |
| `afnn`   | Averaged FNN ([Cao, *Phys D Nonlin Phen*:1997](http://dx.doi.org/10.1016/S0167-2789(97)00118-8))                  |
| `ffnn`   | False First Nearesth Neighbours ([Krakovská et al., *J Complex Sys*:2015](http://dx.doi.org/10.1155/2015/932750)) |

These functions have the optional keyword argument `metric`, with the same meaning and default value as in the calculation of RQA parameters.

**Functions to estimate the optimal threshold (radius) of recurrence matrices**

| Function          | Description                         |
| ----------------- | -----------                         |
| `sorteddistances` | Distances associated to *RR* values |

This function may be used to explore how the *RR* varies as a function of the chosen threshold, in order to define the threshold from a target *RR*, or look for the linear scaling region between those parameters ([Webber & Zbilut, 2005](http://www.nsf.gov/sbe/bcs/pac/nmbs/chap2.pdf), p. 56).

## Windowed RQA

In some cases, specially with very long time series, it may be suitable to perform the analysis at different points, considering only a limited &lsquo;window&rsquo; of data around each observation. The macro `@windowed` modifies the behaviour of the basic functions to calculate RQA parameters in that fashion. For instance, if `rmat` is a 10<sup>4</sup>&times;10<sup>4</sup> recurrence matrix, then
```julia
@windowed determinism(rmat, theiler=2, lmin=3) width=1000 step=100
```
will return a 91-element vector, such that each value is the determinism associated to a 1000-point fragment, starting at every 100 points (i.e. at `1`, `101`, &hellip; `9001`).

The general syntax of that macro is:
```julia
@windowed expr w                 #1
@windowed expr width=w step=s    #2
```
where:

 * `expr` is an expression used to calculate RQA parameters
 * `w` is the width of the window for relevant data around each point.
 * `s` is the step or distance between points where the calculations are done (starting in the first point).

To prevent syntax failures in the expansion of the macro, identify the RQA function (`rqa`, `recurrencerate`, `determinism`,&hellip;) directly by its name (avoid aliases), and use simple variable names (not complex expressions) for the arguments. On the other hand, the windowing options `width` and `step` can be given in any order. If `step` is ommitted, the calculations are done at every point, and the keyword `width` may be ommitted. (However, using `step=1` may be computationally very expensive, and that will provide just overly redundant results around each point, so it is advisable to set `step` a relatively big fraction of the window `width`.)

The value returned by the macro will normally be a vector with the same type of numbers as expected by `expr`. In the case of `@windowed rqa(...) ...`, it will return a dictionary with a similar structure as in the default `rqa` function, but replacing scalar values by vectors.

The macro `@windowed` can also be applied to the functions that calculate recurrence matrices (`recurrencematrix`, `crossrecurrencematrix`, `jointrecurrencematrix`). That creates a sparse matrix with the same size as if the macro was not used, but only containing valid values for pairs of points that belong to the `w` first main diagonals (i.e. the separation in time from one point to the other is `w` or smaller). The &lsquo;step&rsquo; parameter `s` has no effect on those functions. Such &lsquo;windowed&rsquo; matrices can be used as the input arguments to calculate windowed RQA parameters, obtaining the same results as if the complete matrix was used (under certain conditions, see below). For instance, the following calculations are equivalent:

```julia
# Using complete matrix
rmat = recurrencematrix(x, 1.5)
d = @windowed determinism(rmat) width=1000 step=250

# Using windowed matrix
rmatw = @windowed recurrencematrix(x, 1.5) 1000
d = @windowed determinism(rmatw) width=1000 step=250
```

The main difference between the two alternatives is that the second one will be faster and consume less memory. To ensure the equivalence between both approaches, the window width used to create the matrix must be greater than the one used to calculate the RQA parameters. Otherwise, the computation of RQA parameters might involve data points whose value is not well defined. Besides, the threshold to identify recurrences should be referred to a fixed scale. For instance:

```julia
rmat  =           recurrencematrix(x, 0.1, scale=maximum)
rmatw = @windowed recurrencematrix(x, 0.1, scale=maximum) 1000
rmat[1:1000,1:1000] == rmatw[1:1000,1:1000] # FALSE!!!
```
In this example, the `1000×1000` blocks of both matrices differ, because the threshold `0.1` is scaled with respect to the maximum distance between all points of `x` in `rmat`, but in the case of `rmatw` the scale changes between subsets of points. Something similar may happen if the recurrence matrix is calculated for a fixed recurrence rate (with the option `fixedrate=true`).

### Alternative syntax for `@windowed`

The following ways of using the macro `@windowed` are equivalent:

```julia
y = @windowed f(x,...) w
@windowed y=f(x,...) w
y = @windowed(f(x,...), w)
@windowed(y=f(x,...), w)
```

In all four cases, the width parameter `w` might have been qualified with a keyword as `width=w`. If the step parameter is added, the keyword qualification is mandatory.

