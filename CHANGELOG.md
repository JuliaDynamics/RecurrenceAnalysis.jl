# RecurrenceAnalysis.jl News

## v2.0.0

- [**BREAKING** but deprecated] The way of specifying how recurrences are counted in `RecurrenceMatrix` has been completely overhauled and is now based on specifying the type by using actual type instances instead of keywords. See the new docstring of `RecurrenceMatrix`.
- [**BREAKING** but deprecated] Specifying metrics with strings in `RecurrenceMatrix` and co is deprecated and will be removed in the next version. Use proper `Metric` instances like `Euclidean()`.
- [**BREAKING** but deprecated] Giving a number as `scale` in `RecurrenceMatrix` is deprecated (because its pointless, just give `ε = ε*scale` instead).
- New exported function `recurrence_threshold` that calculates the recurrence threshold for a given recurrence specification type.
- RecurrenceAnalysis.jl has been updated to DynamicalSystems.jl (using StateSpaceSets.jl) and hosts its own docs.

## v1.8.1
- Better support for iteration over `AbstractRecurrenceMatrix` with functions from `Base`.

## v1.8.0
- Switch to Graphs.jl from LightGraphs.jl for RNA.

## v1.7.0
- New method of skeletonizing a Recurrence Matrix `skeletonize`, in order to eliminate tangential motion.

## v1.5.2
- New, faster formula for `rna[:averagepath]` allowing disconnected graphs.

## v1.5.1
- Minor documentation fixes.

## v1.5.0
- Integration with `LightGraphs` to create recurrence networks, and `rna` function.

## v1.4.0
- New method of Fixed Amount of Neighbors (FAN) method to create recurrence plots.
- Now `AbstractRecurrenceMatrix` types have a parameter with the possible values:
    * `WithinRange`: when neighboring points are determined by a given distance in the phase space.
    * `NeighborNumber`: when each point has a determined number of neighbors.
- Input arguments of `RecurrenceMatrix`, etc. are generalized to `AbstractDataset`.

## v1.3.2
- Bugfixes of `transitivity`
- Deprecate `AbstractMatrix` inputs to `RecurrenceMatrix` etc., in favor of `Dataset`.

## v1.3.0
- Increased performance of `RecurrenceMatrix`.
- Possible to make `JointRecurrenceMatrix` from input arguments.
- New RQA measure: `transitivity`.
- Dropped support for Julia 1.0, moving to new LTS 1.6.

## v1.2.0
*(no notable changes until so far)*

- Thread-parallelized computation of recurrence plots is now possible.

## v0.12.0
- Doc improvements
- Return type changes in `rqa`
- Breaking changes in implementations of `trend` and denominator of "RR"
- name changes in plotting

## v0.10.0
- Theiler window now has default values that depend on the type: `CrossRecurrenceMatrix` has default 0, otherwise 1.

## v0.9.0
- Recurrence matrices are printed as a plot using `UnicodePlots`. See `textrecurrenceplot`.

## v0.8.0 - v0.8.7
- Big documentation improvements
- Changes to return values for various rqa functions
- Functions that return recurrence structure histograms
- Updated and more tests
- Performance improvements

## v0.7.0
* All `recurrencematrix` have been deprecate to `RecurrenceMatrix` which
  create a dedicated struct for each kind of matrix.

## v0.6.0
* Rework how `recurrencematrix` is computed: now a conversion to `Dataset` is done for small dimensions. In addition, the default metric is now `Euclidean`.

## v0.5.0
* With this update core functionality of delay embedding and estimating embedding dimension is moved (and used from) the Julia package `DelayEmbeddings`.
* `RecurrenceAnalysis` is also joining JuliaDynamics and **DynamicalSystems.jl** from 0.5.0 onwards.
* Extensions have been implemented so that most functionality works also with `Dataset`s.

## 23-09-2018 - v0.4.0

* Function `recurrenceplot` to facilitate the visualization of recurrence plots.
* Support of Manhattan norm for the identification of recurrences.
* Option to create matrices with fixed recurrence rate.
* `rqa` option for obtaining only parameters related to diagonal structures.

## 12-08-2018 - v0.3.0

* Update to support Julia 1.0.0 and 0.7.0
* RQA functions matching outputs of **crqa** (R) and **CRPToolbox** (Matlab)

## 31-07-2018 - v0.2.1

* Bugfixes:
    - `sorteddistances`
    - RQA functions (bugs in parameters of diagonal structures when `theiler==0`)

## 30-07-2018 - v0.3.0-beta

* Break compatibility with 0.6
* Breaking modifications of RQA functions, in order to match the outputs of
  **crqa** (R) and **CRPToolbox** (Matlab):
    - The denominator of RR is fixed to N<sup>2</sup>.
    - Theiler window added to the calculation of LAM and TT.
    - Theiler window defaults to 0 in all functions.
* Bugfix in the calculation of diagonal structures

## 23-06-2018 - v0.2.0

* Update to support Julia 0.6-0.7 (breaking compatibility with 0.5).

## 10-11-2016 - v0.1.0

* Support for windowed RQA
* Change default scale of `recurrencematrix`, etc.
* Improved performance of RQA functions
* Function `rqa` for all RQA parameters
* Improved compatibility with Julia 0.5 and 0.6, backwards compatibility with 0.4
* Bugfixes:
    - `fnn`
    - `maxdiag`
    - `maxvert`

## 25-07-2016 - v0.0.3

* Compatibility with Julia 0.5
* Add `sorteddistances` to calculate optimal thresholds

## 25-05-2016 - v0.0.2

* Add functions to calculate optimal embedding dimensions
* Include documentation and tests

## 18-05-2016 - v0.0.1

First release or RecurrenceAnalysis.jl
