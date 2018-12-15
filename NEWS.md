# RecurrenceAnalysis.jl News

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
