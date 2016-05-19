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
| `ami`             | Average mutual information (Kennel's algorithm)       |
| `gmi`             | Generalised mutual information based on Rényi entropy |
