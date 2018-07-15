from pyunicorn.timeseries import RecurrencePlot
from numpy import loadtxt
from statistics import median
import time

# Measure the times (in ms) of evaluating an expression n times
def measuretime(f,n):
    t = [0]*n
    res = None
    for n in range(n):
        t0 = time.time()
        res = f()
        t[n] = time.time() - t0
    return(1000*t, res)

# Function that will be measured
def fun_rqa(v,metric):
    # Attempt sparse RQA if metric is euclidean
    metric_euc = (metric is "euclidean")
    rp = RecurrencePlot(v, metric=metric, sparse_rqa=metric_euc,
    threshold=1.2, dim=3, tau=6, recurrence_rate=1.2)
    rqa = rp.rqa_summary()
    rqa["Lmax"] = rp.max_diaglength()
    rqa["ENT"] = rp.diag_entropy()
    rqa["TT"] = rp.trapping_time()
    return(rqa)

# Analyse 12 series from 250 to 3000 points 
# (With variable metric)
def benchmark(metric):
    m = loadtxt("rossler.txt")
    for r in range(12):
        x = m[:250*r, 2*r]
        (tt, res) = measuretime(lambda: fun_rqa(x,metric), 5)
        t = median(tt)
        with open("benchmark_rqa_python_%s"%metric, "a") as f:
            f.write("%d\t%f\t"%(r,t))
            for k in ["RR","DET","L","Lmax","ENT","LAM","TT"]:
                f.write("%s\t"%(res[k]))
            f.write("\n")
            

# Do it with max and euclidean norms
benchmark("euclidean")
benchmark("supremum")
