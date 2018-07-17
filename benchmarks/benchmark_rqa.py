from pyunicorn.timeseries import RecurrencePlot
import numpy as np
from statistics import median
import time

# Measure the times (in ms) of evaluating an expression n times
def measuretime(f, n, *args):
    t = [0]*n
    res = f(*args)
    for n in range(n):
        t0 = time.time()
        f(*args)
        t[n] = time.time() - t0
    return(1000*np.array(t), res)

# Function that will be measured
def fun_rqa(v,metric):
    # Attempt sparse RQA if metric is euclidean
    metric_sup = (metric is "supremum")
    rp = RecurrencePlot(v, metric=metric, sparse_rqa=metric_sup,
    threshold=1.2, dim=3, tau=6)
    rqa = rp.rqa_summary()
    rqa["Lmax"] = rp.max_diaglength()
    rqa["ENT"] = rp.diag_entropy()
    rqa["TT"] = rp.trapping_time()
    return(rqa)

# Analyse 12 series from 250 to 3000 points 
# (With variable metric)
def benchmark(metric):
    m = np.loadtxt("rossler.txt")
    for r in range(12):
        x = m[:250*(r+1), 2*r]
        (tt, res) = measuretime(fun_rqa, 5, x, metric)
        t = median(tt)
        with open("benchmark_rqa_python_%s.txt"%metric, "a") as f:
            f.write("%d\t%f\t"%(r,t))
            for k in ["RR","DET","L","Lmax","ENT","LAM","TT"]:
                f.write("%s\t"%(res[k]))
            f.write("\n")
            

# Do it with max and euclidean norms
benchmark("euclidean")
benchmark("supremum")
