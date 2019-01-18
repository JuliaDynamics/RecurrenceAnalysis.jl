using RecurrenceAnalysis
using DynamicalSystemsBase
using SparseArrays
using UnicodePlots

ro = Systems.roessler(a=0.15, b=0.20, c=10.0)
N = 2000; dt = 0.05
tr = trajectory(ro, N*dt; dt = dt, Ttr = 10.0)
R = RecurrenceMatrix(tr, 5.0; metric = "euclidean")


# using PyPlot
# imshow(recurrenceplot(R), cmap = "binary_r")

using UnicodePlots

#prepare scatter data
function scatterdata(R)
   rows = rowvals(R)
   is = zeros(Int, nnz(R))
   js = zeros(Int, nnz(R))
   k = 1;
   m, n = size(R)
   for j = 1:n
     for i in nzrange(R, j)
        is[k] = j
        js[k] = rows[i]
        k += 1
     end
   end
   return is, js
end

is, js = scatterdata(R)
n, m = size(R)


a = UnicodePlots.scatterplot(
is, js, xlim = [1, n], ylim = [1, m], title = summary(R), labels = false,
color = :white, width = 80, height = 40)

# For full ascii method do:
# border = :ascii, canvas = DotCanvas
