using RecurrenceAnalysis, DelimitedFiles, Statistics

# Measure the times (in ms) of evaluating an expression n times
macro measuretime(ex, n)
    quote
        # Train the expression and get the result
        result = $(esc(ex))
        t = zeros($n)
        for i in 1:$n
            t[i] = 1000*(@elapsed $(esc(ex)))
        end
        t, result
    end
end

# Function that will be measured
function fun_rqa(x,metric)
    xe = embed(x,3,6)
    rmat = RecurrenceMatrix(xe,1.2,metric=metric)
    rqa(rmat,theiler=1)
end

# Analyse 12 series from 250 to 3000 points 
# (With variable metric)
function benchmark(metric)
    m = readdlm("rossler.txt")
    for r=1:12
        x=m[1:250r,2r-1]
        tt, res = @measuretime fun_rqa(x,metric) 5
        t = median(tt)
        # Write table of results
        open("benchmark_rqa_julia_$metric.txt","a") do f
            write(f,"$r\t$t\t")
            for k in ["RR","DET","L","Lmax","ENTR","LAM","TT"]
                write(f, "$(res[k])\t")
            end
            write(f,"\n")
        end
    end
end
# Do it with max and euclidean norms
benchmark("euclidean")
benchmark("max")
