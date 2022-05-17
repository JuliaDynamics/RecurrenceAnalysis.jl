const rqa_funs = [
    :recurrencerate,
    :determinism,
    :dl_average,
    :dl_max,
    :divergence,
    :dl_entropy,
    :trend,
    :laminarity,
    :vl_average,
    :vl_max,
    :vl_entropy,
    :trappingtime,
    :rt_average,
    :rt_max,
    :rt_entropy,
    :meanrecurrencetime,
    :nmprt
]

const rqa_types = Dict(
    zip(rqa_funs,
    [eval(:(Base.return_types($f,(AbstractRecurrenceMatrix,))[1])) for f in rqa_funs]
    )
)

# Temporary workaround until `Dataset` can be sliced by single index
# (https://github.com/JuliaDynamics/DelayEmbeddings.jl/pull/73)
# Later on it will suffice with x[i]
_subsetdata(x::AbstractVector, i) = x[i]
_subsetdata(x::AbstractDataset, i) = x[i,:]

"""
    ij_block_rmat(x, y, bsize, dindex, vargs...; kwargs...)

Return the indices of the rows and columns of the nonzero values of a
block-diagonal cross-recurrence matrix.

If `m` is the cross-recurrence matrix of `x` and `y` (created with the
positional and keyword arguments `vargs` and `kwargs`), the indices returned by
this function are limited to the "block-diagonal" indicated by
`dindex ∈ {-1,0,1}`, as in the following graphical representation
(`#` represents the regions that are included, and `O` the excluded regions):

`dindex==-1`    `dindex==0`    `dindex==1`
 OOOOOO          ##OOOO         OO##OO
 OOOOOO          ##OOOO         OO##OO
 ##OOOO          OO##OO         OOOO##
 ##OOOO          OO##OO         OOOO##
 OO##OO          OOOO##         OOOOOO
 OO##OO          OOOO##         OOOOOO

The size of the blocks is `bsize × bsize`. The last block may be smaller if
`bsize` is not a divisor of the size of the whole cross-recurrence matrix.

The returned value is a tuple of two arrays, the first containing the indices
of the rows and columns of the nonzero values within the included regions.
"""
function ij_block_rmat(x, y, bsize, dindex, vargs...; kwargs...)
    n = min(size(x,1), size(y,1))
    brange = 1:bsize
    n_fullblocks = div(n, bsize)
    ix = iy = 0:0
    rws = (Int)[]
    cls = (Int)[]
    for i in abs(dindex):n_fullblocks-1
        ix = i*bsize .+ brange
        iy = i*bsize .+ brange
        if dindex < 0
            iy = iy .+ dindex*bsize
        elseif dindex > 0
            ix = ix .- dindex*bsize
        end
        rmat_b = CrossRecurrenceMatrix(_subsetdata(x, ix), _subsetdata(y, iy), vargs...; kwargs...)
        append!(rws, rowvals(rmat_b) .+ix[1] .- 1)
        append!(cls, colvals(rmat_b) .+iy[1] .- 1)
    end
    ix1 = ix[end]+1
    iy1 = iy[end]+1
    ix2 = min(ix1+bsize-1, n)
    iy2 = min(iy1+bsize-1, n)
    rx = ix1:ix2
    ry = iy1:iy2
    if length(rx) > 0 && length(ry) > 0
        rmat_b = CrossRecurrenceMatrix(_subsetdata(x, ix1:ix2), _subsetdata(y, iy1:iy2), vargs...; kwargs...)
        append!(rws, rowvals(rmat_b) .+ ix1 .- 1)
        append!(cls, colvals(rmat_b) .+ iy1 .- 1)
    end
    rws, cls
end

# check if call is a constructor of a type with given name 
function _check_constructor(call, name)
    call == name && return true
    if (call isa Expr) && call.head == :curly && call.args[1] == name
        return true
    else
        return false
    end
end

"""
    @windowed(f(x,...), width)
    @windowed(f(x,...); width, step=1)

Calculate windowed RQA parameters with a given window width.

`f(x,...)` may be any call to RQA functions (e.g. [`recurrencerate`](@ref),
[`determinism`](@ref), etc.),
with `x` being a named variable that designates the recurrence matrix
(do not use in-place calculations of the recurrence matrix).
The results are returned in a vector with one value for each position of the window.
By default the window moves at one-point intervals, but a longer `step` length
may be specified, together with the window `width`,
by declaring those options as keyword arguments.

This macro may be also used with recurrence matrix constructors
(`RecurrenceMatrix`, `CrossRecurrenceMatrix`, `JointRecurrenceMatrix`),
to create 'incomplete' matrices that are suitable for such windowed RQA.
The values of the resulting matrix in the diagonals within the window width will
be equal to those obtained without the `@windowed` macro, if the distances are
not scaled (using the option `scale=1`, see [`RecurrenceMatrix`](@ref)).
Outside the window width, the values of the recurrence matrix will be undefined
(mostly zero).
"""
macro windowed(ex, options...)
    # Expression can be of type a = f(x...)
    if in(ex.head, [:(=), :kw])
        left, right = (ex.args...,)
        return esc(:($left = @windowed($right,$(options...))))
    end
    # Parse options
    dict_op = Dict{Symbol,Any}(:step=>1)
    if length(options) == 1
        dict_op[:width] = options[1]
    end
    for op in options
        if typeof(op)<:Expr && in(op.head, [:(=), :kw])
            dict_op[op.args[1]] = op.args[2]
        end
    end
    if ex.head == :call
        f = ex.args[1]
        # Iteration of RQA functions
        # fun(x,...) => [fun(x[i+w,i+w],...) for i in 0:s:nw]
        if in(f, rqa_funs)
            x = ex.args[2]
            submat = :(mtype($x[i.+w,i.+w])) # e.g. :(RecurrenceMatrix(x[i.+w,i.+w]))
            ex.args[2] = submat # x is replaced by the windowed matrix
            ret_ex = quote
                local w = 1:$(dict_op[:width])
                local s = $(dict_op[:step])
                # only calculate "complete" blocks - until `nw`
                local nw = size($x,1) - $(dict_op[:width])
                local mtype = typeof($x) # type of the recurrence matrix
                ($(rqa_types[f]))[$ex for i in 0:s:nw]
            end
            return esc(ret_ex)
        end
        # Iteration of all RQA parameters
        if f == :rqa
            x = ex.args[2]
            submat = :(mtype($x[(i-1)*s.+w,(i-1)*s.+w]))
            ex.args[2] = submat
            ret_ex = quote
                local w = 1:$(dict_op[:width])
                local s = $(dict_op[:step])
                local nw = size($x,1) - $(dict_op[:width])
                local ni = div(nw, s)+1 # number of items
                local mtype = typeof($x)
                local rqa_dict = Dict{Symbol, Vector{Float64}}(
                    :RR    => zeros(ni),
                    :TRANS => zeros(ni),
                    :DET   => zeros(ni),
                    :L     => zeros(ni),
                    :Lmax  => zeros(ni),
                    :DIV   => zeros(ni),
                    :ENTR  => zeros(ni),
                    :TREND => zeros(ni),
                    :LAM   => zeros(ni),
                    :TT    => zeros(ni),
                    :Vmax  => zeros(ni),
                    :VENTR => zeros(ni),
                    :MRT   => zeros(ni),
                    :RTE   => zeros(ni),
                    :NMPRT => zeros(ni)
                )
                for i in 1:ni
                    local rqa_i = $ex
                    if i==1 # filter parameters
                        filter!(p->p.first in keys(rqa_i), rqa_dict)
                    end
                    #@show rqa_i
                    for (k,v) in rqa_i
                        rqa_dict[k][i] = v
                    end
                end
                rqa_dict
            end
            return esc(ret_ex)
        end
        # Iteration of matrix construction functions
        if _check_constructor(f, :CrossRecurrenceMatrix)
            # ij_block_rmat(x,y,width,d,...) with d=-1,0,1
            x = ex.args[2]
            y = ex.args[3]
            ex.args[1] = :(RecurrenceAnalysis.ij_block_rmat)
            insert!(ex.args, 4, dict_op[:width])
            insert!(ex.args, 5, -1)
            exd_lower  = :(local i, j = $(parse(string(ex)))) # lower diag block
            ex.args[5] = 0
            exd_center = :(local ii, jj = $(parse(string(ex)))) # central diag block
            ex.args[5] = 1
            exd_upper  = :(local ii, jj = $(parse(string(ex)))) # upper diag block
            ret_ex = quote
                $exd_lower
                $exd_center
                append!(i, ii)
                append!(j, jj)
                $exd_upper
                append!(i, ii)
                append!(j, jj)
                local m = RecurrenceAnalysis.sparse(i,j,true,size($x,1),size($y,1))
                CrossRecurrenceMatrix(m)
            end
            return esc(ret_ex)
        elseif _check_constructor(f, :RecurrenceMatrix)
            # ij_block_rmat(x,x,width,d,...) with d=-1,0
            ex.args[1] = :(RecurrenceAnalysis.ij_block_rmat)
            x = ex.args[2]
            insert!(ex.args, 3, x)
            insert!(ex.args, 4, dict_op[:width])
            insert!(ex.args, 5, -1)
            exd_lower  = :(local i, j = $(parse(string(ex)))) # lower diag block
            ex.args[5] = 0
            exd_center = :(local ii, jj = $(parse(string(ex)))) # central diag block
            ret_ex = quote
                $exd_lower
                i, j = [i;j], [j;i] # the upper diag block is the transpose of the lower
                $exd_center
                append!(i,ii)
                append!(j,jj)
                local n = size($x,1)
                local m = RecurrenceAnalysis.sparse(i,j,true,n,n)
                RecurrenceMatrix(m)
            end
            return esc(ret_ex)
        elseif _check_constructor(f, :JointRecurrenceMatrix)
            x = ex.args[2]
            y = ex.args[3]
            minsz = :(min(size($x,1),size($y,1)))
            subx = :(RecurrenceAnalysis._subsetdata($x, 1:$minsz))
            suby = :(RecurrenceAnalysis._subsetdata($y, 1:$minsz))
            # Call `@windowed RecurrenceMatrix` twice, recycling the argument (`x` and `y`)
            ex.args[1] = :RecurrenceMatrix
            deleteat!(ex.args, 3)
            ex.args[2] = subx
            ex_rmx = :(local rm1 = @windowed($(parse(string(ex))),width=$(dict_op[:width])))
            ex.args[2] = suby
            ex_rmy = :(local rm2 = @windowed($(parse(string(ex))),width=$(dict_op[:width])))
            ret_ex = quote
                $ex_rmx
                $ex_rmy
                JointRecurrenceMatrix(rm1.data .& rm2.data)
            end
            return esc(ret_ex)
        end
        # Throw error if it is not a valid function
        throw(ErrorException("$(string(ex.args[1])) is not a valid function for windowing"))
    end
    # Throw error if it didn't return
    throw(ErrorException("invalid expression for windowing"))
end
