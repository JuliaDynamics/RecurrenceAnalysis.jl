using UnicodePlots
export recurrencetext, recurrencescatter, recurrenceplot
#########################################
# Text-based plotting
#########################################
"""
    recurrencescatter(R) -> xs, ys
Transform the data of a recurrence matrix to scatter data `xs, ys`
(which are the indices of the recurrence points).
"""
function recurrencescatter(R::ARM) # TODO: scan every ÷100 elements depending on size?
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

recurrencetext(R::ARM; kwargs...) = recurrencetext(stdout, R::ARM; kwargs...)

"""
    recurrencetext([io,] R; minh = 25, maxh = 0.5, ascii, kwargs...) -> u

Obtain a text-scatterplot representation of a recurrence matrix `R` to be displayed in
`io` (by default `stdout`). The matrix spans at minimum `minh` rows and at maximum
`maxh*displaysize(io)[1]` (i.e. by default half the display).
As we always try to plot in equal aspect ratio, if the width of the plot is even less,
the minimum height is dictated by the width.

The keyword `ascii::Bool` can ensure that all elements of the plot are ASCII characters
(`true`) or Unicode (`false`).

The rest of the `kwargs` are propagated into `UnicodePlots.scatterplot`.

Internally this function calls `RecurrenceAnalysis.recurrencescatter` to transform
a recurrence matrix into scatter data.
"""
function recurrencetext(io::IO, R::ARM; minh = 25, maxh = 0.5, ascii = nothing, kwargs...)
    @assert maxh ≤ 1
    h, w = displaysize(io)
    h = max(minh, round(Int, maxh * h)) # make matrix as long as half the screen (but not too short)
    s = 2.4 # scale that visually brings width and height to equal aspect ratio
    if w < round(Int, h*s) # ensure equal aspect ratio
        h = round(Int, w/s)
    else
        w = round(Int, h*s)
    end

    is, js = recurrencescatter(R)
    n, m = size(R)

    if ascii == true
        asciidef = (border = :ascii, canvas = DotCanvas)
    elseif ascii == false
        asciidef = (border = :solid, canvas = BrailleCanvas)
    elseif ascii == nothing # default handling
        if isdefined(Main, :IJulia) && Main.IJulia.inited
            # Always use ASCII in IJulia until this issue is fixed:
            # https://github.com/jupyter/notebook/issues/4354
            asciidef = (border = :solid, canvas = DotCanvas)
        elseif isdefined(Main, :Juno) && Main.Juno.isactive()
            asciidef = (border = :solid, canvas = BrailleCanvas)
        else # If not Juno or Jupyter, it is REPL
            if Sys.iswindows()
                asciidef = (border = :solid, canvas = DotCanvas)
            else
                asciidef = (border = :solid, canvas = BrailleCanvas)
            end
        end
    end

    # TODO: Change limits to tuples instead of arrays
    UnicodePlots.scatterplot(
        is, js; xlim = [1, n], ylim = [1, m], title = summary(R), labels = false,
        color = :cyan, width = w, height = h, asciidef..., kwargs...
    )
end

function Base.show(io::IO, ::MIME"text/plain", R::ARM)
    a = recurrencetext(io, R)
    show(io, a)
end

#########################################
# Transform to full plotting
#########################################

# Check the assigned width and height of the plot to ensure that they are
# approximately proportional to the matrix size
function checkgridsize(width::T, height::T, dims::Tuple{T,T}) where T<:Integer
    ratio = dims[2]/dims[1]
    intratios = Dict((height+1,width) => (height+1)/width,
                     (height-1,width) => (height-1)/width,
                     (height,width+1) => height/(width+1),
                     (height,width-1) => height/(width-1))
    distances = Dict(d=>(r-ratio) for (d,r) in intratios)
    distancesigns = sign.(collect(values(distances)))
    if all(x -> x >= 0, distancesigns) || all(x -> x <= 0, distancesigns)
        width_adj = round(Integer, height*dims[1]/dims[2])
        height_adj = round(Integer, width*dims[2]/dims[1])
        width = min(width, width_adj)
        height = min(height, height_adj)
        @warn "The specified dimensions are not proportional to the matrix size. "*
              "The size of the plot will be $width×$height."
    end
    (width, height)
end

# Define an overlapping grid of m blocks through n cells
function overlapgrid(m::T, n::T) where T<:Integer
    blocksize = n/m
    r = 0:blocksize:n
    initial = floor.(T, r[1:end-1].+1)
    final   = ceil.(T, r[2:end])
    collect(zip(initial,final))
end

# Calculate the level of "gray" (0=white, 1=black) corresponding to a matrix block
function block2grayscale(x, rind::Tuple{T,T}, cind::Tuple{T,T}) where T<:Integer
    submat = @view x[rind[1]:rind[2], cind[1]:cind[2]]
    ratio = count(!iszero, submat)/prod(size(submat))
end

"""
    recurrenceplot(x [, bwcode]; width::Int, height::Int, exactsize=false)

Transform the recurrence matrix `x` into a full matrix suitable for plotting as a
grayscale image. By default it returns a matrix with the same size as `x`,
but switched axes, containing "black" values in the cells that represent recurrent points,
and "white" values in the empty cells.

**This function does not do any plotting!** You have to use the return value with
the plotting library of your choice.

The numeric codes for black and white are given in a 2-element tuple as a second
optional argument. Its default value is `(0.0, 1.0)`, i.e. black is coded as `0.0`
(no brightness) and white as `1.0` (full brightness). The type of the elements
in the tuple defines the type of the returned matrix. This must be taken into
account if, for instance, the image is coded as a matrix of integers corresponding
to a grayscale; in such case the black and white codes must be given as numbers
of the required integer type.

The keyword arguments `width` and `height` can be given to define a custom size
of the image. If only one dimension is given, the other is automatically calculated.
If both dimensions are given, by default they are adjusted to keep an aspect
proportional to the original matrix, such that the returned matrix fits into a
matrix of the given dimensions. This automatic adjustment can be disabled by
passing the keyword argument `exactsize=true`.

If the image has different dimensions than `x`, the cells of `x` are distributed
in a grid with the size of the image, and a gray level between white and black
is calculated for each element of the grid, proportional to the number of
recurrent points contained in it. The levels of gray are coded as numbers of the
same type as the black and white codes.
"""
function recurrenceplot(x, bwcode::Tuple{TT,T}=(0.0,1.0);
    exactsize=false, kwargs...) where {TT<:Real, T<:Real}

    dims = size(x)
    kwargs = Dict(kwargs)
    if haskey(kwargs, :width) && !haskey(kwargs, :height)
        width = Integer(kwargs[:width])
        height = round(Integer, width*dims[2]/dims[1])
        return recurrenceplot(x, bwcode; width=width, height=height)
    elseif haskey(kwargs, :height) && !haskey(kwargs, :width)
        height = Integer(kwargs[:height])
        width = round(Integer, height*dims[1]/dims[2])
        return recurrenceplot(x, bwcode; width=width, height=height)
    elseif !haskey(kwargs, :width) || !haskey(kwargs, :height)
        width, height = Integer.(dims)
        return recurrenceplot(x, bwcode; width=width, height=height)
    end
    if exactsize
        width, height = Integer(kwargs[:width]), Integer(kwargs[:height])
    else
        width, height = checkgridsize(Integer(kwargs[:width]),
        Integer(kwargs[:height]), dims)
    end
    # initial and final values of the horizontal and vertical blocks
    rows = overlapgrid(width, dims[1])
    cols = overlapgrid(height, dims[2])
    p = zeros(height, width)
    # Fill img switching dimensions
    for c=1:width, r=1:height
        p[r,c] = block2grayscale(x, rows[c], cols[end-r+1])
    end
    # Normalize values
    p ./= maximum(p)
    # Change to color scale
    p .=  bwcode[1].*p .+ bwcode[2].*(1 .- p)
    pt = (T<:Integer) ? round.(T, p) : T.(p)
end
