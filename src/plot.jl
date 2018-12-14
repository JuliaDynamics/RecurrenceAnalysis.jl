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
    if all(distancesigns .>= 0) || all(distancesigns .<= 0)
        width_adj = round(Integer, height*dims[1]/dims[2])
        height_adj = round(Integer, width*dims[2]/dims[1])
        width = min(width, width_adj)
        height = min(height, height_adj)
        @warn """the specified dimensions are not proportional to the matrix size
              The size of the plot will be $width×$height."""
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
function block2grayscale(x::AbstractMatrix, rind::Tuple{T,T}, cind::Tuple{T,T}) where T<:Integer
    submat = @view x[rind[1]:rind[2], cind[1]:cind[2]]
    ratio = count(!iszero, submat)/prod(size(submat))
end

"""
    recurrenceplot(x [, bwcode]; width::Int, height::Int, exactsize=false)

Transform the recurrence matrix `x` into a full matrix suitable for plotting as a
grayscale image. By default it returns a matrix with the same size as `x`,
but switched axes, containing "black" values in the cells that represent recurrent points,
and "white" values in the empty cells.

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
function recurrenceplot(x::AbstractMatrix, bwcode::Tuple{TT,T}=(0.0,1.0);
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
