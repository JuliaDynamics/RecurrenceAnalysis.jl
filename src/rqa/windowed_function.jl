"""
    windowed(rmat, f, width, step = 1; kwargs...)

A convenience function that applies the RQA function `f`, such as `determinism`,
to windowed views of the given recurrence matrix `rmat` with given window `width`
and `step`. The `kwargs...` are propagated to the call `f(rmat_view; kwargs...)`.
"""
function windowed(rmat::Union{ARM,AbstractMatrix}, f::Function, width::Integer, step=1::Integer; kwargs...)
  (width < 2) && throw(ErrorException(
        "Window width must be must be greater than or equal to 2"))
  (step < 1) && throw(ErrorException(
        "Step size must be greater than or equal to 1"))
  windows = 1:step:(size(rmat, 1)-width+1)
  map(1:length(windows)) do i
    f(rmat[windows[i]:(windows[i]+width-1), windows[i]:(windows[i]+width-1)]; kwargs...)
  end
end