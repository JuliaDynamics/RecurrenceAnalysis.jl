"""
    windowed(rmat, f, width, step = 1; kwargs...)

A convenience function that applies the RQA function `f`, such as `determinism`,
to windowed views of the given recurrence matrix `rmat` with given window `width`
and `step`. The `kwargs...` are propagated to the call `f(rmat_view; kwargs...)`.
"""
function windowed(rmat, f, width, step=1; kwargs...)
  windows = 1:step:(size(rmat, 1)-width)
  map(1:length(windows)) do i
    f(rmat[windows[i]:(windows[i]+width), windows[i]:(windows[i]+width)]; kwargs...)
  end
end