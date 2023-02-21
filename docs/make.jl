cd(@__DIR__)

using RecurrenceAnalysis

pages = [
    "index.md",
    "rplots.md",
    "quantification.md",
    "networks.md",
    "comparison.md",
]

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

build_docs_with_style(pages, RecurrenceAnalysis)