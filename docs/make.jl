cd(@__DIR__)

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/apply_style.jl",
    joinpath(@__DIR__, "apply_style.jl")
)
include("apply_style.jl")

using RecurrenceAnalysis

RQA_PAGES = [
    "index.md",
    "rplots.md",
    "quantification.md",
    "networks.md",
    "windowed.md",
]

makedocs(
    modules = [RecurrenceAnalysis],
    format = Documenter.HTML(
        prettyurls = CI,
        assets = [
            asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
        ],
        collapselevel = 3,
    ),
    sitename = "RecurrenceAnalysis.jl",
    authors = "George Datseris",
    pages = RQA_PAGES,
    doctest = false,
    draft = false,
)

if CI
    deploydocs(
        repo = "github.com/JuliaDynamics/RecurrenceAnalysis.jl.git",
        target = "build",
        push_preview = true
    )
end
