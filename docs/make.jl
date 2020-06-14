using QBase
using Documenter

makedocs(;
    modules=[QBase],
    authors="Brian Doolittle <briand4@illinois.edu> and contributors",
    repo="https://github.com/chitambarlab/QBase.jl/blob/{commit}{path}#L{line}",
    sitename="QBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://chitambarlab.github.io/QBase.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/chitambarlab/QBase.jl",
)
