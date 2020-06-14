using QBase
using Documenter

makedocs(;
    modules=[QBase,QBase.States,QBase.Observables,QBase.Unitaries],
    authors="Brian Doolittle <brian.d.doolittle@gmail.com> and contributors",
    repo="https://github.com/chitambarlab/QBase.jl/blob/{commit}{path}#L{line}",
    sitename="QBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://chitambarlab.github.io/QBase.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Overview" => "overview.md",
        "States" => "States.md",
        "Unitaries" => "Unitaries.md",
        "Observables" => "Observables.md",
        "Information" => "Information.md",
    ],
)

deploydocs(;
    repo="github.com/chitambarlab/QBase.jl",
)
