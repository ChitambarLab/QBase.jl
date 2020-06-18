using QBase
using Documenter

makedocs(;
    modules=[QBase],
    authors="Brian Doolittle <brian.d.doolittle@gmail.com> and contributors",
    repo="https://github.com/ChitambarLab/QBase.jl/blob/{commit}{path}#L{line}",
    sitename="QBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ChitambarLab.github.io/QBase.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "User Guide" => "user_guide.md",
        "Overview" => "overview.md",
        "States" => "States.md",
        "Unitaries" => "Unitaries.md",
        "Observables" => "Observables.md",
        "Information" => "Information.md",
        "QMath" => "QMath.md",
    ],
)

deploydocs(;
    repo="github.com/ChitambarLab/QBase.jl",
)
