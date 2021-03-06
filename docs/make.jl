using QBase
using Documenter

DocMeta.setdocmeta!(QBase, :DocTestSetup, :(using QBase); recursive=true)

makedocs(;
    modules=[QBase],
    authors="Brian Doolittle <brian.d.doolittle@gmail.com> and contributors",
    repo="https://github.com/ChitambarLab/QBase.jl/blob/{commit}{path}#L{line}",
    sitename="QBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ChitambarLab.github.io/QBase.jl",
        assets=String["assets/custom.css"],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Bras and Kets" => "brakets.md",
        "Operators" => "operators.md",
        "States" => "states.md",
        "Evolution" => "evolution.md",
        "Measurements" => "measurements.md",
        "Probabilities" => "probabilities.md",
        "Information Theory" => "information_theory.md",
        "Math Utilities" => "math_utilities.md",
    ],
)

deploydocs(;
    repo="github.com/ChitambarLab/QBase.jl",
    devbranch = "main",
)
