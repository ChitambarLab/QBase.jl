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
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "User Guide" => "user_guide.md",
        "Exports" => "exports.md",
        "Modules" => [
            "States" => "submodules/States.md",
            "Unitaries" => "submodules/Unitaries.md",
            "Observables" => "submodules/Observables.md",
            "Information" => "submodules/Information.md",
            "QMath" => "submodules/QMath.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/ChitambarLab/QBase.jl",
)