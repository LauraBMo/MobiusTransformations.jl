using MobiusTransformations
using Documenter

DocMeta.setdocmeta!(MobiusTransformations, :DocTestSetup, :(using MobiusTransformations); recursive=true)

makedocs(;
    modules=[MobiusTransformations],
    authors="LauBMo <laurea987@gmail.com> and contributors",
    sitename="MobiusTransformations.jl",
    checkdocs=:exports,
    format=Documenter.HTML(;
        canonical="https://LauraBMo.github.io/MobiusTransformations.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Construction" => "construction.md",
        "Evaluation" => "evaluation.md",
        "Operations" => "operations.md",
        "Set Infinity" => "infinity.md",
    ],
)

deploydocs(;
    repo="github.com/LauraBMo/MobiusTransformations.jl",
    devbranch="main",
)
