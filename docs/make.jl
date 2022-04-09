using RetroSignalModel
using RetroSignalModel: RtgMTK
using Documenter

makedocs(;
    modules=[RetroSignalModel],
    authors="stevengogogo <stevengogogo4321@gmail.com> and contributors",
    repo="https://github.com/ntumitolab/RetroSignalModel.jl/blob/{commit}{path}#L{line}",
    sitename="RetroSignalModel.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ntumitolab.github.io/RetroSignalModel.jl",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    branch="gh-pages",
    repo="github.com/stevengogogo/RetroSignalModel.jl.git"
)
