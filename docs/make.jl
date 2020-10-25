using Cardio
using Documenter

makedocs(;
    modules=[Cardio],
    authors="Marius Pille <marius.pille@hu-berlin.de> and contributors",
    repo="https://github.com/mapi1/Cardio.jl/blob/{commit}{path}#L{line}",
    sitename="Cardio.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mapi1.github.io/Cardio.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mapi1/Cardio.jl",
)
