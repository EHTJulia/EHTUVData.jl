using EHTUVData
using Documenter

DocMeta.setdocmeta!(EHTUVData, :DocTestSetup, :(using EHTUVData); recursive=true)

makedocs(;
    modules=[EHTUVData],
    authors="Kazu Akiyama",
    repo="https://github.com/kazuakiyama/EHTUVData.jl/blob/{commit}{path}#{line}",
    sitename="EHTUVData.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kazuakiyama.github.io/EHTUVData.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kazuakiyama/EHTUVData.jl",
    devbranch="main",
)
