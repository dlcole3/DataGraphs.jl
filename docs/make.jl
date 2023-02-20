using Documenter, DataGraphs

const _PAGES = [
    "Introduction" => "index.md",
    "Quick Start"=>"guide.md",
    "API Manual" => "api.md"
]

makedocs(
    sitename = "DataGraphs",
    authors = "David Cole",
    format = Documenter.LaTeX(platform="docker"),
    pages = _PAGES
)

makedocs(
    sitename = "DataGraphs",
    modules = [DataGraphs],
    authors = "David Cole",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        sidebar_sitename = true,
        collapselevel = 1,
    ),
    pages = _PAGES,
    clean = false,
)

deploydocs(
    repo = "github.com/dlcole3/DataGraphs.jl.git",
    target = "build",
    devbranch = "main",
    devurl = "dev",
    push_preview = true,
)
