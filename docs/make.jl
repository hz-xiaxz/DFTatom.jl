using DFTatom
using Documenter

DocMeta.setdocmeta!(DFTatom, :DocTestSetup, :(using DFTatom); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [DFTatom],
    authors = "Xuanzhe Xia, xzxia22@m.fudan.edu.cn",
    repo = "https://github.com/hz-xiaxz/DFTatom.jl/blob/{commit}{path}#{line}",
    sitename = "DFTatom.jl",
    format = Documenter.HTML(; canonical = "https://hz-xiaxz.github.io/DFTatom.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/hz-xiaxz/DFTatom.jl")
