using Documenter
using SchwingerDetectors

makedocs(sitename = "SchwingerDetectors.jl Documentation",
         modules = [SchwingerDetectors],
         pages = [
            "Index" => "index.md",
            "Manual" => ["man/stateprep.md",
                         "man/measure.md",
                         "man/detectors.md",
                         "man/examples.md"],
         ],
         format = Documenter.HTML(prettyurls = false),
         checkdocs = :none)

deploydocs(
    repo = "github.com/srossd/SchwingerDetectors.jl.git",
    devbranch = "main",
)
