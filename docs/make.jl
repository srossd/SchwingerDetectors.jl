using Documenter

# SchwingerDetectors.jl is a loose toolkit of flat files (module `SchwingerToolkit`),
# not a registered package — include it directly so Documenter can resolve the docstrings.
include(joinpath(@__DIR__, "..", "SchwingerToolkit.jl"))
using .SchwingerToolkit

makedocs(sitename = "SchwingerDetectors.jl Documentation",
         modules = [SchwingerToolkit],
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
