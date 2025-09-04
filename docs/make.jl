push!(LOAD_PATH, "../src/")

using AtomicSymmetries
using Documenter

makedocs(sitename="AtomicSymmetries.jl Documentation", format=[Documenter.HTML(), Documenter.LaTeX()],
         pages = ["Home" => "index.md", 
                  "Symmetries in Q space" => "fourier_symmetries.md"])

