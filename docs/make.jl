push!(LOAD_PATH, "../src/")

using AtomicSymmetries
using Documenter

makedocs(sitename="AtomicSymmetries.jl Documentation", format=[Documenter.HTML()],
         pages = ["Home" => "index.md", 
                  "Acoustic Sum Rule" => "acoustic_sum_rule.md",
                  "Symmetries in Q space" => "fourier_symmetries.md"],
         repo = "https://github.com/mesonepigreco/AtomicSymmetries.jl/blob/{commit}{path}#{line}"
        )

if get(ENV, "DOCS_DEPLOY", "false") == "true"
    # deploydocs will detect GitHub Actions and push with GITHUB_TOKEN or use DOCUMENTER_KEY if set
    deploydocs(repo="github.com/mesonepigreco/AtomicSymmetries.jl", devbranch="main")
end
