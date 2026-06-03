# import Pkg # BE CAREFUL TO REMOVE/PUT IN COMMENTS WHEN PUSHING
# Pkg.activate(@__DIR__)

using PeriodicARModels
using Documenter

DocMeta.setdocmeta!(PeriodicARModels, :DocTestSetup, :(using PeriodicARModels); recursive=true)

makedocs(;
    modules=[PeriodicARModels],
    authors="Arnaud Gonin <arnaudcmc@hotmail.com> and David Métivier <david.metivier@inrae.fr>",
    sitename="PeriodicARModels.jl",
    format=Documenter.HTML(;
        canonical="https://dmetivie.github.io/PeriodicARModels.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        # "API" => "API.md",
    ],
    checkdocs=:none
)

deploydocs(;
    repo="github.com/dmetivie/PeriodicARModels.jl",
    devbranch="master",
)

# using LiveServer; # BE CAREFUL TO REMOVE/PUT IN COMMENTS WHEN PUSHING
# serve(dir="docs/build");