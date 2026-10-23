using ArchimedLight
using Documenter
using Bonito

# For interactive rendering:
# using LiveServer
# servedocs()
# or LiveServer.serve(dir=joinpath(@__DIR__, "build"))

DocMeta.setdocmeta!(ArchimedLight, :DocTestSetup, :(using ArchimedLight); recursive=true)

# Update the videos:
# include(joinpath(@__DIR__, "make_video", "make_video.jl"))

# Documenter still evaluates examples, checks references and deploys the site.
# Bonito supplies its static HTML writer and documentation theme.
home = (
    name="ArchimedLight.jl",
    text="Let there be light.",
    tagline="Radiation interception, absorption, and scattering in 3D scenes.",
    image="assets/archimedlight-logo.png",
    actions=[
        (text="Run your first simulation", link="getting_started.html", theme="brand"),
        (text="Explore the model", link="theory_pipeline.html", theme="alt"),
    ],
)

makedocs(;
    modules=[ArchimedLight],
    authors="Rémi Vezy <VEZY@users.noreply.github.com> and contributors",
    repo=Documenter.Remotes.GitHub("VEZY", "ArchimedLight.jl"),
    sitename="ArchimedLight.jl",
    format=Bonito.DocumenterBonito(;
        repo="github.com/VEZY/ArchimedLight.jl",
        devbranch="main",
        devurl="dev",
        version=get(ENV, "GITHUB_REF_TYPE", "") == "tag" ? get(ENV, "GITHUB_REF_NAME", "dev") : "dev",
        home,
        description="Light interception and scattering in three-dimensional plant scenes, computed in Julia.",
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Applications" => [
            "Overview" => "applications/index.md",
            "Wheat–Maize Intercropping" => "applications/intercropping.md",
            "Agrivoltaics" => "applications/agrivoltaics.md",
            "Reconstructed Rice" => "applications/rice.md",
            "Oil-Palm Digital Twins" => "applications/oil-palm.md",
        ],
        "Tutorials" => [
            "Beginner Workflows" => "tutorial_beginner_workflows.md",
            "File-Based Workflow" => "tutorial_files.md",
            "Interactive Workflow" => "tutorial_interactive.md",
            "Copy-Paste Templates" => "onboarding_templates.md",
        ],
        "Reference" => [
            "Configuration And Options" => "reference_config.md",
            "Scene Files And Semantics" => "reference_scene.md",
            "Model Files" => "reference_models.md",
            "Meteo Inputs" => "reference_meteo.md",
        ],
        "Outputs" => "outputs.md",
        "Evaluation" => "evaluation.md",
        "Performance Benchmarks" => "performance_benchmarks.md",
        "Under The Hood" => [
            "Pipeline Overview" => "theory_pipeline.md",
            "First-Order Interception" => "theory_interception.md",
            "Scattering And Assumptions" => "theory_scattering.md",
        ],
        "Advanced Usage" => [
            "PlantSimEngine Coupling" => "plantsimengine.md",
            "File-Based Example" => "full_example.md",
            "Composable Stages" => "stages.md",
            "GPU Backends And Benchmarks" => "gpu_backends.md",
            "Historical ARCHIMED Reference" => "archimed_reference.md",
        ],
        "API Reference" => [
            "Public API" => "api.md",
            "Advanced API" => "advanced_api.md",
        ],
    ],
)

include(joinpath(@__DIR__, "check_static_export.jl"))
finish_static_export()

if get(ENV, "ARCHIMEDLIGHT_DOCS_BUILD_ONLY", "false") != "true"
    deploydocs(;
        repo="github.com/VEZY/ArchimedLight.jl.git",
        devbranch="main",
        push_preview=true,
    )
end
