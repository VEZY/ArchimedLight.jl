function pixel_resolution_check()
    SC = Main.ScatteringComparison
    results = map(SC.cases()) do case
        @info "Pixel refinement check" scene=case.id
        matrices = SC.transfer_matrices(case.objects; equal_area_rings=16, pixel=.02)
        b = zeros(length(case.objects)); b[1] = 100.
        c = [object.scatter for object in case.objects]
        (; id=case.id, before=SC.solve_reference(matrices.before,b,c),
           after=SC.solve_reference(matrices.after,b,c))
    end
    payload = (; sectors=1024, pixel_m=.02, scenes=results)
    open(joinpath(@__DIR__,"pixel-check.json"),"w") do io
        SC.json(io,payload)
    end
    return "Saved 1024-direction, 20 mm comparison for comparison with results-equal1024.json (10 mm)."
end
