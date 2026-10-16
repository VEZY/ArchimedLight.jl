# Local compiler adapter for WasmTarget 0.5.3 on Julia 1.13.1.
# Compiler.compile! now requires external_linkage; the Wasm closed world does not
# emit native external linkage. Patch only this loaded method, never package files.
# This is build-time code; it is not included in the emitted light module.
using WasmTarget
browser_trimcollect_source = read(joinpath(pkgdir(WasmTarget),"src","codegen","trimcollect.jl"),String)
browser_trimcollect_source = replace(browser_trimcollect_source,
    "CC.compile!(codeinfos, workqueue; invokelatest_queue)" =>
        "CC.compile!(codeinfos, workqueue; invokelatest_queue, external_linkage=false)",
    "CC.compile!(codeinfos, invokelatest_queue; invokelatest_queue)" =>
        "CC.compile!(codeinfos, invokelatest_queue; invokelatest_queue, external_linkage=false)",
    "CC.compile!(fresh_ci, fresh_wq; invokelatest_queue=fresh_ilq)" =>
        "CC.compile!(fresh_ci, fresh_wq; invokelatest_queue=fresh_ilq, external_linkage=false)",
    "CC.compile!(fresh_ci, fresh_ilq; invokelatest_queue=fresh_ilq)" =>
        "CC.compile!(fresh_ci, fresh_ilq; invokelatest_queue=fresh_ilq, external_linkage=false)")
for expression in Meta.parseall(browser_trimcollect_source).args
    if expression isa Expr && expression.head == :function &&
       occursin("collect_closed_world(entries",string(expression.args[1]))
        Core.eval(WasmTarget,expression)
    end
end
