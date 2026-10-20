# Append the same plain-JavaScript/CSS shell used for the vendored artifact.
# It only handles responsive presentation and the parent readiness handshake.
function browser_attach_shell!(html_path)
    html = read(html_path, String)
    html = replace(html, r"(?s)\n<!-- ArchimedLight wheat shell -->.*?<!-- /ArchimedLight wheat shell -->\n" => "")
    css = read(joinpath(@__DIR__, "shell.css"), String)
    bridge = read(joinpath(@__DIR__, "ready_bridge.js"), String)
    shell = "\n<!-- ArchimedLight wheat shell -->\n<style>" * css *
        "</style>\n<script>" * bridge * "</script>\n<!-- /ArchimedLight wheat shell -->\n"
    @assert occursin("</body>", html)
    write(html_path, replace(html, "</body>" => shell * "</body>"; count=1))
    return html_path
end
