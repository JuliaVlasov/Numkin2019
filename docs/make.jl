import Remark
using Literate

files =  filter( f -> startswith(f, "0"), readdir("src")) |> collect

for file in files
    Literate.notebook("src/$file", "notebooks",  execute=false)
end

run(pipeline(`cat src/$files`; stdout="slides.jl" ))

slideshowdir = Remark.slideshow("slides.jl", "docs")
