import Remark
using Literate

files =  filter( f -> startswith(f, "0"), readdir("src")) |> collect

for file in files
    Literate.notebook("src/$file", "notebooks",  execute=false)
    slides_path = joinpath("docs",file[1:2])
    mkpath(slides_path)
    s = Remark.slideshow("src/$file", slides_path)
end

# run(pipeline(`cat src/$files`; stdout="slides.jl" ))
