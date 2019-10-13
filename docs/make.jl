import Remark
using Literate

files =  filter( f -> startswith(f, "0"), readdir("src")) |> collect

# run(pipeline(`cat src/$files`; stdout="slides.jl" ))

for file in files
    Literate.notebook("src/$file", "notebooks",  execute=true)
    slides_path = joinpath("docs",file[1:2])
    mkpath(slides_path)
    s = Remark.slideshow("src/$file", slides_path)
end

