import Remark

slides_dir = "slides"
mkpath(slides_dir)

for file in filter( f -> startswith(f, "0"), readdir())
   slideshowdir = Remark.slideshow(file, "slides")
   Remark.open(slideshowdir)
end

