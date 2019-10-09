using Pkg
Pkg.instantiate()
using Conda
Conda.add("jupytext", channel="conda-forge")
Conda.add("r-plot3d", channel="conda-forge")
