pkgs = c("ggplot2", "rmarkdown", "vegan", "RColorBrewer")
ncores = parallel:detectCores()
install.packages(pkgs, Ncpus = ncores)