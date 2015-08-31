# Utility to install packages with convenient local options,
# and using a local library.

# This is basically a simple wrapper for install.packages(),
# which always uses a local library, and uses a particular
# repository (rather than asking each time.)
install.packages.1 = function(pkgs)
  install.packages(pkgs, .libPaths()[1],
    repos = "http://lib.stat.cmu.edu/R/CRAN/")

