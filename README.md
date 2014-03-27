## Reproducing the figures in "Discovery of Brainwide Neural-Behavioral Maps via Multiscale Unsupervised Structure Learning" in Science (March 27, 2014) issue.
   
All the plots in the paper were generated with R.
The latest R can be downloaded from [CRAN](http://cran.us.r-project.org).
At the R prompt, type

```
> source("behaviortypes.R")
```

First, all the required R packages will be installed if
necessary. Depending on your platform, some of system files need to be
installed too, for example, `curl`, `jle`, etc.

Once the package installation is done, each figure will be saved into
a file in the same folder where the above R file is located one at a
time. Some of the images may take a while to be generated, e.g.,
heatmaps, etc.

This code has been tested on Mac and Linux
environments. Unfortunately, it was not tested on Windows machines.
