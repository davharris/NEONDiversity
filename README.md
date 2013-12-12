## NEONDiversity

### About

This package will contain all algorithms used to calculate NEON Data products relating to biodiversity

Currently there are two algorithms implemented

* Hurlbert's PIE
* Individual based rarefaction
* Chao1 estimator

### Quickstart

__Depends__

`devtools`

```coffee
library(devtools)
install_github("NEONDiversity")
```

### Rarefaction


__Matrix example__

```r
library(NEONDiversity)
data(work)
work_rare <- indiv_rare(as.matrix(work))

### Plot results
plot(x = unlist(lapply(work_rare, length)), unlist(lapply(work_rare, max)), 
    ylim = c(0, max(unlist(work_rare))), xlim = c(0, max(unlist(lapply(work_rare, 
        length)))), ylab = "Species", xlab = "Individuals")
for (i in 1:length(work_rare)) {
    lines(work_rare[[i]])
}
```

![plot of chunk unnamed-chunk-1](inst/imagesunnamed-chunk-1.png) 


