# Chromo counting R package

To build and install from scratch do:
- clone repos (git clone URL chromocounting)
- assuming the package directory is chromocounting do in R

```
library(devtools)
document('chromocounting')
check('chromocounting')
build('chromocounting')
install.packages('chromocounting_0.1.tar.gz', repos=NULL)
```
