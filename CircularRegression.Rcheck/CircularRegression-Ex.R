pkgname <- "CircularRegression"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('CircularRegression')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("angular_re")
### * angular_re

flush(stderr()); flush(stdout())

### Name: angular_re
### Title: Angular mixed-effects regression with von Mises random
###   intercepts
### Aliases: angular_re

### ** Examples

# See README example below with Sandhopper data



cleanEx()
nameEx("bison")
### * bison

flush(stderr()); flush(stdout())

### Name: bison
### Title: Bison movement data
### Aliases: bison
### Keywords: datasets

### ** Examples

# View structure
str(bison)





cleanEx()
nameEx("multiplebison")
### * multiplebison

flush(stderr()); flush(stdout())

### Name: multiplebison
### Title: Two plains bison tracked by GPS (Prince Albert National Park,
###   2013)
### Aliases: multiplebison
### Keywords: datasets

### ** Examples

## Basic glimpse
head(multiplebison)


## Example: inter-animal proximity rate (< 500 m)
# mean(multiplebison$Distance < 500, na.rm = TRUE)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
