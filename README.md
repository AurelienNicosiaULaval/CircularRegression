

**CircularRegression** is an R package for performing regression analysis when the **response variable is circular** (e.g., angles, directions, time-of-day measured in radians). It provides estimation methods and visualization tools adapted to this unique data structure.

## ✨ Main Features

- `angular()` — fits a circular regression model
- `consensus()` — fits a circular regression model based on consensus errors
- `summary()`, `plot()`, `coef()`, `residuals()` — S3 methods for fitted models
- `data(bison)` — an example dataset for demonstration

## 📦 Installation

To install the development version from GitHub:

```r
# Install devtools if needed
install.packages("devtools")

# Install from GitHub
devtools::install_github("AurelienNicosiaULaval/CircularRegression")
 
