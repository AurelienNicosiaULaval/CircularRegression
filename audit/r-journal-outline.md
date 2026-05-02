# Suggested R Journal article structure

1. Introduction

   - Motivation: circular response data in movement ecology and related fields.
   - Gap: need for an R implementation of the general angular regression model.
   - Package contribution: homogeneous, consensus, two-step and random-effects
     workflows with a consistent user interface.

2. Statistical background

   - Circular response variables and angle wrapping.
   - Resultant-vector construction of the mean direction.
   - Homogeneous angular model.
   - Consensus model.
   - Random-intercept extension, briefly, with references to the original
     methodological papers.

3. Package design

   - Formula syntax: `x` and `x:z`.
   - Main interface: `circular_regression()`.
   - Preserved expert interfaces: `angular()`, `consensus()`,
     `angular_two_step()` and `angular_re()`.
   - S3 object design and generics.

4. Typical workflow

   - Simulated example with known parameters.
   - Fitting, summarising, predicting and checking residuals.
   - Reference-angle selection.

5. Applied example

   - Bison movement data.
   - Model specification and interpretation.
   - Prediction and residual diagnostics.

6. Implementation details

   - Model-frame construction and NA handling.
   - Numerical stability for Bessel functions.
   - Optimisation and convergence checks.
   - Test suite and reproducibility.

7. Comparison and limitations

   - Relationship to general circular-statistics packages.
   - Scope restrictions: fixed-effect interface and separate random-effects
     function.
   - Computational scaling and practical limits.

8. Discussion

   - Intended users.
   - Future work.
   - Reproducibility statement.
