# CircularRegression 0.1.1

* Improve stability of angular and consensus estimators with QR-based updates and better handling of reference-only models.
* Implement observation weights in `consensus()` and remove the deprecated `model` argument across the API.
* Make `angular_re()` usable without workarounds, add Hessian conditioning checks, and provide more robust predictions.
* Refresh documentation and vignette examples, and add an automated test suite covering key model features.
