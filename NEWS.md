# profileCI 1.1.1

## Bug fixes and minor improvements

* When the argument `epsilon` to `profileCI()` is negative, monotonic cubic spline interpolation is now used instead of quadratic interpolation. Quadratic interpolation could fail in some cases.
* In the documentation of the argument `object` to `profileCI()` it is emphasized that the parameters in the fitted model object must have names.

# profileCI 1.1.0

## New features

* An argument `flat` has been added to `profileCI()` to deal with cases where one or both of the limits is (minus or plus) infinity because the profile log-likelihood fails to drop below the level that defines the confidence interval.
* Arguments `lb` and `ub` have been added to `profileCI()` to place respective lower and upper bounds on the interval over which profiling takes place.
* A `logLikFn` method has been added for objects inheriting from class `nls` and an example added to the documentation of `profileCI()`. 

## Bug fixes and minor improvements

* The argument `optim_args` actually works now! Thank you to Jakob Madsen (#1).
* An explanation of the creation of `coef` and `vcov` methods has been added to the description of `object` in the `profileCI()` documentation.
* Arguments to the user-supplied log-likelihood function are extracted explicitly to simplify the code and avoid potential argument matching issues.

