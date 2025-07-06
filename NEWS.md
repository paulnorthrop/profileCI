# profileCI 1.1.0

## New features

* The argument `optim_args` actually works now! Thank you to Jakob Madsen (#1).
* An argument `flat` has been added to `profileCI()` to deal with cases where one or both of the limits is (minus or plus) infinity because the profile log-likelihood fails to drop below the level that defines the confidence interval.

## Bug fixes and minor improvements

* Arguments to the user-supplied log-likelihood function are extracted explicitly to simplify the code and avoid potential argument matching issues.

