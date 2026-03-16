## Resubmission kuenm2 0.1.1

This is a resubmission. All issues raised by the reviewer have been addressed.

> Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
For more details:
<https://contributor.r-project.org/cran-cookbook/docs_issues.html#missing-value-tags-in-.rd-files>
Missing Rd-tags:
      plot_calibration_hist.Rd: \value
      plot_explore_partition.Rd: \value
      plot_importance.Rd: \value

* Addes \value to .Rd files: plot_calibration_hist.Rd, plot_explore_partition.Rd 
and plot_importance.Rd


## Submission kuenm2 0.1.0
This is the first submission of version 0.1.0.

## Test environments
* Windows 11, R 4.5.1 (local)
* MacOS 15.7.3, R release (GitHub Actions)
* Windows 10.0.26100, R release (GitHub Actions)
* Ubuntu 24.04.3 LTS, R release (GitHub Actions)
* Ubuntu 24.04.3 LTS, R devel (GitHub Actions)
* Ubuntu 24.04.3 LTS, R oldrel-1 (GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Weverton Trindade <wevertonf1993@gmail.com>'
A Note that reminds CRAN maintainers to check that the submission comes actually from his maintainer and not anybody else.

* This is a new release.
