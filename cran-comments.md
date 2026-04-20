## R CMD check results

0 errors | 0 warnings | 1 note

  checking for future file timestamps: unable to verify current time

This note is caused by a local network/firewall issue preventing R from
reaching the IANA time server. It does not appear on win-builder, R-hub,
or CRAN servers and is unrelated to the package itself.

## win-builder (R-devel) previous run fixed issues

- LaTeX Unicode errors (pi, phi, lambda): fixed, replaced with ASCII
- cran-comments.md in tarball: fixed via .Rbuildignore
- Spell-check acronyms (SEIR, SEIRD, SVEIRD, EpiNova): added inst/WORDLIST

## Test suite

76 tests, 0 failures (testthat 3.x).

## New submission

No downstream dependencies.
