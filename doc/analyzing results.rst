Analyzing Results
=================

The output directory contains one file
for each unique forcefield term.
The file has a header line explaining the nature
of the term, and is followed by a list of coefficients.

Parsers which make this output more readable are
available in the https://github.com/frobnitzem/chemparam
package.

Energy-spline files
-------------------

Energy-splines are named as `<type>.rname[.+-]aname.espl`.
They contain one line for each b-spline parameter, and
can be visualized directly to get a rough idea of the
shape of the function, but should be handled with
the pspline utilities to get the "correct" interpolated
function.
answers
