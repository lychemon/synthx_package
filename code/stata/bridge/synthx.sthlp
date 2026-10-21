{smcl}
{* *! version 0.9.1  19oct2026}{...}
{vieweralsosee "synth" "help synth"}{...}
{viewerjumpto "Syntax" "synthx##syntax"}{...}
{viewerjumpto "Description" "synthx##description"}{...}
{viewerjumpto "Options" "synthx##options"}{...}
{viewerjumpto "Methods" "synthx##methods"}{...}
{viewerjumpto "Limitations" "synthx##limits"}{...}
{viewerjumpto "Examples" "synthx##examples"}{...}
{viewerjumpto "Stored results" "synthx##results"}{...}
{viewerjumpto "References" "synthx##refs"}{...}
{title:Title}

{phang}
{bf:synthx} {hline 2} Exact synthetic control: the global optimum over a class of predictor weight matrices


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:synthx} {it:depvar} {it:predictors} [{it:var}{cmd:(}{it:periods}{cmd:)} ...]{cmd:,}
{cmdab:tru:nit(}{it:#}{cmd:)}
{cmdab:trp:eriod(}{it:#}{cmd:)}
{cmdab:pre:window(}{it:# #}{cmd:)}
[{it:options}]

{pstd}
The data must be {helpb xtset} as a panel, or supply {cmd:unit()} and {cmd:time()}.
A predictor written {it:var}{cmd:(}{it:periods}{cmd:)} is averaged over exactly
{it:periods} (a {help numlist}); plain predictors are averaged over {cmd:xperiod()}.

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{p2coldent:* {cmdab:tru:nit(}{it:#}{cmd:)}}panel id of the treated unit{p_end}
{p2coldent:* {cmdab:trp:eriod(}{it:#}{cmd:)}}first treated period{p_end}
{p2coldent:* {cmdab:pre:window(}{it:# #}{cmd:)}}validation window {it:lo hi}: the periods whose
outcomes choose the predictor weights{p_end}
{synopt:{cmdab:xp:eriod(}{it:# #}{cmd:)}}window over which plain predictors are averaged;
default {cmd:prewindow()}{p_end}
{synopt:{cmdab:u:nit(}{it:varname}{cmd:)}}panel id variable; default from {cmd:xtset}{p_end}
{synopt:{cmdab:t:ime(}{it:varname}{cmd:)}}time variable; default from {cmd:xtset}{p_end}
{syntab:Metric}
{synopt:{cmdab:v:class(}{it:class}{cmd:)}}{cmd:diagonal} (default): positive diagonal {it:V};
{cmd:pd}: any positive definite {it:V}; {cmd:invvar}, {cmd:equal}: fixed inverse-variance or identity {it:V}, no optimization{p_end}
{synopt:{cmdab:vw:eights}}report the predictor weight matrices that rationalize the
estimated donor weights: {cmd:e(Vdiag)}, {cmd:e(Vpd)}{p_end}
{synopt:{cmdab:fac:es}}count the faces of the donor hull by dimension (attainable,
visible, all): {cmd:e(faces)}{p_end}
{syntab:Inference and figures}
{synopt:{cmdab:allp:lacebos}}in-space placebo test: refit every unit and rank the treated
unit by its post/pre-treatment MSPE ratio{p_end}
{synopt:{cmdab:post:window(}{it:# #}{cmd:)}}post-treatment window for the MSPE ratio;
default [{cmd:trperiod}, last period]{p_end}
{synopt:{cmd:figure}}draw the treated-vs-synthetic path, or (with {cmd:allplacebos})
the placebo gaps{p_end}
{synopt:{cmdab:treat:line(}{it:#}{cmd:)}}period of the vertical line in the figure; default {cmd:trperiod()}{p_end}
{synopt:{cmdab:maxp:re(}{it:#}{cmd:)}}placebo figure: omit units whose pre-treatment RMSE exceeds {it:#} times the
treated unit's; the test itself uses all units{p_end}
{synopt:{cmd:title(}{it:string}{cmd:)} {cmd:xtitle()} {cmd:ytitle()} {cmd:legend()}}replace the figure's default title, axis
titles and legend{p_end}
{synopt:{cmd:graphopts(}{it:options}{cmd:)}}any other {helpb twoway_options:twoway options} for the figure (scheme, labels,
ranges, {cmd:name()}){p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* required.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:synthx} minimizes validation loss over the donor weights allowed by a
positive diagonal or full positive definite predictor metric. The exterior solver
searches closed regions of visible faces using the same routines as the paper's
Python replication. The interior solver searches attainable cells containing the
treated point. {cmd:e(regime)} records which solver was used.

{pstd}
The command requires Stata 16.1 or later and Python 3.12 or later with NumPy,
SciPy, CVXPY, and Clarabel. Install {cmd:code/python/requirements.txt} into the
interpreter selected by {cmd:python set exec}. Run {cmd:synthx_checkdeps} to check
versions, imports, and a small Clarabel QP. Find the bridge through {cmd:adopath}
or set {cmd:global synthx_bridge "<package>/code/stata/bridge"}.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{cmd:trunit()} identifies the treated unit. {cmd:trperiod()} sets the default start
of the placebo post-period and the figure's vertical line. It does not enter
estimation: weights depend on the predictors and validation outcomes.

{phang}
{cmd:prewindow(lo hi)} supplies the validation outcomes. Plain predictors are
averaged over {cmd:xperiod(lo hi)}, which defaults to {cmd:prewindow()}.
Special predictors, such as {cmd:beer(1984/1988)} or {cmd:schooling(1970 1975)},
use only the listed periods. A leakage-free design keeps predictor and validation
windows separate, as in {cmd:xperiod(1971 1980) prewindow(1981 1990)} for Germany.

{dlgtab:Metric}

{phang}
{cmd:vclass(diagonal)} optimizes over positive diagonal metrics; this is the default.
{cmd:vclass(pd)} uses full positive definite metrics. {cmd:vclass(invvar)} fixes
the inverse-variance diagonal metric and {cmd:vclass(equal)} fixes the identity
metric. Fixed metrics retain their meaning on the original predictor scale.

{phang}
{cmd:vweights} reports predictor matrices that reproduce the fitted donor weights.
Predictors are standardized across treated and donor units, and matrices have
trace one. For each class, the table shows a log-determinant maximizer and an
analytic center. Columns 3 and 4 of {cmd:e(Vdiag)} and {cmd:e(Vpd)} store each
diagonal entry's minimum and maximum. A diagonal range of [0, 0] identifies a
predictor forced to zero in the closed set. The analytic center adds log barriers
only for inactive-donor slacks that can be positive. Exterior fits only.

{phang}
{cmd:faces} counts faces with feasible closed diagonal regions, faces below visible
facets, and all hull faces, by dimension. It checks every visible face and can be
slow. This count is independent of the branch-and-bound optimization. Exterior fits only.

{dlgtab:Inference and figures}

{phang}
{cmd:allplacebos} refits every unit. Each ratio is post-treatment MSPE divided by
pre-treatment MSPE. The treated rank is one plus the number of ratios exceeding
its ratio; {cmd:e(pvalue)} is rank divided by the number of successful fits,
including the treated unit. Failed placebo fits are reported and excluded.

{phang}
{cmd:postwindow(lo hi)} sets the periods used in the ratio's numerator. The default
runs from {cmd:trperiod()} through the last panel period. Set it when the panel
extends beyond the study horizon.

{phang}
{cmd:figure} draws treated and synthetic paths ({cmd:esc_path}), or placebo gaps
with {cmd:allplacebos} ({cmd:esc_placebos}). {cmd:treatline(#)} moves the vertical
line from its default at {cmd:trperiod()}. {cmd:maxpre(#)} hides placebo paths whose
pre-treatment RMSE exceeds this multiple of the treated RMSE; it does not change
the test. Use {cmd:graph export} to save the figure.

{phang}
{cmd:title()}, {cmd:xtitle()}, {cmd:ytitle()}, and {cmd:legend()} replace defaults.
{cmd:graphopts()} passes other {helpb twoway_options:twoway options}, for example
{cmd:graphopts(scheme(s1mono) xlabel(1970(5)2000))}. For custom figures, use
{cmd:e(gaps)} or build outcome paths from {cmd:e(b)} with {cmd:synthx_series}
in {cmd:code/synthx_tools.do}.


{marker methods}{...}
{title:Methods}

{pstd}
For a metric V, donor weights minimize predictor loss over the simplex.
The outer problem chooses weights that minimize validation loss. For exterior
problems, the diagonal solver searches closed regions defined by weak residual
sign constraints; the PD solver uses the visible facets. Lower bounds allow the
search to skip branches that cannot improve the current loss. QPs and
predictor-weight problems use CVXPY/Clarabel; LPs use SciPy/HiGHS.


{marker limits}{...}
{title:Limitations}

{pstd}
The exterior search relies on the paper's general-position assumption.
{cmd:e(exact)} reports checks of visible facets and the selected face: 1 means
these local checks passed and 0 means they failed. Passing does not verify the
full assumption for the data. The result is absent when the checks do not apply.

{pstd}
Optimized classes require more donors than predictors. The diagonal class drops
constant donor predictors with a note; PD retains them. Near collinearity produces
a note. An exactly lower-dimensional donor hull stops
the diagonal solver. For PD, the solver reduces to the donor span when the
treated point is in that span; outside the span, it minimizes outcome loss over
the full simplex. A near-collinearity fallback to the inverse-variance metric is reported.
Treated points exactly on the donor-hull boundary are unsupported by the current
interior/exterior routing.

{pstd}
Missing predictors or validation outcomes exclude a donor. Small donor weights
are set to zero and renormalized. The interior search raises an error after
200,000 cells. {cmd:vweights} and {cmd:faces} are available for exterior fits;
reports that cannot be computed are identified in the Results window.


{marker examples}{...}
{title:Examples}

{pstd}California, Proposition 99 (Abadie, Diamond, and Hainmueller 2010): all
pre-treatment years as the validation sample{p_end}
{phang2}{cmd:. use data/smoking.dta}{p_end}
{phang2}{cmd:. synthx cigsale lnincome age15to24 retprice beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), trunit(3) trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) vweights}{p_end}

{pstd}The same under the full positive definite class{p_end}
{phang2}{cmd:. synthx cigsale lnincome age15to24 retprice beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), trunit(3) trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) vclass(pd)}{p_end}

{pstd}West Germany, leakage-free design: predictors from 1971-1980, validation
outcomes 1981-1990{p_end}
{phang2}{cmd:. use data/germany.dta}{p_end}
{phang2}{cmd:. synthx gdp gdp trade infrate industry schooling(1970 1975) invest70(1980), trunit(17) trperiod(1990) prewindow(1981 1990) xperiod(1971 1980) vweights}{p_end}

{pstd}Placebo test with a capped post-treatment window and a figure whose vertical line
marks the passage of Proposition 99 in 1988{p_end}
{phang2}{cmd:. synthx cigsale lnincome age15to24 retprice beer(1984/1988) cigsale(1988) cigsale(1980) cigsale(1975), trunit(3) trperiod(1989) prewindow(1970 1988) xperiod(1980 1988) allplacebos postwindow(1989 2000) figure treatline(1988) maxpre(1.5)}{p_end}

{pstd}A complete worked example is {cmd:code/example_california.do} in the package.{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:synthx} is an {help estimation command:e-class} command. It stores:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(loss)}}minimized validation loss{p_end}
{synopt:{cmd:e(slack)}}hull-test slack for a single fit; positive in the strict interior{p_end}
{synopt:{cmd:e(exact)}}1 if local geometry checks pass, 0 if not; absent when not applicable (see Limitations){p_end}
{synopt:{cmd:e(pd_mineig)}}smallest eigenvalue of the log-det maximizer ({cmd:vweights}){p_end}
{synopt:{cmd:e(pvalue)}}placebo p-value = rank/{it:N} ({cmd:allplacebos}){p_end}
{synopt:{cmd:e(rank)}}rank of the treated unit by the MSPE ratio ({cmd:allplacebos}){p_end}
{synopt:{cmd:e(n)}}number of units in the placebo test ({cmd:allplacebos}){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:synthx}{p_end}
{synopt:{cmd:e(regime)}}{cmd:exterior} or {cmd:interior}{p_end}
{synopt:{cmd:e(vclass)}}{cmd:diagonal}, {cmd:pd}, {cmd:invvar} or {cmd:equal}{p_end}
{synopt:{cmd:e(plabels)}}predictor labels, semicolon-separated, in row order of
{cmd:e(balance)}, {cmd:e(Vdiag)}, {cmd:e(Vpd)}{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}donor weights (1 x {it:J}, columns named by donor id){p_end}
{synopt:{cmd:e(balance)}}predictor balance, {it:K} x 3: treated / synthetic / donor average{p_end}
{synopt:{cmd:e(Vdiag)}}{it:K} x 4: log-det maximizer, analytic center, min, max over the diagonal rationalizers ({cmd:vweights}){p_end}
{synopt:{cmd:e(Vpd)}}{it:K} x 4: log-det maximizer, analytic center, min, max over the
positive definite rationalizers ({cmd:vweights}){p_end}
{synopt:{cmd:e(faces)}}{it:K} x 3, row {it:d}+1 = dimension {it:d}: attainable / visible / all ({cmd:faces}){p_end}
{synopt:{cmd:e(placebos)}}per unit: id / pre-RMSE / post-RMSE / ratio ({cmd:allplacebos}){p_end}
{synopt:{cmd:e(gaps)}}period column followed by each successful unit's gap path ({cmd:allplacebos}){p_end}
{p2colreset}{...}


{marker refs}{...}
{title:References}

{phang}
Abadie, A., and Y. Shang. 2026. Synthetic control: An exact geometric formulation.
Working paper, MIT.

{phang}
Abadie, A., A. Diamond, and J. Hainmueller. 2010. Synthetic control methods for
comparative case studies: Estimating the effect of California's tobacco control
program. {it:Journal of the American Statistical Association} 105: 493-505.

{phang}
Abadie, A., A. Diamond, and J. Hainmueller. 2015. Comparative politics and the
synthetic control method. {it:American Journal of Political Science} 59: 495-510.


{title:Also see}

{psee}
{helpb synth} (if installed), {helpb xtset}, {helpb python}
{p_end}
