# Centered-axis fixed-shell method fixture

This is a separate, pre-run method-validation input. It reuses the exact outer,
editable-domain, and fixed-shell meshes from
`cherchi-rotated-box-fixed-shell-v1`; only the target frame differs. The target
axis is a principal symmetry axis through the mesh origin. This choice checks
whether the generic request path correctly incorporates a nonzero fixed shell
and unequal material densities without introducing an off-center balance
condition. The previously frozen translated-axis request remains unchanged and
is retained as a failed numerical diagnostic; this fixture does not replace or
upgrade that result.

The solver start is generated solely from the request and deterministic seed
0. No terminal coefficients, previous result, or target solution are supplied.
The run is method validation only, not an author asset or paper-case result.
Before optimization, the accepted gates are: current source/build/probe
freshness; exact outer/design/fixed partition audit; fixed-shell moments
independently recomputed from its mesh; and run cap of eight restarts, 60
function evaluations per restart, eight fixed hyperboloid radius doublings,
24 GiB memory, zero swap, and 60 minutes wall time. No setting may change during
the run. A successful outcome must still pass the terminal complement,
independent ten-moment, exact represented-domain, and artifact-sealing checks.

If it fails, retain the full failure ledger and do not search alternative
origins, axes, densities, or seeds within this fixture.
