# Independent Task5h acceptance derivation (before implementation report)

The new fine grid is the old cage union subdivided, not a new fine-level source-only envelope. For original [.25,.75]^3, unit root/coarse level3, closed-contact coarse indices1..6: fine level4 selected indices2..13 inclusive per axis. Thus12^3=1728 cells,13^3=2197nodes and3*12*13^2=6084edges. Exactly218 original coarse boundary nodes are handles, all map to even fine coordinates without nearest-node lookup. The envelope volume is unchanged27/64, NOT material volume.

For query(.3,.45,.6), fine cell key(4,7,9), local fractions(4/5,1/5,3/5). The eight x-bit-first trilinear coefficients independently are [.064,.256,.016,.064,.096,.384,.024,.096]. At(.5,.45,.6), closed lex-first cell is(7,7,9); at(.5,.5,.5), cell(7,7,7) and only corner7 is1. At(.28125,.34375,.40625), key(4,5,6), each coefficient1/8. No lookup or interpolation helper is used to derive these literals.

Analytic fields[1,x,y,z,xy,xz,yz,xyz] at(.3,.45,.6) equal[1,.3,.45,.6,.135,.18,.27,.081]. Fields are arbitrary values, not solved BBW. In particular a constant field2 must remain2; silently normalizing it would be wrong. Testing polynomial reconstruction validates trilinear field transfer only, NOT linear precision of BBW handle-coordinate weights (a separate mathematical property).

Sparse coarse level9 singlecell(100,100,100), fine10 has8cells/27nodes/54edges, original8handles map to corners200/202 eachaxis. A deeper level with bounded local cell count likewise must not allocate a full root lattice. Tests label this geometry/interpolation, not optimization.

Bounded hollow-shell test fine16-grid: the coarse cage envelope is unitroot minus(.375,.625)^3. Refinement retains the identical open cavity; its center(.5,.5,.5) must reject. Through-hole center likewise rejects; disjoint-island gap rejects. These tests must not falsely treat source cavities narrower than the chosen coarse grid as guaranteed preserved: only the actual validated cage union is inherited.

Outside nextafter query and exact boundary query must distinguish strictly outside from inside without tolerance expansion. Non-dyadic represented planes require FMA and local actual widths; tests compare affine/trilinear values rather than assuming exact decimal arithmetic.
