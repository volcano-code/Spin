# Independent mathematical fixtures for Task5e

Computed from rectangular-domain integrals and exact Python Fraction arithmetic, without native source/moment/Boolean helpers. These are SYNTHETIC_REGRESSION fixtures, not original paper cases or recovered parameters.

## Actual source-solver/export fixture

Original outer O=[-1,1]^2 x[0,2]; original design D=[-.5,.5]^2 x[1,1.5]; root side2, depth4, cell width1/8. Inclusive boundary exclusion removes cells touching design planes, so accepted carveable domain is A=[-3/8,3/8]^2 x[9/8,11/8]. Its lengths(3/4,3/4,1/4), center(0,0,5/4), volume9/64 and firstz45/256 follow directly. Source finest FULL coverage has72cells, though stored source cuboids may merge them.

For a box of volumeV, centeru and side lengths h, first moment=Vu, second moment Q=V(uu^T+diag(h^2)/12). Thus A has Qxx=Qyy=27/4096,Qzz=903/4096,cross terms0. O has volume8,Lz8,Qxx=Qyy=8/3,Qzz32/3. At full removal the final material is O\A:

| Quantity | Exact value |
|---|---|
| Volume / mass at density1 | 503/64 |
| Firstx,firsty,firstz | 0,0,2003/256 |
| Qxx,Qyy | 32687/12288 |
| Qzz | 128363/12288 |
| Qxy,Qxz,Qyz | 0,0,0 |
| Centerz | 2003/2012 |
| Central Ixx,Iyy | 16432021/3090432 |
| Central Izz | 32687/6144 |
| Minimum original-boundary/cavity separation | 5/8 |

Static objective is c_z^2 plus nonnegative graph regularizer. Every removable point has z>=9/8>original center1. For any beta in[0,1], the remaining center is positive and <=1; its derivative w.r.t. a cell's removed fraction is v(c_z-z_cell)/m<0. Full removal is feasible for x/y moment constraints by symmetry, minimizes center over the whole box and has zero Laplacian energy at common beta1. Consequently full removal is the unique global optimum for the fixed positive center weight, independent of which admissible partition represents A. This mathematical proof is not a substitute for actual audited native convergence.

Exact clearance derives from original faces: horizontal side gaps1-3/8=5/8, top gap2-11/8=5/8, bottom9/8. Therefore threshold5/8 passes and nextafter(5/8,+infinity) fails with unchanged geometry. Requested1/2 passes. This is outer-to-cavity clearance only, not all internal webs/accessibility/manufacturability.

## Nonzero-cross-moment translation fixture

Translate the same final geometry by t=(1/8,-1/4,3/8). With m=503/64, use L'=L+mt and Q'=Q+tL^T+Lt^T+mtt^T, derived by expanding(x+t)(x+t)^T. Exact translated values:

- L'=(503/512,-503/256,5515/512).
- Q'xx=8549/3072,Q'yy=38723/12288,Q'zz=53513/3072.
- Q'xy=-503/2048,Q'xz=5515/4096,Q'yz=-5515/2048.
- Central inertia is unchanged by translation; density2 doubles m,L,Q,I but not center or clearance.

This translation is a direct geometry/integration test, not automatic target-axis rotation or a new static-success claim.

## Hollow26cube union

Union26unitcubes in[0,3]^3 with[1,2]^3 omitted is one connected shell enclosing one void, not a filled27cube block. Subtract the central unit cube from the3cube box: V26, center(3/2,3/2,3/2), each first moment39, each diagonal raw Q236/3, each cross Q117/2. Each central Q=27*(3^2/12)-1*(1/12)=121/6; central inertia diag121/3 and offdiag0. The inner boundary must have negative signed orientation in the union; a source material formed by removing this union retains an island inside, so orientation must follow signed-set nesting rather than treating every shell as an independent outward solid.
