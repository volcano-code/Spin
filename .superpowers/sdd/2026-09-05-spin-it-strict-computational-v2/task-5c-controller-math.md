# Independent Task5c oracle notes

These are controller hand-integral/integer-grid checks, not a production solver or original-paper experiment. No production function is called. Fixed fixture and boundary convention were specified before implementation.

## Conservative volume enclosure

A connected closed cube not intersecting any source boundary has constant material membership. An interior terminal is therefore entirely material. An exterior terminal contributes nothing. All missing material lies in the union of unresolved terminal cubes, whose interiors are disjoint. Thus `0 <= exact_design_volume - approximate_volume <= unresolved_volume`, subject to represented-number arithmetic error in the final floating sum. This does not bound all ten moments by the same scalar relative tolerance and is not an infeasibility certificate. A boundary test based only on corners/center cannot establish the premise.

## Exact dyadic box moments at level9

Root [0,1]^3; design [1/4,3/4]^3; h=1/512. Because touching the boundary counts as intersection, accepted material is the box [1/4+h,3/4-h]^3. Let a,b be these limits and s=b-a.

`V=s^3`, `L_i=s^2*(b^2-a^2)/2`, `Q_ii=s^2*(b^3-a^3)/3`, `Q_ij=s*((b^2-a^2)/2)^2`.

Exact rational results calculated with Python fractions independently:

| Quantity | Exact value |
|---|---|
| V | 2048383/16777216 |
| each L | 2048383/33554432 |
| each diagonal Q | 435766854271/13194139533312 |
| each cross Q | 2048383/67108864 |
| exact volume minus approximate | 48769/16777216 |
| unresolved volume | 49153/8388608 |

The unresolved region is the width-h closed grid layer on either side of the six box faces: `(1/2+2h)^3-(1/2-2h)^3`. Its computed volume is greater than the actual deficit; it is intentionally conservative.

For a reversed internal cubic cavity [7/16,9/16]^3, its excluded grid neighborhood expands to [7/16-h,9/16+h]^3. Accepted volume becomes 1006223/8388608 (at the same depth). Disconnected boundaries must be detected before membership testing.

## Independent hierarchy count (integer arithmetic, no CGAL)

Direct recursive cell-box classification on integer grid [0,512]^3, faces at 128/384, early full/empty termination, all eight children per boundary internal:

- internal: 262249
- full terminal: 506360
- empty terminal: 542936
- unresolved boundary terminal: 786448
- total nodes: 2097993
- per-level actual visits: 1, 8, 64, 512, 1664, 6272, 24704, 98432, 393344, 1572992.

This was run independently in 1.55s under the controller's Python process. It verifies a specific hierarchy convention, not production runtime or resource consumption. A production representation may compact terminal records, but cannot substitute requested depth for actual traversal. This geometry test alone is not level9 optimization.

## Pinned predicate inspection

Concrete risk reported by implementer: pinned CGAL has Triangle/Bbox rather than the expected Triangle/Iso_cuboid overload. Controller independently read the complete triangle and plane Bbox implementations at CGAL HEAD f7c3c8212b56c0d6dae63787efc99093f4383415. Separation rejects strict interval separation; extrema become kernel points, and projection/plane predicates use kernel FT. With Simple_cartesian<Gmpq>, these predicates operate on exact represented binary64 endpoints, including contact. This source inspection does not by itself verify the new caller/BVH or cell-coordinate construction; task review and actual boundary fixtures still required.

- `Bbox_3_Triangle_3_do_intersect.h` SHA256 6453b54da90e02b320551d8980dd33c061d95bb0291c640872f74112642e7073.
- `Bbox_3_Plane_3_do_intersect.h` SHA256 5226550e602f682eeae7d0283f6168158f24f264e1345e6e9310cf11f7148045.

Both under `libigl/external/cgal/Intersections_3/include/CGAL/internal/Intersections_3/`. No dependency modification.
