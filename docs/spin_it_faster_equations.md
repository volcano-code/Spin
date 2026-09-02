# Spin-It Faster equation contract

This document is an implementation contract derived from Hafner, Ly, and
Wojtan (2024), DOI `10.1145/3658194`, and its supplemental material.  It is a
formula reconstruction, not a claim of equivalence with the unreleased Matlab
code or author CAD assets.

## Coordinates and stored moments

Optimization occurs in an approximately unit-scale right-handed frame.  The
target spin/gravity axis is `e3`; `e1,e2` span its orthogonal plane.  For a
material characteristic function `chi`, the raw-moment order is

```text
(v, l1, l2, l3, Q11, Q22, Q33, Q23, Q13, Q12)
```

with `v=integral chi`, `l=integral chi*x`, and
`Q=integral chi*x*x^T`.  For constant density `rho`,

```text
m = rho*v
c = l/v
Io = rho*(trace(Q)*Identity - Q)
Ic = Io - m*((c dot c)*Identity - c*c^T)
```

Fixed-shell raw moments are computed once and added to the variable-domain
moments before evaluating objectives and constraints.

## Quadric material distribution

The degree-two function and its coefficient order are

```text
beta(x) = a + b dot x + x^T*A*x
(a, b1, b2, b3, A11, A22, A33, A23, A13, A12)
```

`A` is symmetric; stored cross coefficients are matrix entries, so evaluating
`x^T*A*x` produces the corresponding factor of two.  Material is present
exactly where `beta(x) > 0`.  Positive rescaling does not change this set, so
the coefficient vector `p` is normalized to `||p||=1` on `S^9`.

## Application problems

- Static, maximum mass: minimize `-v` subject to `l1=0,l2=0`.  Its optimal
  interface is a vertical plane.
- Static, lower center of mass: minimize `c3` subject to `l1=0,l2=0`.  Its
  optimal interface is a non-vertical plane through the center of mass, with
  the lower halfspace solid.
- Fidget spinner in 2D: minimize `-(Q11+Q22)` subject to `l1=0,l2=0`.  Its
  empty region is a disk bounded by a circle through the target origin.
- Yo-yo: minimize `((Ic11)^2+(Ic22)^2+2*(Ic12)^2)/(Ic33)^2`, subject to
  `l1=l2=l3=0` and `Ic13=Ic23=0`.
- Spinning top: minimize the yo-yo objective plus `wc*c3`, subject to
  `l1=l2=0` and `Ic13=Ic23=0`.

All constraint coordinates are evaluated after transforming the combined
fixed and variable moments into the target frame.

## Optimality system and update

For equality constraints `g(r)=0`, solve the paper's Eq. 14:

```text
grad_r(f) + mu*p - Jac_r(g)^T*lambda = 0
g(r(p)) = 0
||p|| = 1
```

At every iterate, recompute `(lambda,mu)` by the least-squares projection of
`grad_r(f)` onto the span of `{-p, grad_r(g1),...,grad_r(gk)}`.  The Newton
direction is constrained to the tangent space of the coefficient sphere.  For
tangent direction `Delta`, use the geodesic update

```text
p_next = cos(t)*p + sin(t)*Delta/||Delta||.
```

Initialize `t=min(||Delta||,t_max)` and halve it until the complete residual
decreases.  The locked paper values are `alpha = 1/2` and `t_max = 1/2`.
Terminate only when the normalized Eq. 14 L2 residual is at most `1e-12`.
Restart after at most 60 function evaluations, reject a starting sample whose
moment Jacobian vanishes, and accept a converged point as a strict local
minimum only when `mu>0`.

## Quadric discretization and moment evaluation

For a nonsingular quadric, the supplement maps a fixed unit-icosphere normal
`nu` through the inverse Gauss map.  With
`beta(x)=x^T*A*x+b^T*x+a`, compute

```text
w^2 = 4*(nu^T*A^-1*nu)/(b^T*A^-1*b - 4*a)
x = A^-1*(nu/w - b/2).
```

A negative discriminant has no Euclidean/ideal preimage.  `w=0` is an ideal
point.  For hyperboloids, intersect icosphere edges with
`nu^T*A^-1*nu=0`, stitch the valid-normal region, initialize truncation beyond
the design-domain diameter, double the radius until the required two closed
loops exist, and cap beyond the domain.

The material domain is the exact mesh Boolean
`Omega intersection {beta > 0}`.  Raw moments are integrated from the closed
boundary after one divergence-theorem reduction using a triangle rule exact
for the required polynomial degree.  Voxel sampling is not an allowed
fallback.  Independently recomputed physical constraints must have relative
residual at most `1e-6`.
