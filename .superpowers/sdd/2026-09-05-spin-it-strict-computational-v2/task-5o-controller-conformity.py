"""Exact independent warning fixture for the NEXT deformed-material gate.

This is not a project solver run or evidence that SourceCageBinding deforms
anything. Even an injective trilinear map need not preserve consistency between
a coarse triangulated face and a finer neighboring face evaluated at new points.
"""
from fractions import Fraction as F

alpha = F(1, 16)
# Map (x,y,z) -> (x,y,z+alpha*x*y), det Jacobian identically 1.
# Shared rest face z=0, [0,1]^2. Coarse diagonal endpoints (0,0),(1,1).
coarse_center_z = (F(0) + alpha) / 2
fine_center_z = alpha * F(1,2) * F(1,2)
gap = coarse_center_z - fine_center_z
assert gap == F(1,64) > 0
print("INDEPENDENT_GEOMETRY_COUNTEREXAMPLE (not a solver result)")
print(f"injective_map_determinant=1 coarse_face_center={coarse_center_z} fine_face_center={fine_center_z} mismatch={gap}")
print("Required next gate: conforming shared-face topology and actual deformed partition audit; rest binding alone grants neither.")
