# Task5f independent coordinate/density checks

Hand-derived rectangular integrals and Python Fraction arithmetic only; no native or Python production moment/normalization helper used. SYNTHETIC_REGRESSION, not a paper input.

Fixed Task5e outer O=[-1,1]^2x[0,2], design[-.5,.5]^2x[1,1.5], conservative full carve A=[-3/8,3/8]^2x[9/8,11/8]. The native solver remains in original model coordinates. SpinArtifactProblem normalization is origin o=(0,0,1), length L=2 and reference density R=rho for uniform material.

For density-weighted model moments (m,b,Q) and x'=(x-o)/L, rho'=rho/R:

- m'=m/(R L^3)
- b'=(b-mo)/(R L^4)
- Q'=(Q-ob^T-bo^T+m oo^T)/(R L^5)
- c'=(c-o)/L; central inertia I_c'=I_c/(R L^5).

The geometry is **not** transformed for native solve or exported outer; these formulas adapt only the moment claim into the independently checked normalized frame. Applying translation after density/length scaling must use consistent units. Cross terms need both linear terms, not only m oo^T.

Moment order M,LX,LY,LZ,QXX,QYY,QZZ,QYZ,QXZ,QXY. At rho1, exact native-world tuple is:

`(503/64, 0, 0, 2003/256, 32687/12288, 32687/12288, 128363/12288, 0, 0, 0)`.

Exact normalized tuple:

`(503/512, 0, 0, -9/4096, 32687/393216, 32687/393216, 32651/393216, 0, 0, 0)`.

Normalized center z=-9/4024. Model central inertia diag=(16432021/3090432,16432021/3090432,32687/6144). Normalized central inertia diag=(16432021/98893824,16432021/98893824,32687/196608). All off-diagonal central terms zero for this fixture. Density2 doubles world m,b,Q,I_c, but R also doubles; normalized tuple/inertia and worldcenter/clearance remain identical.

For a separate nonzero-cross-term transformation correctness test, translate same geometry by t=(1/8,-1/4,3/8). This is not an automatically reoriented successful solve. Exact world moments (rho1):

`(503/64,503/512,-503/256,5515/512,8549/3072,38723/12288,53513/3072,-5515/2048,5515/4096,-503/2048)`.

Problem normalization now o=(1/8,-1/4,11/8),L2,R1. The normalized tuple must equal the fixed tuple above, testing all mixed terms and centering. Native target restriction remains unchanged; do not pass this as evidence of arbitrary-target optimization.

Observation-policy note: existing PaperRunEvidenceV2 includes actual absolute output file paths in its determinism field, so two directories can have distinct observation digests. Compare bridge scientific_digest, sealed scientific_digest, complete native trace/content hashes and scientific payload bytes for path-independent scientific determinism. Keep differing provenance/evidence records; do not rewrite generic evidence hashing or delete paths to force equality.

For independent decomposition checks at density1, fixed shell O-minus-D world10 is `(15/2,0,0,59/8,21/8,21/8,79/8,0,0,0)`. Retained design D-minus-A world10 is `(23/64,0,0,115/256,431/12288,431/12288,7019/12288,0,0,0)`. Their sum equals full material tuple. Fixed normalized LZ=-1/128; retained normalized LZ=23/4096, sum=-9/4096. These literals detect omitting fixed shell or multiplying density twice without deriving expectations from the seal itself.
