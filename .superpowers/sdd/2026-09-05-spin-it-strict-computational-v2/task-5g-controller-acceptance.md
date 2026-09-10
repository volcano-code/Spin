# Independent Task5g analytical expectations

These expectations are derived from cuboid volumes and integer-grid topology, not SourceSurfaceCage outputs. They do not assume any original paper input is available.

For source[1/4,3/4]^3 in unit root, level3 CLOSED contact selects1..6 in each axis. Envelope[1/8,7/8]^3. Cells216; handles7^3-5^3=218; surface quads6*6^2=216; triangles432. Triangulated indexed boundary has648 undirected edges (Euler218-648+432=2). Interior corners count5^3=125 excluded. Total selected-grid corners343; handles218. These independent counts distinguish volume-grid vertices from boundary handles.

Geometric ten moments in x/y/z order are V=27/64, first each27/128, diagonal second each513/4096, cross each27/256. Central inertia diagonal each81/2048, off-diagonal0. These follow center(1/2,1/2,1/2), side3/4, Qii=V*(1/4+(3/4)^2/12). They apply to the CAGE ENVELOPE, not original material, and must never replace original object's mass.

Level2 envelope is unit cube:64cells,98boundaryhandles,192triangles,288undirectededges, Euler2. Level0 onecell8handles12triangles.

Potential bounded-cavity independent fixture: original outer[1/8,7/8]^3 minus open inner cube(1/4,3/4)^3, level3 unit root. Closed-contact envelope outerunitcube minus cavity(3/8,5/8)^3. Exactly512-8=504 selectedcells; boundaryhandles386+26=412; triangles768+48=816. Two separate boundary components with signs+/-; materialvolume1-1/64=63/64. This fixture uses original shell triangles, not a caller-filled cell list. It is an independent option, not a demand to change worker's own negative/positive fixture.

Potential through-hole variant: original square frame with outerxy[1/8,7/8]^2, holexy(1/4,3/4)^2 through z[1/8,7/8], level3. Coarse envelope leaves2x2 through-hole across8z cells:480cells, volume15/16; one genus1 boundary (Euler0), not capped or filled. Its triangles are768-16+128=880 (remove two end2x2face patches, add four2x8innerwalls). Handles386-2+8*7=440 (remove two endface interiorholepoints, add eightinnerwall perimeterpoints at seveninteriorzplanes). Same envelope for worker's source outer(.1,.9) and hole(.3,.7). Controller flagged worker's initial448/896counts before production GREEN; preserve this independent derivation rather than fitting counts to implementation. Corner/edge-contact singular negative: two closed source boxes strictly contained in diagonally touching coarsecells are separate valid source components but produce a nonmanifold grid union; reject before usablecage.

Real level9 extraction of tiny source strictly inside coarsecell(100,100,100) must keep exactlyonecell with 8corners/12triangles, no512^3 root scan/allocation. A budget failure must report RESOURCE_LIMIT, never fewerhandles/lowerlevel/infeasible.

Final acceptance combines actual focused/native tests, allregistered regression results, feature-isolation test, owned immutable input/output/hashes and independent audit. No screenshot, testcount or cage-envelope moments grant original-paper replication, physical validity, BBW or joint solve.
