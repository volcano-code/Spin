# Task7v1 independent controller acceptance expectations

Do not use renderer-generated expected data or a view screenshot as geometry proof.

## Literal independent fixture facts

Actual Task5f source fixture outer[-1,1]^2x[0,2],design[-1/2,1/2]^2x[1,3/2],carved[-3/8,3/8]^2x[9/8,11/8]. Outer8vertices/12triangles; signed fixed shell16vertices/24triangles; design8/12; cavity8/12; finalmaterial16/24. These follow box topology; retained design16/24. Carved volume9/64; material503/64; center(0,0,2003/2012); outer volume8; fixedshell15/2; retaineddesign23/64. Layer records may split rather than concatenate, but sums/in-memory source bytes must correspond.

Material central inertia diag(16432021/3090432,16432021/3090432,32687/6144). First two values repeated: only z principal-axis line is unique; xy is a2D eigenspace, not two identifiable physical axes. Cube with isotropic inertia has zero unique axes. Non-diagonal analytic unit-density2x3x4box, centered then rotated45degrees aroundz: mass24, central inertia [[45,5,0],[5,45,0],[0,0,26]], unique eigenvalues26,40,50 with axis lines z,(1,-1,0),(1,1,0); signs conventional. No numerical eigensolver is used to derive these independent expectations.

Target native+z/origin0, support_points[(0,0,0)] is only one target point, not a finite-area support polygon. Viewer must not present a bbox-square plane as certified support. Apparent box stability is not solver target feasibility. Geometry seal remains UNKNOWN for support/rotation target and cannot grant numerical convergence from native-run.json.

## Admission checks

- actual selected extension/path; independently load current regenerated bundle with external locked input digest; invaliddigest/oldnativepayload/fakelegacy flag fails beforewindow;
- originalpose/faceordering/signedzero preservation; fixedshell distinct fromouter, finalwholematerial distinct fromretaineddesign;
- display clipping/explode changes rendered copies only; sourcebundle SHAmap unchangedbefore/after;
- mode/status always visible, including degenerateaxes, geometry-only andnotphysical, no title misleadingpaperoriginal;
- actual C++ scene conversion and validation-only command success/failure at frozencompiledHEAD;
- same-camera original/final viewports, visibly inspectable cavity, no OpaqueAllLayers hidings;
- hostactualDISPLAY=:1 physicalDP-0 screen1920x1080, NVIDIA directrenderRTX5090D; noXvfb/windowManagerconfigurationchanges.

NativeGUI observation is separate from automatedgeometry/CLItests. Controller will capture only the task's viewer window, not unrelated desktop. CurrentPillowConda lacksXCB, ImageMagick/xdotool/wmctrl/mss absent; xwd exists. System /usr/bin/python3 Pillow XCB=True is available for targeted taskwindow screenshot. Prefer native renderer screenshot option if implemented, otherwise systemPillow actualscreenclip after observedxwininfo bounds. Do not install systempackages or invent a screenshot.
