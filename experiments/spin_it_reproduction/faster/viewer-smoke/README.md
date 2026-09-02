# Spin-It Faster real-desktop OpenGL smoke

The viewer was built with OpenGL/GLFW/ImGui and opened on the host GNOME Xorg
display `:1`; no Xvfb or other virtual display was used.  All four solver modes
created a 1280×800 interactive window from sealed artifacts.  The static case
also received a real X11 button event that disabled the “Cavity / signed
quadric” layer, and the resulting host-window image is preserved in
`static-cavity-off.png`.

The UI independently exposes fixed outer shell, design domain, retained
material, cavity/signed quadric region, CoM/target/principal axes, and the
target/support plane.  The source diagnostics and viewer executable hashes are
recorded in `report.json`.
