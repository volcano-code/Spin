#include "viewer/SealedSpinView.h"
#include <cstdint>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <GLFW/glfw3.h>
#include <pybind11/embed.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <charconv>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace sealed_spin_view {
namespace {
namespace py = pybind11;
using Vertices = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using Faces = Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>;
struct Mesh { std::string key, kind; Vertices vertices; Faces faces; };
struct Scene {
    std::vector<Mesh> meshes;
    Eigen::Vector3d center, targetOrigin, targetAxis;
    std::vector<Eigen::Vector3d> axes, support;
    std::vector<int> multiplicities;
    std::vector<std::string> labels;
    std::string summary, artifact, geometryStatus;
    double mass;
};
void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
std::string modulePath(const py::handle& module) {
    return std::filesystem::canonical(py::str(module.attr("__file__")).cast<std::string>()).string();
}

Scene admit(const std::string& bundle, const std::string& digest) {
    PyConfig config;
    PyConfig_InitIsolatedConfig(&config);
    config.install_signal_handlers = 0;
    config.write_bytecode = 0;
    const auto status = PyConfig_SetBytesString(&config, &config.program_name, SEALED_SPIN_PYTHON);
    if (PyStatus_Exception(status)) {
        const std::string error = status.err_msg ? status.err_msg : "CPython configuration failed";
        PyConfig_Clear(&config);
        throw std::runtime_error(error);
    }
    py::scoped_interpreter interpreter(&config, 0, nullptr, false);
    try {
        const auto source = std::filesystem::canonical(SEALED_SPIN_SOURCE_PYTHON);
        const auto extension = std::filesystem::canonical(SEALED_SPIN_EXTENSION);
        auto sys = py::module_::import("sys");
        require(std::filesystem::canonical(sys.attr("executable").cast<std::string>()) ==
                std::filesystem::canonical(SEALED_SPIN_PYTHON), "unexpected embedded interpreter path");
        sys.attr("path").attr("insert")(0, extension.parent_path().string());
        sys.attr("path").attr("insert")(0, source.string());
        auto package = py::module_::import("make_it_stand");
        require(modulePath(package) == (source / "make_it_stand/__init__.py").string(),
                "unexpected make_it_stand package path");
        auto native = py::module_::import("_make_it_stand");
        require(modulePath(native) == extension.string(), "unexpected exact extension path");
        auto adapter = py::module_::import("make_it_stand.verification.spin_view");
        require(modulePath(adapter) == (source / "make_it_stand/verification/spin_view.py").string(),
                "unexpected scene adapter path");
        auto value = adapter.attr("load_spin_scene")(bundle, py::arg("expected_problem_digest") = digest).cast<py::dict>();
        py::dict loaded;
        loaded["package"] = modulePath(package); loaded["native"] = modulePath(native);
        loaded["adapter"] = modulePath(adapter);
        loaded["interpreter"] = sys.attr("executable");
        for (auto item : sys.attr("modules").cast<py::dict>()) {
            const auto name = py::str(item.first).cast<std::string>();
            if (name.rfind("make_it_stand.", 0) == 0 && py::hasattr(item.second, "__file__")) {
                const auto actual = modulePath(item.second);
                require(actual.rfind((source / "make_it_stand").string() + "/", 0) == 0,
                        "unexpected project module path: " + name);
                loaded[py::str(name)] = actual;
            }
        }
        Scene scene;
        py::dict summary, layers;
        auto hash = py::module_::import("hashlib").attr("sha256");
        for (auto item : value["meshes"].cast<py::list>()) {
            auto row = item.cast<py::dict>();
            Mesh mesh{row["key"].cast<std::string>(), row["kind"].cast<std::string>(),
                      row["vertices"].cast<Vertices>(), row["faces"].cast<Faces>()};
            require(mesh.vertices.allFinite(), "nonfinite scene coordinates");
            require(mesh.faces.size() == 0 || (mesh.faces.minCoeff() >= 0 &&
                    mesh.faces.maxCoeff() < mesh.vertices.rows()), "scene index out of range");
            // Observe actual copy-owned C++ bytes after conversion, not only Python descriptors.
            auto vh = hash(py::bytes(reinterpret_cast<const char*>(mesh.vertices.data()),
                           mesh.vertices.size()*sizeof(double))).attr("hexdigest")();
            auto fh = hash(py::bytes(reinterpret_cast<const char*>(mesh.faces.data()),
                           mesh.faces.size()*sizeof(int))).attr("hexdigest")();
            auto identity = row["identity"].cast<py::dict>();
            require(vh.equal(identity["vertices"]["sha256"]) && fh.equal(identity["faces"]["sha256"]),
                    "scene conversion changed geometry bytes");
            py::dict counts;
            counts["vertices"] = mesh.vertices.rows(); counts["faces"] = mesh.faces.rows();
            counts["vertices_sha256"] = vh; counts["faces_sha256"] = fh; counts["density"] = row["density"];
            layers[py::str(mesh.key)] = counts;
            scene.meshes.push_back(std::move(mesh));
        }
        scene.center = value["physical"]["center_of_mass"].cast<Eigen::Vector3d>();
        scene.mass = value["physical"]["mass"].cast<double>();
        scene.targetOrigin = value["target"]["origin"].cast<Eigen::Vector3d>();
        scene.targetAxis = value["target"]["axis"].cast<Eigen::Vector3d>();
        for (auto point : value["target"]["support_points"].cast<py::list>())
            scene.support.push_back(point.cast<Eigen::Vector3d>());
        for (auto axis : value["principal_axes"]["unique_axes"].cast<py::list>())
            scene.axes.push_back(axis["direction"].cast<Eigen::Vector3d>());
        scene.multiplicities = value["principal_axes"]["multiplicities"].cast<std::vector<int>>();
        scene.geometryStatus = value["geometry_status"].cast<std::string>();
        scene.artifact = value["artifact_digest"].cast<std::string>();
        scene.labels.push_back("HOLLOW_ONLY | Original model posture retained | Display is not simulation");
        scene.labels.push_back("Mode: " + value["target"]["mode"].cast<std::string>() + " | Units: " +
            value["length_unit"].cast<std::string>() + ", " + value["mass_unit"].cast<std::string>());
        for (auto item : value["disclosures"].cast<py::dict>())
            scene.labels.push_back(py::str(item.first).cast<std::string>() + ": " + py::str(item.second).cast<std::string>());
        scene.labels.push_back("Support: target points only; no certified hull. No waterline.");
        for (const char* name : {"artifact_digest", "evidence_digest", "problem_digest", "geometry_status",
             "disclosures", "physical", "principal_axes", "support_display", "target", "validation_scope"})
            summary[py::str(name)] = value[py::str(name)];
        summary["layers"] = layers; summary["loaded_paths"] = loaded;
        scene.summary = py::module_::import("json").attr("dumps")(summary, py::arg("sort_keys") = true).cast<std::string>();
        return scene;
    } catch (const py::error_already_set& error) {
        // Format the error before the interpreter and all Python references die.
        throw std::runtime_error(error.what());
    }
}

Eigen::RowVector3d color(const std::string& kind) {
    if (kind == "cavity") return {230/255., 26/255., 20/255.};
    if (kind == "retained" || kind == "material") return {242/255.,173/255.,31/255.};
    if (kind == "design") return {46/255.,122/255.,224/255.};
    return {199/255.,214/255.,242/255.};
}

int display(const Scene& scene, int frames) {
    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    menu.callback_draw_viewer_menu = [] {};
    const auto& original = scene.meshes.front();
    const Eigen::RowVector3d lo = original.vertices.colwise().minCoeff(), hi = original.vertices.colwise().maxCoeff();
    const double span = (hi-lo).maxCoeff(), axisLength = span*.3;
    std::map<std::string, bool> visible{{"outer",false},{"fixed",false},{"material",true},
                                      {"retained",false},{"design",false},{"cavity",true}};
    bool cutaway = true, isolate = false, explode = false, axes = true;
    float cutFraction = .5f;
    unsigned int left = viewer.core().id, right = 0;
    std::vector<int> ids;
    for (std::size_t i=0; i<scene.meshes.size(); ++i) {
        if (i) viewer.append_mesh();
        ids.push_back(viewer.data().id);
    }
    viewer.append_mesh(); const int marks = viewer.data().id;
    auto refresh = [&] {
        for (std::size_t i=0; i<scene.meshes.size(); ++i) {
            const auto& mesh = scene.meshes[i]; auto& data = viewer.data(ids[i]); data.clear();
            if (mesh.faces.rows()) {
                Vertices vertices = mesh.vertices; Faces faces = mesh.faces;
                if (mesh.kind != "original" && mesh.kind != "cavity" && cutaway) {
                    const double plane = lo.x()+cutFraction*(hi.x()-lo.x());
                    int count = 0;
                    for (int row=0; row<mesh.faces.rows(); ++row) {
                        bool keep = true;
                        for (int col=0; col<3; ++col) keep = keep && mesh.vertices(mesh.faces(row,col),0) <= plane;
                        if (keep) faces.row(count++) = mesh.faces.row(row);
                    }
                    faces.conservativeResize(count,3);
                }
                if (mesh.kind == "cavity" && explode) vertices.col(0).array() += span*.65;
                if (faces.rows()) { data.set_mesh(vertices,faces); data.set_colors(color(mesh.kind)); }
            }
            data.show_lines = false; data.double_sided = true;
            data.set_visible(mesh.kind == "original", left);
            data.set_visible(mesh.kind != "original" && visible[mesh.kind] && (!isolate || mesh.kind == "cavity"), right);
        }
        auto& data = viewer.data(marks); data.clear();
        if (axes) {
            data.add_points(scene.center.transpose(), Eigen::RowVector3d(.95,.15,.2));
            for (const auto& axis : scene.axes)
                data.add_edges((scene.center-axisLength*axis).transpose(), (scene.center+axisLength*axis).transpose(),
                               Eigen::RowVector3d(.2,.85,.45));
            const Eigen::Vector3d direction = scene.targetAxis.stableNormalized();
            data.add_edges(scene.targetOrigin.transpose(), (scene.targetOrigin+axisLength*direction).transpose(),
                           Eigen::RowVector3d(.85,.35,.85));
            for (const auto& point : scene.support) data.add_points(point.transpose(), Eigen::RowVector3d(.9,.85,.1));
        }
        data.point_size = 9; data.line_width = 3;
        data.set_visible(false,left); data.set_visible(true,right);
    };
    auto resize = [&](int w,int h) {
        int windowWidth, windowHeight;
        glfwGetWindowSize(viewer.window,&windowWidth,&windowHeight);
        const int reserved = 330*h/std::max(1,windowHeight);
        const int displayHeight = std::max(1,h-reserved);
        viewer.core(left).viewport = Eigen::Vector4f(0,0,w/2,displayHeight);
        viewer.core(right).viewport = Eigen::Vector4f(w/2,0,w-w/2,displayHeight);
    };
    viewer.callback_init = [&](auto&) {
        right = viewer.append_core(Eigen::Vector4f(750,0,750,670));
        int w,h; glfwGetFramebufferSize(viewer.window,&w,&h); resize(w,h);
        viewer.selected_core_index = 0;
        for (auto& core : viewer.core_list) core.background_color = Eigen::Vector4f(.08,.10,.14,1);
        refresh(); return false;
    };
    viewer.callback_post_resize = [&](auto&,int w,int h) { resize(w,h); return true; };
    viewer.callback_pre_draw = [&](auto&) {
        const auto& from = viewer.core(); auto& to = viewer.core(from.id == left ? right : left);
        to.trackball_angle = from.trackball_angle; to.camera_zoom = from.camera_zoom;
        to.camera_translation = from.camera_translation; to.camera_base_zoom = from.camera_base_zoom;
        to.camera_base_translation = from.camera_base_translation;
        to.orthographic = from.orthographic;
        return false;
    };
    menu.callback_draw_custom_window = [&] {
        ImGui::SetNextWindowPos(ImVec2(0,0),ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x,330),ImGuiCond_Always);
        ImGui::Begin("Sealed spin geometry observation",nullptr,
            ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
        ImGui::Text("%s | %s",scene.geometryStatus.c_str(),scene.artifact.c_str());
        ImGui::Columns(2,nullptr,false);
        for (const auto& label : scene.labels) ImGui::TextWrapped("%s",label.c_str());
        ImGui::NextColumn();
        bool changed = false; int index = 0;
        for (auto& entry : visible) {
            const std::string label = entry.first == "fixed" ? "Fixed shell material" :
                                      entry.first == "outer" ? "Outer surface (not shell)" : entry.first;
            changed |= ImGui::Checkbox(label.c_str(),&entry.second);
            if (++index % 3) ImGui::SameLine();
        }
        changed |= ImGui::Checkbox("Display cutaway",&cutaway); ImGui::SameLine();
        ImGui::PushItemWidth(170);
        changed |= ImGui::SliderFloat("x fraction",&cutFraction,0.f,1.f); ImGui::PopItemWidth();
        changed |= ImGui::Checkbox("Isolate cavity",&isolate); ImGui::SameLine();
        changed |= ImGui::Checkbox("Explode cavity (display only)",&explode);
        changed |= ImGui::Checkbox("Material CoM / unique axes / target points",&axes);
        ImGui::TextWrapped("Cutaway filters triangles; NO section caps. Explode translates a display copy. Neither is exportable material.");
        ImGui::Text("Complete material mass %.9g; CoM (%.6g, %.6g, %.6g)",scene.mass,scene.center.x(),scene.center.y(),scene.center.z());
        ImGui::Text("Unique principal axes: %zu (sign conventional)",scene.axes.size());
        for (int n : scene.multiplicities) if (n>1) ImGui::Text("Degenerate eigenspace: dimension %d; arbitrary axes suppressed",n);
        ImGui::TextWrapped("Axes: full world central inertia; near-equal eigenvalues grouped at relative 1e-10 display tolerance.");
        ImGui::TextWrapped("Red: CoM | Green: unique axes | Purple: target axis | Yellow: support points");
        ImGui::Columns(1); ImGui::Separator();
        ImGui::Text("LEFT: original model                         RIGHT: sealed final material (gold), cavity (red)");
        ImGui::End(); if (changed) refresh();
    };
    int rendered = 0;
    viewer.callback_post_draw = [&](auto&) {
        if (frames>0 && ++rendered>=frames) glfwSetWindowShouldClose(viewer.window,GLFW_TRUE);
        return false;
    };
    viewer.core().is_animating = frames>0;
    require(viewer.launch_init(true,false,"Sealed spin | HOLLOW_ONLY | PARTIAL_BLOCKED",1500,1000) == EXIT_SUCCESS,
            "native GLFW window initialization failed");
    glfwSetWindowSizeLimits(viewer.window,1200,850,GLFW_DONT_CARE,GLFW_DONT_CARE);
    // launch_init automatically frames every visible mesh after callback_init;
    // apply the single original-model camera only after that upstream behavior.
    viewer.selected_core_index = 0;
    viewer.core(left).align_camera_center(original.vertices,original.faces);
    viewer.core(left).camera_zoom *= .72f;
    viewer.core(left).trackball_angle = Eigen::Quaternionf(Eigen::AngleAxisf(.28f,Eigen::Vector3f::UnitX()) *
                                                         Eigen::AngleAxisf(-.5f,Eigen::Vector3f::UnitY()));
    viewer.launch_rendering(true); viewer.launch_shut(); return 0;
}
} // namespace

int run(int argc,char* argv[]) {
    if (argc<2) return -1;
    const std::string option = argv[1];
    if (option != "--sealed-spin-bundle" && option != "--check-sealed-spin-bundle") return -1;
    try {
        const bool check = option == "--check-sealed-spin-bundle";
        require(argc==4 || (!check && argc==6 && std::string(argv[4])=="--frames"),
            "usage: make-it-stand --[check-]sealed-spin-bundle BUNDLE EXPECTED_PROBLEM_SHA256 [--frames N]");
        int frames = 0;
        if (argc==6) {
            const std::string count = argv[5];
            const auto parsed = std::from_chars(count.data(),count.data()+count.size(),frames);
            require(parsed.ec==std::errc() && parsed.ptr==count.data()+count.size() && frames>0 && frames<=10000,
                    "--frames requires an integer in [1,10000]");
        }
        const Scene scene = admit(argv[2],argv[3]);
        if (check) { std::cout << scene.summary << '\n'; return 0; }
        return display(scene,frames);
    } catch (const std::exception& error) {
        std::cerr << "sealed-spin admission/view: " << error.what() << '\n'; return 1;
    }
}
} // namespace sealed_spin_view
