"""Regression for the adapter's real index mapping, not a Boolean predicate.

This compiles the production mapping function unchanged with a minimal
non-injective approximate-coordinate provider. Geometry/provenance integration
is separately exercised by the actual-operand Boolean replays.
"""
from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.fixture(scope="module")
def binding_program(tmp_path_factory):
    compiler = shutil.which("c++")
    if compiler is None:
        pytest.fail("C++ compiler required for actual Cherchi index-binding regression")
    source = (Path(__file__).resolve().parents[2] /
              "experiments/cherchi_reference_probe/probe.cpp").read_text()
    begin = source.index("void alignExactPointsWithApproximateCoordinates(")
    end = source.index("\nvoid replaceWithNearestExactCoordinates(", begin)
    harness = r'''
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
struct genericPoint {
    double coordinate=1.0;
    bool getApproxXYZCoordinates(double& x,double& y,double& z) const {
        x=coordinate; y=0.0; z=0.0; return true;
    }
};
std::string canonicalCoordinate(double value) {
    char buffer[64]; std::snprintf(buffer,sizeof(buffer),"%.17g",value==0.0?0.0:value);
    return buffer;
}
''' + source[begin:end] + r'''
int main(int argc,char** argv) {
    if(argc!=2) return 9;
    const std::string mode=argv[1];
    genericPoint a,b;
    if(mode!="collision") b.coordinate=2.0;
    std::vector<const genericPoint*> points{&a,&b};
    if(mode=="mismatch") {
        try {
            alignExactPointsWithApproximateCoordinates(points,1.0,{2.,0.,0.,1.,0.,0.});
        } catch(const std::runtime_error&) {
            return points[0]==&a && points[1]==&b ? 0 : 8;
        }
        return 7;
    }
    if(mode!="collision" && mode!="unique") return 9;
    alignExactPointsWithApproximateCoordinates(points,1.0,
        {1.,0.,0.,b.coordinate,0.,0.});
    return points[0]==&a && points[1]==&b ? 0 : 7;
}
'''
    directory = tmp_path_factory.mktemp("cherchi-output-binding")
    cpp = directory / "regression.cpp"
    binary = directory / "regression"
    cpp.write_text(harness)
    built = subprocess.run([compiler,"-std=c++17","-O0",str(cpp),"-o",str(binary)],
                           capture_output=True,text=True,timeout=30)
    assert built.returncode == 0, built.stderr
    return binary


def test_duplicate_display_coordinates_preserve_bound_point_identity(binding_program):
    result = subprocess.run([str(binding_program),"collision"],capture_output=True,timeout=5)
    assert result.returncode == 0, "display collisions must not permute exact point identities"


def test_distinct_already_bound_vertices_preserve_identity(binding_program):
    result = subprocess.run([str(binding_program),"unique"],capture_output=True,timeout=5)
    assert result.returncode == 0


def test_index_coordinate_disagreement_rejects_without_rematching(binding_program):
    result = subprocess.run([str(binding_program),"mismatch"],capture_output=True,timeout=5)
    assert result.returncode == 0, "index disagreement must fail closed without changing binding"
