#include "spin2014/SurfaceRawMomentDifferential.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <fstream>
#include <iostream>
#include <limits>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
using namespace surface_test;
void checkOverflowCategories() {
    // Catches changing either legacy API's normal-overflow category or
    // accidentally translating later polynomial overflow to invalid_argument.
    bool categoriesPreserved=true;
    for(bool normalOverflow:{true,false}) {
        auto surface=cube(0.,1.);
        surface.v*=normalOverflow?1e200:1e100;
        const Eigen::Vector3d u=(surface.v.row(2)-surface.v.row(0)).transpose();
        const Eigen::Vector3d v=(surface.v.row(1)-surface.v.row(0)).transpose();
        require(surface.v.allFinite(),"overflow fixture coordinates must be finite");
        require(u.cross(v).allFinite()!=normalOverflow,"distinguish normal from later polynomial overflow");
        for(bool differential:{false,true}) {
            bool expected=false;
            const char* category="none";
            try {
                if(differential) spin2014::differentiateSurfaceRawMoments(surface.v,surface.f);
                else spin2014::integrateSurfaceRawMoments(surface.v,surface.f);
            } catch(const std::invalid_argument&) {
                category="invalid_argument"; expected=normalOverflow;
            } catch(const std::overflow_error&) {
                category="overflow_error"; expected=!normalOverflow;
            }
            std::cout<<(differential?"differential":"value")<<' '
                     <<(normalOverflow?"normal":"polynomial")<<" overflow category="<<category<<'\n';
            categoriesPreserved &= expected;
        }
    }
    require(categoriesPreserved,"both legacy APIs must preserve distinct normal and polynomial overflow categories");
}
int main(int argc,char**) {
    try {
        auto b=box({-2,1,3},{1,5,7});
        if(argc>1) {
            // Allocate caller storage before imposing an actual process address-space
            // cap. A 10x3N Jacobian exceeds the remaining allowance by >200 MB.
            Eigen::MatrixXd vertices=Eigen::MatrixXd::Zero(1000000,3);
            vertices.topRows(8)=b.v;
            auto child=[&](bool differential) {
                const pid_t pid=fork(); require(pid>=0,"fork failed");
                if(pid==0) {
                    std::size_t pages=0; std::ifstream("/proc/self/statm")>>pages;
                    rlimit limit{}; if(getrlimit(RLIMIT_AS,&limit)!=0) _exit(7);
                    limit.rlim_cur=pages*static_cast<std::size_t>(sysconf(_SC_PAGESIZE))+32*1024*1024;
                    if(setrlimit(RLIMIT_AS,&limit)!=0) _exit(3);
                    try {
                        if(differential) spin2014::differentiateSurfaceRawMoments(vertices,b.f);
                        else close(spin2014::integrateSurfaceRawMoments(vertices,b.f).values(),boxMoments({-2,1,3},{1,5,7}),"allocation value");
                        _exit(differential?4:0);
                    } catch(const std::bad_alloc&) { _exit(differential?0:5); }
                    catch(...) { _exit(6); }
                }
                int status=0; require(waitpid(pid,&status,0)==pid,"allocation subprocess wait failed");
                if(!WIFEXITED(status) || WEXITSTATUS(status)!=0) {
                    std::cerr<<"allocation child differential="<<differential<<" wait_status="<<status<<'\n';
                    throw std::runtime_error("allocation-sensitive subprocess failed");
                }
            };
            child(false); child(true);
            std::cout<<"value-only succeeds with 32 MiB headroom; dense differential allocation fails as expected, N=1000000\n"; return 0;
        }
        close(spin2014::integrateSurfaceRawMoments(b.v,b.f,2.5).values(),2.5*boxMoments({-2,1,3},{1,5,7}),"translated box all ten/density");
        Mesh t; t.v.resize(4,3); t.v<<1,2,3,3,2.5,3,1.25,5,3.5,1.5,2.25,7;
        t.f.resize(4,3); t.f<<1,2,3,0,2,1,0,1,3,0,3,2;
        const auto moments=spin2014::integrateSurfaceRawMoments(t.v,t.f);
        close(moments.values(),tetraMoments(t.v),"asymmetric tetra cross moments");
        Eigen::Vector3d shift(-3,4,1); t.v.rowwise()+=shift.transpose();
        close(spin2014::integrateSurfaceRawMoments(t.v,t.f,3).values(),3*tetraMoments(t.v),"translated tetra density");
        close(spin2014::integrateSurfaceRawMoments(t.v,t.f).values(),moments.transformed(Eigen::Matrix3d::Identity(),shift).values(),"translation covariance");
        auto cavity=join(b,box({-1,2,4},{0,3,5}),true);
        close(spin2014::integrateSurfaceRawMoments(cavity.v,cavity.f).values(),boxMoments({-2,1,3},{1,5,7})-boxMoments({-1,2,4},{0,3,5}),"signed cavity");
        close(spin2014::integrateSurfaceRawMoments(t.v,t.f).values(),spin2014::differentiateSurfaceRawMoments(t.v,t.f).moments.values(),"differential compatibility");
        for(double rho:{0.,-1.,std::numeric_limits<double>::infinity()}) rejects([&]{spin2014::integrateSurfaceRawMoments(b.v,b.f,rho);},"invalid density");
        auto bad=b; bad.v(0,0)=std::numeric_limits<double>::quiet_NaN();
        rejects([&]{spin2014::integrateSurfaceRawMoments(bad.v,bad.f);},"nonfinite vertex");
        rejects([&]{spin2014::integrateSurfaceRawMoments(b.v,b.f.topRows(11));},"open chain");
        bad=b; bad.f(0,0)=99; rejects([&]{spin2014::integrateSurfaceRawMoments(bad.v,bad.f);},"bad index");
        bad=b; bad.v*=1e100; rejects([&]{spin2014::integrateSurfaceRawMoments(bad.v,bad.f);},"moment overflow");
        checkOverflowCategories();
        std::cout<<"surface value analytic, signed, covariance, invalid-input checks passed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
