#include "spin/detail/FixedNormalSphere.h"
#include <iostream>
#include <iomanip>
int main() {
    const auto sphere=spin::detail::fixedNormalSphere(2);
    std::cout<<std::setprecision(17)<<'[';
    for(int k=0;k<2;++k) {
        if(k)std::cout<<',';
        const auto& n=sphere.vertices[k==0?40:135];
        std::cout<<'['<<n.x()<<','<<n.y()<<','<<n.z()<<']';
    }
    std::cout<<"]\n";
}
