#include "viewer/SealedSpinView.h"

// The existing strict viewer, isolated from optimization and legacy loaders.
int main(int argc, char** argv) {
    const int result = sealed_spin_view::run(argc, argv);
    return result < 0 ? 2 : result;
}
