#pragma once
#include <CGAL/Gmpq.h>
#include <string>

namespace verification {

// Accept only signed decimal integer / positive decimal integer. This is a
// transport for exact construction coordinates, not a floating conversion.
CGAL::Gmpq parseExactRationalCoordinate(const std::string& value);

} // namespace verification
