#ifndef UTILS_MATH_H
#define UTILS_MATH_H

#include <vector>
namespace utils {
/* Test if the product of two numbers is bounded by a third number.
   Safe against overflow. The caller must guarantee
   0 <= factor1, factor2 <= limit; failing this is an error. */
extern bool is_product_within_limit(int factor1, int factor2, int limit);
extern bool is_product_within_limit(double factor1, double factor2, double limit);
extern std::pair<double,double> avg_and_standard_deviation(std::vector<double> data);
}

#endif
