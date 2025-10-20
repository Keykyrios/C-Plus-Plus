/**
 * @file
 * @brief Implementation of Neville's Algorithm for polynomial interpolation.
 * @details Neville's algorithm is an iterative method to find the value of a
 * polynomial that passes through a given set of n+1 points. It evaluates the
 * polynomial at a specific point `x` without needing to compute the
 * polynomial's coefficients. For more details, see:
 * https://en.wikipedia.org/wiki/Neville%27s_algorithm
 * @author [Mitrajit Ghorui](https://github.com/KeyKyrios)
 */

#include <cassert>        /// for assert
#include <cmath>          /// for fabs
#include <iostream>       /// for IO operations
#include <stdexcept>      /// for std::invalid_argument
#include <unordered_set>  /// for std::unordered_set
#include <vector>         /// for std::vector

/**
 * @namespace numerical_methods
 * @brief Numerical Methods
 */
namespace numerical_methods {
/**
 * @namespace nevilles_algorithm
 * @brief Functions for Neville's Algorithm
 */
namespace nevilles_algorithm {
/**
 * @brief Evaluates the interpolating polynomial at a specific point.
 * @param x_coords Vector of x-coordinates of the given points.
 * @param y_coords Vector of y-coordinates of the given points.
 * @param target The x-coordinate at which to evaluate the polynomial.
 * @returns The interpolated y-value at the target x-coordinate.
 * @throws std::invalid_argument if vectors are empty, have mismatched
 * sizes, or if x-coordinates are not unique.
 */
double interpolate(const std::vector<double>& x_coords,
                   const std::vector<double>& y_coords, double target) {
    if (x_coords.empty() || y_coords.empty()) {
        throw std::invalid_argument("Input vectors cannot be empty.");
    }
    if (x_coords.size() != y_coords.size()) {
        throw std::invalid_argument(
            "x and y coordinate vectors must be the same size.");
    }

    std::unordered_set<double> seenX;
    for (double val : x_coords) {
        if (!seenX.insert(val).second) {
            throw std::invalid_argument("Input x-coordinates must be unique.");
        }
    }

    size_t n = x_coords.size();
    std::vector<double> p = y_coords;  // Initialize p with y values

    for (size_t k = 1; k < n; ++k) {
        for (size_t i = 0; i < n - k; ++i) {
            p[i] = ((target - x_coords[i + k]) * p[i] +
                    (x_coords[i] - target) * p[i + 1]) /
                   (x_coords[i] - x_coords[i + k]);
        }
    }

    return p[0];
}
}  // namespace nevilles_algorithm
}  // namespace numerical_methods

/**
 * @brief Self-test implementations
 * @returns void
 */
static void test() {
    // Test 1: Linear interpolation y = 2x + 1
    // Points (0, 1), (2, 5). Target x=1, expected y=3.
    std::vector<double> x1 = {0, 2};
    std::vector<double> y1 = {1, 5};
    double target1 = 1;
    double expected1 = 3.0;
    assert(fabs(numerical_methods::nevilles_algorithm::interpolate(x1, y1, target1) -
                expected1) < 1e-9);
    std::cout << "Linear test passed." << std::endl;

    // Test 2: Quadratic interpolation y = x^2
    // Points (0, 0), (1, 1), (3, 9). Target x=2, expected y=4.
    std::vector<double> x2 = {0, 1, 3};
    std::vector<double> y2 = {0, 1, 9};
    double target2 = 2;
    double expected2 = 4.0;
    assert(fabs(numerical_methods::nevilles_algorithm::interpolate(x2, y2, target2) -
                expected2) < 1e-9);
    std::cout << "Quadratic test passed." << std::endl;

    // Test 3: Negative numbers y = x^2 - 2x + 1
    // Points (-1, 4), (0, 1), (2, 1). Target x=1, expected y=0.
    std::vector<double> x3 = {-1, 0, 2};
    std::vector<double> y3 = {4, 1, 1};
    double target3 = 1;
    double expected3 = 0.0;
    assert(fabs(numerical_methods::nevilles_algorithm::interpolate(x3, y3, target3) -
                expected3) < 1e-9);
    std::cout << "Negative numbers test passed." << std::endl;

    // Test 4: Exception for duplicate x-coordinates
    std::vector<double> x4 = {1, 2, 1};
    std::vector<double> y4 = {5, 8, 3};
    bool exception_thrown = false;
    try {
        numerical_methods::nevilles_algorithm::interpolate(x4, y4, 1.5);
    } catch (const std::invalid_argument& e) {
        exception_thrown = true;
    }
    assert(exception_thrown);
    std::cout << "Duplicate X coordinate test passed." << std::endl;

    std::cout << "All tests have successfully passed!" << std::endl;
}

/**
 * @brief Main function
 * @returns 0 on exit
 */
int main() {
    test();  // run self-test implementations
    return 0;
}
