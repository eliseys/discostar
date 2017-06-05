//
// Created by Konstantin Malanchev on 05/06/2017.
//

#define BOOST_TEST_MODULE Roche_lobe_star
#include <boost/test/included/unit_test.hpp>

#include "commons.hpp"
#include "BinaryParameters.hpp"
#include "LightCurve.hpp"
#include "Renderer.hpp"

using namespace discostar;

BOOST_AUTO_TEST_CASE(Roche_lobe_star)
{
	const size_t phases = 100;

    try{
        const BinaryParameters bp(
                static_cast<float>(1 * SOLAR_RADIUS), // a
                1 * 2e33f, // Mx
                1 * 2e33f, // Mstar
                static_cast<float>(90.0 / 180.0 * M_PI), // i
                8400, // Tstar_pole
                static_cast<float>(SOLAR_RADIUS), // Rstar
				0.25f, // grav_darkness
                0 * 5000, // Tdisk
                0 * static_cast<float>(SOLAR_RADIUS), // Rdisk
                0 * 0.05f, // z0Rdisk
                0 * 1e37f // Lx
        );
        const LightCurve lc(bp);

        const auto start = std::chrono::high_resolution_clock::now();
        const auto dots = lc.calc(phases);
        const auto end = std::chrono::high_resolution_clock::now();

		const double maximum = *std::max_element(dots.begin(), dots.end());
		const double minimum = 0.5 * (*std::min_element(dots.begin(), dots.begin()+dots.size()/2)
									 +*std::min_element(dots.begin()+dots.size()/2, dots.end()));
		const double amplitude = 2.5 * log10(maximum / minimum);

		BOOST_TEST(static_cast<int>(amplitude * 1000) == 194);
    } catch ( GlfwException e ){
		std::cerr << e.what() << std::endl;
    }
}