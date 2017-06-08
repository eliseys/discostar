//
// Created by Konstantin Malanchev on 05/06/2017.
//

#define BOOST_TEST_MODULE single_star

#include <array>
#include <boost/test/included/unit_test.hpp>

#include "commons.hpp"
#include "BinaryParameters.hpp"
#include "LightCurve.hpp"


using namespace discostar;


struct Amplitudes{
	const float mass_ratio;
	const float limb_darkness;
	const float grav_darkness;
	const float inclination_degrees;
	const BinaryParameters bp;
	const std::vector<double> dots;
	const double avermin_max;
	const double min_min;

	Amplitudes(float mass_ratio, float limb_darkness, float grav_darkness, float inclination_degrees):
			mass_ratio(mass_ratio),
			limb_darkness(limb_darkness),
			grav_darkness(grav_darkness),
			inclination_degrees(inclination_degrees),
			bp(static_cast<float>(1 * SOLAR_RADIUS), // a
			   mass_ratio * 2e33f, // Mx
			   1 * 2e33f, // Mstar
			   static_cast<float>(inclination_degrees / 180.0 * M_PI), // i
			   limb_darkness, // limb_darkness
			   8400, // Tstar_pole
			   static_cast<float>(SOLAR_RADIUS), // Rstar
			   grav_darkness, // grav_darkness
			   0 * 5000, // Tdisk
			   0 * static_cast<float>(SOLAR_RADIUS), // Rdisk
			   0 * 0.05f, // z0Rdisk
			   0 * 1e37f // Lx
			),
			dots( LightCurve(bp, false).calc(4) ),
			avermin_max( delta_mag((dots[1] + dots[3])/2, dots[0]) ),
			min_min( delta_mag(dots[3], dots[1]) ){
		BOOST_TEST( std::round(delta_mag(dots[0], dots[2]) * 1e5) == 0 );
	}
};


BOOST_AUTO_TEST_CASE(Roche_lobe_star){
	BOOST_TEST( std::round(Amplitudes(1.0f, 0.0f, 0.00f,  0.0f).avermin_max * 1000) ==   0 );
	BOOST_TEST( std::round(Amplitudes(1.0f, 0.0f, 0.00f, 22.5f).avermin_max * 1000) ==  27 );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(1.0f, 0.00f, 45.0f) * 1000)) == 88  );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(1.0f, 0.00f, 67.5f) * 1000)) == 137 );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(1.0f, 0.00f, 90.0f) * 1000)) == 156 );
//
//	BOOST_TEST( static_cast<int>(std::round(amplitude(3.2f, 0.00f,  0.0f) * 1000)) == 0   );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(3.2f, 0.00f, 22.5f) * 1000)) == 32  );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(3.2f, 0.00f, 45.0f) * 1000)) == 109 );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(3.2f, 0.00f, 67.5f) * 1000)) == 177 );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(3.2f, 0.00f, 90.0f) * 1000)) == 203 );

//	BOOST_TEST( static_cast<int>(std::round(amplitude(1.0f, 1.00f,  0.0f) * 1000)) == 0   );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(1.0f, 1.00f, 22.5f) * 1000)) == 36  );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(1.0f, 1.00f, 45.0f) * 1000)) == 134 );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(1.0f, 1.00f, 67.5f) * 1000)) == 241 );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(1.0f, 1.00f, 90.0f) * 1000)) == 289 );
//
//	BOOST_TEST( static_cast<int>(std::round(amplitude(0.1f, 0.25f,  0.0f) * 1000)) == 0   );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(0.1f, 0.25f, 22.5f) * 1000)) == 15  );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(0.1f, 0.25f, 45.0f) * 1000)) == 48  );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(0.1f, 0.25f, 67.5f) * 1000)) == 75  );
//	BOOST_TEST( static_cast<int>(std::round(amplitude(0.1f, 0.25f, 90.0f) * 1000)) == 86  );


	BOOST_TEST( std::round(Amplitudes(1.0f, 0.4f, 0.00f,  0.0f).min_min * 1000 ) ==   0 );
	BOOST_TEST( std::round(Amplitudes(1.0f, 0.4f, 0.00f, 22.5f).min_min * 1000 ) == -15 );
	BOOST_TEST( std::round(Amplitudes(1.0f, 0.4f, 0.00f, 45.0f).min_min * 10000) == -50 );
	BOOST_TEST( std::round(Amplitudes(1.0f, 0.4f, 0.00f, 67.5f).min_min * 1000 ) ==  21 );
	BOOST_TEST( std::round(Amplitudes(1.0f, 0.4f, 0.00f, 90.0f).min_min * 1000 ) ==  35 );
}