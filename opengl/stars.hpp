//
// Created by Konstantin Malanchev on 05/04/2017.
//

#ifndef HERX1_STARS_HPP
#define HERX1_STARS_HPP


#include "commons.hpp"
#include "DescriteToGeometric.hpp"
#include "software_render.hpp"


namespace discostar {
namespace geometry {

class SphericalStar: public Sphere{
public:
    const float r;
    const float Q;

    SphericalStar(size_t bin_splits, float r, float Q):
            Sphere(bin_splits),
            r(r),
            Q(Q)
    {};

    virtual ObjectModel get_object_model() const{
        auto om = Sphere::get_object_model();

        for( auto &vert : om.vertices ){
            vert *= r;
        }

        return om;
    }

    virtual TextureImage get_texture_image() const { return TextureImage(1, Q, 0.0f, 0.0f); }
};


class RocheLobeStar: public Sphere{
public:
	const float mass_ratio;
	const float degree_of_feeling;
	const float Q_pole;
	const float gravitation_darkness;
	const float omega;

	RocheLobeStar(size_t bin_splits, float mass_ratio, float degree_of_feeling, float Q_pole, float gravitation_darkness):
			Sphere(bin_splits),
			mass_ratio(mass_ratio),
			Q_pole(Q_pole),
			gravitation_darkness(gravitation_darkness),
			degree_of_feeling(degree_of_feeling),
			omega( static_cast<float>( discostar::softrend::omg(mass_ratio, degree_of_feeling) ) )
	{};

	virtual ObjectModel get_object_model() const{
		auto om = Sphere::get_object_model();

		for ( auto &vert : om.vertices ){
			const auto vert_polar = glm::polar(vert);
			const float r = static_cast<float>( discostar::softrend::radius_star(
					vert_polar.y + M_PI_2,
					M_PI_2 - vert_polar.x, // softrend uses polar angle, glm::polar outputs latitude
					mass_ratio,
					omega
			) );
			vert *= r; // initial length of vert is 1
		}

		for ( auto &norm : om.normals ){
			const auto vert_polar = glm::polar(norm); // For sphere normal is equal to vertices coordinate
			double *gradient = discostar::softrend::gradient(
					vert_polar.y + M_PI_2,
					M_PI_2 - vert_polar.x, // softrend uses polar angle, glm::polar outputs latitude
					mass_ratio,
					omega
			);

			// Convert from softrend coordinates to Open GL coordinates
			norm.x = static_cast<float>(-gradient[1]);
			norm.y = static_cast<float>( gradient[3]);
			norm.z = static_cast<float>( gradient[2]);

			free(gradient); // gradient is allocated by malloc
		}

		return om;
	}

	virtual TextureImage get_texture_image() const{
		size_t texture_size =
				exp2_int(splits + 1) > GL_MAX_TEXTURE_SIZE ? GL_MAX_TEXTURE_SIZE : exp2_int(splits + 1);
		TextureImage ti(texture_size);

		double *gradient_pole = discostar::softrend::gradient(0, 0, mass_ratio, omega);
		const double g_pole = gradient_pole[0];
		free(gradient_pole);  // gradient_pole is allocated by malloc

		for ( size_t i_theta = 0; i_theta < texture_size; ++i_theta ){
			const float theta = static_cast<float>(M_PI) * (static_cast<float>(i_theta) / static_cast<float>(texture_size) - 0.5f);
			for ( size_t i_phi = 0; i_phi < texture_size; ++i_phi ){
				const auto phi = 2 * static_cast<float>(M_PI) * (static_cast<float>(i_phi) / static_cast<float>(texture_size) + 1);

				double *gradient = discostar::softrend::gradient(
						phi + M_PI_2,
						M_PI_2 - theta, // softrend uses polar angle, Sphere uses latitude
						mass_ratio,
						omega
				);
				const float darkness = static_cast<float>( pow(gradient[0] / g_pole, gravitation_darkness) );
				free(gradient); // gradient is allocated by malloc

				ti[i_phi * texture_size * 3 + i_theta * 3 + 0] = Q_pole * darkness;
				ti[i_phi * texture_size * 3 + i_theta * 3 + 1] = 0.f;
				ti[i_phi * texture_size * 3 + i_theta * 3 + 2] = 0.f;
			}
		}

		return ti;
	}
};


}} // namespace discostar::geometry

#endif //HERX1_STARS_HPP
