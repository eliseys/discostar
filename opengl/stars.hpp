//
// Created by Konstantin Malanchev on 05/04/2017.
//

#ifndef HERX1_STARS_HPP
#define HERX1_STARS_HPP


#include "commons.hpp"
#include "DescriteToGeometric.hpp"


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


}} // namespace discostar::geometry

#endif //HERX1_STARS_HPP
