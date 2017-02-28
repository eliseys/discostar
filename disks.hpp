//
// Created by Konstantin Malanchev on 28/02/2017.
//

#ifndef HERX1_DISKS_HPP
#define HERX1_DISKS_HPP


#include <cmath>
#include <glm/glm.hpp>

#include "DescriteToGeometric.hpp"


class StandardDisk: public Circle{
public:
    const float r_out, H_r_out, H_out, n;
    const size_t full_size;

    StandardDisk(size_t bin_splits, float r_out, float H_r_out, float dlogH_dlogr):
            Circle(bin_splits),
            r_out(r_out),
            H_r_out(H_r_out),
            n(dlogH_dlogr),
            H_out(r_out * H_r_out),
            full_size(2 * size){}

    float semiheight(const glm::vec2 &cyl) const{
        return H_out * powf(cyl.r, n);
    }

    float tangent(const glm::vec2 &cyl) const{
        return H_r_out * n * powf(cyl.r, n - 1);
    }

    glm::vec3 upper_normal(const glm::vec2 &cyl) const{
        if ( cyl.r == 0 ){
            return glm::vec3(0, 1, 0);
        }
        const glm::vec2 pol( M_PI_2 - atanf(tangent(cyl)), cyl.g - M_PI );
        return glm::euclidean(pol);
    }

    const unsigned short lower_index(size_t rho, size_t psi) const{
        return static_cast<unsigned short>(size) + index(rho, psi);
    }

    const unsigned short lower_index(const DiscreteCoordinate &dc) const{
        return lower_index(dc.rho, dc.psi);
    }

    const ObjectModel get_object_model() const{
        ObjectModel om;

        om.vertices.resize(full_size);
        om.uvs     .resize(full_size);
        om.normals .resize(full_size);

        for ( size_t t = 0; t < this->tr; ++t ) {
            for (auto it = this->triangle_begin(t); it != this->triangle_end(t); ++it) {
                const auto cyl = cylindrical(*it);
                const auto r   = static_cast<value_type>( cyl.r * r_out );
                const auto x   = static_cast<value_type>( r * sinf(cyl.g) );
                const auto y   = static_cast<value_type>( semiheight(cyl) );
                const auto z   = static_cast<value_type>( r * cosf(cyl.g) );

                const auto i1 = index(*it);
                const auto i2 = lower_index(*it);

                om.vertices[i1] = glm::vec3(x, y, z);
                om.normals [i1] = upper_normal(cyl);
                om.uvs     [i1] = glm::vec2(cyl.r, static_cast<value_type>(cyl.g / (2 * M_PI)));

                om.vertices[i2] = glm::vec3(x, -y, z);
                om.normals [i2] = glm::vec3(om.normals[i1].x, -om.normals[i1].y, om.normals[i1].z);
                om.uvs     [i2] = om.uvs[i1];
            }
        }

        om.elements = this->get_elements();
        auto lower_half = om.elements;
        for (auto &i : lower_half ){
            i = lower_index(coordinate(i));
        }
        om.elements.insert(om.elements.end(), lower_half.begin(), lower_half.end());
        {
            const size_t rho = this->rho_size - 1;
            for ( size_t psi = 0; psi < this->psi_size(rho); ++psi ){
                om.elements.push_back(      index(rho, psi    ));
                om.elements.push_back(      index(rho, psi + 1));
                om.elements.push_back(lower_index(rho, psi + 1));

                om.elements.push_back(      index(rho, psi    ));
                om.elements.push_back(lower_index(rho, psi    ));
                om.elements.push_back(lower_index(rho, psi + 1));
            }
        }

        return om;
    }
};


#endif //HERX1_DISKS_HPP
