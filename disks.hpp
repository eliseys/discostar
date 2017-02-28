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

    StandardDisk(size_t bin_splits, float r_out, float H_r_out, float dlogH_dlogr):
            Circle(bin_splits),
            r_out(r_out),
            H_r_out(H_r_out),
            n(dlogH_dlogr),
            H_out(r_out * H_r_out){}

    float semiheight(const glm::vec2 &cyl) const{
        return H_out * powf(cyl.r, n);
    }

    float tangent(const glm::vec2 &cyl) const{
        return H_r_out * n * powf(cyl.r, n - 1);
    }

    glm::vec3 upper_normal(const glm::vec2 &cyl) const{
        const glm::vec2 pol( M_PI_2 - atanf(tangent(cyl)), cyl.g - M_PI );
        std::cout << pol.x << "\t" << pol.y << "\t|\t"
                  << glm::euclidean(pol).x << "\t" << glm::euclidean(pol).y << "\t" << glm::euclidean(pol).z << "\t|\t"
                  << cyl.x << "\t" << cyl.y << std::endl;
        return glm::euclidean(pol);
    }

    const ObjectModel get_object_model() const{
        ObjectModel om;

        size_t full_size = 2 * this->size - 1;

        om.vertices.resize(full_size);
        om.uvs     .resize(full_size);
        om.normals .resize(full_size);

        for ( size_t t = 0; t < this->tr; ++t ) {
            for (auto it = this->triangle_begin(t); it != this->triangle_end(t); ++it) {
                const auto cyl = cylindrical(*it);
                const auto r   = static_cast<value_type>( cyl.r * r_out );
                const auto x   = static_cast<value_type>( r * sin(cyl.g) );
                const auto y   = static_cast<value_type>( semiheight(cyl) );
                const auto z   = static_cast<value_type>( r * cos(cyl.g) );

                const auto i = index(*it);

                    om.vertices[i] = glm::vec3(x, y, z);
                    om.normals [i] = upper_normal(cyl);
                    om.uvs     [i] = glm::vec2(cyl.r, static_cast<value_type>(cyl.g / (2 * M_PI)));

                if ( i != 0 ) {
                    om.vertices[full_size - i] = glm::vec3(x, -y, z);
                    om.normals [full_size - i] = glm::vec3(om.normals[i].x, -om.normals[i].y, om.normals[i].z);
                    om.uvs     [full_size - i] = om.uvs[i];
                }
            }
        }

        om.elements = this->get_elements();
        auto lower_half = om.elements;
        for (auto &i : lower_half ){
            if ( i == this->size ){
                i = 0;
            } else{
                i = static_cast<unsigned short>(full_size - i);
            }
        }
        om.elements.insert(om.elements.end(), lower_half.begin(), lower_half.end());

        return om;
    }
};


#endif //HERX1_DISKS_HPP
