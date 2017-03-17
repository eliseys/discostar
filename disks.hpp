//
// Created by Konstantin Malanchev on 28/02/2017.
//

#ifndef HERX1_DISKS_HPP
#define HERX1_DISKS_HPP


#include <cmath>
#include <map>

#include <GL/glew.h>
#include <glm/glm.hpp>

#include "DescriteToGeometric.hpp"


class StandardDisk: public Circle{
public:
    const float r_out, H_r_out, H_out, n, Q_out;
    const size_t full_size;

    StandardDisk(size_t bin_splits, float r_out, float H_r_out, float dlogH_dlogr, float Q_out):
            Circle(bin_splits),
            r_out(r_out),
            H_r_out(H_r_out),
            n(dlogH_dlogr),
            Q_out(Q_out),
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
//        { // TODO: move cylinder to separate class
//            const size_t rho = this->rho_size - 1;
//            for ( size_t psi = 0; psi < this->psi_size(rho); ++psi ){
//                om.elements.push_back(      index(rho, psi    ));
//                om.elements.push_back(      index(rho, psi + 1));
//                om.elements.push_back(lower_index(rho, psi + 1));
//
//                om.elements.push_back(      index(rho, psi    ));
//                om.elements.push_back(lower_index(rho, psi    ));
//                om.elements.push_back(lower_index(rho, psi + 1));
//            }
//        }

        return om;
    }

    const TextureImage get_texture_image() const{
        size_t texture_size = exp2_int(splits + 1) > GL_MAX_TEXTURE_SIZE ? GL_MAX_TEXTURE_SIZE : exp2_int(splits + 1);
        TextureImage ti(texture_size);
        for (size_t i_phi = 0; i_phi < texture_size; ++i_phi ){
//                const float phi = 2 * static_cast<float>(M_PI) * static_cast<float>(i_phi) / static_cast<float>(texture_size);
            for ( size_t i_r = 1; i_r < texture_size; ++i_r ){
                const float r = static_cast<float>(i_r) / static_cast<float>(texture_size-1);
                ti[ i_phi * texture_size * 3 + i_r * 3 + 0 ] = Q_out * powf(r, -3.0f);
                ti[ i_phi * texture_size * 3 + i_r * 3 + 1 ] = 0.0;
                ti[ i_phi * texture_size * 3 + i_r * 3 + 2 ] = 0.0;
            }
            // for i_r = 0 use value from i_r = 1
            ti[i_phi * texture_size * 3 + 0] =  ti[ i_phi * texture_size * 3 + 3 + 0 ];
            ti[i_phi * texture_size * 3 + 1] =  ti[ i_phi * texture_size * 3 + 3 + 1 ];
            ti[i_phi * texture_size * 3 + 2] =  ti[ i_phi * texture_size * 3 + 3 + 2 ];
        }
        return ti;
    }
};


class Belt: public Basic3DObject{
public:
    const StandardDisk disk;
    const size_t rho;
    const size_t psi_size;

    Belt(const StandardDisk &disk):
            disk(disk),
            rho( disk.rho_size - 1 ),
            psi_size( disk.psi_size(rho) )
    {};

    unsigned short index(size_t psi) const{
        return static_cast<unsigned short>((psi % psi_size) * 2);
    }

    unsigned short lower_index(size_t psi) const{
        return static_cast<unsigned short>(index(psi) + 1);
    }

    const ObjectModel get_object_model() const{
        ObjectModel om;
        const auto disk_om = disk.get_object_model();

        for ( size_t psi = 0; psi < disk.psi_size(rho); ++psi ){
            const auto disk_i1 = disk.index(rho, psi);
            const auto disk_i2 = disk.lower_index(rho, psi);

            om.vertices.push_back(disk_om.vertices[disk_i1]);
            om.normals .push_back(disk_om.vertices[disk_i1]);
            om.uvs     .push_back(glm::vec2(0, disk_om.uvs[disk_i1].g));

            om.vertices.push_back(disk_om.vertices[disk_i2]);
            om.normals .push_back(disk_om.vertices[disk_i2]);
            om.uvs     .push_back(glm::vec2(1, disk_om.uvs[disk_i2].g));

            om.elements.push_back(      index(psi    ));
            om.elements.push_back(      index(psi + 1));
            om.elements.push_back(lower_index(psi + 1));

            om.elements.push_back(      index(psi    ));
            om.elements.push_back(lower_index(psi    ));
            om.elements.push_back(lower_index(psi + 1));
        }
        return om;
    }

    const TextureImage get_texture_image() const{
        TextureImage ti(1);
        ti[0] = disk.Q_out;
        ti[1] = 0.0f;
        ti[2] = 0.0f;
        return ti;
    }
};


#endif //HERX1_DISKS_HPP
