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


namespace discostar {
namespace geometry {


template <typename T_INDEX>
class Disk : public Circle<T_INDEX> {
public:
    Disk<T_INDEX>(T_INDEX bin_splits) : Circle<T_INDEX>(bin_splits) {};

    T_INDEX lower_index(T_INDEX rho, T_INDEX psi) const {
        return this->size + this->index(rho, psi);
    }

    T_INDEX lower_index(const DiscreteCoordinate<T_INDEX> &dc) const {
        return lower_index(dc.rho, dc.psi);
    }
};


template <typename T_INDEX>
class StandardDisk : public Disk<T_INDEX> {
public:
    const float r_out, H_r_out, H_out, n, Q_out;
    const T_INDEX full_size;

    StandardDisk<T_INDEX>(T_INDEX bin_splits, float r_out, float H_r_out, float dlogH_dlogr, float Q_out) :
            Disk<T_INDEX>(bin_splits),
            r_out(r_out),
            H_r_out(H_r_out),
            n(dlogH_dlogr),
            Q_out(Q_out),
            H_out(r_out * H_r_out),
            full_size(2 * this->size) {}

    float semiheight(const glm::vec2 &cyl) const {
        return H_out * powf(cyl.r, n);
    }

    float tangent(const glm::vec2 &cyl) const {
        return H_r_out * n * powf(cyl.r, n - 1);
    }

    glm::vec3 upper_normal(const glm::vec2 &cyl) const {
        if (cyl.r == 0) {
            return glm::vec3(0, 1, 0);
        }
        const glm::vec2 pol(M_PI_2 - atanf(tangent(cyl)), cyl.g - M_PI);
        return glm::euclidean(pol);
    }

    virtual ObjectModel<T_INDEX> get_object_model() const {
        ObjectModel<T_INDEX> om;

        om.vertices.resize(full_size);
        om.uvs.resize(full_size);
        om.normals.resize(full_size);

        for (T_INDEX t = 0; t < this->tr; ++t) {
            for (auto it = this->triangle_begin(t); it != this->triangle_end(t); ++it) {
                const auto cyl = this->cylindrical(*it);
                const auto r = static_cast<value_type>( cyl.r * r_out );
                const auto x = static_cast<value_type>( r * sinf(cyl.g));
                const auto y = static_cast<value_type>( semiheight(cyl));
                const auto z = static_cast<value_type>( r * cosf(cyl.g));

                const auto i1 = this->index(*it);
                const auto i2 = this->lower_index(*it);

                om.vertices[i1] = glm::vec3(x, y, z);
                om.normals[i1] = upper_normal(cyl);
                om.uvs[i1] = glm::vec2(cyl.r, (cyl.g / static_cast<value_type>(2 * M_PI)));

                om.vertices[i2] = glm::vec3(x, -y, z);
                om.normals[i2] = glm::vec3(om.normals[i1].x, -om.normals[i1].y, om.normals[i1].z);
                om.uvs[i2] = om.uvs[i1];
            }
        }

        om.elements = this->get_elements();
        auto lower_half = om.elements;
        for (auto &i : lower_half) {
            i = this->lower_index(this->coordinate(i));
        }
        om.elements.insert(om.elements.end(), lower_half.begin(), lower_half.end());

        return om;
    }

    virtual TextureImage get_texture_image() const {
        unsigned int texture_size =
                exp2_int(static_cast<unsigned int>(this->splits)  + 1) > GL_MAX_TEXTURE_SIZE
                ? GL_MAX_TEXTURE_SIZE
                : exp2_int(static_cast<unsigned int>(this->splits) + 1);
        TextureImage ti(texture_size);
        for (unsigned int i_phi = 0; i_phi < texture_size; ++i_phi) {
//                const float phi = 2 * static_cast<float>(M_PI) * static_cast<float>(i_phi) / static_cast<float>(texture_size);
            for (unsigned int i_r = 1; i_r < texture_size; ++i_r) {
                const float r = static_cast<float>(i_r) / static_cast<float>(texture_size - 1);
                ti[i_phi * texture_size * 3 + i_r * 3 + 0] = Q_out * powf(r, -3.0f);
                ti[i_phi * texture_size * 3 + i_r * 3 + 1] = 0.0;
                ti[i_phi * texture_size * 3 + i_r * 3 + 2] = 0.0;
            }
            // for i_r = 0 use value from i_r = 1
            ti[i_phi * texture_size * 3 + 0] = ti[i_phi * texture_size * 3 + 3 + 0];
            ti[i_phi * texture_size * 3 + 1] = ti[i_phi * texture_size * 3 + 3 + 1];
            ti[i_phi * texture_size * 3 + 2] = ti[i_phi * texture_size * 3 + 3 + 2];
        }
        return ti;
    }
};


template <typename T_INDEX>
class DiskBelt : public Basic3DObject<T_INDEX> {
public:
    const StandardDisk<T_INDEX> disk;
    const T_INDEX rho;
    const T_INDEX psi_size;

    DiskBelt<T_INDEX>(const StandardDisk<T_INDEX> &disk) :
            disk(disk),
            rho(disk.rho_size - 1),
            psi_size(disk.psi_size(rho)) {};

    T_INDEX index(T_INDEX psi) const {
        return (psi % psi_size) * 2;
    }

    T_INDEX lower_index(T_INDEX psi) const {
        return index(psi) + 1;
    }

    virtual ObjectModel<T_INDEX> get_object_model() const {
        ObjectModel<T_INDEX> om;
        const auto disk_om = disk.get_object_model();

        for (T_INDEX psi = 0; psi < disk.psi_size(rho); ++psi) {
            const auto disk_i1 = disk.index(rho, psi);
            const auto disk_i2 = disk.lower_index(rho, psi);

            om.vertices.push_back(disk_om.vertices[disk_i1]);
            om.normals.push_back(disk_om.vertices[disk_i1]);
            om.uvs.push_back(glm::vec2(0, disk_om.uvs[disk_i1].g));

            om.vertices.push_back(disk_om.vertices[disk_i2]);
            om.normals.push_back(disk_om.vertices[disk_i2]);
            om.uvs.push_back(glm::vec2(1, disk_om.uvs[disk_i2].g));

            om.elements.push_back(index(psi));
            om.elements.push_back(index(psi + 1));
            om.elements.push_back(lower_index(psi + 1));

            om.elements.push_back(index(psi));
            om.elements.push_back(lower_index(psi));
            om.elements.push_back(lower_index(psi + 1));
        }
        return om;
    }

    virtual TextureImage get_texture_image() const {
        const auto ti_disk = disk.get_texture_image();

        TextureImage ti(1);
        ti[0] = ti_disk[ti_disk.size() - 3];
        ti[1] = ti_disk[ti_disk.size() - 2];
        ti[2] = ti_disk[ti_disk.size() - 1];

        return ti;
    }
};


}} // namespace discostar::geometry

#endif //HERX1_DISKS_HPP
