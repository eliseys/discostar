//
// Created by Konstantin Malanchev on 21/02/2017.
//

#ifndef TUTORIALS_DESCRITETOGEOMETRIC_HPP
#define TUTORIALS_DESCRITETOGEOMETRIC_HPP

#include <cmath>
#include <vector>
#include <glm/gtx/polar_coordinates.hpp>

#include "ObjectModel.hpp"
#include "TriangleDiscreteCoordinates.hpp"


typedef typename glm::vec3::value_type value_type;


class Sphere: public TriangleDiscreteCoordinates{
protected:
    std::vector<glm::vec3> northern_hemisphere() const{
        std::vector<glm::vec3> vertices(this->size);

        const glm::vec3 A = glm::euclidean(glm::vec2(M_PI_2, 0));
        for ( size_t t = 0; t < this->tr; ++t ){
            const glm::vec3 B = glm::euclidean(glm::vec2(0, M_PI_2 * t));
            const glm::vec3 C = glm::euclidean(glm::vec2(0, M_PI_2 * (t + 1)));
            for ( auto it = this->triangle_begin(t); it != this->triangle_end(t); ++it ){
                const auto rho        = static_cast<value_type>( it->rho );
                const auto rho_length = static_cast<value_type>( this->rho_size ) - 1;
                const auto psi_length = static_cast<value_type>( this->psi_triangle_size(it->rho) );
                const auto psi        = static_cast<value_type>( it->psi % this->psi_triangle_size(it->rho) );

                value_type alpha ( 1 - rho / rho_length );
                value_type beta  ( (1 - alpha) * ( 1 - psi / psi_length ) );
                value_type gamma ( 1 - alpha - beta );
                if ( it->rho == 0 ) {
                    alpha = 1;
                    beta  = 0;
                    gamma = 0;
                }
                if ( it->rho == 1 ){
                    beta  = 1 - alpha;
                    gamma = 0;
                }
                vertices[this->index(*it)] = glm::normalize( alpha * A + beta * B + gamma * C );
            }
        }

        return vertices;
    }

public:
    Sphere(size_t bin_splits): TriangleDiscreteCoordinates(4, bin_splits){}

    const ObjectModel get_object_model() const{
        ObjectModel om;

        om.vertices = northern_hemisphere();
        // Add southern hemisphere vertices
        const size_t last_line = this->tr * (this->rho_size - 2) * (this->rho_size - 1) / 2 + 1;
        for ( size_t i = 0; i < last_line; ++i ){
            auto southern_vert = om.vertices[i];
            southern_vert.y *= -1;
            om.vertices.push_back(southern_vert);
        }

        for ( auto &vert : om.vertices ){
            const auto polar  = glm::polar(vert);
            const auto theta  = static_cast<value_type>(polar.x);
            const auto phi    = static_cast<value_type>(polar.y);
            const auto v      = static_cast<value_type>( (theta + M_PI_2) / M_PI );
            const auto u      = static_cast<value_type>( (phi + M_PI) / ( 2 * M_PI ) );
            om.uvs.emplace_back( u, v );
        }

        om.normals.insert(om.normals.begin(), om.vertices.begin(), om.vertices.end());

        om.elements = this->get_elements();
        // generate southern hemisphere triangles
        auto southern_ind = om.elements;
        for (auto &i : southern_ind){
            if (i < last_line){
                i += this->size;
            }
        }
        om.elements.insert(om.elements.end(), southern_ind.begin(), southern_ind.end());

        return om;
    }
};


class Circle: public TriangleDiscreteCoordinates{
public:
    Circle(size_t bin_splits): TriangleDiscreteCoordinates(6, bin_splits){}

    glm::vec2 cylindrical(const DiscreteCoordinate &dc) const{
        const auto rho        = static_cast<value_type>( dc.rho );
        const auto rho_length = static_cast<value_type>( this->rho_size ) - 1;
        const auto psi_length = static_cast<value_type>( this->psi_size(dc.rho) );
        const auto psi        = static_cast<value_type>( dc.psi );
        const auto r          = static_cast<value_type>( rho / rho_length );
        const auto phi        = static_cast<value_type>( 2 * M_PI * psi / psi_length );
        return glm::vec2(r, phi);
    }

    const ObjectModel get_object_model() const{
        ObjectModel om;

        om.vertices.resize(this->size);
        om.uvs     .resize(this->size);
        om.normals .resize(this->size);
        for ( size_t t = 0; t < this->tr; ++t ) {
            for (auto it = this->triangle_begin(t); it != this->triangle_end(t); ++it) {
                const auto cyl = cylindrical(*it);
                const auto x   = static_cast<value_type>( cyl.r * sinf(cyl.g) );
                const auto y   = static_cast<value_type>( 0 );
                const auto z   = static_cast<value_type>( cyl.r * cosf(cyl.g) );

                om.vertices[this->index(*it)] = glm::vec3(x, y, z);
                om.normals [this->index(*it)] = glm::vec3(0, 1, 0);
                om.uvs     [this->index(*it)] = glm::vec2(cyl.r, static_cast<value_type>(cyl.g / (2 * M_PI)));
            }
        }

        om.elements = this->get_elements();

        return om;
    }
};



#endif //TUTORIALS_DESCRITETOGEOMETRIC_HPP
