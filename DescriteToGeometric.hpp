//
// Created by Konstantin Malanchev on 21/02/2017.
//

#ifndef TUTORIALS_DESCRITETOGEOMETRIC_HPP
#define TUTORIALS_DESCRITETOGEOMETRIC_HPP

#include <cmath>
#include <vector>
#include <glm/gtx/polar_coordinates.hpp>

#include "TriangleDiscreteCoordinates.hpp"


typedef typename glm::vec3::value_type value_type;


template<typename T_INDEX>
class Sphere: public TriangleDiscreteCoordinates<T_INDEX>{
protected:
    std::vector<glm::vec3> northern_hemisphere() const{
        std::vector<glm::vec3> vertices(this->size);

        const glm::vec3 A = glm::euclidean(glm::vec2(M_PI_2, 0));
        for ( size_t t = 0; t < this->tr; ++t ){
            const glm::vec3 B = glm::euclidean(glm::vec2(0, M_PI_2 * t));
            const glm::vec3 C = glm::euclidean(glm::vec2(0, M_PI_2 * (t + 1)));
            for ( auto it = this->triangle_begin(t); it != this->triangle_end(t); ++it ){
                const auto rho        = static_cast<value_type>( it->rho );
                const auto rho_length = static_cast<value_type>( this->rho_size );
                const auto psi_length = static_cast<value_type>( this->psi_triangle_size(it->rho) );
                const auto psi        = static_cast<value_type>( it->psi % this->psi_triangle_size(it->rho) );

                value_type alpha ( 1 - rho / (rho_length - 1) );
                value_type beta  ( (1 - alpha) * ( 1 - psi / (psi_length - 1) ) );
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
    Sphere(size_t bin_splits): TriangleDiscreteCoordinates<T_INDEX>(4, bin_splits){}

    void get_viun(std::vector<glm::vec3> &vertices,
                  std::vector<T_INDEX>   &indeces,
                  std::vector<glm::vec2> &uvs,
                  std::vector<glm::vec3> &normals){
        auto verts = northern_hemisphere();
        // Add southern hemisphere vertices
        const size_t last_line = this->tr * (this->rho_size - 2) * (this->rho_size - 1) / 2 + 1;
        for ( size_t i = 0; i < last_line; ++i ){
            auto southern_vert = verts[i];
            southern_vert.y *= -1;
            verts.push_back(southern_vert);
        }

        vertices.insert(vertices.begin(), verts.begin(), verts.end());

        auto ind = this->get_small_triangles();
        // generate southern hemisphere triangles
        auto southern_ind = ind;
        for (auto &i : southern_ind){
            if (i < last_line){
                i += this->size;
            }
        }
        ind.insert(ind.end(), southern_ind.begin(), southern_ind.end());

        indeces.insert(indeces.begin(), ind.begin(), ind.end());

        uvs.clear();
        for ( auto &vert : verts ){
            const auto polar  = glm::polar(vert);
            const auto theta  = static_cast<value_type>(polar.x);
            const auto phi    = static_cast<value_type>(polar.y);
            const auto v      = static_cast<value_type>( (theta + M_PI_2) / M_PI );
            const auto u      = static_cast<value_type>( (phi + M_PI) / ( 2 * M_PI ) );
            uvs.emplace_back( u, v );
        }

        normals.clear();
        normals.insert(normals.begin(), vertices.begin(), vertices.end());
    }
};


template<typename T_INDEX>
class Circle: public TriangleDiscreteCoordinates<T_INDEX> {

public:
    Circle(size_t bin_splits): TriangleDiscreteCoordinates<T_INDEX>(6, bin_splits){}

    void get_viun(std::vector<glm::vec3> &vertices,
                  std::vector<T_INDEX>   &indeces,
                  std::vector<glm::vec2> &uvs,
                  std::vector<glm::vec3> &normals) {
        const auto ind = this->get_small_triangles();
        indeces.clear();
        indeces.insert(indeces.begin(), ind.begin(), ind.end());

        vertices.resize(this->size);
        uvs     .resize(this->size);
        normals .resize(this->size);
        for ( size_t t = 0; t < this->tr; ++t ) {
            for (auto it = this->triangle_begin(t); it != this->triangle_end(t); ++it) {
                const auto rho        = static_cast<value_type>( it->rho );
                const auto rho_length = static_cast<value_type>( this->rho_size );
                const auto psi_length = static_cast<value_type>( this->psi_size(it->rho));
                const auto psi        = static_cast<value_type>( it->psi % this->psi_size(it->rho));
                const auto r          = static_cast<value_type>( rho / (rho_length - 1) );
                const auto phi        = static_cast<value_type>( 2 * M_PI * psi / (psi_length - 1) );
                const auto x          = static_cast<value_type>( r * sin(phi) );
                const auto y          = static_cast<value_type>( 0 );
                const auto z          = static_cast<value_type>( r * cos(phi) );

                vertices[this->index(*it)] = glm::vec3(x, y, z);
                normals [this->index(*it)] = glm::vec3(0, 1, 0);
                uvs     [this->index(*it)] = glm::vec2(r, static_cast<value_type>(psi / (psi_length - 1)));
            }
        }
    }
};



#endif //TUTORIALS_DESCRITETOGEOMETRIC_HPP
