//
// Created by Konstantin Malanchev on 21/02/2017.
//


#ifndef TUTORIALS_TRIANGLEDESCRETECOORDINATES_HPP
#define TUTORIALS_TRIANGLEDESCRETECOORDINATES_HPP


#include <boost/iterator/filter_iterator.hpp>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "commons.hpp"
#include "ObjectModel.hpp"
#include "TextureImage.hpp"


namespace discostar {
namespace geometry {


template <typename T_INDEX>
class Basic3DObject {
public:
    virtual ObjectModel<T_INDEX> get_object_model() const { return ObjectModel<T_INDEX>(); }

    virtual TextureImage get_texture_image() const { return TextureImage(1, 0.1f); }
};


template <typename T_INDEX>
struct DiscreteCoordinate {
	const T_INDEX rho;
	const T_INDEX psi;

	DiscreteCoordinate<T_INDEX>(T_INDEX rho, T_INDEX psi) : rho(rho), psi(psi) {};
};


template <typename T_INDEX>
class TriangleDiscreteCoordinates : public Basic3DObject<T_INDEX> {
public:
    typedef typename std::vector< DiscreteCoordinate<T_INDEX> > dcvec;
    typedef typename dcvec::const_iterator dcvec_const_it;

protected:
    dcvec idx;

public:
    const T_INDEX tr;
    const T_INDEX splits;
    const T_INDEX rho_size;
    const T_INDEX size;

    TriangleDiscreteCoordinates<T_INDEX>(T_INDEX triangles, T_INDEX bin_splits) :
            tr(triangles),
            splits(bin_splits),
            rho_size(discostar::exp2_int(splits) + 1),
            size(tr * rho_size * (rho_size - 1) / 2 + 1) {
        for (T_INDEX rho = 0; rho < rho_size; ++rho) {
            for (T_INDEX psi = 0; psi < psi_size(rho); ++psi) {
                idx.emplace_back(rho, psi);
            }
        }
    }

    dcvec_const_it begin() const {
        return idx.cbegin();
    }

    dcvec_const_it end() const {
        return idx.cend();
    }

    boost::filter_iterator<std::function<bool(DiscreteCoordinate<T_INDEX>)>, dcvec_const_it>
    triangle_begin(T_INDEX triangle) const {
        std::function<bool(DiscreteCoordinate<T_INDEX>)> predicator = [this, triangle](DiscreteCoordinate<T_INDEX> x) {
            return this->what_triangle(x) == triangle;
        };
        return boost::make_filter_iterator(
                predicator,
                begin(), end()
        );
    }

    boost::filter_iterator<std::function<bool(DiscreteCoordinate<T_INDEX>)>, dcvec_const_it>
    triangle_end(T_INDEX triangle) const {
        std::function<bool(DiscreteCoordinate<T_INDEX>)> predicator = [this, triangle](DiscreteCoordinate<T_INDEX> x) {
            return this->what_triangle(x) == triangle;
        };
        return boost::make_filter_iterator(
                predicator,
                end(), end()
        );
    }

    DiscreteCoordinate<T_INDEX> coordinate(T_INDEX index) const {
        return idx[index];
    }

    T_INDEX what_triangle(T_INDEX rho, T_INDEX psi) const {
        psi %= psi_size(rho);
        return psi / psi_triangle_size(rho);
    }

    T_INDEX what_triangle(const DiscreteCoordinate<T_INDEX> &dc) const {
        return what_triangle(dc.rho, dc.psi);
    }

    T_INDEX psi_size(T_INDEX rho) const {
        if (rho == 0) {
            return 1;
        }
        return tr * rho;
    }

    T_INDEX psi_triangle_size(T_INDEX rho) const {
        if (rho == 0) {
            return 1;
        }
        return rho;
    }

    T_INDEX index(T_INDEX rho, T_INDEX psi) const {
        if (rho == 0) {
            return 0;
        }
        psi %= psi_size(rho);
        return tr * rho * (rho - 1) / 2 + psi + 1;
    }

    T_INDEX index(const DiscreteCoordinate<T_INDEX> &dc) const {
        return index(dc.rho, dc.psi);
    }

    std::vector<T_INDEX> get_elements() const {
        std::vector<T_INDEX> triangles;
        for (T_INDEX rho = 1; rho < rho_size; ++rho) {
            for (T_INDEX psi = 0; psi < psi_size(rho); ++psi) {
                triangles.push_back(index(rho, psi));
                triangles.push_back(index(rho, psi + 1));
                triangles.push_back(index(rho - 1, psi - what_triangle(rho, psi)));
                if (rho < rho_size - 1) {
                    triangles.push_back(index(rho, psi + 1));
                    triangles.push_back(index(rho, psi));
                    triangles.push_back(index(rho + 1, psi + what_triangle(rho, psi) + 1));
                }
            }
        }
        return triangles;
    }
};

}} // namespace discostar::geometry

#endif //TUTORIALS_TRIANGLEDESCRETECOORDINATES_HPP
