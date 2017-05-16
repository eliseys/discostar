//
// Created by Konstantin Malanchev on 15/03/2017.
//

#ifndef HERX1_TEXTUREIMAGE_HPP
#define HERX1_TEXTUREIMAGE_HPP


#include <vector>

#include "commons.hpp"


//typedef typename std::vector<float> TextureImage;


namespace discostar {
namespace geometry {

class TextureImage : public std::vector<float> {
public:
    const size_t n;

    TextureImage(size_t n) : n(n), vector(n * n * 3) {};

    TextureImage(size_t n, float v) : n(n), vector(n * n * 3, v) {};

    TextureImage(size_t n, float r, float g, float b) :
            TextureImage(n) {
        for (size_t i = 0; i < n; ++i) {
            this->at(3 * i + 0) = r;
            this->at(3 * i + 1) = g;
            this->at(3 * i + 2) = b;
        }
    };
};

}} // namespace discostar::geometry

#endif //HERX1_TEXTUREIMAGE_HPP
