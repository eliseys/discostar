//
// Created by Konstantin Malanchev on 15/03/2017.
//

#ifndef HERX1_TEXTUREIMAGE_HPP
#define HERX1_TEXTUREIMAGE_HPP


#include <vector>

#include "commons.hpp"


//typedef typename std::vector<float> TextureImage;


class TextureImage: public std::vector<float>{
public:
    const size_t n;

    TextureImage(size_t n):          n(n), vector(n * n * 3)   {};
    TextureImage(size_t n, float v): n(n), vector(n * n * 3, v){};
};


#endif //HERX1_TEXTUREIMAGE_HPP
