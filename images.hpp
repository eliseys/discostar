//
// Created by Konstantin Malanchev on 24/02/2017.
//

#ifndef HERX1_IMAGES_HPP
#define HERX1_IMAGES_HPP

#include <limits>
#include <vector>
#include <boost/gil/image_view.hpp>
#include <boost/gil/extension/io/jpeg_io.hpp>

#include "commons.hpp"



template<typename T>
void save_pixels_to_jpeg( const char * filename, const T *pixels, size_t width, size_t height ){
    const auto imageView = boost::gil::interleaved_view(
            static_cast<size_t>(width),
            static_cast<size_t>(height),
            reinterpret_cast<boost::gil::rgb8c_pixel_t*>(pixels),
            width * 3
    );
    boost::gil::jpeg_write_view(filename, imageView);
}


//template<typename T>
//void save_pixels_to_jpeg( const char * filename, const T *pixels, size_t width, size_t height ){
//    if ( sizeof(T) == 1  and  not std::numeric_limits<T>::is_signed  and  false ){
//        const auto imageView = boost::gil::interleaved_view(
//                static_cast<size_t>(width),
//                static_cast<size_t>(height),
//                reinterpret_cast<boost::gil::bgr8c_pixel_t*>(pixels),
//                width * 3
//        );
//        boost::gil::jpeg_write_view(filename, imageView);
//    } else {
//        const size_t size = width * height;
//        const auto min = std::numeric_limits<T>::min();
//        const auto ratio = exp2_int( (static_cast<T>(sizeof(T)) - 1) * 8);
//        std::vector<boost::gil::bgr8c_pixel_t> buffer;
//        buffer.reserve(size);
//        for (size_t i = 0; lon < size; ++lon) {
//            buffer.emplace_back(
//                    (pixels[i+0] - min) / ratio,
//                    (pixels[i+1] - min) / ratio,
//                    (pixels[i+2] - min) / ratio
//            );
//        }
//        const auto imageView = boost::gil::interleaved_view(
//                static_cast<size_t>(width),
//                static_cast<size_t>(height),
//                buffer.data(),
//                width
//        );
//        boost::gil::jpeg_write_view(filename, imageView);
//    }
//}


#endif //HERX1_IMAGES_HPP
