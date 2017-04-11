#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>

#include "BinaryParameters.hpp"
#include "DescriteToGeometric.hpp"
#include "disks.hpp"
#include "LightCurve.hpp"
#include "LimbDarking.hpp"
#include "MvpFromInput.hpp"
#include "Renderer.hpp"

int main() {
    int return_code = 0;



//    const float star_position_x =  1.5f;
//    const float disk_position_x = -1.5f;
//    const LightSource light(disk_position_x, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 10.0f);
//    const auto limb_darking = limbDarking(-0.3f, 0.0f, 0.0f, 0.0f);
//    const Sphere star(4);
//    const StandardDisk disk(5, 1.0f, 0.1f, 1.125f, 1e-3);
//    const DiskBelt belt(disk);
//    const std::vector<ObjectModel>  oms { star.get_object_model(),  disk.get_object_model(),  belt.get_object_model() };
//    const std::vector<TextureImage> tis { star.get_texture_image(), disk.get_texture_image(), belt.get_texture_image() };
//    const MVP base_mvp_star(3, 0, 0, star_position_x);
//    MVP base_mvp_disk(3, 0, 0, disk_position_x);
//    base_mvp_disk.change_rotation(static_cast<float>(M_PI / 36), static_cast<float>(M_PI_4));
//
//    try{
//        Renderer r(512, 512, oms, tis, light, limb_darking, true);
//        MvpFromInput mvp_input_star( r.window, base_mvp_star );
//        MvpFromInput mvp_input_disk( r.window, base_mvp_disk );
//        while( true ){
//            const auto start = std::chrono::high_resolution_clock::now();
//            const auto rgb = r.get_rgb_fluxes( {mvp_input_star.get(), mvp_input_disk.get(), mvp_input_disk.get()} );
//            const auto end = std::chrono::high_resolution_clock::now();
//            const float distance = mvp_input_star.get_distance();
//            std::cout << rgb[0] * distance * distance << "\t"
//                      << rgb[1] * distance * distance << "\t"
//                      << rgb[2] * distance * distance << "\t"
//                      << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/1000. << "ms"
//                      << std::endl;
//        }
//    } catch ( GlfwException ){
//        return_code = -1;
//    }



    const size_t phases = 100;
    try{
        const BinaryParameters bp(
                static_cast<float>(3 * SOLAR_RADIUS), // a
                0 * 2e33f, // Mx
                1 * 2e33f, // Mstar
                static_cast<float>(90.0 / 180.0 * M_PI), // i
                6000, // Tstar
                static_cast<float>(SOLAR_RADIUS), // Rstar
                0 * 5000, // Tdisk
                0 * 3e11f, // Rdisk
                0 * 0.05f, // z0Rdisk
                0 * 1e37f // Lx
        );
        const LightCurve lc(bp);

        const auto start = std::chrono::high_resolution_clock::now();
        const auto dots = lc.calc(phases);
        const auto end = std::chrono::high_resolution_clock::now();

        std::ofstream output("dots.dat");
        for ( auto &x : dots ) {
            output << x << "\n";
        }

        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/1000.
                  << "\t"
                  << *std::max_element(dots.begin(), dots.end()) / *std::min_element(dots.begin(), dots.end()) - 1
                  << "\t"
                  << fabs(M_PI * bp.Rstar*bp.Rstar * SIGMA_SB * powf(bp.Tstar,4) / dots.front() - 1 ) * 100
                  << "%"
                  << std::endl;
    } catch ( GlfwException ){
        return_code = -1;
    }



    return return_code;
}