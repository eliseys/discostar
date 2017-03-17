#include <chrono>
#include <iostream>

#include "DescriteToGeometric.hpp"
#include "disks.hpp"
#include "LimbDarking.hpp"
#include "MvpFromInput.hpp"
#include "Renderer.hpp"

int main() {
    int return_code = 0;

    const float star_position_x =  1.5f;
    const float disk_position_x = -1.5f;
    const LightSource light(disk_position_x, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 10.0f);
    const auto limb_darking = limbDarking(-0.3f, 0.0f, 0.0f, 0.0f);
    const Sphere star(4);
    const StandardDisk disk(5, 1.0f, 0.1f, 1.125f, 1e-3);
    const DiskBelt belt(disk);
    const std::vector<ObjectModel>  oms { star.get_object_model(),  disk.get_object_model(),  belt.get_object_model() };
    const std::vector<TextureImage> tis { star.get_texture_image(), disk.get_texture_image(), belt.get_texture_image() };
    const MVP base_mvp_star(3, 0, 0, star_position_x);
    const MVP base_mvp_disk(3, 0, 0, disk_position_x);

    try{
        Renderer r(512, 512, oms, tis, light, limb_darking);
        MvpFromInput mvp_input_star( r.window, base_mvp_star );
        MvpFromInput mvp_input_disk( r.window, base_mvp_disk );
        while( true ){
            auto start = std::chrono::high_resolution_clock::now();
            auto rgb = r.get_rgb_fluxes( {mvp_input_star.get(), mvp_input_disk.get(), mvp_input_disk.get()} );
            auto end = std::chrono::high_resolution_clock::now();
            std::cout << rgb[0] << "\t"
                      << rgb[1] << "\t"
                      << rgb[2] << "\t"
                      << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/1000. << "ms"
                      << std::endl;
        }
    } catch ( GlfwException ){
        return_code = -1;
    }

    return return_code;
}