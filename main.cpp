#include <iostream>

#include "DescriteToGeometric.hpp"
#include "MvpFromInput.hpp"
#include "Renderer.hpp"

int main() {
    int return_code = 0;

    const float star_position_x =  1.5f;
    const float disk_position_x = -1.5f;
    const LightSource light(disk_position_x, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 10.0f);
    const Sphere star(1);
    const Circle disk(1);
    const std::vector<ObjectModel> oms { star.get_object_model(), disk.get_object_model() };
    const MVP base_mvp_star(3, 0, 0, star_position_x);
    const MVP base_mvp_disk(3, 0, 0, disk_position_x);

    try{
        Renderer r(512, 512, oms, light);
        MvpFromInput mvp_input_star( r.window, base_mvp_star );
        MvpFromInput mvp_input_disk( r.window, base_mvp_disk );
        while( true ){
            auto rgb = r.get_rgb_fluxes( {mvp_input_star.get(), mvp_input_disk.get()} );
            std::cout << rgb[0] << "\t" << rgb[1] << "\t" << rgb[2] << std::endl;
        }
    } catch ( GlfwException ){
        return_code = -1;
    }

    return return_code;
}