#include <iostream>

#include "DescriteToGeometric.hpp"
#include "Renderer.hpp"

int main() {
    int return_code = 0;

    LightSource light(-1.5f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 10.0f);
    Sphere star(4);
    Circle disk(4);
    std::vector<ObjectModel> oms { star.get_object_model(), disk.get_object_model() };

    try{
        Renderer r(512, 512, oms, light);
    } catch ( GlfwException ){
        return_code = 1;
    }

    return return_code;
}