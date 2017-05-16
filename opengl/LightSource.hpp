//
// Created by Konstantin Malanchev on 25/02/2017.
//

#ifndef HERX1_LIGHTSOURCE_HPP
#define HERX1_LIGHTSOURCE_HPP


#include <glm/glm.hpp>


struct LightSource{
    const glm::vec3 position;
    const glm::vec3 color;
    LightSource(const glm::vec3 &position, const glm::vec3 &color, float power=1):
            position( position ),
            color( power * color ){}
    LightSource(float x, float y, float z, float r, float g, float b, float power=1):
            LightSource( glm::vec3(x,y,z), glm::vec3(r,g,b), power ){}
};


#endif //HERX1_LIGHTSOURCE_HPP
