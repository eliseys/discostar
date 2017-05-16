//
// Created by Konstantin Malanchev on 25/02/2017.
//

#ifndef HERX1_OBJECTMODEL_HPP
#define HERX1_OBJECTMODEL_HPP


#include <vector>

#include <glm/glm.hpp>


struct ObjectModel{
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec2> uvs;
    std::vector<glm::vec3> normals;
    std::vector<unsigned short> elements;
};


#endif //HERX1_OBJECTMODEL_HPP
