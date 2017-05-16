//
// Created by Konstantin Malanchev on 25/02/2017.
//

#ifndef HERX1_MVP_HPP
#define HERX1_MVP_HPP


#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/polar_coordinates.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/transform.hpp>


namespace discostar {
namespace mvp {


class MVP {
public:
    typedef typename glm::mat4 mat_t;

    mat_t scale, rotation, translation, view, projection;

    mat_t get_model() const {
        return translation * rotation * scale;
    }

    mat_t mv() const {
        return view * get_model();
    }

    mat_t mvp() const {
        return projection * view * get_model();
    }

    MVP() :
            scale(1),
            rotation(1),
            translation(glm::translate(glm::vec3(0, 0, 0))),
            view(glm::lookAt(glm::vec3(0, 0, 1), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0))),
            projection(glm::ortho(-1.0f, 1.0f, -1.0f, 1.0f, 0.0f, 100.0f)) {}

    MVP(float distance, const glm::vec2 &camera_spherical_coordinates, const glm::vec3 &translate_vector) :
            MVP() {
        change_view(distance, camera_spherical_coordinates);
        change_translation(translate_vector);
    }

    MVP(float distance, float lat, float lon, float translate_x) :
            MVP(distance, glm::vec2(lat, lon), glm::vec3(translate_x, 0, 0)) {}

    void change_rotation(float theta, float phi) {
        const auto q_rotate_around_y = glm::angleAxis(phi, glm::vec3(0, 1, 0));
        const auto q_rotate_around_x = glm::angleAxis(theta, glm::vec3(1, 0, 0));
        rotation = glm::toMat4(q_rotate_around_y * q_rotate_around_x);
    }

    void change_translation(const glm::vec3 &translate_vector) {
        translation = glm::translate(translate_vector);
    }

    void change_translation(float x, float y = 0, float z = 0) {
        change_translation(glm::vec3(x, y, z));
    }

    void change_view(float distance, const glm::vec2 &camera_spherical_coordinates) {
        view = glm::lookAt(
                distance * euclidean(camera_spherical_coordinates),
                glm::vec3(0, 0, 0),
                glm::vec3(0, 1, 0)
        );
        projection = glm::ortho(
                -distance, distance,
                -distance, distance,
                0.0f, distance * 100.0f
        );
    }

    void change_view(float distance, float lat, float lon) {
        change_view(distance, glm::vec2(lat, lon));
    }
};

}} // namespace discostar::mvp

#endif //HERX1_MVP_HPP
