//
// Created by Konstantin Malanchev on 25/02/2017.
//

#ifndef HERX1_MVPFROMINPUT_HPP
#define HERX1_MVPFROMINPUT_HPP


#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "MVP.hpp"


double ScrollXOffset = 0.;
double ScrollYOffset = 0.;
void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset){
    ScrollXOffset += xoffset;
    ScrollYOffset += yoffset;
}


class MvpFromInput{
protected:
    const float mouseSpeed = 0.005f;
    GLFWwindow *window;
    int window_width, window_height;
    const float base_distance;
    const MVP base_mvp;

public:
    MvpFromInput(GLFWwindow *window, const MVP &base_mvp):
            window(window),
            base_distance( 1 / base_mvp.projection[0][0] ),
            base_mvp( base_mvp ){
        glfwGetWindowSize(window, &window_width, &window_height);
        glfwSetScrollCallback(window, ScrollCallback);
    }

    MVP get() const{
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        float lon = mouseSpeed * ( static_cast<float>(window_width/2)  - static_cast<float>(xpos) );
        float lat = mouseSpeed * ( static_cast<float>(window_height/2) - static_cast<float>(ypos) );
        float distance = (float) exp(ScrollYOffset) * base_distance;

        auto mvp = base_mvp;
        mvp.change_view(distance, lat, lon);
        return mvp;
    }

};


#endif //HERX1_MVPFROMINPUT_HPP
