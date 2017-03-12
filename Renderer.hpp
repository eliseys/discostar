//
// Created by Konstantin Malanchev on 25/02/2017.
//

#ifndef HERX1_RENDERER_HPP
#define HERX1_RENDERER_HPP


#include <array>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "LightSource.hpp"
#include "MVP.hpp"
#include "ObjectModel.hpp"
#include "shader.hpp"


struct GlfwException: public std::runtime_error{
    GlfwException(const char *what_arg):        std::runtime_error(what_arg){}
    GlfwException(const std::string &what_arg): std::runtime_error(what_arg){}
};


class Renderer{
protected:
    int color_frame_width, color_frame_height;
    int shadow_frame_size;

    GLuint vertex_arrayID;

    // IDs for shadow map
    GLuint shadow_programID;
    GLint  shadow_MVP_ID;
    GLuint shadow_framebuffer_name = 0;
    GLuint rendered_shadow_texture;

    // IDs for render into buffer
    GLuint programID;
    GLint MVP_ID, MV_ID, View_ID, Model_ID, LightPosID, LightColID, ShadowBiasID, ShadowMapID;

    // IDs for color frame buffer
    GLuint quad_programID;
    GLint  quad_textureID;
    GLuint quad_vertex_buffer;
    GLuint quad_uv_buffer;
    // The framebuffer, which regroups 0, 1, or more textures, and 0 or 1 depth buffer.
    GLuint color_framebuffer_name = 0;
    GLuint rendered_color_texture;
    GLuint depth_render_buffer;
    GLenum draw_buffers[1];

    // Model object attributes
    std::vector<GLuint> vertex_buffers;
    std::vector<GLuint> uv_buffers;
    std::vector<GLuint> normal_buffers;
    std::vector<GLuint> element_buffers;

    const glm::mat4 biasMatrix = glm::mat4(
        0.5f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.5f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.5f, 0.0f,
        0.5f, 0.5f, 0.5f, 1.0f
    );

    // Screen model attributes
    static constexpr GLfloat g_quad_vertex_buffer_data[] = {
            -1.0f, -1.0f, 0.0f,
            1.0f, -1.0f, 0.0f,
            -1.0f,  1.0f, 0.0f,
            -1.0f,  1.0f, 0.0f,
            1.0f, -1.0f, 0.0f,
            1.0f,  1.0f, 0.0f,
    };
    static constexpr GLfloat g_quad_uv_buffer_data[] = {
            0.0f, 0.0f,
            1.0f, 0.0f,
            0.0f, 1.0f,
            0.0f, 1.0f,
            1.0f, 0.0f,
            1.0f, 1.0f,
    };

    void initialize_glfw(){
        // Initialise GLFW
        if( not glfwInit() ){
            throw GlfwException("Failed to initialize GLFW");
        }

        glfwWindowHint(GLFW_SAMPLES, 4);
        glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
//    	glfwWindowHint(GLFW_DECORATED, GL_TRUE);
//    	glfwWindowHint(GLFW_VISIBLE, GL_TRUE);
//    	glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);

        // TODO: Can we work on headless systems?
        window = glfwCreateWindow( window_width, window_height, "Binary", NULL, NULL);
        if( window == NULL ){
            throw GlfwException("Failed to open GLFW window");
        }
        glfwGetFramebufferSize( window, &color_frame_height, &color_frame_width );
        shadow_frame_size = (color_frame_height > color_frame_width ? color_frame_height : color_frame_width) * shadow_to_color_size;
        glfwMakeContextCurrent(window);

        // Initialize GLEW
        glewExperimental = static_cast<GLboolean>(true); // Needed for core profile
        if (glewInit() != GLEW_OK) {
            throw GlfwException("Failed to initialize GLEW");
        }

        // Enable depth test
        glEnable(GL_DEPTH_TEST);
        // Accept fragment if it closer to the camera than the former one
        glDepthFunc(GL_LESS);
        // Cull triangles which normal is not towards the camera
//        glEnable(GL_CULL_FACE);

        // TODO: WTF?
        GLuint VertexArrayID;
        glGenVertexArrays(1, &VertexArrayID);
        glBindVertexArray(VertexArrayID);
    }

    void initialize_shadow_renderer(){
        // Depth shadow shadow frame buffer
        shadow_programID = LoadShaders(
                "shaders/ShadowVertexShader.glsl",
                "shaders/ShadowFragmentShader.glsl"
        );
        shadow_MVP_ID = glGetUniformLocation(shadow_programID, "shadowMVP");
        glGenFramebuffers(1, &shadow_framebuffer_name);
        glBindFramebuffer(GL_FRAMEBUFFER, shadow_framebuffer_name);
        glGenTextures(1, &rendered_shadow_texture);
        glBindTexture(GL_TEXTURE_2D, rendered_shadow_texture);
        glTexImage2D(
                GL_TEXTURE_2D,
                0,
                GL_DEPTH_COMPONENT,
                shadow_frame_size,
                shadow_frame_size,
                0,
                GL_DEPTH_COMPONENT,
                GL_FLOAT,
                0
        );
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, rendered_shadow_texture, 0);
        glDrawBuffer(GL_NONE);
        glReadBuffer(GL_NONE);
        if ( glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE ) {
            GlfwException("Cannot create shadow frame buffer");
        }
    }

    void initialize_object_buffers(){
        for (size_t i = 0; i < object_models.size(); ++i) {
            glGenBuffers(1, &vertex_buffers[i]);
            glBindBuffer(GL_ARRAY_BUFFER, vertex_buffers[i]);
            glBufferData(
                    GL_ARRAY_BUFFER,
                    object_models[i].vertices.size() * sizeof(glm::vec3),
                    &object_models[i].vertices[0],
                    GL_STATIC_DRAW
            );
            glGenBuffers(1, &uv_buffers[i]);
            glBindBuffer(GL_ARRAY_BUFFER, uv_buffers[i]);
            glBufferData(
                    GL_ARRAY_BUFFER,
                    object_models[i].uvs.size() * sizeof(glm::vec2),
                    &object_models[i].uvs[0],
                    GL_STATIC_DRAW
            );
            glGenBuffers(1, &normal_buffers[i]);
            glBindBuffer(GL_ARRAY_BUFFER, normal_buffers[i]);
            glBufferData(
                    GL_ARRAY_BUFFER,
                    object_models[i].normals.size() * sizeof(glm::vec3),
                    &object_models[i].normals[0],
                    GL_STATIC_DRAW
            );
            glGenBuffers(1, &element_buffers[i]);
            glBindBuffer(GL_ARRAY_BUFFER, element_buffers[i]);
            glBufferData(
                    GL_ARRAY_BUFFER,
                    object_models[i].elements.size() * sizeof(unsigned short),
                    &object_models[i].elements[0],
                    GL_STATIC_DRAW
            );
        }
    }

    void initialize_color_renderer() {
        // TODO: rewrite shader.?pp
        programID = LoadShaders(
                "shaders/ColorVertexShader.glsl",
                "shaders/ColorFragmentShader.glsl"
        );
        MVP_ID = glGetUniformLocation(programID, "MVP");
        MV_ID = glGetUniformLocation(programID, "MV");
        View_ID = glGetUniformLocation(programID, "V");
        Model_ID = glGetUniformLocation(programID, "M");
        LightPosID = glGetUniformLocation(programID, "LightPosition_worldspace");
        LightColID = glGetUniformLocation(programID, "LightColor");
        ShadowBiasID = glGetUniformLocation(programID, "DepthBiasMVP");
        ShadowMapID = glGetUniformLocation(programID, "shadowMap");
        glGenFramebuffers(1, &color_framebuffer_name);
        glBindFramebuffer(GL_FRAMEBUFFER, color_framebuffer_name);
        glGenTextures(1, &rendered_color_texture);
        glBindTexture(GL_TEXTURE_2D, rendered_color_texture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, color_frame_width, color_frame_height, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glGenRenderbuffers(1, &depth_render_buffer);
        glBindRenderbuffer(GL_RENDERBUFFER, depth_render_buffer);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, color_frame_width, color_frame_height);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth_render_buffer);
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, rendered_color_texture, 0);
        draw_buffers[0] = GL_COLOR_ATTACHMENT0;
        glDrawBuffers(1, draw_buffers); // "1" is the size of DrawBuffers
        if ( glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE ) {
            GlfwException("Cannot create color frame buffer");
        }
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    void initialize_window_renderer(){
        quad_programID = LoadShaders(
                "shaders/TextureVertexShader.glsl",
                "shaders/TextureFragmentShader.glsl"
        );
        quad_textureID = glGetUniformLocation(quad_programID, "renderedTexture");
        glGenBuffers(1, &quad_vertex_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, quad_vertex_buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_vertex_buffer_data), g_quad_vertex_buffer_data, GL_STATIC_DRAW);
        glGenBuffers(1, &quad_uv_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, quad_uv_buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_uv_buffer_data), g_quad_uv_buffer_data, GL_STATIC_DRAW);
    }

public:
    const unsigned short window_width;
    const unsigned short window_height;
    const std::vector<ObjectModel> object_models;
    const LightSource light_source;

    const int shadow_to_color_size = 1;

    GLFWwindow *window;

    Renderer(
            unsigned short width,
            unsigned short height,
            const std::vector<ObjectModel> &object_models,
            const LightSource &light_source
    ) throw(GlfwException):
            window_width(width),
            window_height(height),
            object_models(object_models),
            light_source(light_source),
            vertex_buffers (object_models.size()),
            uv_buffers     (object_models.size()),
            normal_buffers (object_models.size()),
            element_buffers(object_models.size())
    {
        initialize_glfw();

        // Black background
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

        // TODO: ??
        // Create Vertex Array Object (VAO)
        glGenVertexArrays(1, &vertex_arrayID);
        glBindVertexArray(vertex_arrayID);

        initialize_object_buffers();
        initialize_shadow_renderer();
        initialize_color_renderer();
        initialize_window_renderer();
    }


    std::array<double, 3> get_rgb_fluxes(const std::vector<MVP> &mvps){
        std::array<double, 3> fluxes{0,0,0};

        // TODO: move somewhere
        const glm::mat4 shadow_projection = glm::perspective( 90.0f, 1.0f, 0.1f, 10.0f );
//        const glm::mat4 shadow_projection = glm::ortho( -10.0f, 10.0f, -10.0f, 10.0f, -10.0f, 20.0f );
        const glm::mat4 shadow_view = glm::lookAt(light_source.position, glm::vec3(0,0,0), glm::vec3(0,1,0));

        // Render shadow depth map
        glBindFramebuffer(GL_FRAMEBUFFER, shadow_framebuffer_name);
        glViewport(0, 0, shadow_frame_size, shadow_frame_size);
        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        // Use our shader
        glUseProgram(shadow_programID);
        std::vector<glm::mat4> shadow_mvps;
        for ( size_t i = 0; i < object_models.size(); ++i ) {
            shadow_mvps.push_back( shadow_projection * shadow_view * mvps[i].get_model() );
            glUniformMatrix4fv(shadow_MVP_ID, 1, GL_FALSE, &shadow_mvps.back()[0][0]);

            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, vertex_buffers[i]);
            glVertexAttribPointer(
                    0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                    3,                  // size
                    GL_FLOAT,           // type
                    GL_FALSE,           // normalized?
                    0,                  // stride
                    (void *) 0          // array buffer offset
            );

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffers[i]);
            glDrawElements(
                    GL_TRIANGLES,                               // mode
                    (GLsizei) object_models[i].elements.size(), // count
                    GL_UNSIGNED_SHORT,                          // type
                    (void *) 0                                  // element array buffer offset
            );

            glDisableVertexAttribArray(0);
        }

//        {
//            // Find nonzero pixels in the shadow texture
//            const auto pixels_size = static_cast<size_t>(shadow_frame_size * shadow_frame_size);
//            float pixels[pixels_size];
//            glGetTexImage(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT16, GL_FLOAT, pixels);
//            for ( size_t i = 0; i < pixels_size; ++i ) {
//                if ( pixels[i] != 0.0f ){
//                    std::cout << pixels[i] << std::endl;
//                }
//            }
//        }

        // Render color map
        glBindFramebuffer(GL_FRAMEBUFFER, color_framebuffer_name);
        glViewport(0, 0, color_frame_width, color_frame_height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(programID);
        glUniform3f(LightPosID, light_source.position.x, light_source.position.y, light_source.position.z);
        glUniform3f(LightColID, light_source.color   .x, light_source.color   .y, light_source.color   .z);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, rendered_shadow_texture);
        glUniform1i(ShadowMapID, 0);
        for ( size_t i = 0; i < object_models.size(); ++i ) {
            const auto shadow_bias_mvp = biasMatrix * shadow_mvps[i];

            glUniformMatrix4fv(MVP_ID,       1, GL_FALSE, &mvps[i].mvp()      [0][0]);
            glUniformMatrix4fv(MV_ID,        1, GL_FALSE, &mvps[i].mv()       [0][0]);
            glUniformMatrix4fv(View_ID,      1, GL_FALSE, &mvps[i].view       [0][0]);
            glUniformMatrix4fv(Model_ID,     1, GL_FALSE, &mvps[i].get_model()[0][0]);
            glUniformMatrix4fv(ShadowBiasID, 1, GL_FALSE, &shadow_bias_mvp    [0][0]);

            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, vertex_buffers[i]);
            glVertexAttribPointer(
                    0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                    3,                  // size
                    GL_FLOAT,           // type
                    GL_FALSE,           // normalized?
                    0,                  // stride
                    (void *) 0          // array buffer offset
            );

            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, uv_buffers[i]);
            glVertexAttribPointer(
                    1,                  // attribute
                    2,                  // size
                    GL_FLOAT,           // type
                    GL_FALSE,           // normalized?
                    0,                  // stride
                    (void *) 0          // array buffer offset
            );

            // 3rd attribute buffer : normals
            glEnableVertexAttribArray(2);
            glBindBuffer(GL_ARRAY_BUFFER, normal_buffers[i]);
            glVertexAttribPointer(
                    2,                  // attribute
                    3,                  // size
                    GL_FLOAT,           // type
                    GL_FALSE,           // normalized?
                    0,                  // stride
                    (void *) 0          // array buffer offset
            );

            // Index buffer
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffers[i]);
            // Draw the triangles!
            glDrawElements(
                    GL_TRIANGLES,                               // mode
                    (GLsizei) object_models[i].elements.size(), // count
                    GL_UNSIGNED_SHORT,                          // type
                    (void *) 0                                  // element array buffer offset
            );

            glDisableVertexAttribArray(0);
            glDisableVertexAttribArray(1);
            glDisableVertexAttribArray(2);
        }


        // Render to the screen
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, color_frame_width, color_frame_height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(quad_programID);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, rendered_color_texture);
        glUniform1i(quad_textureID, 0);

        // Calculate fluxes
        {
            // TODO: initialize once
            const auto pixels_size = static_cast<size_t>(color_frame_width * color_frame_height * fluxes.size());
            unsigned short pixels[pixels_size];
            glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_UNSIGNED_SHORT, pixels);
            const auto max = static_cast<double>(std::numeric_limits<unsigned short>::max());
//#pragma omp parallel for private(i) shared(pixels, max) reduction(+:flux_red)
            for (size_t i = 0; i < pixels_size; i += fluxes.size()) {
                for (size_t j = 0; j < fluxes.size(); ++j) {
                    fluxes[j] += static_cast<double>(pixels[i + j]) / max;
                }
            }
        }

        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, quad_vertex_buffer);
        glVertexAttribPointer(
                0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                3,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void *) 0          // array buffer offset
        );

        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, quad_uv_buffer);
        glVertexAttribPointer(
                1,                  // attribute 1
                2,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void *) 0          // array buffer offset
        );
        glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);


        glViewport(0, 0, color_frame_width/2, color_frame_height/2);
        glUseProgram(quad_programID);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, rendered_shadow_texture);
        glUniform1i(quad_textureID, 0);

        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, quad_vertex_buffer);
        glVertexAttribPointer(
                0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                3,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void *) 0          // array buffer offset
        );

        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, quad_uv_buffer);
        glVertexAttribPointer(
                1,                  // attribute 1.
                2,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void *) 0          // array buffer offset
        );
        glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);


        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();

        return fluxes;
    }


    ~Renderer(){
        for ( size_t i = 0; i < object_models.size(); ++i ){
            glDeleteBuffers(1, &vertex_buffers[i]);
            glDeleteBuffers(1, &uv_buffers[i]);
            glDeleteBuffers(1, &normal_buffers[i]);
            glDeleteBuffers(1, &element_buffers[i]);
        }
        glDeleteProgram(shadow_programID);
        glDeleteProgram(programID);
        glDeleteProgram(quad_programID);
        glDeleteVertexArrays(1, &vertex_arrayID);
        glDeleteFramebuffers(1, &color_framebuffer_name);
        glDeleteFramebuffers(1, &shadow_framebuffer_name);
        glDeleteTextures(1, &rendered_color_texture);
        glDeleteTextures(1, &rendered_shadow_texture);
        glDeleteBuffers(1, &quad_vertex_buffer);
        glDeleteBuffers(1, &quad_uv_buffer);

        glfwTerminate();
    }
};

constexpr GLfloat Renderer::g_quad_vertex_buffer_data[];
constexpr GLfloat Renderer::g_quad_uv_buffer_data[];


#endif //HERX1_RENDERER_HPP
