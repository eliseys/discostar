//
// Created by Konstantin Malanchev on 25/02/2017.
//

#ifndef HERX1_RENDERER_HPP
#define HERX1_RENDERER_HPP


#include <algorithm>
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

#ifdef ENABLE_OPENMP
#include "omp.h"
#endif // ENABLE_OPENMP


namespace discostar {

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
    GLint MVP_ID, MV_ID, View_ID, Model_ID, LightPosID, LightColID, ShadowBiasID, ShadowMapID, TextureID, DarkingID;

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
    std::vector<GLuint> textures;

    std::vector<float> pixels;

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
        window = glfwCreateWindow(window_width, window_height, "Binary", NULL, NULL);
        if (window == NULL) {
            throw GlfwException("Failed to open GLFW window");
        }
        glfwGetFramebufferSize(window, &color_frame_height, &color_frame_width);
        glfwMakeContextCurrent(window);
        shadow_frame_size = std::min(color_frame_height, color_frame_width) * shadow_to_color_size;
        const auto pixels_size = static_cast<size_t>(color_frame_width * color_frame_height * 3);
        pixels.resize(pixels_size);

        // Initialize GLEW
        glewExperimental = static_cast<GLboolean>(true); // Needed for core profile
        if (glewInit() != GLEW_OK) {
            throw GlfwException("Failed to initialize GLEW");
        }

        // Enable depth test
        glEnable(GL_DEPTH_TEST);
        // Enable Multisampling
        glEnable(GL_MULTISAMPLE);
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
        for ( size_t i = 0; i < object_models.size(); ++i ) {
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

    void texture_image_buffers(){
        for ( size_t i = 0; i < texture_images.size(); ++i ){
            glGenTextures(1, &textures[i]);
            glBindTexture(GL_TEXTURE_2D, textures[i]);
            glTexImage2D(
                    GL_TEXTURE_2D,
                    0,
                    GL_RGB32F,
                    static_cast<GLsizei>(texture_images[i].n),
                    static_cast<GLsizei>(texture_images[i].n),
                    0,
                    GL_RGB,
                    GL_FLOAT,
                    texture_images[i].data()
            );
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);  // u coordinate
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);  // v coordinate
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); // It is good!
            glGenerateMipmap(GL_TEXTURE_2D);
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
        TextureID = glGetUniformLocation(programID, "textureSampler");
        DarkingID = glGetUniformLocation(programID, "limbDarkingCoeffs");
        glGenFramebuffers(1, &color_framebuffer_name);
        glBindFramebuffer(GL_FRAMEBUFFER, color_framebuffer_name);
        glGenTextures(1, &rendered_color_texture);
        glBindTexture(GL_TEXTURE_2D, rendered_color_texture);
        glTexImage2D(
				GL_TEXTURE_2D,
				0,
				GL_RGB32F,
				color_frame_width,
				color_frame_height,
				0,
				GL_RGB,
				GL_FLOAT,
				0
		);
//        glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, 4, GL_RGB32F, color_frame_width, color_frame_height, GL_TRUE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glGenRenderbuffers(1, &depth_render_buffer);
        glBindRenderbuffer(GL_RENDERBUFFER, depth_render_buffer);
//        glRenderbufferStorageMultisample(GL_RENDERBUFFER, 4, GL_DEPTH_COMPONENT, color_frame_width, color_frame_height);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, color_frame_width, color_frame_height);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth_render_buffer);
//        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, rendered_color_texture, 0);
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, rendered_color_texture, 0);
        if ( glGetError() != GL_NO_ERROR ){
            throw GlfwException("Framebuffer error");
        }
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

    std::array<double, 3> calculate_rgb_fluxes(){
        std::array<double, 3> fluxes{0,0,0};
        glReadPixels(0, 0, color_frame_width, color_frame_height, GL_RGB, GL_FLOAT, pixels.data());
        size_t i;
        for( size_t j = 0; j < fluxes.size(); ++j ) {
            double flux = fluxes[j];
#ifdef ENABLE_OPENMP
#pragma omp parallel for private(i) shared(pixels, j) reduction(+:flux)
#endif // ENABLE_OPENMP
            for ( i = 0; i < pixels.size(); i += fluxes.size() ) {
                const double df = static_cast<double>(pixels[i + j]);
                flux += pow(df, 4);
            }
            fluxes[j] = flux / static_cast<double>( pixels.size() / fluxes.size() );
        }
        return fluxes;
    }

public:
    const unsigned short window_width;
    const unsigned short window_height;
    const std::vector<geometry::ObjectModel> object_models;
    const std::vector<geometry::TextureImage> texture_images;
    const light::LightSource light_source;
    const glm::vec4 limb_darking;

    const int shadow_to_color_size = 1;

    GLFWwindow *window;
    const bool show_in_window;

    Renderer(
            unsigned short width,
            unsigned short height,
            const std::vector<geometry::ObjectModel> &object_models,
            const std::vector<geometry::TextureImage> &texture_images,
            const light::LightSource &light_source,
            const glm::vec4 &limb_darking,
            bool show_in_windpw = true
    ) throw(GlfwException):
            window_width(width),
            window_height(height),
            object_models(object_models),
            texture_images(texture_images),
            light_source(light_source),
            limb_darking(limb_darking),
            show_in_window(show_in_windpw),
            vertex_buffers (object_models.size()),
            uv_buffers     (object_models.size()),
            normal_buffers (object_models.size()),
            element_buffers(object_models.size()),
            textures(texture_images.size())
    {
        initialize_glfw();

        // Black background
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

        // TODO: ??
        // Create Vertex Array Object (VAO)
        glGenVertexArrays(1, &vertex_arrayID);
        glBindVertexArray(vertex_arrayID);

        initialize_object_buffers();
        texture_image_buffers();
        initialize_shadow_renderer();
        initialize_color_renderer();
        initialize_window_renderer();
    }


    std::array<double, 3> get_rgb_fluxes(const std::vector<mvp::MVP> &mvps){
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
        if ( glGetError() != GL_NO_ERROR ){
            throw GlfwException("Shading error");
        }
        // Render color map
        glBindFramebuffer(GL_FRAMEBUFFER, color_framebuffer_name);
        glViewport(0, 0, color_frame_width, color_frame_height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(programID);
        glUniform3f(LightPosID, light_source.position.x, light_source.position.y, light_source.position.z);
        glUniform3f(LightColID, light_source.color   .x, light_source.color   .y, light_source.color   .z);
        glUniform4f(DarkingID, limb_darking[0], limb_darking[1], limb_darking[2], limb_darking[3]);
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

            glActiveTexture(GL_TEXTURE1);
            glBindTexture(GL_TEXTURE_2D, textures[i]);
            glUniform1i(TextureID, 1);

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

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, rendered_color_texture);
//        glBindFramebuffer(GL_READ_FRAMEBUFFER, color_framebuffer_name);
//        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, rendered_color_texture);
//        glBlitFramebuffer(0, 0, color_frame_width, color_frame_height, 0, 0, color_frame_width, color_frame_height, GL_COLOR_BUFFER_BIT, GL_LINEAR);

        auto fluxes = calculate_rgb_fluxes();

        // Render to the screen
        if ( show_in_window ) {
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            glViewport(0, 0, color_frame_width, color_frame_height);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glUseProgram(quad_programID);
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
        }

        // Swap buffers
        if( show_in_window ) {
            glfwSwapBuffers(window);
            glfwPollEvents();
        }

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

} // namespace discostar::renderer_opengl

#endif //HERX1_RENDERER_HPP
