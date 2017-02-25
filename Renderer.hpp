//
// Created by Konstantin Malanchev on 25/02/2017.
//

#ifndef HERX1_RENDERER_HPP
#define HERX1_RENDERER_HPP


#include <exception>
#include <iostream>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "ObjectModel.hpp"
#include "LightSource.hpp"
#include "shader.hpp"


struct GlfwException: public std::runtime_error{
    GlfwException(const char *what_arg):        std::runtime_error(what_arg){}
    GlfwException(const std::string &what_arg): std::runtime_error(what_arg){}
};


class Renderer{
protected:
    int frame_width;
    int frame_height;

    GLuint vertex_arrayID;

    // IDs for render into buffer
    GLuint programID;
    GLint MVP_ID, MV_ID, View_ID, Model_ID, LightID;

    // IDs for showing on screen
    GLuint quad_programID;
    GLint rendered_textureID;
    GLuint quad_vertex_buffer;
    GLuint quad_uv_buffer;
    // The framebuffer, which regroups 0, 1, or more textures, and 0 or 1 depth buffer.
    GLuint framebuffer_name = 0;
    GLuint rendered_texture;
    GLenum draw_buffers[1];

    // Model object attributes
    std::vector<GLuint> vertex_buffers;
    std::vector<GLuint> uv_buffers;
    std::vector<GLuint> normal_buffers;
    std::vector<GLuint> element_buffers;

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

public:
    const unsigned short window_width;
    const unsigned short window_height;
    const std::vector<ObjectModel> object_models;
    const LightSource light_source;

    const float scale0 = 3.0f;
    const float disk_position_x = -1.5f;
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
            light_source(light_source){
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

        // TODO: Can we work without window?
        window = glfwCreateWindow( window_width, window_height, "Playground", NULL, NULL);
        if( window == NULL ){
            throw GlfwException("Failed to open GLFW window");
        }
        glfwGetFramebufferSize( window, &frame_height, &frame_width );
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

        // Black background
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

        // TODO: ??
        // Create Vertex Array Object (VAO)
        glGenVertexArrays(1, &vertex_arrayID);
        glBindVertexArray(vertex_arrayID);

        // TODO: rewrite shader.?pp
        programID = LoadShaders(
                "VertexShader.glsl",
                "FragmentShader.glsl"
        );
        MVP_ID   = glGetUniformLocation(programID, "MVP");
        MV_ID    = glGetUniformLocation(programID, "MV");
        View_ID  = glGetUniformLocation(programID, "V");
        Model_ID = glGetUniformLocation(programID, "M");
        LightID  = glGetUniformLocation(programID, "LightPosition_worldspace");
        vertex_buffers .resize(object_models.size());
        uv_buffers     .resize(object_models.size());
        normal_buffers .resize(object_models.size());
        element_buffers.resize(object_models.size());
        for ( size_t i = 0; i < object_models.size(); ++i ){
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
                    &object_models[i].elements[i],
                    GL_STATIC_DRAW
            );
        }

        quad_programID = LoadShaders(
                "TextureVertexShader.glsl",
                "TextureFragmentShader.glsl"
        );
        rendered_textureID = glGetUniformLocation(quad_programID, "renderedTexture");
        glGenBuffers(1, &quad_vertex_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, quad_vertex_buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_vertex_buffer_data), g_quad_vertex_buffer_data, GL_STATIC_DRAW);
        glGenBuffers(1, &quad_uv_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, quad_uv_buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_uv_buffer_data), g_quad_uv_buffer_data, GL_STATIC_DRAW);
        glGenFramebuffers(1, &framebuffer_name);
        glBindFramebuffer(GL_FRAMEBUFFER, framebuffer_name);
        glGenTextures(1, &rendered_texture);
        glBindTexture(GL_TEXTURE_2D, rendered_texture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, rendered_texture, 0);
        draw_buffers[0] = GL_COLOR_ATTACHMENT0;
        glDrawBuffers(1, draw_buffers); // "1" is the size of DrawBuffers
        if ( glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE ) {
            GlfwException("Cannot create fram buffer");
        }
    }

    ~Renderer(){
        for ( size_t i = 0; i < object_models.size(); ++i ){
            glDeleteBuffers(1, &vertex_buffers[i]);
            glDeleteBuffers(1, &uv_buffers[i]);
            glDeleteBuffers(1, &normal_buffers[i]);
            glDeleteBuffers(1, &element_buffers[i]);
        }
        glDeleteProgram(programID);
        glDeleteProgram(quad_programID);
        glDeleteVertexArrays(1, &vertex_arrayID);
        glDeleteFramebuffers(1, &framebuffer_name);
        glDeleteTextures(1, &rendered_texture);
        glDeleteBuffers(1, &quad_vertex_buffer);

        glfwTerminate();
//        delete window;
    }
};

constexpr GLfloat Renderer::g_quad_vertex_buffer_data[];
constexpr GLfloat Renderer::g_quad_uv_buffer_data[];


#endif //HERX1_RENDERER_HPP
