#include <stdio.h>
#include <iostream>
#include <array>
#include <chrono>
#include <vector>

#include <GL/glew.h>

#include <GLFW/glfw3.h>
GLFWwindow* window;

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/polar_coordinates.hpp>
using namespace glm;

#include "shader.hpp"
#include "TriangleDiscreteCoordinates.hpp"
#include "DescriteToGeometric.hpp"
#include "images.hpp"




double ScrollXOffset = 0.;
double ScrollYOffset = 0.;
void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset){
    ScrollXOffset += xoffset;
    ScrollYOffset += yoffset;
}


class MVP{
protected:
	const float scale0;
	const float translate_x;
	float scale;
    const float mouseSpeed = 0.005f;
    float lat = 0.0f;
    float lon = 0.0f;

    void positionFromInput(){
        // Get coursor position
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        // Get window size
        int width, height;
        glfwGetWindowSize(window, &width, &height);

        lon += mouseSpeed * float(width/2  - xpos);
        lat += mouseSpeed * float(height/2 - ypos);

		// Set scale via scroll
        scale = (float) exp(ScrollYOffset) * scale0;
    }

public:
    MVP(float scale, float translate_x=0): scale0(scale), scale(scale), translate_x(translate_x){
		positionFromInput();
    };

	mat4 getProjectionMatrix(){
		return ortho(
				-scale,scale,
				-scale,scale,
				0.0f,100.0f
		);
	}

	mat4 getViewMatrix(vec2 spherical_coordinates){
		return lookAt(
				scale * euclidean(spherical_coordinates),
				vec3(0,0,0), // and looks at the origin
				vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
		);
	}

	mat4 getViewMatrix(){
		return getViewMatrix(vec2(lat, lon));
	}

	mat4 getTranslationMatrix(){
		return translate(vec3(translate_x,0,0));
	}

	mat4 getRotationMatrix(){
		return mat4(1);
	}

	mat4 getScaleMatrix(){
		return mat4(1);
	}

	mat4 getModelMatrix(){
		return getTranslationMatrix() * getRotationMatrix() * getScaleMatrix();
	}

	mat4 mv(){
		return getViewMatrix() * getModelMatrix();
	}

	mat4 mvp(){
		return getProjectionMatrix() * getViewMatrix() * getModelMatrix();
	}
};


int playground()
{
    const int widthWindow = 512;
    const int heightWindow = 512;
    const float scale0 = 3.0f;
    const float disk_position_x = -1.5f;
    const bool saveImage = false;

	// Initialise GLFW
	if( !glfwInit() )
	{
		fprintf( stderr, "Failed to initialize GLFW\n" );
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
//	glfwWindowHint(GLFW_DECORATED, GL_TRUE);
//	glfwWindowHint(GLFW_VISIBLE, GL_TRUE);
//	glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow( widthWindow, heightWindow, "Playground", NULL, NULL);
    // But on MacOS X with a retina screen it'll be larger
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
	if( window == NULL ){
		fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
    glewExperimental = static_cast<GLboolean>(true); // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	// Cull triangles which normal is not towards the camera
//	glEnable(GL_CULL_FACE);

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	glfwSetScrollCallback(window, ScrollCallback);

	// Black background
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    // Create Vertex Array Object (VAO)
    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    GLuint programID = LoadShaders(
            "VertexShader.glsl",
            "FragmentShader.glsl"
    );

	GLuint quad_programID = LoadShaders(
			"TextureVertexShader.glsl",
			"TextureFragmentShader.glsl"
	);
	GLint texID = glGetUniformLocation(quad_programID, "renderedTexture");
	// The fullscreen quad's FBO
	static const GLfloat g_quad_vertex_buffer_data[] = {
			-1.0f, -1.0f, 0.0f,
			1.0f, -1.0f, 0.0f,
			-1.0f,  1.0f, 0.0f,
			-1.0f,  1.0f, 0.0f,
			1.0f, -1.0f, 0.0f,
			1.0f,  1.0f, 0.0f,
	};
	GLuint quad_vertexbuffer;
	glGenBuffers(1, &quad_vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_vertex_buffer_data), g_quad_vertex_buffer_data, GL_STATIC_DRAW);
	static const GLfloat g_quad_uv_buffer_data[] = {
			0.0f, 0.0f,
			1.0f, 0.0f,
			0.0f, 1.0f,
			0.0f, 1.0f,
			1.0f, 0.0f,
			1.0f, 1.0f,
	};
	GLuint quad_uvbuffer;
	glGenBuffers(1, &quad_uvbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, quad_uvbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_uv_buffer_data), g_quad_uv_buffer_data, GL_STATIC_DRAW);


	GLint MVP_ID =   glGetUniformLocation(programID, "MVP");
	GLint MV_ID =    glGetUniformLocation(programID, "MV");
	GLint View_ID =  glGetUniformLocation(programID, "V");
	GLint Model_ID = glGetUniformLocation(programID, "M");
	GLint LightID =  glGetUniformLocation(programID, "LightPosition_worldspace");

	glm::vec3 lightPos = glm::vec3(disk_position_x,0,0);

	std::vector<vec3>           star_vertices;
	std::vector<unsigned short> star_indeces;
	std::vector<vec2>           star_uvs;
	std::vector<vec3>           star_normals;
    auto star = Sphere<unsigned  short>(4);
	star.get_viun(star_vertices, star_indeces, star_uvs, star_normals);
	GLuint star_vertex_buffer;
	glGenBuffers(1, &star_vertex_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, star_vertex_buffer);
	glBufferData(GL_ARRAY_BUFFER, star_vertices.size() * sizeof(vec3), &star_vertices[0], GL_STATIC_DRAW);
	GLuint star_index_buffer;
	glGenBuffers(1, &star_index_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, star_index_buffer);
	glBufferData(GL_ARRAY_BUFFER, star_indeces.size() * sizeof(unsigned short), &star_indeces[0], GL_STATIC_DRAW);
	GLuint star_uv_buffer;
	glGenBuffers(1, &star_uv_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, star_uv_buffer);
	glBufferData(GL_ARRAY_BUFFER, star_uvs.size() * sizeof(vec2), &star_uvs[0], GL_STATIC_DRAW);
	GLuint star_normal_buffer;
	glGenBuffers(1, &star_normal_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, star_normal_buffer);
	glBufferData(GL_ARRAY_BUFFER, star_normals.size() * sizeof(vec3), &star_normals[0], GL_STATIC_DRAW);

	std::vector<vec3>           disk_vertices;
	std::vector<unsigned short> disk_indeces;
	std::vector<vec2>           disk_uvs;
	std::vector<vec3>           disk_normals;
	auto disk = Circle<unsigned  short>(4);
	disk.get_viun(disk_vertices, disk_indeces, disk_uvs, disk_normals);
	GLuint disk_vertex_buffer;
	glGenBuffers(1, &disk_vertex_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, disk_vertex_buffer);
	glBufferData(GL_ARRAY_BUFFER, disk_vertices.size() * sizeof(vec3), &disk_vertices[0], GL_STATIC_DRAW);
	GLuint disk_index_buffer;
	glGenBuffers(1, &disk_index_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, disk_index_buffer);
	glBufferData(GL_ARRAY_BUFFER, disk_indeces.size() * sizeof(unsigned short), &disk_indeces[0], GL_STATIC_DRAW);
	GLuint disk_uv_buffer;
	glGenBuffers(1, &disk_uv_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, disk_uv_buffer);
	glBufferData(GL_ARRAY_BUFFER, disk_uvs.size() * sizeof(vec2), &disk_uvs[0], GL_STATIC_DRAW);
	GLuint disk_normal_buffer;
	glGenBuffers(1, &disk_normal_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, disk_normal_buffer);
	glBufferData(GL_ARRAY_BUFFER, disk_normals.size() * sizeof(vec3), &disk_normals[0], GL_STATIC_DRAW);

	// The framebuffer, which regroups 0, 1, or more textures, and 0 or 1 depth buffer.
	GLuint FramebufferName = 0;
	glGenFramebuffers(1, &FramebufferName);
	glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
	// The texture we're going to render to
	GLuint renderedTexture;
	glGenTextures(1, &renderedTexture);
	// "Bind" the newly created texture : all future texture functions will modify this texture
	glBindTexture(GL_TEXTURE_2D, renderedTexture);
	// Give an empty image to OpenGL ( the last "0" means "empty" )
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
	// Poor filtering
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	// From tutorial 14 source file
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	// The depth buffer - optional
	GLuint depthrenderbuffer;
	glGenRenderbuffers(1, &depthrenderbuffer);
	glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, width, height);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);
	// Set "renderedTexture" as our colour attachement #0
	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, renderedTexture, 0);
	// Set the list of draw buffers.
	GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
	glDrawBuffers(1, DrawBuffers); // "1" is the size of DrawBuffers
	// Always check that our framebuffer is ok
	if ( glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE ) {
		return false;
	}

	do{
		// Render to our framebuffer
		glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
		glViewport(0, 0, width, height); // Render on the whole framebuffer, complete from the lower left corner to the upper right

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Use our shader
		glUseProgram(programID);

		glUniform3f(LightID, lightPos.x, lightPos.y, lightPos.z);

		{
			auto mvp = MVP(scale0, 1.5f);

			// Send our transformation to the currently bound shader, in the "MVP" uniform
			// This is done in the main loop since each model will have a different MVP matrix (At least for the M part)
			glUniformMatrix4fv(MVP_ID, 1, GL_FALSE, &mvp.mvp()[0][0]);
			glUniformMatrix4fv(MV_ID, 1, GL_FALSE, &mvp.mv()[0][0]);
			glUniformMatrix4fv(View_ID, 1, GL_FALSE, &mvp.getViewMatrix()[0][0]);
			glUniformMatrix4fv(Model_ID, 1, GL_FALSE, &mvp.getModelMatrix()[0][0]);

			// 1rst attribute buffer : vertices
			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, star_vertex_buffer);
			glVertexAttribPointer(
					0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
					3,                  // size
					GL_FLOAT,           // type
					GL_FALSE,           // normalized?
					0,                  // stride
					(void *) 0            // array buffer offset
			);

			// 2nd attribute buffer : UVs
			glEnableVertexAttribArray(1);
			glBindBuffer(GL_ARRAY_BUFFER, star_uv_buffer);
			glVertexAttribPointer(
					1,                                // attribute
					2,                                // size
					GL_FLOAT,                         // type
					GL_FALSE,                         // normalized?
					0,                                // stride
					(void *) 0                          // array buffer offset
			);

			// 3rd attribute buffer : normals
			glEnableVertexAttribArray(2);
			glBindBuffer(GL_ARRAY_BUFFER, star_normal_buffer);
			glVertexAttribPointer(
					2,                                // attribute
					3,                                // size
					GL_FLOAT,                         // type
					GL_FALSE,                         // normalized?
					0,                                // stride
					(void *) 0                          // array buffer offset
			);

			// Index buffer
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, star_index_buffer);
			// Draw the triangles!
			glDrawElements(
					GL_TRIANGLES,             // mode
					(GLsizei) star_indeces.size(), // count
					GL_UNSIGNED_SHORT,        // type
					(void *) 0                  // element array buffer offset
			);

			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
			glDisableVertexAttribArray(2);
		}

		{
			auto mvp = MVP(scale0, disk_position_x);

			// Send our transformation to the currently bound shader, in the "MVP" uniform
			// This is done in the main loop since each model will have a different MVP matrix (At least for the M part)
			glUniformMatrix4fv(MVP_ID, 1, GL_FALSE, &mvp.mvp()[0][0]);
			glUniformMatrix4fv(MV_ID, 1, GL_FALSE, &mvp.mv()[0][0]);
			glUniformMatrix4fv(View_ID, 1, GL_FALSE, &mvp.getViewMatrix()[0][0]);
			glUniformMatrix4fv(Model_ID, 1, GL_FALSE, &mvp.getModelMatrix()[0][0]);

//			// Bind our texture in Texture Unit 0
//			glActiveTexture(GL_TEXTURE0);
//			glBindTexture(GL_TEXTURE_2D, Texture);
//			// Set our "myTextureSampler" sampler to user Texture Unit 0
//			glUniform1i(TextureID, 0);

			// 1rst attribute buffer : vertices
			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, disk_vertex_buffer);
			glVertexAttribPointer(
					0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
					3,                  // size
					GL_FLOAT,           // type
					GL_FALSE,           // normalized?
					0,                  // stride
					(void *) 0            // array buffer offset
			);

			// 2nd attribute buffer : UVs
			glEnableVertexAttribArray(1);
			glBindBuffer(GL_ARRAY_BUFFER, disk_uv_buffer);
			glVertexAttribPointer(
					1,                                // attribute
					2,                                // size
					GL_FLOAT,                         // type
					GL_FALSE,                         // normalized?
					0,                                // stride
					(void *) 0                          // array buffer offset
			);

			// 3rd attribute buffer : normals
			glEnableVertexAttribArray(2);
			glBindBuffer(GL_ARRAY_BUFFER, disk_normal_buffer);
			glVertexAttribPointer(
					2,                                // attribute
					3,                                // size
					GL_FLOAT,                         // type
					GL_FALSE,                         // normalized?
					0,                                // stride
					(void *) 0                          // array buffer offset
			);

			// Index buffer
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, disk_index_buffer);
			// Draw the triangles!
			glDrawElements(
					GL_TRIANGLES,             // mode
					(GLsizei) disk_indeces.size(), // count
					GL_UNSIGNED_SHORT,        // type
					(void *) 0                  // element array buffer offset
			);

			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
			glDisableVertexAttribArray(2);
		}

        {
            const auto pixels_size = static_cast<size_t>(width * height * 3);
            unsigned short pixels[pixels_size];
            glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_UNSIGNED_SHORT, pixels);
            float flux_red = 0.0f;
            const float max = static_cast<float>(std::numeric_limits<unsigned short>::max()) + 1.0f;
            const auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for private(i) shared(pixels, max) reduction(+:flux_red)
            for( size_t i = 0; i < pixels_size; i += 3 ){
                flux_red += static_cast<float>(pixels[i+0]) / max;
            }
            const auto end = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " mcs\t"
                      << flux_red << std::endl;
        }

        if ( saveImage ){
            unsigned char pixels[width * height * 3];
            glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
            save_pixels_to_jpeg("image.jpeg", pixels, static_cast<size_t>(width), static_cast<size_t>(height));
        }

		{
			// Render to the screen
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			// Render on the whole framebuffer, complete from the lower left corner to the upper right
			glViewport(0, 0, width, height);
			// Clear the screen
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			// Use our shader
			glUseProgram(quad_programID);
			// Bind our texture in Texture Unit 0
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, renderedTexture);
			// Set our "renderedTexture" sampler to user Texture Unit 0
			glUniform1i(texID, 0);
			// 1rst attribute buffer : vertices
			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
			glVertexAttribPointer(
					0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
					3,                  // size
					GL_FLOAT,           // type
					GL_FALSE,           // normalized?
					0,                  // stride
					(void *) 0          // array buffer offset
			);
			glEnableVertexAttribArray(1);
			glBindBuffer(GL_ARRAY_BUFFER, quad_uvbuffer);
			glVertexAttribPointer(
					1,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
					2,                  // size
					GL_FLOAT,           // type
					GL_FALSE,           // normalized?
					0,                  // stride
					(void *) 0          // array buffer offset
			);
			// Draw the triangles !
			glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles

			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
		}

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();
	} // Check if the Space key was pressed or the window was closed
	while( glfwGetKey(window, GLFW_KEY_SPACE ) != GLFW_PRESS &&
		   glfwWindowShouldClose(window) == 0 );

	// Cleanup VBO and shader
	glDeleteBuffers(1, &star_vertex_buffer);
	glDeleteBuffers(1, &star_uv_buffer);
	glDeleteBuffers(1, &star_normal_buffer);
	glDeleteBuffers(1, &star_index_buffer);
	glDeleteBuffers(1, &disk_vertex_buffer);
	glDeleteBuffers(1, &disk_uv_buffer);
	glDeleteBuffers(1, &disk_normal_buffer);
	glDeleteBuffers(1, &disk_index_buffer);
	glDeleteProgram(programID);
	glDeleteProgram(quad_programID);
	glDeleteVertexArrays(1, &VertexArrayID);
	glDeleteFramebuffers(1, &FramebufferName);
	glDeleteTextures(1, &renderedTexture);
	glDeleteRenderbuffers(1, &depthrenderbuffer);
	glDeleteBuffers(1, &quad_vertexbuffer);
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}

