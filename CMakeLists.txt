cmake_minimum_required(VERSION 3.6)
project(HerX1)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DENABLE_OPENMP)
endif()

find_package(Boost REQUIRED)
find_package(jpeg REQUIRED)
if (BOOST_FOUND)
    include_directories(${BOOST_INCLUDE_DIRS})
    link_libraries(${BOOST_LIBRARIES})
endif()

if (jpeg_FOUND)
    include_directories(${jpeg_INCLUDE_DIRS})
    link_libraries(${jpeg_LIBRARIES})
endif()

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})
link_libraries(${OPENGL_gl_LIBRARY})
#link_libraries(${OPENGL_glu_LIBRARY})
find_package(GLEW REQUIRED)
if (GLEW_FOUND)
    include_directories(${GLEW_INCLUDE_DIRS})
    link_libraries(${GLEW_LIBRARIES})
endif()
find_package(glfw3 REQUIRED)
if (glfw3_FOUND)
    include_directories(${glfw3_INCLUDE_DIRS})
    link_libraries(${glfw3_LIBRARIES})
endif()
find_package(glm REQUIRED)
if (glm_FOUND)
    include_directories(${glm_INCLUDE_DIRS})
    link_libraries(${glm_LIBRARIES})
endif()

set(CMAKE_CXX_STANDARD 11)

set(ALL_LIBS
        ${OPENGL_LIBRARY}
        glfw
        jpeg)
link_libraries(${ALL_LIBS})


set(SOURCE_FILES
        main.cpp
        shader.cpp shader.hpp
        DescriteToGeometric.hpp disks.hpp
        TriangleDiscreteCoordinates.hpp
        images.hpp
        commons.hpp
        Renderer.hpp
        LightSource.hpp
        ObjectModel.hpp
        MVP.hpp MvpFromInput.hpp TextureImage.hpp LimbDarking.hpp)
add_executable(main ${SOURCE_FILES})