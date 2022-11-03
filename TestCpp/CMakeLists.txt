cmake_minimum_required(VERSION 3.15)

# 项目信息
project (TestVtkCpp
         VERSION 1.0.0.0
         LANGUAGES CXX
)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})

vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
    )