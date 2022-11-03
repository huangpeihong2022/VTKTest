cmake_minimum_required (VERSION 3.15)

# 项目信息
project (TestCyLinder
         VERSION 1.0.0.0
         LANGUAGES CXX
)

## 源文件列表
file(GLOB_RECURSE QML_FILES ${CMAKE_CURRENT_SOURCE_DIR}/qml/*.qml)
foreach(QML_FILE ${QML_FILES})
   file(RELATIVE_PATH QML_RELATIVE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${QML_FILE})
   list(APPEND QML_RELATIVE_PATHS ${QML_RELATIVE_PATH})
endforeach()

# 启动编译，链接依赖
add_executable(${PROJECT_NAME} main.cpp qml.qrc)
target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES}
                                                Qt6::Core Qt6::Gui Qt6::Qml Qt6::Quick)

add_compile_definitions(QuickVtkImportPath="D:/VTK_/master/debug/lib/qml/Debug")

qt6_add_qml_module(${PROJECT_NAME}
    URI "${PROJECT_NAME}App"
    VERSION 1.0
    QML_FILES ${QML_RELATIVE_PATHS}
    NO_RESOURCE_TARGET_PATH
    NO_GENERATE_PLUGIN_SOURCE
)

vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
)
