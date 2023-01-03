cmake_minimum_required(VERSION 3.17)

set(PROJECT "GradientExplorer")

PROJECT(${PROJECT})

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOMOC ON)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /DWIN32 /EHsc /MP")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
endif(MSVC)

file(TO_CMAKE_PATH $ENV{HDPS_INSTALL_DIR} INSTALL_DIR)

# -----------------------------------------------------------------------------
# Install libraries
# -----------------------------------------------------------------------------
# Set the directory where CMake include files are to be found
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(flann_VERSION 1.9.1)
set(lz4_VERSION 1.9.3)

set(LIBRARY_INSTALL_DIR ${PROJECT_BINARY_DIR})
set(ARTIFACTORY_LIBS_INSTALLED TRUE)

include(InstallArtifactoryPackage)
message(STATUS "Beep: ${ARTIFACTORY_LIBS_INSTALLED}")
if (NOT ARTIFACTORY_LIBS_INSTALLED) 
    message(STATUS "Installing artifactory packages to: ${LIBRARY_INSTALL_DIR}")
    # Both HDILib and flann are available prebuilt in the lkeb-artifactory as combined Debug/Release packages
    # lz4 is also available in the lkb-artifactory in separate Debug and |Release packages
    install_artifactory_package(flann ${flann_VERSION} lkeb TRUE)
    install_artifactory_package(lz4 ${lz4_VERSION} _ FALSE)


    message(STATUS "module path ${CMAKE_MODULE_PATH}")

    set(ARTIFACTORY_LIBS_INSTALLED TRUE CACHE BOOL "Use the prebuilt libraries from artifactory" FORCE)
endif()

# Include flann includes via target_include_directories to the current project
macro(set_flann_project_includes target)
    if(1)
        MESSAGE( STATUS "Including PREBUILT flann ${LIBRARY_INSTALL_DIR}/flann")
        target_include_directories("${target}" PRIVATE "${LIBRARY_INSTALL_DIR}/flann/include")
    else()
        if(MSVC)
            target_include_directories ("${target}" PRIVATE "${LIBRARY_INSTALL_DIR}/tSNE/lib/Flann/Win/include")
        endif()

        if(APPLE)
            target_include_directories ("${target}" PRIVATE "${LIBRARY_INSTALL_DIR}/tSNE/lib/Flann/OSX/include")
        endif()

        if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
            target_include_directories ("${target}" PRIVATE "${LIBRARY_INSTALL_DIR}/tSNE/lib/Flann/Linux/include")
        endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    endif()
endmacro()

# Link the flann libraries via target_link_library to the current project
macro(set_flann_project_link_libraries target)
    if(1)
        MESSAGE( STATUS "Linking PREBUILT flann libraries...")
        if(MSVC)
            target_link_libraries("${target}" debug "${LIBRARY_INSTALL_DIR}/flann/lib/flann_cpp_s.lib")
            target_link_libraries("${target}" optimized "${LIBRARY_INSTALL_DIR}/flann/lib/flann_cpp_s.lib")
        endif()
    else()
        if(MSVC)
            MESSAGE( STATUS "Linking Windows flann libraries...")
            target_link_libraries("${target}" debug "${LIBRARY_INSTALL_DIR}/tSNE/lib/Flann/Win/Debug/flann_cpp_s.lib")
            target_link_libraries("${target}" optimized "${LIBRARY_INSTALL_DIR}/tSNE/lib/Flann/Win/Release/flann_cpp_s.lib")
        endif()

        if(APPLE)
            MESSAGE( STATUS "Linking Mac OS X libraries...")
            target_link_libraries("${target}" "${LIBRARY_INSTALL_DIR}/tSNE/lib/Flann/OSX/libflann_s.a")
        endif()
    endif()
endmacro()

# lz4 currently only with prebuilt libs - is used in flann 1.8.4 and greater.
macro(set_lz4_project_includes target)
    if(1)
        MESSAGE( STATUS "Including PREBUILT lz4 ${LIBRARY_INSTALL_DIR}/lz4")
        target_include_directories("${target}" PRIVATE "${LIBRARY_INSTALL_DIR}/lz4/Release/include")
    else()
    endif()
endmacro()

macro(set_lz4_project_link_libraries target)
    if(1)
        MESSAGE( STATUS "Linking PREBUILT lz4 libraries...")
        if(MSVC)
            target_link_libraries("${target}" debug "${LIBRARY_INSTALL_DIR}/lz4/Release/lib/lz4.lib")
            target_link_libraries("${target}" optimized "${LIBRARY_INSTALL_DIR}/lz4/Release/lib/lz4.lib")
        endif()
        if(APPLE)
            target_link_libraries("${target}" debug "${LIBRARY_INSTALL_DIR}/lz4/Debug/lib/liblz4.a")
            target_link_libraries("${target}" optimized "${LIBRARY_INSTALL_DIR}/lz4/Release/lib/liblz4.a")
        endif()
    else()
    endif()
endmacro()

# -----------------------------------------------------------------------------
# Rest
# -----------------------------------------------------------------------------
find_package(Qt6 COMPONENTS Widgets WebEngineWidgets OpenGL OpenGLWidgets Charts REQUIRED)

set(PLUGIN
    src/Common.h
    src/ScatterplotPlugin.h
    src/ScatterplotPlugin.cpp
)

set(UI
    src/ScatterplotWidget.h
    src/ScatterplotWidget.cpp
    src/ProjectionView.h
    src/ProjectionView.cpp
    src/GradientGraph.h
    src/GradientGraph.cpp
    src/ExportImageDialog.h
    src/ExportImageDialog.cpp
)

set(Actions
    src/Actions/ColoringAction.h
    src/Actions/ColoringAction.cpp
    src/Actions/DensityPlotAction.h
    src/Actions/DensityPlotAction.cpp
    src/Actions/ManualClusteringAction.h
    src/Actions/ManualClusteringAction.cpp
    src/Actions/MiscellaneousAction.h
    src/Actions/MiscellaneousAction.cpp
    src/Actions/PlotAction.h
    src/Actions/PlotAction.cpp
    src/Actions/PluginAction.h
    src/Actions/PluginAction.cpp
    src/Actions/PointPlotAction.h
    src/Actions/PointPlotAction.cpp
    src/Actions/PositionAction.h
    src/Actions/PositionAction.cpp
    src/Actions/RenderModeAction.h
    src/Actions/RenderModeAction.cpp
    src/Actions/ScalarAction.h
    src/Actions/ScalarAction.cpp
    src/Actions/ScalarSourceAction.h
    src/Actions/ScalarSourceAction.cpp
    #src/Actions/SelectionAction.h
    #src/Actions/SelectionAction.cpp
    src/Actions/SettingsAction.h
    src/Actions/SettingsAction.cpp
    src/Actions/SubsetAction.h
    src/Actions/SubsetAction.cpp
    src/Actions/ExportImageAction.h
    src/Actions/ExportImageAction.cpp
)

set(Models
    src/Models/ColorSourceModel.h
    src/Models/ColorSourceModel.cpp
    src/Models/ScalarSourceModel.h
    src/Models/ScalarSourceModel.cpp
)

set(Compute
    src/Compute/LocalDimensionality.h
    src/Compute/LocalDimensionality.cpp
    src/Compute/RandomWalks.h
    src/Compute/RandomWalks.cpp
    src/Compute/KnnGraph.h
    src/Compute/KnnGraph.cpp
    src/Compute/Filters.h
    src/Compute/Filters.cpp
)

set(SHADERS
    res/shaders/SelectionTool.frag
    res/shaders/SelectionTool.vert
)

set(AUX
    res/Resources.qrc
    src/GradientExplorerPlugin.json
)

set(SOURCES ${PLUGIN} ${UI} ${Actions} ${Models} ${Compute})

source_group(Plugin FILES ${PLUGIN})
source_group(UI FILES ${UI})
source_group(Actions FILES ${Actions})
source_group(Models FILES ${Models})
source_group(Compute FILES ${Compute})
source_group(Shaders FILES ${SHADERS})
source_group(Aux FILES ${AUX})

include_directories("${INSTALL_DIR}/$<CONFIGURATION>/include/")
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/Eigen/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/Actions)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/Models)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/Compute)

add_library(${PROJECT} SHARED ${SOURCES} ${SHADERS} ${AUX})

set_lz4_project_includes(${PROJECT})
set_flann_project_includes(${PROJECT})

# Set CXX_STANDARD to avoid osx 10.14.4/xcode10.3 build error "no member named 'make_unique'..."
# (C++14 should be sufficient to avoid the error, but C++17 is even better!)
set_property(TARGET ${PROJECT} PROPERTY CXX_STANDARD 17)

if(MSVC)
    set(LIB_LINK_SUFFIX ".lib") 
else()
    set(LIB_LINK_SUFFIX "${CMAKE_SHARED_LIBRARY_SUFFIX}")
endif(MSVC)
set(HDPS_LINK_LIBRARY "${INSTALL_DIR}/$<CONFIGURATION>/lib/${CMAKE_SHARED_LIBRARY_PREFIX}HDPS_Public${LIB_LINK_SUFFIX}")
set(POINTDATA_LINK_LIBRARY "${INSTALL_DIR}/$<CONFIGURATION>/lib/${CMAKE_SHARED_LIBRARY_PREFIX}PointData${LIB_LINK_SUFFIX}") 
set(CLUSTERDATA_LINK_LIBRARY "${INSTALL_DIR}/$<CONFIGURATION>/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ClusterData${LIB_LINK_SUFFIX}") 

target_link_libraries(${PROJECT} Qt6::Widgets)
target_link_libraries(${PROJECT} Qt6::WebEngineWidgets)
target_link_libraries(${PROJECT} Qt6::OpenGL)
target_link_libraries(${PROJECT} Qt6::OpenGLWidgets)
target_link_libraries(${PROJECT} Qt6::Charts)
target_link_libraries(${PROJECT} "${HDPS_LINK_LIBRARY}")
target_link_libraries(${PROJECT} "${POINTDATA_LINK_LIBRARY}")
target_link_libraries(${PROJECT} "${CLUSTERDATA_LINK_LIBRARY}")
set_lz4_project_link_libraries(${PROJECT})
set_flann_project_link_libraries(${PROJECT})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT} OpenMP::OpenMP_CXX)
endif()

install(TARGETS ${PROJECT}
   RUNTIME DESTINATION Plugins COMPONENT SHAREDLIB
   LIBRARY DESTINATION Plugins COMPONENT LINKLIB
   ARCHIVE DESTINATION lib COMPONENT LINKLIB
)

add_custom_command(TARGET ${PROJECT} POST_BUILD
    COMMAND "${CMAKE_COMMAND}"
        --install ${CMAKE_CURRENT_BINARY_DIR}
        --config $<CONFIGURATION>
        --prefix ${INSTALL_DIR}/$<CONFIGURATION>
)

# Automatically set the debug environment (command + working directory) for MSVC in debug mode
set_property(TARGET ${PROJECT} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY $<IF:$<CONFIG:DEBUG>,${INSTALL_DIR}/debug,${INSTALL_DIR}/release>)
set_property(TARGET ${PROJECT} PROPERTY VS_DEBUGGER_COMMAND $<IF:$<CONFIG:DEBUG>,${INSTALL_DIR}/debug/HDPS.exe,${INSTALL_DIR}/release/HDPS.exe>)