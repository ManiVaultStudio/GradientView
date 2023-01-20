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

find_package(Qt6 COMPONENTS Widgets WebEngineWidgets OpenGL OpenGLWidgets Charts REQUIRED)

set(PLUGIN
    src/Common.h
    src/ScatterplotPlugin.h
    src/ScatterplotPlugin.cpp
    src/CudaCheck.h
    src/Timer.h
    src/Timer.cpp
)

set(UI
    src/ScatterplotWidget.h
    src/ScatterplotWidget.cpp
    src/CellRenderer.h
    src/CellRenderer.cpp
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
    src/Compute/KnnIndex.h
    src/Compute/KnnIndex.cpp
    src/Compute/KnnGraph.h
    src/Compute/KnnGraph.cpp
    src/Compute/Filters.h
    src/Compute/Filters.cpp
    src/Compute/DataTransformations.h
    src/Compute/DataTransformations.cpp
    src/Compute/SecondaryDistanceMeasures.h
    src/Compute/SecondaryDistanceMeasures.cpp
)

set(SHADERS
    res/shaders/SelectionTool.frag
    res/shaders/SelectionTool.vert
    res/shaders/TrianglePlot.frag
    res/shaders/TrianglePlot.vert
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
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/faiss/include)
include_directories("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.0/include")
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/jcv/include)
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/annoy/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/Actions)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/Models)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/Compute)

add_library(${PROJECT} SHARED ${SOURCES} ${SHADERS} ${AUX})

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

find_package(CUDA 12.0 REQUIRED)

message(STATUS "${MKL_LIBRARIES}")
target_link_libraries(${PROJECT} Qt6::Widgets)
target_link_libraries(${PROJECT} Qt6::WebEngineWidgets)
target_link_libraries(${PROJECT} Qt6::OpenGL)
target_link_libraries(${PROJECT} Qt6::OpenGLWidgets)
target_link_libraries(${PROJECT} Qt6::Charts)
target_link_libraries(${PROJECT} ${PROJECT_SOURCE_DIR}/thirdparty/faiss/lib/$<CONFIG>/faiss.lib)
target_link_libraries(${PROJECT} ${CUDA_LIBRARIES})
target_link_libraries(${PROJECT} ${CUDA_CUBLAS_LIBRARIES})
target_link_libraries(${PROJECT} "${HDPS_LINK_LIBRARY}")
target_link_libraries(${PROJECT} "${POINTDATA_LINK_LIBRARY}")
target_link_libraries(${PROJECT} "${CLUSTERDATA_LINK_LIBRARY}")

# find_package(MKL CONFIG REQUIRED)
# if(MKL_FOUND)
  # message(STATUS "${MKL_LIBRARIES}")
  # target_link_libraries(${PROJECT} MKL::MKL)
# endif()

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