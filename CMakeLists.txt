
PROJECT(SPK)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
#target_link_libraries(SPK ${PYTHON_LIBRARIES})

add_subdirectory(goal_implementations)
add_executable (spk spkmeans.c spkmeans.h
        utils.c utils.h
        kmeans.c kmeans.h
        goal_implementations/wam.c goal_implementations/wam.h
        goal_implementations/ddg.c goal_implementations/ddg.h
        goal_implementations/lnorm.c goal_implementations/lnorm.h
        goal_implementations/spk.c goal_implementations/spk.h
        )
