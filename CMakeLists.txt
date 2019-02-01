Set(INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/AraDisplay ${CMAKE_SOURCE_DIR}/AraCorrelator ${CMAKE_SOURCE_DIR}/AraEvent ${LIBROOTFFTWWRAPPER_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/analysis/a23_analysis_tools )

Set(ANALYSIS_DIRECTORY ${CMAKE_SOURCE_DIR}/analysis/a23_analysis_programs)

include_directories( ${INCLUDE_DIRECTORIES} )
link_directories( ${CMAKE_SOURCE_DIR}/AraSim_1.0 ${CMAKE_SOURCE_DIR}/RayTraceCorrelator )
add_library(filterEvent a23_analysis_tools/tools_filterEvent.cxx)

add_executable(CompareFilters ${ANALYSIS_DIRECTORY}/daily/2018.12.13-Compare-Filters.cxx)
target_link_libraries(CompareFilters AraEvent filterEvent AraSim ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(Plot_CompareFilters ${ANALYSIS_DIRECTORY}/daily/2018.12.13-Compare-Filters-Plot.cxx)
target_link_libraries(Plot_CompareFilters AraEvent filterEvent AraSim ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

install(TARGETS CompareFilters Plot_CompareFilters DESTINATION ${ARAROOT_INSTALL_PATH}/bin)