Set(INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/AraDisplay ${CMAKE_SOURCE_DIR}/AraCorrelator ${CMAKE_SOURCE_DIR}/AraEvent ${LIBROOTFFTWWRAPPER_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/analysis/a23_analysis_tools )

Set(ANALYSIS_DIRECTORY ${CMAKE_SOURCE_DIR}/analysis/a23_analysis_programs)

include_directories( ${INCLUDE_DIRECTORIES} )
link_directories( ${CMAKE_SOURCE_DIR}/AraSim_1.0 ${CMAKE_SOURCE_DIR}/RayTraceCorrelator )
add_library(filterEvent a23_analysis_tools/tools_filterEvent.cxx)

#primary analysis code

add_executable(v2_Baselines ${ANALYSIS_DIRECTORY}/diffuse/v2_analysis_baselines.cxx)
target_link_libraries(v2_Baselines AraEvent AraSim ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(v2_CWID ${ANALYSIS_DIRECTORY}/diffuse/v2_analysis_CWID.cxx)
target_link_libraries(v2_CWID AraEvent AraSim ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(v2_Filter ${ANALYSIS_DIRECTORY}/diffuse/v2_analysis_filter.cxx)
target_link_libraries(v2_Filter AraEvent filterEvent AraSim ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(v2_FinalPlot ${ANALYSIS_DIRECTORY}/diffuse/v2_final_plots.cxx)
target_link_libraries(v2_FinalPlot ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(v2_Join ${ANALYSIS_DIRECTORY}/diffuse/v2_analysis_join_filter_41_300.cxx)
target_link_libraries(v2_Join AraEvent filterEvent AraSim ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(v2_Reco ${ANALYSIS_DIRECTORY}/diffuse/v2_analysis_reco.cxx)
target_link_libraries(v2_Reco AraEvent filterEvent AraSim RayTraceCorrelator ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(v2_RunSummary ${ANALYSIS_DIRECTORY}/diffuse/v2_analysis_run_summary.cxx)
target_link_libraries(v2_RunSummary AraEvent filterEvent AraSim ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(v2_SaveVals ${ANALYSIS_DIRECTORY}/diffuse/v2_save_vals.cxx)
target_link_libraries(v2_SaveVals AraEvent filterEvent AraSim RayTraceCorrelator ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

#other debug stuff

add_executable(CompareFilters ${ANALYSIS_DIRECTORY}/daily/2018.12.13-Compare-Filters.cxx)
target_link_libraries(CompareFilters AraEvent filterEvent AraSim ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(Plot_CompareFilters ${ANALYSIS_DIRECTORY}/daily/2018.12.13-Compare-Filters-Plot.cxx)
target_link_libraries(Plot_CompareFilters AraEvent filterEvent AraSim ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(Debug_PrintEvent ${ANALYSIS_DIRECTORY}/diffuse/debug_print_event.cxx)
target_link_libraries(Debug_PrintEvent AraEvent AraSim RayTraceCorrelator ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES})


install(TARGETS v2_Baselines v2_CWID v2_Filter v2_FinalPlot v2_Join v2_Reco v2_RunSummary v2_SaveVals CompareFilters Plot_CompareFilters Debug_PrintEvent DESTINATION ${ARAROOT_INSTALL_PATH}/bin)