

add_executable(movecells movecells.cc)
target_link_libraries(movecells libtftt util)

add_executable(demoline demoline.cc)
target_link_libraries(demoline libtftt util)

add_executable(distrib distrib.cc)
target_link_libraries(distrib libtftt util)

add_executable(borderanalysis borderanalysis.cc)
target_link_libraries(borderanalysis libtftt util)
