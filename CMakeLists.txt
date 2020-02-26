project(TH_reconstruction)

cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)


add_ngsolve_python_module(reclib_TH
	reconstruction_utils.cpp
	reconstruction_vertex.cpp
	reconstruction_element.cpp
	reconstruction_vertex_3D.cpp
  )

# check if CMAKE_INSTALL_PREFIX is set by user, otherwise install to user-local python modules directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import site; print(site.USER_SITE)" OUTPUT_VARIABLE install_dir OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_INSTALL_PREFIX ${install_dir} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

message("The python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS reclib_TH DESTINATION .)
