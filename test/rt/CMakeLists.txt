
set(rt_sources
    main.cc
    boundarycond.cc 
    boundarycond.h 
    compute.cc 
    compute.h 
    computePoisCoef.cc 
    computePoisCoef.h 
    init.cc 
    init.h 
    pars.cc 
    pars.h 
    phys.cc 
    phys.h 
    relax.cc 
    slvsml.cc 
    solve.h 
    readData.cc
    readData.h 
    resid.cc 
)
add_executable(rt ${rt_sources})
target_link_libraries(rt libtftt)
