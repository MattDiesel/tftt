
set(tfttstruct_SRCS
    treegroup.h 
    treegroup.cc
    treecell.h 
    treecell.cc
    tree.h 
    tree.cc
    treevertex.cc 
    treevertex.h 
    treeface.h 
    treeface.cc 
)

add_library(tfttstruct STATIC
    ${tfttstruct_SRCS}
)
