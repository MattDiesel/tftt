
set(tfttstruct_SRCS
    treegroup.cc
    treecell.cc
    tree.cc
    treevertex.cc 
    treeface.cc 
)

add_library(tfttstruct STATIC
    ${tfttstruct_SRCS}
)
