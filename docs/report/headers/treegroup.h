struct TreeGroup {
    ident_t id;
    double origin[DIM];
    std::array<TreeCell, 1<<DIM> cells;

    // FTT
    std::array<CellRef,DIM*2> neighbours;
    CellRef parent;

    // For SFC
    int orientation;
    CellRef next;
    CellRef prev;
};
