struct Tree {
    TreeGroup* root;
    size_t ccells;

    double size[DIM];

    // For thread
    CellRef first;
    CellRef last;

    std::vector<std::set<CellRef, CellRef::parless>> ghosts;
    std::vector<std::set<CellRef, CellRef::parless>> borders;

    std::vector<data_t*> ghostData;
    std::vector<data_t*> borderData;
};
