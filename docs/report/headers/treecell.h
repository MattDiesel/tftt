struct TreeCell {
    // For Calculating the offset to the group
    int8_t index;

    data_t data;
    TreeGroup* children;
    node_t rank;

    int8_t poisNgbC;
    CellRef poisNgb[8];
    double poisCoef[8];
    double cenCoef;

    ADAPTFLAGS adaptFlags;
};
