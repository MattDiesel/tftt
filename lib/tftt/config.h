
#ifndef TFTT_CONFIG_H
#define TFTT_CONFIG_H


#include <cstdint>


#ifndef DIM
    #define DIM 2
#endif

// #define TFTT_DEBUG
// #define TFTT_FACES
// #define TFTT_VERTICES


struct rt_data {
    double cc;
    // double vof1;
    // double vof2;
    // double vof3;
    double P;
    // double Pstar;
    double rhs;
    double res;

    // double vof;
    // double rho;

    // double dive;
};

struct rt_facedata {
    double poisCoef;
};

struct rt_vertexdata {
    int placeholder;
    // double V[DIM];
    // double v[DIM];
    // double F1[DIM];
    // double F2[DIM];
    // double D1[DIM];
    // double D2[DIM];
};



namespace tftt {


template<typename T> struct TreeId;

typedef rt_data data_t;
typedef rt_facedata facedata_t;
typedef rt_vertexdata vertexdata_t;
typedef TreeId<uint64_t> ident_t;
typedef int8_t node_t;

// Functions to be supplied by user
typedef double (*fnDataNorm)(data_t& d, int max);
typedef double (*fnData)(data_t& d);
typedef double& (*fnDataRef)(data_t& d);

extern struct TFTTOPTIONS {
    int two2oneFlag; // 0 - None, 1 - strict, 2 - incl. corners, 3 - 3-2-1
} options;

struct CellRef;
typedef CellRef cell_t;

#ifdef TFTT_VERTICES
    struct VertexRef;
    typedef VertexRef vertex_t;
#endif

#ifdef TFTT_FACES
    struct FaceRef;
    typedef FaceRef face_t;
#endif


typedef bool (*fnCheckCell)(cell_t& cl);
typedef double (*fnCell)(cell_t& cl);


} // namespace tftt


#endif
