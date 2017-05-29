
#ifndef TFTT_IO_SERIALIZE_H
#define TFTT_IO_SERIALIZE_H


#include <string>
#include <ostream>


namespace tftt {


struct TrHeader {
    static constexpr uint16_t Magic = (uint16_t('L') << 8) | uint16_t('T');
    static constexpr int VersionMajor = 4;
    static constexpr int VersionMinor = 0;

    // File metadata
    uint16_t magic;
    uint8_t versionMajor;
    uint8_t versionMinor;

    // Config Requirements
    uint8_t dimensions;
    uint8_t identSize;
    uint16_t dataSize;

    // Problem metadata
    double domainOrigin[DIM];
    double domainSize[DIM];

    // World Info
    uint8_t worldSize;
    uint8_t worldRank;

    TrHeader();
    TrHeader(int rank, int world);

    void check();
};


void saveTree(std::string fname);
void saveTree(std::ostream& os);

// Todo: implement
// void loadTree(std::string fname);
// void loadTree(std::ostream& os);

void saveParTree(std::string fnameFmt, int world);
void loadParTree(std::string fname);


} // namespace tftt


#endif
