
#ifndef TFTT_IO_SERIALIZE_H
#define TFTT_IO_SERIALIZE_H


namespace tftt {


void saveTree(std::string fname);
void saveTree(std::ostream& os);
void splitToDisk(std::string fnameFmt);


void loadTree(std::string fname, int n = -1);


} // namespace tftt


#endif
