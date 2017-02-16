
#include "gray.h"


namespace tftt {
namespace utils {


int toGray(int n) {
    return (n ^ (n >> 1));
}
    

}} // namespace tftt::utils
