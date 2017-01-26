
#include "gray.h"


namespace tftt::utils {


int toGray(int n) {
    return (n ^ (n >> 1));
}
    

} // namespace tftt::utils
