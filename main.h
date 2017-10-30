#include <divsufsort64.h> // include header for suffix sort
#include <sdsl/rmq_support.hpp> //include header for range minimum queries
#include <sdsl/bit_vectors.hpp> // include header for bit vectors

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

using namespace sdsl;
using namespace std;
