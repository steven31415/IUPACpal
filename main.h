#include <divsufsort64.h> // Header for suffix sort
#include <sdsl/rmq_support.hpp> // Header for Range Minimum Queries
#include <sdsl/bit_vectors.hpp> // Header for bit vectors

// Enforce use of 64-bit integers if intructed by compiler to do so
#ifdef _USE_64
typedef int64_t INT;
#endif

// Enforce use of 32-bit integers if intructed by compiler to do so
#ifdef _USE_32
typedef int32_t INT;
#endif

void rmq_preprocess(INT *m, INT *v, INT n);

static __inline INT flog2(INT v) {
  return (INT) floor (log2((double)(v)));
}

using namespace sdsl;
using namespace std;
