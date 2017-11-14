#include <divsufsort64.h> // include header for suffix sort
#include <sdsl/rmq_support.hpp> //include header for range minimum queries
#include <sdsl/bit_vectors.hpp> // include header for bit vectors

#define ALLOC_SIZE              1048576
#define DEL                     '$'
#define DEL_STR                 "$"

#define DNA                     "ACGTN"                         //DNA alphabet
#define PROT                    "ARNDCQEGHILKMFPSTWYVOUBZJX*"          //Proteins alphabet
#define IUPAC                   "ACGTUWSMKRYBDHVN"          	//IUPAC alphabet
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

using namespace sdsl;
using namespace std;


struct TSwitch
 {

   char               * input_filename;
   char               * output_filename;
   unsigned int         k,m,g,x;
   unsigned int         total_length;
 };


double gettime( void );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
void usage ( void );





