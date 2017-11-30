#include <iostream>
#include <string>
#include <map>
#include <set>
#include <list>
#include <tuple>
#include <sys/stat.h>
#include <limits.h>
#include "main.h" // include header for main

typedef map< char, set<char> > I_map;

////////////////////////
//  HELPER FUNCTIONS  //
////////////////////////

// check file exists
int exist(const char *name) {
  struct stat   buffer;
  return (stat (name, &buffer) == 0);
}

// print usage of IUPACpal
void usage() {
	fprintf ( stdout, "\n" );
	fprintf ( stdout, "      PARAMETER       TYPE      DEFAULT         DESCRIPTION\n" );
	fprintf ( stdout, "  -f, input_file      <str>     -               Input filename of plaintext DNA data (IUPAC).\n" );
	fprintf ( stdout, "  -m, min_len         <int>     10              Minimum length.\n");
	fprintf ( stdout, "  -M, max_len         <int>     100             Maximum length.\n");
	fprintf ( stdout, "  -g, max_gap         <int>     100             Maximum permissible gap.\n");
	fprintf ( stdout, "  -x, mismatches      <int>     0               Maximum permissible mismatches.\n");
	fprintf ( stdout, "  -o, output_file     <str>     IUPACpal.out    Output filename.\n" );
	fprintf ( stdout, "\n" );
	}

// nlogn RMQ preprocssing
void rmq_preprocess( INT * m, INT * v, INT n ) 
{
  INT i, j;
  INT lgn = flog2(n); 
  for (i = 0; i < n; i++)
    m[i*lgn] = i;
  for (j = 1; ((( INT ) 1) << j) <= n; j++)
    for (i = 0; i + (1 << j) - 1 < n; i++)
      if (v[m[i*lgn + j - 1]] < v[m[(i + (1 << (j - 1)))*lgn + j - 1]])
        m[i*lgn + j] = m[i*lgn + j - 1];
      else
        m[i*lgn + j] = m[(i + (1 << (j - 1)))*lgn + j - 1];
}

// nlogn RMQ query
static __inline INT rmq(INT *m, INT *v, INT n, INT i, INT j) { //query call
  INT lgn = flog2(n); 
  if (i > j) {INT tmp = j; j = i; i = tmp;}
  i++; //for LCP
  if (i == j) return i;
  INT k = flog2(j-i+1); 
  INT a = m[i * lgn + k]; 
  INT shift = ((( INT ) 1) << k);
  INT b = m[(j - shift + 1) * lgn + k];  
  return v[a]>v[b]?b:a;
}

// Convert a string to any type, if possible
template <typename T>
T ConvertString( const string &data, int* return_code) {
  if ( !data.empty())
  {
    T ret;
    istringstream iss( data );
    if ( data.find( "0x" ) != std::string::npos )
    {
      iss >> std::hex >> ret;
    }
    else
    {
      iss >> std::dec >> ret;
    }

    if ( iss.fail( ))
    {
      *return_code = -1; // Error
      return T ( ); // Empty Type
    }
    *return_code = 0; // Valid input
    return ret; // Value of converted input
  }
  *return_code = 1; // Empty input
  return T( ); // Empty Type
}

// Get character length of integer
int getDigitCount(int x) {
    stringstream iss;
    iss << x;
    return iss.str().size();
}

// Print an array of any type
template<typename T>
void print_array(string title, T seq, long int n, bool print_indices=false) {
    if (print_indices) {
        for(int i = 0 ; i < title.size() + 2 ; i++) {
            cout << " ";
        }

        for(int i = 0 ; i < n ; i++) {
            stringstream data;
            data << seq[i];
            int data_size = data.str().size();

            stringstream index;
            index << i;
            int index_size = index.str().size();

            cout << i;

            while (data_size > index_size - 2) {
                cout << " ";
                data_size--;
            }
        }
    }

    cout << endl;

    cout << title << ": ";
    for(int i=0 ; i<n ; i++) {
        cout << seq[i] << "  ";
    }

    cout << endl << endl;
}

// Inserts IUPAC characters mapping information into an STSL map
void IUPAC_map_insert(I_map* IUPAC_map, int IUPAC_to_value[], char IUPAC_char, set<char> mapped_chars) {
    static int index = 0;
    IUPAC_map->insert( pair< char, set<char> > (IUPAC_char, mapped_chars) );
    IUPAC_to_value[IUPAC_char] = index++;
}

// Get a value from a matrix stored as a 1D contiguous sequence of memory
template<typename T>
T get_value_from_matrix(T* mtx, uint32_t ncols, uint32_t i, uint32_t j)
{
    return mtx[i + i * (ncols - 1) + j];
}

// Set a value in a matrix stored as a 1D contiguous sequence of memory
template<typename T>
void set_value_for_matrix(T* mtx_ptr, uint32_t ncols, uint32_t i, uint32_t j, T val) 
{
    mtx_ptr[i + i * (ncols - 1) + j] = val;
}

class MatchMatrix {
    public:
        static bool* match_matrix;
        static int IUPAC_map_count;
        static int* IUPAC_to_value;

    static bool match(char a, char b) {
        return get_value_from_matrix(match_matrix, IUPAC_map_count, IUPAC_to_value[a], IUPAC_to_value[b]);
    }
};

bool* MatchMatrix::match_matrix;
int MatchMatrix::IUPAC_map_count;
int* MatchMatrix::IUPAC_to_value;

////////////////////////
//  STRING FUNCTIONS  //
////////////////////////

// Calculates the Longest Common Prefix array of a text and stores in given variable LCP
// Requires Text and Text Length, Suffix Array and Longest Common Prefix
unsigned int LCParray(unsigned char *text, INT n, INT * SA, INT * invSA, INT * LCP)
{                                                                               
        INT i = 0, j = 0;

        LCP[0] = 0;
        for ( i = 0; i < n; i++ ) // compute LCP[invSA[i]]
                if ( invSA[i] != 0 ) 
                {
                        if ( i == 0) j = 0;
                        else j = (LCP[invSA[i-1]] >= 2) ? LCP[invSA[i-1]]-1 : 0;
                        while ( text[i+j] == text[SA[invSA[i]-1]+j] )
                                j++;
                        LCP[invSA[i]] = j;
                }

        return ( 1 );
}

// Returns the Longest Common Extension between position i and j (i, j order doesn't matter)
// Requires Text Length, Inverse Suffix Array,
// Longest Common Prefix Array (LCP) and RMQ function over LCP


#ifdef _USE_NLOGN_RMQ
unsigned int LCE(INT i, INT j, INT n, INT * invSA, INT * LCP, INT * A) {
#else
unsigned int LCE(INT i, INT j, INT n, INT * invSA, INT * LCP, rmq_succinct_sct<> rmq) {
#endif
    if (i == j) {
            return n - i;
    }

    INT a = invSA[i];
    INT b = invSA[j];
    INT c = 0;

    if (a > b) {
            c = a;
            a = b;
            b = c;
    }

    #ifdef _USE_NLOGN_RMQ
        return LCP[rmq(A, LCP, n, a, b)]; // rmq(a, b) does not include "a" value
    #else
        return LCP[rmq(a + 1, b)]; // rmq(a, b) does include "a" value
    #endif
}

// Calculates a list of longest common extensions, corresponding to 0, 1, 2, ... allowed mismatches
// Only considers "real" mismatches (degenerate string mismatching)
// Takes into account the matching possibility of non A, C, G, T/U characters
// Longest Common Extension calculated from position i and j (i, j order doesn't matter)
// Requires Text, Text Length, Inverse Suffix Array,
// Longest Common Prefix Array (LCP), RMQ function over LCP, Match Matrix, mismatches (max allowed) and initial_gap
// Only starts counting mismatches that occur after initial gap, however earlier mismatches are still stored
// Should only be used after MatchMatrix has been instantiated with necessary data
#ifdef _USE_NLOGN_RMQ
void realLCE_mismatches(unsigned char* text, INT i, INT j, INT n, INT * invSA, INT * LCP, INT * A, int mismatches, int initial_gap, list<double>* mismatch_locs) {
#else
void realLCE_mismatches(unsigned char* text, INT i, INT j, INT n, INT * invSA, INT * LCP, rmq_succinct_sct<> rmq, int mismatches, int initial_gap, list<double>* mismatch_locs) {
#endif
    if ( i == j ) {
        mismatch_locs->push_back( n - i );
    }
    else {
        int real_lce = 0;

        while (mismatches >= 0) {
            #ifdef _USE_NLOGN_RMQ
            real_lce = real_lce + LCE(i + real_lce, j + real_lce, n, invSA, LCP, A);
            #else
            real_lce = real_lce + LCE(i + real_lce, j + real_lce, n, invSA, LCP, rmq);
            #endif

            if ( i + real_lce >= (n / 2) or j + real_lce >= n ) {
                break;
            }

            char s1 = text[i + real_lce];
            char s2 = text[j + real_lce];

            if ( !MatchMatrix::match(s1, s2) ) {

                mismatch_locs->push_back( real_lce );
                if (real_lce >= initial_gap) {
                    mismatches--;
                }
            }

            real_lce++;
        }
    }
}

#ifdef _USE_NLOGN_RMQ
void addPalindromes(set<tuple<int, int, int>>* palindromes, unsigned char* S, int S_n, int n, INT * invSA, INT * LCP, INT * A, tuple<int, int, int, int> params) {
#else
void addPalindromes(set<tuple<int, int, int>>* palindromes, unsigned char* S, int S_n, int n, INT * invSA, INT * LCP, rmq_succinct_sct<> rmq, tuple<int, int, int, int> params) {
#endif
    int min_len = get<0>(params);
    int max_len = get<1>(params);
    int mismatches = get<2>(params);
    int max_gap = get<3>(params);

    for (double c = 0; c <= (n - 1); c += 0.5 ) {

        int i, j;

        // Determine if centre corresponds to odd or even palindrome
        bool isOdd = (trunc(c) == c);

        // Strategically choose i and j to determine maximum extension
        if (isOdd) {
            i = int(c + 1.0);
            j = int(2.0 * n + 1.0 - c);
        } else {
            i = int(c + 0.5);
            j = int(2.0 * n + 1.0 - (c + 0.5));
        }

        int initial_gap;

        if (max_gap % 2 == 1) {
            initial_gap = (max_gap - 1) / 2;
        }
        else {
            if (isOdd) {
                initial_gap = (max_gap - 2) / 2;
            } else {
                initial_gap = max_gap / 2;
            }
        }

        list<double> mismatch_locs;

        #ifdef _USE_NLOGN_RMQ
            realLCE_mismatches(S, i, j, S_n, invSA, LCP, A, mismatches, initial_gap, &mismatch_locs);
        #else
            realLCE_mismatches(S, i, j, S_n, invSA, LCP, rmq, mismatches, initial_gap, &mismatch_locs);
        #endif

        mismatch_locs.push_front(-1.0);

        list<pair<double, int>> valid_start_locs; // (mismatch_location, mismatch_index)
        list<pair<double, int>> valid_end_locs; // (mismatch_location, mismatch_index)

        // Determine list of valid start and end locations
        int mismatch_index = 0;
        for (list<double>::iterator it = mismatch_locs.begin(); it != mismatch_locs.end(); ++it){
            if (next(it) != mismatch_locs.end() and *next(it) != *it + 1.0) {
                valid_start_locs.push_back(pair<double, int> (*it, mismatch_index));
            }

            if (it != mismatch_locs.begin() and *prev(it) != *it - 1.0) {
                valid_end_locs.push_back(pair<double, int> (*it, mismatch_index));
            }

            mismatch_index++;
        }

        // Optional printing of mismatch locations relative to centre
        if (false) {
            cout << "c=" << c << "\t";
            cout << "[ ";
            for (list<double>::iterator it = mismatch_locs.begin(); it != mismatch_locs.end(); ++it){
                cout << *it << " ";
            }
            cout << "]" << endl;
        }

        // Optional printing of valid start and end locations
        if (false) {
            cout << "start" << "\t";
            cout << "[ ";
            for (list<pair<double, int>>::iterator it = valid_start_locs.begin(); it != valid_start_locs.end(); ++it){
                cout << "(" << it->first << ", " << it->second << ") ";
            }
            cout << "]" << endl;

            cout << "end" << "\t";
            cout << "[ ";
            for (list<pair<double, int>>::iterator it = valid_end_locs.begin(); it != valid_end_locs.end(); ++it){
                cout << "(" << it->first << ", " << it->second << ") ";
            }
            cout << "]" << endl;

            cout << endl;
        }

        if ( !valid_start_locs.empty() and !valid_end_locs.empty() ) {
            int mismatch_diff;
            list<pair<double, int>>::iterator start_it = valid_start_locs.begin();
            list<pair<double, int>>::iterator end_it = valid_end_locs.begin();
            int left, right, gap;
        	double start_mismatch, end_mismatch;

            while( start_it != valid_start_locs.end() and end_it != valid_end_locs.end()) {
                mismatch_diff = end_it->second - start_it->second - 1;

                while (mismatch_diff > mismatches) {
                    start_it = next(start_it);
                    mismatch_diff = end_it->second - start_it->second - 1;
                }

                while (mismatch_diff <= mismatches and end_it != valid_end_locs.end()) {
                    end_it = next(end_it);
                    mismatch_diff = end_it->second - start_it->second - 1;
                }

                start_mismatch = start_it->first;
                end_mismatch = prev(end_it)->first;

                if (start_mismatch >= initial_gap ) {
                	break;
                }

                if (isOdd) {
	                left = int(c - end_mismatch);
	                right = int(c + end_mismatch);
	                gap = int(2.0 * (start_mismatch + 1.0) + 1.0);
	            }
	            else {
	                left = int(c - 0.5 - (end_mismatch - 1.0));
	                right = int(c + 0.5 + (end_mismatch - 1.0));
	                gap = int(2.0 * (start_mismatch + 1.0));
	            }

	            int outer_left = left + 1;
			    int outer_right = right + 1;

	            if ( (right - left + 1 - gap) / 2 >= min_len and (right - left + 1 - gap) / 2 <= max_len ) {
	            	palindromes->insert(tuple<int, int, int>(left, right, gap));
	            }

                start_it = next(start_it);
            }
        }
    }
}

////////////////////
//  MAIN PROGRAM  //
////////////////////

int main(int argc, char* argv[]) {

        ///////////////////////
        //  REQUEST OPTIONS  //
        ///////////////////////

        string input_file = "";
        int min_len = 10;
        int max_len = 100;
        int max_gap = 100;
        int mismatches = 0;
        string output_file = "IUPACpal.out";

	    int c;
	    while( ( c = getopt (argc, argv, "f:m:M:g:x:o:") ) != -1 ) 
	    {
	        switch(c)
	        {
	            case 'f':
	                if(optarg) input_file = optarg;
	                break;
	            case 'm':
	                if(optarg) min_len = std::atoi(optarg);
	                break;
	            case 'M':
	                if(optarg) max_len = std::atoi(optarg);
	                break;
	            case 'g':
	                if(optarg) max_gap = std::atoi(optarg);
	                break;
	            case 'x':
	                if(optarg) mismatches = std::atoi(optarg);
	                break;
	            case 'o':
	                if(optarg) output_file = optarg;
	                break;
	        }
	    }

	    ///////////////////////
        //  VERIFY OPTIONS   //
        ///////////////////////

        if (!exist(input_file.c_str())) {  usage(); cout << "Error: input_file not found" << endl; return -1; }

	    ifstream in(input_file);
        string contents((std::istreambuf_iterator<char>(in)), 
        istreambuf_iterator<char>());
        long int n = contents.length();
        unsigned char * seq = ( unsigned char* ) malloc( ( n ) * sizeof( unsigned char ) );

        for (int i = 0; i < n; ++i) {
            seq[i] = contents[i];
        }

        if (min_len < 2) { usage(); cout << "Error: min_len must not be less than 2" << endl; return -1; }
        if (min_len > INT_MAX) { usage(); cout << "Error: min_len must not greater than " << INT_MAX << endl; return -1; }
        if (max_len < 0) { usage(); cout << "Error: max_len must not be a negative value" << endl; return -1; }
        if (max_len > INT_MAX) { usage(); cout << "Error: max_len must not greater than " << INT_MAX << endl; return -1; }
        if (max_gap < 0) { usage(); cout << "Error: max_gap must not be a negative value" << endl; return -1; }
        if (max_gap > INT_MAX) { usage(); cout << "Error: max_gap must not greater than " << INT_MAX << endl; return -1; }
        if (mismatches < 0) { usage(); cout << "Error: mismatches must not be a negative value" << endl; return -1; }
        if (mismatches > INT_MAX) { usage(); cout << "Error: mismatches must not greater than " << INT_MAX << endl; return -1; }

        if (min_len >= n) { usage(); cout << "Error: min_len must be less than sequence length" << endl; return -1; }
        if (max_len < min_len) { usage(); cout << "Error: max_len must not be less than min_len" << endl; return -1; }
        if (max_gap >= n) { usage(); cout << "Error: max_gap must be less than sequence length" << endl; return -1; }
        if (min_len >= n) { usage(); cout << "Error: min_len must be less than sequence length" << endl; return -1; }
        if (mismatches >= n) { usage(); cout << "Error: mismatches must be less than sequence length" << endl; return -1; }
        if (mismatches >= min_len) { usage(); cout << "Error: mismatches must be less than min_len" << endl; return -1; }
        
        // Display user given options
        if (true) {
        	cout << endl;
            cout << "input_file: " << input_file << endl;

            // Optional printing of sequence data
            if (false) {
            	cout << "sequence: " << contents << endl;
            }

            cout << "min_len: " << min_len << endl;
            cout << "max_len: " << max_len << endl;
            cout << "max_gap: " << max_gap << endl;
            cout << "mismatches: " << mismatches << endl;
            cout << "output_file: " << output_file << endl;
            cout << endl;
        }

        //////////////////////////////
        //  Determine Match Matrix  //
        //////////////////////////////

        I_map IUPAC_map;
        int IUPAC_map_count;
        int IUPAC_to_value[128];
        int complement[128];

        // Initialize IUPAC_to_value and base_complement
        for (int i = 0; i < 128; ++i) {
            IUPAC_to_value[i] = -1;
            complement[i] = -1;
        }

        // Build IUPAC_map and IUPAC_to_value
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'a', {'a'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'c', {'c'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'g', {'g'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 't', {'t'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'u', {'t'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'r', {'a', 'g'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'y', {'c', 't'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 's', {'g', 'c'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'w', {'a', 't'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'k', {'g', 't'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'm', {'a', 'c'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'b', {'c', 'g', 't'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'd', {'a', 'g', 't'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'h', {'a', 'c', 't'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'v', {'a', 'c', 'g'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, 'n', {'a', 'c', 'g', 't'});

        // Note: 'u' and 't' are identical

        // Non-IUPAC characters $ and # to be used within suffix tree
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, '$', {'$'});
        IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, '#', {'#'});

        IUPAC_map_count = IUPAC_map.size();

        char value_to_IUPAC[IUPAC_map_count];
        bool* match_matrix = (bool*) malloc(IUPAC_map_count * IUPAC_map_count * sizeof(bool));

        // Build value_to_IUPAC
        for (int i = 0; i < 128; ++i) {
            if (IUPAC_to_value[i] != -1) {
                value_to_IUPAC[IUPAC_to_value[i]] = i;
            }
        }

        // IUPAC_to_value: Given a character from IUPAC_map, returns a 0-based index
        // according to the order of insertion into IUPAC_map
        //
        // value_to_IUPAC: Given a value from 0 to n - 1 (n = IUPAC_map size)
        // returns the corresponding character, according to IUPAC_to_value
        //
        // i.e. IUPAC_to_value and value_to_IUPAC are inverse functions

        // Build match_matrix
        for(I_map::iterator it1 = IUPAC_map.begin(); it1 != IUPAC_map.end(); it1++) {
            for(I_map::iterator it2 = IUPAC_map.begin(); it2 != IUPAC_map.end(); it2++) {
                char char1 = it1->first;
                char char2 = it2->first;
                int i = IUPAC_to_value[char1];
                int j = IUPAC_to_value[char2];

                bool match = false;
                for(set<char>::iterator it_set1 = it1->second.begin(); it_set1 != it1->second.end(); it_set1++) {
                    for(set<char>::iterator it_set2 = it2->second.begin(); it_set2 != it2->second.end(); it_set2++) {
                        if (*it_set1 == *it_set2) {
                            match = true;
                            break;
                        }
                    }
                }

                if (match) {
                    set_value_for_matrix(match_matrix, IUPAC_map_count, i, j, true);
                }
                else {
                    set_value_for_matrix(match_matrix, IUPAC_map_count, i, j, false);
                }
            }
        }

        // Assign data structures to MatchMatrix class to enable use of MatchMatrix::match function
        MatchMatrix::match_matrix = match_matrix;
        MatchMatrix::IUPAC_map_count = IUPAC_map_count;
        MatchMatrix::IUPAC_to_value = IUPAC_to_value;

        // Optionally print match_matrix
        if (false) {
            string letters = "acgturyswkmbdhvn$#";

            cout << "IUPAC_map_count " << IUPAC_map_count << endl;

            cout << "  ";
            for (int i = 0; i < IUPAC_map_count; ++i) {
                cout << letters[i] << " ";
            }
            cout << endl;

            for (int i = 0; i < IUPAC_map_count; ++i) {
                cout << letters[i] << " ";
                for (int j = 0; j < IUPAC_map_count; ++j) {
                    cout << MatchMatrix::match(letters[i], letters[j]) << " ";
                }
                cout << endl;
            }

            cout << endl;
        }

        // Build complement array
        complement['a'] = 't';
        complement['c'] = 'g';
        complement['g'] = 'c';
        complement['t'] = 'a';
        complement['u'] = 'a';
        complement['r'] = 'y';
        complement['y'] = 'r';
        complement['s'] = 's';
        complement['w'] = 'w';
        complement['k'] = 'm';
        complement['m'] = 'k';
        complement['b'] = 'v';
        complement['d'] = 'h';
        complement['h'] = 'd';
        complement['v'] = 'b';
        complement['n'] = 'n';

        ////////////////////////////////////////////////////
        //  Construct S = seq $ complement(reverse(seq) # //
        ////////////////////////////////////////////////////

        int S_n = 2 * n + 2;
        unsigned char S[S_n];

        for (int i = 0; i < n; ++i) {
            S[i] = seq[i];
        }

        S[n] = '$';

        for (int i = 0; i < n; ++i) {
            S[n + 1 + i] = complement[seq[n - 1 - i]];
        }

        S[2 * n + 1] = '#';

        ////////////////////
        //  Calculate SA  //
        ////////////////////

        INT * SA;
        SA = ( INT * ) malloc( ( S_n ) * sizeof( INT ) );

        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }

        #ifdef _USE_64
                if( divsufsort64( S, SA,  S_n ) != 0 )
                {
                        fprintf(stderr, " Error: SA computation failed.\n" );
                        exit( EXIT_FAILURE );
                }
        #endif

        #ifdef _USE_32
                if( divsufsort( S, SA,  S_n ) != 0 )
                {
                        fprintf(stderr, " Error: SA computation failed.\n" );
                        exit( EXIT_FAILURE );
                }
        #endif

        /////////////////////
        // Calculate invSA //
        /////////////////////

        INT * invSA;
        invSA = ( INT * ) malloc( S_n * sizeof( INT ) );

        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < S_n; i ++ )
        {
                invSA [SA[i]] = i;
        }

        /////////////////////
        //  Calculate LCP  //
        /////////////////////

        INT * LCP;
        LCP = ( INT * ) malloc  ( S_n * sizeof( INT ) );

        if( ( LCP == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        if( LCParray( S, S_n, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }

        ////////////////////////////
        //  Calculate RMQ of LCP  //
        ////////////////////////////

        #ifdef _USE_NLOGN_RMQ
            INT l = S_n;
            INT lgl = flog2( l );
            INT * A = ( INT * ) calloc( ( INT ) l * lgl, sizeof(INT) );
            rmq_preprocess(A, LCP, l);
        #else
            int_vector<> v(S_n , 0); // Create a vector of length n and initialize it with 0s
            for ( INT i = 0; i < S_n; i ++ )
            {
                    v[i] = LCP[i];
            }

            rmq_succinct_sct<> rmq(&v);
        #endif

        // Optional printing of data structures
        if (false) {
            cout << endl << endl;
            print_array("  seq", seq, n);
            print_array("    S", S, S_n, true);
            print_array("   SA", SA, S_n, true);
            print_array("invSA", invSA, S_n, true);
            print_array("  LCP", LCP, S_n, true);
            cout << endl << endl;
        }

        /////////////////////////////
        //  Calculate Palindromes  //
        /////////////////////////////

        set<tuple<int, int, int>> palindromes;

        // All palindromes calculate and stored
        #ifdef _USE_NLOGN_RMQ
            addPalindromes(&palindromes, S, S_n, n, invSA, LCP, A, tuple<int, int, int, int>(min_len, max_len, mismatches, max_gap));
        #else
            addPalindromes(&palindromes, S, S_n, n, invSA, LCP, rmq, tuple<int, int, int, int>(min_len, max_len, mismatches, max_gap));
        #endif

        /////////////////////////
        //  Print Palindromes  //
        /////////////////////////

        ofstream file;
        file.open(output_file);

        file << "Palindromes of:  " << input_file << endl;
        file << "Sequence length is: " << n << endl;
        file << "Start at position: " << 1 << endl;
        file << "End at position: " << n << endl;
        file << "Minimum length of Palindromes is: "  << min_len << endl;
        file << "Maximum length of Palindromes is: "  << max_len << endl;
        file << "Maximum gap between elements is: "  << max_gap << endl;
        file << "Number of mismatches allowed in Palindrome: " << mismatches << endl;
        file << endl << endl << endl;
        file << "Palindromes:" << endl;

        // Dummy entry to ensure all previous palindromes are encountered during sorting
        palindromes.insert(tuple<int, int, int>(S_n, S_n, 0));

        if (!palindromes.empty()) {
	        int prev_left = get<0>(*palindromes.begin());

	        set<tuple<int, int, int>>::iterator it_left_to_right = palindromes.begin();
	        set<tuple<int, int, int>>::iterator it;

	        while ( it_left_to_right != palindromes.end() ) {

	        	int left = get<0>(*it_left_to_right);

	            if (prev_left != left) {
	            	it = prev(it_left_to_right);
	            	while (it != prev(palindromes.begin()) and get<0>(*it) == prev_left) {
	            		
	            		int left = get<0>(*it);
			            int right = get<1>(*it);
			            int gap = get<2>(*it);

			            int outer_left = left + 1;
			            int outer_right = right + 1;
			            int inner_left = (outer_left + outer_right - 1 - gap) / 2;
			            int inner_right = (outer_right + outer_left + 1 + gap) / 2;

			            string pad = "         ";
			            int pad_length = pad.size();

			            file << outer_left;
			            for (int i = 0; i < pad_length - getDigitCount(outer_left); ++i) { file << " "; }
			            for (int i = outer_left; i <= inner_left; ++i) { file << seq[i - 1]; }
			            for (int i = 0; i < pad_length - getDigitCount(inner_left); ++i) { file << " "; }
			            file << inner_left;

			            file << "\n";

			            file << pad;
			            for (int i = 0; i < (inner_left - outer_left + 1); ++i) {
			                file << ( (MatchMatrix::match(seq[ outer_left - 1 + i ], complement[ seq[ outer_right - 1 - i ] ])) ? "|" : " " );
			            }
			            
			            file << "\n";

			            file << outer_right;
			            for (int i = 0; i < pad_length - getDigitCount(outer_right); ++i) { file << " "; }
			            for (int i = outer_right; i >= inner_right; --i) { file << seq[i - 1]; }
			            for (int i = 0; i < pad_length - getDigitCount(inner_right); ++i) { file << " "; }
			            file << inner_right;

			            file << "\n" << "\n";

	            		it--;
	            	}
	            }

	            prev_left = left;
	            it_left_to_right++;
	        }
	    }

        file << endl << endl << endl;

        file.close();

        return 0;
}