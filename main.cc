#include <iostream>
#include <string>
#include <map>
#include <set>
#include <list>
#include <tuple>
#include <sys/stat.h>
#include <limits.h>
#include "main.h"

typedef map< char, set<char> > I_map;


///////////////////////////////////////////
//  RANGE MINIMUM QUERY (RMQ) FUNCTIONS  //
///////////////////////////////////////////

// Range Minimum Query (Type 1)
static __inline INT rmq(INT *m, INT *v, INT n, INT i, INT j) {
    INT lgn = flog2(n);

    if (i > j) {INT tmp = j; j = i; i = tmp;}
    i++;
    if (i == j) return i;

    INT k = flog2(j-i+1);
    INT a = m[i * lgn + k];
    INT shift = ((( INT ) 1) << k);
    INT b = m[(j - shift + 1) * lgn + k];

    return v[a]>v[b]?b:a;
}

// O(nlogn)-time preprocessing function for Type 1 Range Minimum Queries
void rmq_preprocess(INT * m, INT * v, INT n)
{
    INT i, j;
    INT lgn = flog2(n);

    for (i = 0; i < n; i++) {
        m[i*lgn] = i;
    }

    for (j = 1; ((( INT ) 1) << j) <= n; j++) {
        for (i = 0; i + (1 << j) - 1 < n; i++) {
            if (v[m[i*lgn + j - 1]] < v[m[(i + (1 << (j - 1)))*lgn + j - 1]]) {
                m[i*lgn + j] = m[i*lgn + j - 1];
            }
            else {
                m[i*lgn + j] = m[(i + (1 << (j - 1)))*lgn + j - 1];
            }
        }
    }
}


//////////////////////////////
//  BASIC HELPER FUNCTIONS  //
//////////////////////////////

// Checks if a filename exists
int exist(const char *name) {
    struct stat buffer;
    return (stat (name, &buffer) == 0);
}

// Prints documentation on usage of IUPACpal arguments
void usage() {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  FLAG  PARAMETER       TYPE      DEFAULT         DESCRIPTION\n" );
    fprintf ( stdout, "  -f    input_file      <str>     input.fasta     Input filename (FASTA).\n" );
    fprintf ( stdout, "  -s    seq_name        <str>     seq0            Input sequence name.\n");
    fprintf ( stdout, "  -m    min_len         <int>     10              Minimum length.\n");
    fprintf ( stdout, "  -M    max_len         <int>     100             Maximum length.\n");
    fprintf ( stdout, "  -g    max_gap         <int>     100             Maximum permissible gap.\n");
    fprintf ( stdout, "  -x    mismatches      <int>     0               Maximum permissible mismatches.\n");
    fprintf ( stdout, "  -o    output_file     <str>     IUPACpal.out    Output filename.\n" );
    fprintf ( stdout, "\n" );
}

// Convert a string to any type, if possible
template <typename T>
T ConvertString(const string &data, int* return_code) {
    if (!data.empty()) {
        T ret;
        istringstream iss( data );

        if (data.find( "0x" ) != std::string::npos) {
            iss >> std::hex >> ret;
        }
        else {
            iss >> std::dec >> ret;
        }

        if (iss.fail()) {
            *return_code = -1; // Error
            return T ( ); // Empty type
        }

        *return_code = 0; // Valid input
        return ret; // Value of converted input
    }

    *return_code = 1; // Empty input

    return T( ); // Empty type
}

// Get character length of decimal representation of an integer
int getDigitCount(int x) {
    stringstream iss;
    iss << x;
    return iss.str().size();
}

// Print an array of any type
template<typename T>
void print_array(string title, T seq, long int n, bool print_indices=false) {
    if (print_indices) {
        for(int i = 0; i < title.size() + 2; i++) {
            cout << " ";
        }

        for(int i = 0; i < n; i++) {
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
    for(int i=0; i<n; i++) {
        cout << seq[i] << "  ";
    }

    cout << endl << endl;
}


////////////////////////////////////////////////
//  IUPAC CHARACTER MATRIX CLASS & FUNCTIONS  //
////////////////////////////////////////////////

// Inserts IUPAC characters mapping information into an STSL map
void IUPAC_map_insert(I_map* IUPAC_map, int IUPAC_to_value[], char IUPAC_char, set<char> mapped_chars) {
    static int index = 0;
    IUPAC_map->insert( pair< char, set<char> > (IUPAC_char, mapped_chars) );
    IUPAC_to_value[IUPAC_char] = index++;
}

// Set a value in a matrix stored as a 1D contiguous sequence of memory
template<typename T>
void set_value_for_matrix(T* mtx_ptr, uint32_t ncols, uint32_t i, uint32_t j, T val) {
    mtx_ptr[i + i * (ncols - 1) + j] = val;
}

// Get a value from a matrix stored as a 1D contiguous sequence of memory
template<typename T>
T get_value_from_matrix(T* mtx, uint32_t ncols, uint32_t i, uint32_t j) {
    return mtx[i + i * (ncols - 1) + j];
}

// Class to store preprocessed information relating to the match status of IUPAC characters
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

// Calculates the Longest Common Prefix array of a text and stores value in given variable LCP
//
// INPUT:
// - Text
// - Text length
// - Suffix Array
// - Longest Common Prefix data structure (empty)
unsigned int LCParray(unsigned char *text, INT n, INT * SA, INT * invSA, INT * LCP)
{
    INT i = 0, j = 0;
    LCP[0] = 0;

    for (i = 0; i < n; i++) {
        if (invSA[i] != 0) {
            if ( i == 0) {
                j = 0;
            }
            else {
                j = (LCP[invSA[i-1]] >= 2) ? LCP[invSA[i-1]]-1 : 0;
            }

            while ( text[i+j] == text[SA[invSA[i]-1]+j] ) {
                j++;
            }

            LCP[invSA[i]] = j;
        }
    }

    return 1;
}

// Returns the Longest Common Extension between position i and j (order of i, j input does not matter)
//
// INPUT:
// - Indexes i and j
// - Text length
// - Inverse Suffix Array
// - Longest Common Prefix Array data structure (filled)
// - Data structure (filled) with preprocessed values to perform Range Minimum Queries (Type 1: 'A', Type 2: 'rmq')
#ifdef _USE_NLOGN_RMQ
// Using Type 1 RMQs
unsigned int LCE(INT i, INT j, INT n, INT * invSA, INT * LCP, INT * A) {
#else
// Using Type 2 RMQs
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
        return LCP[rmq(A, LCP, n, a, b)]; // Using Type 1 RMQs: rmq(a, b) does not include 'a' value in range
    #else
        return LCP[rmq(a + 1, b)]; // Using Type 2 RMQs: rmq(a, b) does include 'a' value in range
    #endif
}

// Calculates a list of Longest Common Extensions, corresponding to 0, 1, 2, etc. allowed mismatches, up to maximum number of allowed mismatches
//
// EXTRA INFO:
// - Only considers "real" mismatches (degenerate string mismatching according to IUPAC character matrix)
// - Takes into account the matching possibility of non A, C, G, T/U characters
// - Longest Common Extension calculated from positions i and j (order of i, j input does not matter)
// - Only starts counting number of allowed mismatches that occur after the given initial gap, however earlier mismatches are still stored
// - Should only be used after MatchMatrix has been instantiated with necessary data
//
// INPUT:
// - Text
// - Indexes i and j
// - Text length
// - Inverse Suffix Array
// - Longest Common Prefix Array (LCP)
// - Data structure (filled) with preprocessed values to perform Range Minimum Queries (Type 1: 'A', Type 2: 'rmq')
// - Maximum number of allowed mismatches
// - Initial gap
// - Data structure to store resulting mismatch locations
#ifdef _USE_NLOGN_RMQ
// Using Type 1 RMQs
void realLCE_mismatches(unsigned char* text, INT i, INT j, INT n, INT * invSA, INT * LCP, INT * A, int mismatches, int initial_gap, list<int>* mismatch_locs) {
#else
// Using Type 2 RMQs
void realLCE_mismatches(unsigned char* text, INT i, INT j, INT n, INT * invSA, INT * LCP, rmq_succinct_sct<> rmq, int mismatches, int initial_gap, list<int>* mismatch_locs) {
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

// Finds all inverted repeats (palindromes) with given parameters and adds them to an output set
//
// INPUT:
// - Data structure (set of integer 3-tuples) to store palindromes in form (left_index, right_index, gap)
// - S = text + '$' + complement(reverse(text) + '#'
// - Length of S
// - Text length
// - Inverse Suffix Array
// - Longest Common Prefix Array (LCP)
// - Data structure (filled) with preprocessed values to perform Range Minimum Queries (Type 1: 'A', Type 2: 'rmq')
// - Tuple of parameters for palindromes to be found (minimum_length, maximum_length, maximum_allowed_number_of_mismatches, maximum_gap)
#ifdef _USE_NLOGN_RMQ
// Using Type 1 RMQs
void addPalindromes(set<tuple<int, int, int>>* palindromes, unsigned char* S, int S_n, int n, INT * invSA, INT * LCP, INT * A, tuple<int, int, int, int> params) {
#else
// Using Type 2 RMQs
void addPalindromes(set<tuple<int, int, int>>* palindromes, unsigned char* S, int S_n, int n, INT * invSA, INT * LCP, rmq_succinct_sct<> rmq, tuple<int, int, int, int> params) {
#endif
    // Retrieve parameters from tuple
    int min_len = get<0>(params);
    int max_len = get<1>(params);
    int mismatches = get<2>(params);
    int max_gap = get<3>(params);

    // Cycle through possible centres of text from left to right
    for (double c = 0; c <= (n - 1); c += 0.5 ) {
        // Determine if value of centre corresponds to an odd or even palindrome
        bool isOdd = (trunc(c) == c);

        // Strategically choose i and j to determine maximum extension with text
        int i, j;

        if (isOdd) {
            i = int(c + 1.0);
            j = int(2.0 * n + 1.0 - c);
        } else {
            i = int(c + 0.5);
            j = int(2.0 * n + 1.0 - (c + 0.5));
        }

        // Calculate initial number of characters ignored when performing Longest Common Extensions to determine palindromes
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

        // Calculate list of relevant mismatch locations when performing Longest Common Extensions in 'kangaroo method' style
        list<int> mismatch_locs;

        #ifdef _USE_NLOGN_RMQ
            realLCE_mismatches(S, i, j, S_n, invSA, LCP, A, mismatches, initial_gap, &mismatch_locs);
        #else
            realLCE_mismatches(S, i, j, S_n, invSA, LCP, rmq, mismatches, initial_gap, &mismatch_locs);
        #endif

        // Always set -1 as a mismatch location
        mismatch_locs.push_front(-1);

        // Create data structures to store valid start and end mismatch locations alongside an ID value (mismatch_location, mismatch_id)
        list<pair<int, int>> valid_start_locs, valid_end_locs;

        // Determine list of valid start and end mismatch locations (that could mark the potential start or end of a palindrome)
        int mismatch_id = 0;
        for (list<int>::iterator it = mismatch_locs.begin(); it != mismatch_locs.end(); ++it){
            if (next(it) != mismatch_locs.end() and *next(it) != *it + 1) {
                valid_start_locs.push_back(pair<int, int> (*it, mismatch_id));
            }

            if (it != mismatch_locs.begin() and *prev(it) != *it - 1) {
                valid_end_locs.push_back(pair<int, int> (*it, mismatch_id));
            }

            mismatch_id++;
        }

        // Optional printing of mismatch locations relative to centre, valid start locations and valid end locations
        #ifdef _DIAGNOSTICS
            cout << "centre = " << c << "\n";
            cout << "mismatches: " << "\t";
            cout << "[ ";
            for (list<int>::iterator it = mismatch_locs.begin(); it != mismatch_locs.end(); ++it){
                cout << *it << " ";
            }
            cout << "]" << endl;
            cout << "starts: " << "\t";
            cout << "[ ";
            for (list<pair<int, int>>::iterator it = valid_start_locs.begin(); it != valid_start_locs.end(); ++it){
                cout << "(" << it->first << ", " << it->second << ") ";
            }
            cout << "]" << endl;

            cout << "ends: " << "\t\t";
            cout << "[ ";
            for (list<pair<int, int>>::iterator it = valid_end_locs.begin(); it != valid_end_locs.end(); ++it){
                cout << "(" << it->first << ", " << it->second << ") ";
            }
            cout << "]" << endl;

            cout << endl;
        #endif

        // Check if valid start and end mismatch locations have been found
        if ( !valid_start_locs.empty() and !valid_end_locs.empty() ) {
            list<pair<int, int>>::iterator start_it = valid_start_locs.begin();
            list<pair<int, int>>::iterator end_it = valid_end_locs.begin();

            int mismatch_diff, left, right, gap;
        	int start_mismatch, end_mismatch;

            // Loop while both start and end mismatch locations have not reached the end of their respective lists
            while( start_it != valid_start_locs.end() and end_it != valid_end_locs.end()) {
                // Count the difference in mismatches between the start and end location
                mismatch_diff = end_it->second - start_it->second - 1;

                // While mismatch difference is too large, move start location to the right until mismatch difference is within acceptable bound
                while (mismatch_diff > mismatches) {
                    start_it = next(start_it);
                    mismatch_diff = end_it->second - start_it->second - 1;
                }

                // While mismatch difference is within acceptable bound, move end location to the right until mismatch difference becomes unacceptable
                while (mismatch_diff <= mismatches and end_it != valid_end_locs.end()) {
                    end_it = next(end_it);
                    mismatch_diff = end_it->second - start_it->second - 1;
                }

                start_mismatch = start_it->first; // Pick the current start mismatch
                end_mismatch = prev(end_it)->first; // Pick the end mismatch directly after the current end mismatch

                // Skip this iteration if the start mismatch chosen is such that the gap is not within the acceptable bound
                if (start_mismatch >= initial_gap ) {
                	break;
                }

                // Optionally view diagnostics information
                #ifdef _DIAGNOSTICS
                    cout << "(start_mismatch, end_mismatch) = " << start_mismatch << " " << end_mismatch << endl;
                #endif

                // Set left, right indexes and gap of potential palindrome, according to chosen start and end mismatch
                if (isOdd) {
	                left = int(c - end_mismatch);
	                right = int(c + end_mismatch);
	                gap = 2 * (start_mismatch + 1) + 1;
	            }
	            else {
	                left = int(c - 0.5 - (end_mismatch - 1.0));
	                right = int(c + 0.5 + (end_mismatch - 1.0));
	                gap = 2 * (start_mismatch + 1);
	            }

                // Optionally view diagnostics information
                #ifdef _DIAGNOSTICS
                    cout << "(left, gap, right) = " << left << " " << right << " " << gap << endl << endl;
                #endif

                // Check that potential palindrome is not too short
	            if ((right - left + 1 - gap) / 2 >= min_len) {
                    // Check that potentialinput_file palindrome is not too long
                    if ((right - left + 1 - gap) / 2 <= max_len) {
                        // Palindrome is not too long, so add to output
                        palindromes->insert(tuple<int, int, int>(left, right, gap));
                    }
                    else {
                        // Palindrome is too long, so attempt truncation
                        int prev_end_mismatch = prev(prev(end_it))->first;
                        int mismatch_gap = end_mismatch - prev_end_mismatch - 1;
                        int overshoot = ( (right - left + 1 - gap) / 2 ) - max_len;

                        // Check if truncation results in the potential palindrome ending in a mismatch
                        if (overshoot != mismatch_gap) {
                            // Potential palindrome does not end in a mismatch, so add to output
                            palindromes->insert(tuple<int, int, int>(left + overshoot, right - overshoot, gap));
                        }
                        else {
                            // Potential palindrome does end in a mismatch, so truncate an additional 1 character either side then add to output
                            palindromes->insert(tuple<int, int, int>(left + overshoot + 1, right - overshoot - 1, gap));
                        }
                    }
	            }

                // Go to next start mismatch in list and loop
                start_it = next(start_it);
            }
        }

        #ifdef _DIAGNOSTICS
            cout << "----------------" << endl << endl;
        #endif
    }
}


//////////////////////
//  MAIN EXECUTION  //
//////////////////////

int main(int argc, char* argv[]) {

    ///////////////////////
    //  REQUEST OPTIONS  //
    ///////////////////////

    // Default parameters
    string input_file = "input.fasta";
    string seq_name = "seq0";
    int min_len = 10;
    int max_len = 100;
    int max_gap = 100;
    int mismatches = 0;
    string output_file = "IUPACpal.out";

    // Parse command line arguments
    int c;
    while( ( c = getopt (argc, argv, "f:s:m:M:g:x:o:") ) != -1 )
    {
        switch(c)
        {
            case 'f':
                if(optarg) input_file = optarg;
                break;
            case 's':
                if(optarg) seq_name = optarg;
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
;
    ///////////////////////
    //  VERIFY OPTIONS   //
    ///////////////////////

    // Check input file exists, exit if it does not
    if (!exist(input_file.c_str())) {  usage(); cout << "Error: File '" + input_file + "' not found." << endl; return -1; }

    // Read input file
    ifstream input( input_file );
    string name = "";
    string contents = "";
    bool found_seq = false;

    // Parse file line by line
    for( string line; getline( input, line ); )
    {
        int line_length = line.length();

        // Look for sequence name
        if (found_seq == false) {
            if (line_length > 0 && line[0] == '>') {
                int i = 1;
                while ( i < line_length && line[i] == ' ') {
                    i++;
                }
                while ( i < line_length && line[i] != ' ' ) {
                    name += line[i];
                    i++;
                }
            }

            // Sequence name found
            if (name == seq_name) {
                found_seq = true;
            }

            name = "";
        }
        // Once sequence name is found, extract sequence
        else {
            if (line[0] != ' ' && line[0] != '>' && line[0] != ';') {
                contents += line;
            }
            else {
                // End of sequence
                break;
            }
        }
    }

    // Check if sequence name was found, exit if not
    if (!found_seq) {  usage(); cout << "Error: Sequence '" + seq_name + "' not found in file '" + input_file + "'." << endl; return -1; }

    long int n = contents.length();
    unsigned char * seq = ( unsigned char* ) malloc( ( n ) * sizeof( unsigned char ) );

    // Convert extracted sequence to character array, all lowercase
    for (int i = 0; i < n; ++i) {
        seq[i] = contents[i];
        seq[i] = tolower(seq[i]);
    }

    // Verify arguments are valid with respect to individual limits
    if (min_len < 2) { usage(); cout << "Error: min_len must not be less than 2." << endl; return -1; }
    if (min_len > INT_MAX) { usage(); cout << "Error: min_len must not greater than " << INT_MAX << "." << endl; return -1; }
    if (max_len < 0) { usage(); cout << "Error: max_len must not be a negative value." << endl; return -1; }
    if (max_len > INT_MAX) { usage(); cout << "Error: max_len must not greater than " << INT_MAX << "." << endl; return -1; }
    if (max_gap < 0) { usage(); cout << "Error: max_gap must not be a negative value." << endl; return -1; }
    if (max_gap > INT_MAX) { usage(); cout << "Error: max_gap must not greater than " << INT_MAX << "." << endl; return -1; }
    if (mismatches < 0) { usage(); cout << "Error: mismatches must not be a negative value." << endl; return -1; }
    if (mismatches > INT_MAX) { usage(); cout << "Error: mismatches must not greater than " << INT_MAX << "." << endl; return -1; }

    // Verify arguments are valid with respect to each other
    if (min_len >= n) { usage(); cout << "Error: min_len must be less than sequence length." << endl; return -1; }
    if (max_len < min_len) { usage(); cout << "Error: max_len must not be less than min_len." << endl; return -1; }
    if (max_gap >= n) { usage(); cout << "Error: max_gap must be less than sequence length." << endl; return -1; }
    if (min_len >= n) { usage(); cout << "Error: min_len must be less than sequence length." << endl; return -1; }
    if (mismatches >= n) { usage(); cout << "Error: mismatches must be less than sequence length." << endl; return -1; }
    if (mismatches >= min_len) { usage(); cout << "Error: mismatches must be less than min_len." << endl; return -1; }

    // Optionally display user given options
    if (true) {
    	cout << endl;
        cout << "input_file: " << input_file << endl;
        cout << "seq_name: " << seq_name << endl;

        // Optionally print sequence data
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

    //////////////////////////
    //  BUILD MATCH MATRIX  //
    //////////////////////////

    I_map IUPAC_map;
    int IUPAC_map_count;
    int IUPAC_to_value[128];
    int complement[128];

    // Initialize arrays
    for (int i = 0; i < 128; ++i) {
        IUPAC_to_value[i] = -1;
        complement[i] = -1;
    }

    // Build IUPAC_map and IUPAC_to_value
    // Note: 'u' and 't' are considered identical
    // Note: 'n' and '*' and '-' are considered identical
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
    IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, '*', {'a', 'c', 'g', 't'});
    IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, '-', {'a', 'c', 'g', 't'});

    // Non-IUPAC characters $ and # will be used within suffix tree
    IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, '$', {'$'});
    IUPAC_map_insert(&IUPAC_map, IUPAC_to_value, '#', {'#'});

    IUPAC_map_count = IUPAC_map.size();

    char value_to_IUPAC[IUPAC_map_count];
    bool* match_matrix = (bool*) malloc(IUPAC_map_count * IUPAC_map_count * sizeof(bool));

    // Build array
    for (int i = 0; i < 128; ++i) {
        if (IUPAC_to_value[i] != -1) {
            value_to_IUPAC[IUPAC_to_value[i]] = i;
        }
    }

    // INFO:
    //
    // IUPAC_to_value:
    // Given a character from IUPAC_map, returns a 0-based index
    // according to the order of insertion into IUPAC_map
    //
    // value_to_IUPAC:
    // Given a value from 0 to n - 1 (n = IUPAC_map size)
    // returns the corresponding character, according to IUPAC_to_value
    //
    // i.e. IUPAC_to_value and value_to_IUPAC are inverse functions

    // Build match matrix
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

    // Assign static data structures to MatchMatrix class to enable use of MatchMatrix::match function
    MatchMatrix::match_matrix = match_matrix;
    MatchMatrix::IUPAC_map_count = IUPAC_map_count;
    MatchMatrix::IUPAC_to_value = IUPAC_to_value;

    // Optionally print match matrix
    #ifdef _DIAGNOSTICS
        string letters = "acgturyswkmbdhvn*-$#";

        cout << "Match Matrix:" << endl;

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

        cout << endl << endl;
    #endif

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
    complement['*'] = 'n';
    complement['-'] = 'n';

    //////////////////////////////////////////////////////////////
    //  CONSTRUCT S = seq + '$' + complement(reverse(seq) + '#' //
    //////////////////////////////////////////////////////////////

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

    ///////////////////////////////////
    //  CALCULATE Suffix Array (SA)  //
    ///////////////////////////////////

    INT * SA;
    SA = ( INT * ) malloc( ( S_n ) * sizeof( INT ) );

    if( ( SA == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
        return 0;
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

    //////////////////////////////////////////////
    //  CALCULATE Inverse Suffix Array (invSA)  //
    //////////////////////////////////////////////

    INT * invSA;
    invSA = ( INT * ) malloc( S_n * sizeof( INT ) );

    if( ( invSA == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        return 0;
    }

    for ( INT i = 0; i < S_n; i ++ )
    {
        invSA [SA[i]] = i;
    }

    ///////////////////////////////////////////////////
    //  CALCULATE Longest Common Prefix Array (LCP)  //
    ///////////////////////////////////////////////////

    INT * LCP;
    LCP = ( INT * ) malloc  ( S_n * sizeof( INT ) );

    if( ( LCP == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        return 0;
    }

    if( LCParray( S, S_n, SA, invSA, LCP ) != 1 )
    {
        fprintf(stderr, " Error: LCP computation failed.\n" );
        exit( EXIT_FAILURE );
    }

    ////////////////////////////
    //  CALCULATE RMQ of LCP  //
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
    //  CALCULATE PALINDROMES  //
    /////////////////////////////

    set<tuple<int, int, int>> palindromes;

    // All palindromes calculate and stored
    #ifdef _USE_NLOGN_RMQ
        addPalindromes(&palindromes, S, S_n, n, invSA, LCP, A, tuple<int, int, int, int>(min_len, max_len, mismatches, max_gap));
    #else
        addPalindromes(&palindromes, S, S_n, n, invSA, LCP, rmq, tuple<int, int, int, int>(min_len, max_len, mismatches, max_gap));
    #endif

    /////////////////////////
    //  PRINT PALINDROMES  //
    /////////////////////////

    ofstream file;
    file.open(output_file);

    file << "Palindromes of: " << input_file << endl;
    file << "Sequence name: " << seq_name << endl;
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

    cout << "Search complete!" << endl;

    free(match_matrix);
    free(seq);
    free(SA);
    free(invSA);
    free(LCP);

    return 0;
}
