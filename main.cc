#include <iostream>
#include <string>
#include <map>
#include <set>
#include <list>
#include <tuple>
#include "main.h" // include header for main

typedef map< char, set<char> > I_map;

////////////////////////
//  HELPER FUNCTIONS  //
////////////////////////

// Convert a string to any type, if possible
template <typename T>
T ConvertString( const string &data, int* return_code)
{
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

// Request an integer and perform error checks on input, default value required
void RequestInt(string request_message, int* input, int min_value, int max_value, int default_value) {
    int return_code = -1;
    string user_input;

    do {
        cout << request_message;
        cout <<  " [" << default_value << "]" << ": ";

        cin.clear();
        getline(cin, user_input);
        
        *input = ConvertString<int>(user_input, &return_code);
        if ( return_code == -1 ) {
            cout << "Error: Invalid integer value '" << user_input << "'" << endl;
        }
        if (return_code == 1 ) {
            *input = default_value;
        }
        if (return_code == 0 ) {
            if (*input < min_value) {
                cout << "Warning: integer value out of range " << *input << " less than (reset to) " << min_value << endl;
                *input = min_value;
            }
            if (*input > max_value) {
                cout << "Warning: integer value out of range " << *input << " less than (reset to) " << max_value << endl;
                *input = max_value;
            }
        }
    }
    while ( return_code == -1 );
}

// Request a string and perform error checks on input, default value optional
void RequestString(string request_message, string* input, string default_value = "") {
    int return_code = -1;
    string user_input;

    cout << request_message;

    if (default_value.compare("")) {
        cout <<  " [" << default_value << "]" << ": ";
    }
    else {
        cout << ": ";
    }

    cin.clear();
    getline(cin, user_input);

    if (user_input.empty()) {
        *input = default_value;
    }
    else {
        *input = user_input;
    }
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
unsigned int LCE(INT i, INT j, INT n, INT * invSA, INT * LCP, rmq_succinct_sct<> rmq) {
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

        return LCP[rmq(a + 1, b)];
}

// Calculates a list of longest common extensions, corresponding to 0, 1, 2, ... allowed mismatches
// Only considers "real" mismatches (degenerate string mismatching)
// Takes into account the matching possibility of non A, C, G, T/U characters
// Longest Common Extension calculated from position i and j (i, j order doesn't matter)
// Requires Text, Text Length, Inverse Suffix Array,
// Longest Common Prefix Array (LCP), RMQ function over LCP, Match Matrix, mismatches (max allowed) and initial_gap
// Only starts counting mismatches that occur after initial gap, however earlier mismatches are still stored
// Should only be used after MatchMatrix has been instantiated with necessary data
void realLCE_mismatches(unsigned char* text, INT i, INT j, INT n, INT * invSA, INT * LCP, rmq_succinct_sct<> rmq, int mismatches, int initial_gap, list<double>* mismatch_locs) {
    if ( i == j ) {
        mismatch_locs->push_back( n - i );
    }
    else {
        int real_lce = 0;

        while (mismatches >= 0) {
            real_lce = real_lce + LCE(i + real_lce, j + real_lce, n, invSA, LCP, rmq);

            if ( i + real_lce >= (n / 2) or j + real_lce >= n) {
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

void addPalindromes(set<tuple<int, int, int>>* palindromes, unsigned char* S, int S_n, int n, INT * invSA, INT * LCP, rmq_succinct_sct<> rmq, int mismatches, int max_gap) {
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
        realLCE_mismatches(S, i, j, S_n, invSA, LCP, rmq, mismatches, initial_gap, &mismatch_locs);
        mismatch_locs.push_front(-1.0);

        // Optional printing of mismatch locations relative to centre
        if (false) {
            cout << "[ ";
            for (list<double>::iterator it = mismatch_locs.begin(); it != mismatch_locs.end(); ++it){
                cout << *it << " ";
            }
            cout << "]" << endl;
        }

        list<double>::iterator it = mismatch_locs.begin();
        list<double>::iterator it_offset;

        if (mismatch_locs.size() > mismatches + 1) {
            it_offset = next(mismatch_locs.begin(), mismatches + 1);
        } else {
            it_offset = next(mismatch_locs.begin(), mismatch_locs.size() - 1);
        }

        int left, right, gap;
        double start_mismatch, end_mismatch;

        for ( ; it_offset != mismatch_locs.end(); ++it, ++it_offset){

            start_mismatch = *it;
            end_mismatch = *it_offset;

            if (*it_offset == *prev(it_offset) + 1.0) {
                if (it == mismatch_locs.begin()) {
                    list<double>::iterator temp_it_offset = it_offset;

                    while (temp_it_offset != mismatch_locs.begin() and *temp_it_offset == *prev(temp_it_offset) + 1.0) {
                        temp_it_offset = prev(temp_it_offset);
                    }

                    if (temp_it_offset == mismatch_locs.begin()) {
                        continue;
                    } else {
                        end_mismatch = *temp_it_offset;
                    }
                } else {
                    continue;
                }
            }

             // MAKE THIS WORK IN A MORE EFFICIENT WAY
            /*
            if (*it == *next(it) - 1.0) {
                list<double>::iterator temp_it = it;

                while (temp_it != mismatch_locs.end() and *temp_it == *next(temp_it) - 1.0) {
                    temp_it = next(temp_it);
                }

                start_mismatch = *temp_it;
            }
            */

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

            palindromes->insert(tuple<int, int, int>(left, right, gap));
        }
    }
}

////////////////////
//  MAIN PROGRAM  //
////////////////////

int main() {

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

        ///////////////////////
        //  REQUEST OPTIONS  //
        ///////////////////////

        int min_len, max_len, max_gap, mismatches;
        string input_file, output_file;

        cout << "Find inverted repeats in IUPAC nucleotide sequence(s)" << endl;

        RequestString("Input nucleotide sequence(s)", &input_file, "input.txt");

        ifstream in(input_file);
        string contents((std::istreambuf_iterator<char>(in)), 
        istreambuf_iterator<char>());
        long int n = contents.length();
        unsigned char * seq = ( unsigned char* ) malloc( ( n ) * sizeof( unsigned char ) );

        for (int i = 0; i < n; ++i) {
            seq[i] = contents[i];
        }

        RequestInt("Enter minimum length of palindrome", &min_len, 1, 1000000, 10);
        RequestInt("Enter maximum length of palindrome", &min_len, 1, 1000000, 100);
        RequestInt("Enter maximum gap between repeated regions", &max_gap, 0, 1000000, 100);
        RequestInt("Number of mismatches allowed", &mismatches, 0, 1000000, 0);
        RequestString("Input nucleotide sequence(s)", &output_file, "output.txt");

        // Display user given options
        if (false) {
            cout << "input_file: " << input_file << endl;
            cout << "sequence: " << seq << endl;
            cout << "min_len: " << min_len << endl;
            cout << "min_len: " << min_len << endl;
            cout << "max_gap: " << max_gap << endl;
            cout << "mismatches: " << mismatches << endl;
            cout << "output_file: " << output_file << endl;
        }

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

        int_vector<> v(S_n , 0); // Create a vector of length n and initialize it with 0s
        for ( INT i = 0; i < S_n; i ++ )
        {
                v[i] = LCP[i];
        }

        rmq_succinct_sct<> rmq(&v);

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
        addPalindromes(&palindromes, S, S_n, n, invSA, LCP, rmq, mismatches, max_gap); 

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

        for (set<tuple<int, int, int>>::iterator it = palindromes.begin(); it != palindromes.end(); it++) {
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

            file << endl;

            file << pad;
            for (int i = 0; i < (inner_left - outer_left + 1); ++i) {
                file << ( (MatchMatrix::match(seq[ outer_left - 1 + i ], complement[ seq[ outer_right - 1 - i ] ])) ? "|" : " " );
            }
            file << pad;
            
            file << endl;

            file << outer_right;
            for (int i = 0; i < pad_length - getDigitCount(outer_right); ++i) { file << " "; }
            for (int i = outer_right; i >= inner_right; --i) { file << seq[i - 1]; }
            for (int i = 0; i < pad_length - getDigitCount(inner_right); ++i) { file << " "; }
            file << inner_right;

            file << endl << endl;
        }

        file << endl << endl << endl;

        file.close();

        ///////////////////
        //  Free Memory  //
        ///////////////////

        free(match_matrix);
        free(seq);
        free(SA);
        free(invSA);
        free(LCP);

        return 0;
}