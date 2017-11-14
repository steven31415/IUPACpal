

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "main.h"
using namespace std;
int main(int argc, char **argv)
{

	struct TSwitch  sw;

	FILE *          in_fd;                  // the input file descriptor
        char *          input_filename;         // the input file name
        char *          output_filename;        // the output file name
        unsigned char * seq    = NULL;         	// the sequence in memory

        unsigned char * seq_id = NULL;         	// the sequence id in memory
        unsigned int    num_seqs = 0;           // the total number of sequences considered
	char *          alphabet;               // the alphabet
	unsigned int    i, j;
    	unsigned int    k, K, r;
	double   * C ;
    	INT m, M,G,X;

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	/* Check the arguments */
        if ( i < 6 )
        {
                usage ();
                return ( 1 );
        }
        else
        {
            

		K      = sw . k;
	        M      = sw . m; 
                G      = sw.g;
                X      =sw.x;

                input_filename          = sw . input_filename;
                output_filename         = sw . output_filename;

        }

	//double start = gettime();

	/* Read the (Multi)FASTA file in memory */
	if ( ! ( in_fd = fopen ( input_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
		return ( 1 );
	}

	char c;
	c = fgetc( in_fd );
	do
	{
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
			return ( 1 );
		}
		else
		{
			unsigned int max_alloc_seq_id = 0;
			unsigned int seq_id_len = 0;
			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id )
				{
					seq_id = ( unsigned char * ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id += ALLOC_SIZE;
				}
				seq_id[ seq_id_len++ ] = c;
			}
			seq_id[ seq_id_len ] = '\0';
			
		}

		unsigned int max_alloc_seq = 0;
		INT seq_len = 0;

		while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
		{
			if( seq_len == 0 && c == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file %s!\n", input_filename );
				c = fgetc( in_fd );
				break;
			}
			if( c == '\n' ) continue;

			c = toupper( c );

			if ( seq_len >= max_alloc_seq )
			{
				seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) );
         
				max_alloc_seq += ALLOC_SIZE;
			}

			if( strchr ( alphabet, c ) )
			{
				seq[ seq_len++ ] = c;
			}
			else
			{
				if ( strchr ( IUPAC, c ) )
				{
					seq[ seq_len++ ] = 'N';
				}
				else
				{
					fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", input_filename, c );
					return ( 1 );
				}
			}

        }
        
        
        
		if( seq_len != 0 )
		{

			num_seqs++;
			if ( seq_len >= max_alloc_seq )
			{
				seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) );
              
				max_alloc_seq += ALLOC_SIZE;
			}
			seq[ seq_len ] = '\0';

			C = ( double   *  ) calloc ( ( seq_len ) , sizeof(  double ) );
   


            
		

		}
		
		free ( C );
        	C = NULL;
		free ( seq );
		seq = NULL;
		free ( seq_id );
		seq_id = NULL;

	} while( c != EOF );

	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	//double end = gettime();

      //fprintf( stderr, "Elapsed time for processing %d sequence(s): %lf secs\n", num_seqs, ( end - start ) );
	
        free ( sw . input_filename );
        free ( sw . output_filename );
       /* free ( sw . g );
        free ( sw . m );

        free ( sw . x );
        free ( sw . k );*/


	return ( 0 );
}
