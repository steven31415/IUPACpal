

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "main.h"

static struct option long_options[] =
 {
  
   { "input-file",              required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "nummismatches ",              required_argument, NULL, 'k' },
   { "minpallen ",              required_argument, NULL, 'm' },
 { "maxpallen  ",              required_argument, NULL, 'x' },
 { "gaplimit ",              required_argument, NULL, 'g' },
   { "help",                    no_argument,       NULL, 'h' },
   { NULL,                      0,                 NULL, 0   }
 };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */

   sw -> input_filename                 = NULL;
   sw -> output_filename                = NULL;
   sw -> k                              = 0;
   sw -> m                              = 10;
  sw -> x                             = 100;
  sw -> g                              = 100;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "i:o:k:m:x:g:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
       

         case 'i':
           sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> input_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
           break;

         case 'k':
           val = strtol ( optarg, &ep, 0);
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> k = val;
           args ++;
           break;
           case 'm':
               val = strtol ( optarg, &ep, 10 );
               if ( optarg == ep )
               {
                   return ( 0 );
               }
               sw -> m = val;
               args ++;
               break;

    case 'x':
               val = strtol ( optarg, &ep, 100 );
               if ( optarg == ep )
               {
                   return ( 0 );
               }
               sw -> m = val;
               args ++;
               break;
    case 'g':
               val = strtol ( optarg, &ep, 100 );
               if ( optarg == ep )
               {
                   return ( 0 );
               }
               sw -> m = val;
               args ++;
               break;


         case 'h':
           return ( 0 );
       }
    }

   if ( args < 6 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }


/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, " Usage: map <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );

   fprintf ( stdout, "  -i, --input-file          <str>     (Multi)FASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file         <str>     Output filename.\n" );
   fprintf ( stdout, "  -k, --nummismatches          <int>     The max Hamming distance k.\n");
   fprintf ( stdout, "  -m, --minpallen          <int>     The len of the substring m.\n");
  fprintf ( stdout, "  -x, --maxpallen         <int>     The len of the substring m.\n");
  fprintf ( stdout, "  -g, --gaplimit          <int>     The len of the substring m.\n");
 }

