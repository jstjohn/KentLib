/**
 * helloWorld
 *  Author: John St. John
 *  Date: 2/8/2012
 *
 *  This program is an example shell of a C program
 *  for working with the kent source library.
 *
 */



/**
 * Includes
 *
 */

#include "common.h"
#include "options.h"
#include "linefile.h"




/**
 * Global Variables and definitions
 *
 */ 








/**
 * Command line options
 */ 

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "helloWorld -- Writes a greeting to the screen\n"
      "\tUsage: qseqToFastq [options] \n"
      "\noptions:\n"
      "\t-hello\tWrites hello world to the screen.\n"
      "\t-goodBye\tWrites good bye world to the screen.\n"
      "\t-help\tWrites this help to the screen.\n"
      );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"goodBye",OPTION_BOOLEAN},
    {"hello",OPTION_BOOLEAN},
    {"help",OPTION_BOOLEAN},
    {NULL, 0}
}; //end options()










/**
 * Program's functions
 *
 */




int main(int argc, char *argv[])
/* Process command line. */
{

  optionInit(&argc, argv, options);
  boolean help = optionExists("help");
  if(help) usage();

  boolean goodbye = optionExists("goodBye");
  boolean hello = optionExists("hello");
  if(!goodbye && !hello) usage();
  if(hello) printf("Hello World!\n");
  if(goodbye) printf("Good Bye World!\n");

  return 0;
} //end main()


