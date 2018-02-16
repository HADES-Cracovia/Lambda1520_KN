#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <getopt.h>
#include <string>

#include <hloop.h>

//#include "fwdet_tests.h"
#include "fwdet_res.h"

int main(int argc, char **argv)
{
	TROOT Analysis("Analysis","compiled analysis macro");
	gStyle->SetOptStat(0);

	// argc is the number of arguments in char* array argv
	// CAUTION: argv[0] contains the progname
	// argc has to be nargs+1


	/* Flag set by ‘--verbose’. */
	static int verbose_flag;
	static int flag_sim;

	int c;
	long int events = -1;
	TString output = "output.root";
	float beam_momentum = 690.;

	while (1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"verbose", no_argument,       &verbose_flag, 1},
			{"brief",   no_argument,       &verbose_flag, 0},
			{"sim", 	no_argument,	   &flag_sim, 1},
			{"exp",		no_argument, 	   &flag_sim, 0},
			/* These options don’t set a flag.
			 *              We distinguish them by their indices. */
//			{"add",     no_argument,       0, 'a'},
//			{"append",  no_argument,       0, 'b'},
//			{"delete",  required_argument, 0, 'd'},
//			{"create",  required_argument, 0, 'c'},
//			{"file",    required_argument, 0, 'f'},
			{"output",	required_argument,	0,	'o'},
			{"events",	required_argument,	0,	'e'},
			{"mom",		required_argument,	0,	'm'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "abc:d:f:e:o:",
				long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
			case 0:
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag != 0)
					break;
				printf ("option %s", long_options[option_index].name);
				if (optarg)
					printf (" with arg %s", optarg);
				printf ("\n");
				break;

			case 'a':
				puts ("option -a\n");
				break;

			case 'b':
				puts ("option -b\n");
				break;

			case 'c':
				printf ("option -c with value `%s'\n", optarg);
				break;

			case 'd':
				printf ("option -d with value `%s'\n", optarg);
				break;

			case 'f':
				printf ("option -f with value `%s'\n", optarg);
				break;

			case 'e':
				events = atol(optarg);
				break;

			case 'o':
				output = optarg;
				break;

			case 'm':
				beam_momentum = atof(optarg);
				break;

			case '?':
				/* getopt_long already printed an error message. */
				break;

			default:
				abort ();
		}
	}

	/* Instead of reporting ‘--verbose’
	 *      and ‘--brief’ as they are encountered,
	 *           we report the final status resulting from them. */
	if (verbose_flag)
		puts ("verbose flag is set");

	HLoop * loop = new HLoop(kTRUE);

	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		Bool_t ret;
// 		printf ("non-option ARGV-elements: ");
		while (optind < argc)
		{
			TString infile = argv[optind++];
			if		(infile.Contains(","))		ret = loop->addMultFiles(infile);
			else if (infile.Contains(".root"))	ret = loop->addFiles(infile);
			else								ret = loop->addFilesList(infile);

			if (!ret)
			{
				std::cerr << "READBACK: ERROR : cannot find inputfiles : " << infile.Data() << endl;
				std::exit(EXIT_FAILURE);
			}
		}
	}

	AnaParameters anapars;
	anapars.outfile = output;
	anapars.events = events;
	anapars.beam_momentum = beam_momentum;

	if (flag_sim==1) {
			anapars.sim = true;}
	if (flag_sim==0){
			anapars.sim = false;}

	fwdet_tests(loop, anapars);

	exit(0);
}


