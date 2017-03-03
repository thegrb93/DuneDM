// main program
#include <iostream>
#include <vector>
#include <map>
#include "DMAnalysis.h"

const option::Descriptor usage[] =
{
	{OPT_UNKNOWN, 0,"" , ""    ,option::Arg::None, "USAGE: DuneDM [options] files\n\n"
		                                     "Options:" },
	{OPT_HELP,    0,"h" , "help", option::Arg::None, "  --help, -h  \tPrint usage and exit." },
	{OPT_MODE,    0,"", "mode", option::Arg::Optional, "  --mode, -p  \tIncrement count." },
	{OPT_PARTICLE, 0,"" ,  "particle"   ,option::Arg::Optional, "Examples:\n"
		                                     "  example --unknown -- --this_is_no_option\n"
		                                     "  example -unk --plus -ppp file1 file2\n" },
	{OPT_PARTICLEATTRIBUTE, 0, "", "attribute", option::Arg::Optional, ""},
	{0,0,0,0,0,0}
};
 
option::Option* gOptions;

int main (int argc, char** argv)
{
	if(argc<2){
		option::printUsage(std::cout, usage);
		return 0;
	}
	
	option::Stats  stats(usage, argc, argv);
	option::Option options[stats.options_max], buffer[stats.buffer_max];
	option::Parser parse(usage, argc, argv, options, buffer);
	gOptions = options;

	if (parse.error())
		return 1;

	if (options[OPT_HELP]) {
		option::printUsage(std::cout, usage);
		return 0;
	}

	for (option::Option* opt = options[OPT_UNKNOWN]; opt; opt = opt->next())
		std::cout << "Unknown option: " << opt->name << "\n";

	if(parse.nonOptionsCount()>0)
	{
		std::string mode;
		if(options[OPT_MODE])
			mode = std::string(options[OPT_MODE].last()->arg);
		else
			mode = "detector";
			
		std::map<std::string,ROOTAnalysis*(*)()> modes = {
			{"statisics", &StatisticsAnalysis::create},
			{"histograms", &DarkMatterDistribution::create},
			{"detector", &DetectorAnalysis::create}
		};

		auto constructor = modes.find(mode);
		if(constructor != modes.end())
		{
			ROOTAnalysis* analysis = (constructor->second)();
			
			for (int i = 0; i < parse.nonOptionsCount(); ++i)
				analysis->files.push_back(parse.nonOption(i));
			
			analysis->Process();
			delete analysis;
		}
		else
			std::cout << "Invalid mode: " << mode << std::endl;
	}
	else
		std::cout << "Expected input root files. Got none.\n";

	return 0;
}

