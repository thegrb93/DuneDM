// main program
#include <iostream>
#include <vector>
#include <map>
#include <TApplication.h>
#include "cxxopts.hpp"
#include "DMAnalysis.h"

TApplication* gApp;

int main (int argc, char** argv) {
	cxxopts::Options options("DuneDM", "USAGE: DuneDM [options] files");
	options.add_options()
            ("mode", "Which analysis mode to use. Can be 'statistics', 'histograms', or 'detector'.", cxxopts::value<std::string>(), "")
            ("particle", "Particle pdgcode to analyze", cxxopts::value<int>(), "")
            ("attribute", "The particle attribute to analyze", cxxopts::value<std::string>(), "")
            ("detector", "Which detector type to use. Can be 'DUNE'.", cxxopts::value<std::string>(), "")
            ("files","The root files to analyze", cxxopts::value<std::vector<std::string>>(), "");
    options.parse_positional("files");
	options.parse(argc, argv);

    //TApplication app("tapp", &argc, argv);
    //gApp = &app;

    const std::vector<std::string>& files = options["files"].as<std::vector<std::string>>();
	if(files.size()>0) {
		std::string mode;
        try { mode = options["mode"].as<std::string>(); }
        catch(...){ mode = "detector"; }
        if(mode.empty())
            mode = "detector";
			
		std::map<std::string,DMAnalysis*(*)()> modes = {
			{"statistics", &StatisticsAnalysis::create},
			{"histograms", &DarkMatterDistribution::create},
			{"detector", &DetectorAnalysis::create}
		};

		auto constructor = modes.find(mode);
		if(constructor != modes.end()) {
			DMAnalysis* analysis = (constructor->second)();
            analysis->files.assign(files.begin(), files.end());
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

