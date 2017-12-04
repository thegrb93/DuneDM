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
            ("mode", "Which analysis mode to use. Can be 'statistics', 'sensitivity', or 'detector'.", cxxopts::value<std::string>(), "")
            ("particle", "Particle pdgcode to analyze", cxxopts::value<int>(), "")
            ("attribute", "The particle attribute to analyze", cxxopts::value<std::string>(), "")
            ("detector", "Which detector type to use. Can be 'DUNE'.", cxxopts::value<std::string>(), "")
            ("files","The folder containing root files to analyze", cxxopts::value<std::vector<std::string>>(), "");
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
		
		auto statistics = [](const std::vector<std::string>& folders){
            StatisticsAnalysis analysis;
            return analysis.Process(folders[0]);
		};
		auto detector = [](const std::vector<std::string>& folders){
            DetectorAnalysis analysis;
            return analysis.Process(folders[0]);
		};
		auto sensitivity = [](const std::vector<std::string>& folders){
            SensitivityScan analysis;
            return analysis.Process(folders);
		};
		
		std::map<std::string,int(*)(const std::vector<std::string>&)> modes = {
			{"statistics", statistics},
			{"detector", detector},
			{"sensitivity", sensitivity}
		};

		auto func = modes.find(mode);
		if(func != modes.end())
		    func->second(files);
		else
			std::cout << "Invalid mode: " << mode << std::endl;
	}
	else
		std::cout << "Missing input directory.\n";

	return 0;
}

