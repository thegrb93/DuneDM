#pragma once
#include <string>

class TBranch;
class TGraph2D;

class DarkMatterXMomentum
{
	TGraph2D* graph;
	int nfiles;
	int index;
	
	public:
	DarkMatterXMomentum(int nfiles);
	~DarkMatterXMomentum();
	
	void Fill(const std::string& file, TBranch* branch, int nentries);
	void Save();
};

