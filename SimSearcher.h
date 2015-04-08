#pragma once
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <map>
using namespace std;
const int SUCCESS = 0;
const int FAILURE = 1;

class SimSearcher
{
public:
	SimSearcher();
	~SimSearcher();

	vector<string> context;
	int** context_hash;
	int* context_word_len;
	vector<int>** jaccard_list;
//	map<string,int> jaccard_hash;
//	vector< vector<int> > jaccard_list;
//	map<string,int> ed_hash;
//	vector< vector<int> > ed_list;	
	vector<int> ql;
//	string** jaccard_hash;
	long long* jaccard_hash;
	int	min_context_len;
	string sep;
	int q_gram;
	vector<int>** ed_list;
	long long* ed_hash;
	void clear();
	int createIndex(const char *filename, unsigned q);
	int searchJaccard(const char *query, double threshold, std::vector<std::pair<unsigned, double> > &result);
	int searchED(const char *query, unsigned threshold, std::vector<std::pair<unsigned, unsigned> > &result);
};

