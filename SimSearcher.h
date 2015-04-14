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

	vector<string> context_pre;
	vector<const char*> context;
	int** context_hash;
	int* context_word_len;
	int* context_len;
	vector<int>** jaccard_list;
	int* candid;
	int* candid_num;
	pair<int,int>* candid_list;
	int* candid_set;
	int* candid_ys;
//	map<string,int> jaccard_hash;
//	vector< vector<int> > jaccard_list;
//	map<string,int> ed_hash;
//	vector< vector<int> > ed_list;	
	vector<int> ql;
//	string** jaccard_hash;
	int* jaccard_hash;
	int	min_context_len;
	string sep;
	int q_gram;
	vector<int>** ed_list;
	int* ed_hash;
	void clear();
	void convert(const char*, vector<int>&, int );
	int createIndex(const char *filename, unsigned q);
	int searchJaccard(const char *query, double threshold, std::vector<std::pair<unsigned, double> > &result);
	int searchED(const char *query, unsigned threshold, std::vector<std::pair<unsigned, unsigned> > &result);
};

