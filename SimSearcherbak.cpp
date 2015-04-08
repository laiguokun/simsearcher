#include "SimSearcher.h"
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstring>
const int mode = 999997;
const int hash_size = 1000000;
const int max_int = 1000000;
using namespace std;

SimSearcher::SimSearcher()
{
	sep = " ";
}

SimSearcher::~SimSearcher()
{
}

inline int myabs(int a)
{
	if (a > 0) return a; else return -a;
}

inline int getmax(int a, int b)
{
	if (a > b) return a; else return b;
}

inline int getmin(int a, int b)
{
	if (a < b) return a; else return b;
}
void split(string& s, string& delim, vector<string>& ret)  
{  
	ret.clear();
    size_t last = 0;  
    size_t index=s.find_first_of(delim,last);  
    while (index!=string::npos)  
    {  
        ret.push_back(s.substr(last,index-last));  
        last=index+1;  
        index=s.find_first_of(delim,last);  
    }  
    if (index-last>0)  
    {  
        ret.push_back(s.substr(last,index-last));  
    }  
}  

int myhash(string& str, string** dataset)
{
	int res = 0;
	for (int i = 0; i < str.length(); i++)
		res = (res * 29 + str[i]) % mode;
	while (dataset[res] != NULL && str.compare(*dataset[res]) != 0)
	{
		res = (res + 1) % mode;
	}
	if (dataset[res] == NULL)
		dataset[res] = new string(str);
	return res;
}

int search(string& str, string** dataset)
{
	int res = 0;
	for (int i = 0; i < str.length(); i++)
		res = (res * 29 + str[i]) % mode;
	while (dataset[res] != NULL && str.compare(*dataset[res]) != 0)
	{
		res = (res + 1) % mode;
	}
	if (dataset[res] == NULL)
		return -1;
	return res;
}

double verify_jaccard(int* query, int s1, int* dataset, int s2, double threshold)
{
	int h1 = 0, h2 = 0;
	int cnt = 0;
	while (h1 < s1 && h2 < s2)
	{
		if (query[h1] == dataset[h2])
		{
			h1 ++;
			h2 ++;
			cnt ++;
			continue;
		}
		if (query[h1] > dataset[h2])
			h2 ++;
		else
			h1 ++;
	}
	return ((double)cnt / (double)(s1+s2-cnt));
}
 
void SimSearcher::clear()
{

}

int SimSearcher::createIndex(const char *filename, unsigned q)
{
	ifstream fin(filename);
	min_context_len = max_int;
	string str;
	context.clear();
	while (getline(fin,str))
		context.push_back(str);
	//build jaccard index
	//build invert list of jaccard
	//clear
	clear();
	jaccard_hash = new string*[hash_size];
	jaccard_list = new vector<int>*[hash_size];
	context_hash = new int*[context.size()];
	context_word_len = new int[context.size()];
	for (int i = 0; i < hash_size; i++)
	{
		jaccard_hash[i] = NULL;
		jaccard_list[i] = NULL;
	}
	vector<string> word;
	for (int i = 0; i < context.size(); i++)
	{
		word.clear();
		split(context[i], sep, word);
		context_hash[i] = new int[word.size()];
		for (int j = 0; j < word.size(); j++)
		{
			int mark = myhash(word[j], jaccard_hash);
			if (jaccard_list[mark] == NULL)
			{
				jaccard_list[mark] = new vector<int>;
				jaccard_list[mark]->clear();
			}
			jaccard_list[mark]->push_back(i);
			context_hash[i][j] = mark;
			context_word_len[i] = word.size();
			if (min_context_len > context_word_len[i])
				min_context_len = context_word_len[i];
		}
		sort(context_hash[i],context_hash[i]+word.size());
//		cout << i << endl;
	}
	
	//build ed index---qgram
	//clear
	ql.clear();
	q_gram = q;
	ed_hash = new string*[hash_size];
	ed_list = new vector<int>*[hash_size];
	for (int i = 0; i < hash_size; i++)
	{
		ed_hash[i] = NULL;
		ed_list[i] = NULL;
	}
	for (int i = 0; i < context.size(); i++)
	{
		if (context[i].length() >= q)
			for (int j = 0; j < context[i].length() - q + 1; j++)
			{
				string str = context[i].substr(j, q);
				int mark = myhash(str, ed_hash);
				if (ed_list[mark] == NULL)
				{
					ed_list[mark] = new vector<int>;
					ed_list[mark]->clear();
				}
				ed_list[mark]->push_back(i);
			}
		else
		{
			ql.push_back(i);
		}
	}
	
	return SUCCESS;
}

int SimSearcher::searchJaccard(const char *query, double threshold, vector<pair<unsigned, double> > &result)
{
	result.clear();
	string str(query);
	//find candidate
	vector<string> word;
	vector< pair<int, int> > candid_list;
	split(str, sep, word);
	int word_num = word.size();
	int* word_list = new int[word_num];
	int tmp = 0;
	//calc threshold
	double t1 = (double) ((word_num + min_context_len) * threshold) / (1.0 + threshold);
	int thres = getmax(ceil(threshold * word_num), ceil(t1));
	vector<int> candid;
	candid.clear();
	int* str_hash = new int[word_num];
	for (int i = 0; i < word_num; i++)
	{
		word_list[i] = search(word[i], jaccard_hash);
		str_hash[i] = word_list[i];
	}
	sort(str_hash, str_hash + word_num);
	int* candid_set = new int[context.size()];
	memset(candid_set, 0, sizeof(int) * context.size());
	for (int i = 0; i < word_num; i++)
	{
		if (word_list[i] != -1)
			for (int j = 0; j < jaccard_list[word_list[i]]->size(); j++)
			{
				candid_set[(*jaccard_list[word_list[i]])[j]] ++;
			}
	}
	for (int i = 0; i < context.size(); i++)
	{
		if (candid_set[i] >= thres)
			candid.push_back(i);
	}
	for (int i = 0; i < candid.size(); i++)
	{
		double jaccard = (verify_jaccard(str_hash, word_num, context_hash[candid[i]], context_word_len[candid[i]], threshold));
		if (jaccard >= threshold)
			result.push_back(make_pair(candid[i], jaccard));
	}
	delete[] str_hash;
	delete[] word_list;
	delete[] candid_set;
	return SUCCESS;
}

inline bool okf(int x, int y, int threshold)
{
	if (x >=0 && y >= 0 && myabs(y - x) <= threshold)
		return true;
	else
		return false;
}

int verify_ed(string s1, string s2, int threshold)
{
//	cout << s2.length() << " " << s1.length() << endl;
//	cout << myabs(s2.length() - s1.length()) << endl;
//	return 0;
	if (myabs(s2.length() - s1.length()) > threshold) return max_int;
	if (s2.length() < s1.length())
	{
		string tmp = s1;
		s1 = s2;
		s2 = tmp;
	}
	int** dp = new int*[2];
	dp[0] = new int[s2.length()];
	dp[1] = new int[s2.length()]; 
	for (int i = 0; i < s2.length(); i++)
	{
		if (s1[0] == s2[i]) dp[0][i] = i; else dp[0][i] = i + 1;
		if (i > 0) dp[0][i] = getmin(dp[0][i], dp[0][i-1] + 1);
	}
	int now = 0;
//	cout << s1 << " " << s2  <<endl;
	if (s1.length() > 1)
		for (int i = 1; i < s1.length(); i++)
		{
			now = now ^ 1;
			int l,r;
			l = getmax(i - threshold, 0);
			r = getmin(i + threshold + 1, s2.length());
			bool ok = true;
			for (int j = l; j < r; j++)
			{
				dp[now][j] = max_int;
				if (okf(i-1, j-1, threshold)) 
				{
					if (s1[i] == s2[j]) dp[now][j] = getmin(dp[now][j], dp[now ^ 1][j-1]);
						else dp[now][j] = getmin(dp[now][j], dp[now ^ 1][j-1] + 1);
				}
				if (okf(i-1, j, threshold))
					dp[now][j] = getmin(dp[now][j], dp[now ^ 1][j] + 1);
				if (okf(i, j-1, threshold))
					dp[now][j] = getmin(dp[now][j], dp[now][j-1] + 1);
				if (dp[now][j] <= threshold) ok = false;
//				cout << i << " " << j << " " << dp[now][j] << endl;
			}	
			if (ok) 
			{
				for (int j = 0; j < 2; j++)
					delete[] dp[j];
				delete[] dp;
				return max_int;
			}
		}
	int result = dp[now][s2.length() - 1];
	for (int j = 0; j < 2; j++)
		delete[] dp[j];
	delete[] dp;
	return result;
}


int SimSearcher::searchED(const char *query, unsigned threshold, vector<pair<unsigned, unsigned> > &result)
{
	result.clear(); 
	string str(query);
	//calc threshold
	int thres = getmax(str.length() - q_gram + 1 - threshold * q_gram, 0);
//	thres = 0;
//	cout << thres << endl;
	int tmp;
	vector<int> candid;
	candid.clear();
	if (thres == 0)
		for (int i = 0; i < context.size(); i++)
			candid.push_back(i);
	else
	{
		for (int i = 0; i < ql.size(); i++)
			candid.push_back(ql[i]);
		int word_num = str.length() - q_gram + 1;
		int* word_list = new int[word_num];
		for (int i = 0; i < word_num; i++)
		{
			string word = str.substr(i, q_gram);
			word_list[i] = search(word, ed_hash);
		}
		int* candid_set = new int[context.size()];
		memset(candid_set, 0, sizeof(int) * context.size());
		for (int i = 0; i < word_num; i++)
		{
			if (word_list[i] != -1)
				for (int j = 0; j < ed_list[word_list[i]]->size(); j++)
					candid_set[(*ed_list[word_list[i]])[j]] ++;
		}
		for (int i = 0; i < context.size(); i++)
		{
			if (candid_set[i] >= thres)
			{
				candid.push_back(i);
			}
		}
		delete[] word_list;
		delete[] candid_set;
	}
//	sort(candid.begin(),candid.end());
	for (int i = 0; i < candid.size(); i++)
	{
		int num = verify_ed(str, context[candid[i]], threshold);
//		cout << candid[i] << " " << num <<endl;
		if (num <= threshold)
			result.push_back(make_pair(candid[i], num));
	}
	return SUCCESS;
}

