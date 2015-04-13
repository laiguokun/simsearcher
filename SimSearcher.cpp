#include "SimSearcher.h"
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstring>
const int mode = 1000003;
const int mode2 = 10000003;
const int hash_size = 1000005;
const int max_int = 1000000;
const double u = 0.0085;
const double u1 = 0.0085;
double M = 0;
double M1 = 0;
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
	for (int i = 0; i < s.length(); i++)
		if (s[i] == ' ') 
		{
			ret.push_back(s.substr(last, (i - last)));
			last = i + 1;
		}
	ret.push_back(s.substr(last, (s.length() - last)));
}  

int myhash(string& str, int* dataset)
{
	int res = 0;
	int res2 = 0;
	for (int i = 0; i < str.length(); i++)
	{
		res = (res * 29 + str[i]);
		if (res >= mode) res = res % mode;
		res2 = (res2 * 39 + str[i]);
		if (res2 >= mode2) res2 = res % mode2;
	}
	res2 ++;
	while (dataset[res] != 0 && dataset[res] != res2)
	{
		res = (res + 1);
		if (res >= mode) res -= mode;
	}
	if (dataset[res] == 0)
		dataset[res] = res2;
	return res;
}

int search(string& str, int* dataset)
{
	int res = 0;
	int res2 = 0;
	for (int i = 0; i < str.length(); i++)
	{
		res = (res * 29 + str[i]);
		if (res >= mode) res = res % mode;
		res2 = (res2 * 39 + str[i]);
		if (res2 >= mode2) res2 = res % mode2;
	}
	res2 ++;
	while (dataset[res] != 0 && dataset[res] != res2)
	{
		res = (res + 1);
		if (res >= mode) res -= mode;
	}
	if (dataset[res] == 0)
		return -1;
	return res;
}

double verify_jaccard(int* query, int& s1, int* dataset, int& s2, double& threshold)
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
vector<string> w;
vector<int> tmp;	
void SimSearcher::convert(string &str, string &sep, vector<int>& word)
{
	w.clear();
	tmp.clear();
	split(str, sep, w);
	for (int i = 0; i < w.size(); i++)
		tmp.push_back(myhash(w[i], jaccard_hash));
	sort(tmp.begin(), tmp.end());
	word.push_back(tmp[0]);
	for (int i = 1; i < tmp.size(); i++)
		if (tmp[i] != tmp[i-1])
			word.push_back(tmp[i]);
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
	jaccard_hash = new int[hash_size];
	jaccard_list = new vector<int>*[hash_size];
	context_hash = new int*[context.size()];
	context_word_len = new int[context.size()];
	memset(jaccard_hash, 0, sizeof(int) * hash_size);
	for (int i = 0; i < hash_size; i++)
	{
		jaccard_list[i] = NULL;
	}
	vector<int> word;
	for (int i = 0; i < context.size(); i++)
	{
		word.clear();
		convert(context[i], sep, word); 
		context_hash[i] = new int[word.size()];
		for (int j = 0; j < word.size(); j++)
		{
			int mark = word[j];
			if (jaccard_list[mark] == NULL)
			{
				jaccard_list[mark] = new vector<int>;
				jaccard_list[mark]->clear();
			}
			jaccard_list[mark]->push_back(i);
			if (jaccard_list[mark]->size() > M1) M1 = jaccard_list[mark]->size();
			context_hash[i][j] = mark;
		}
		context_word_len[i] = word.size();
		if (min_context_len > context_word_len[i])
			min_context_len = context_word_len[i];
		sort(context_hash[i],context_hash[i]+word.size());
//		cout << i << endl;
	}
	//build ed index---qgram
	//clear
	ql.clear();
	q_gram = q;
	ed_hash = new int[hash_size];
	ed_list = new vector<int>*[hash_size];
	memset(ed_hash, 0, sizeof(int) * hash_size);
	for (int i = 0; i < hash_size; i++)
	{
		ed_list[i] = NULL;
	}
	for (int i = 0; i < context.size(); i++)
	{
		int r = context[i].length() - q + 1;
		for (int j = 0; j < r; j++)
		{
			string str = context[i].substr(j, q);
			int mark = myhash(str, ed_hash);
			if (ed_list[mark] == NULL)
			{
				ed_list[mark] = new vector<int>;
				ed_list[mark]->clear();
			}
			ed_list[mark]->push_back(i);
			if (ed_list[mark]->size() > M) M = ed_list[mark]->size();
		}
	}
//	cout << "fuck" <<endl;		
	return SUCCESS;
}

inline int cmp(pair<int, int> a, pair<int, int> b)
{
	return a.first < b.first;
}

vector<int> candid_num;
vector<string> word;
vector<int> candid;
vector< pair<int, int> > candid_list;
int SimSearcher::searchJaccard(const char *query, double threshold, vector<pair<unsigned, double> > &result)
{
	result.clear();
	string str(query);
	//find candidate
	candid_list.clear();
	vector<int> word;
	word.clear();
	convert(str, sep, word); 
	int word_num = word.size();
	int* word_list = new int[word_num];
	int tmp = 0;
	//calc threshold
	double t1 = (double) ((word_num + min_context_len) * threshold) / (1.0 + threshold);
	int thres;
	if (ceil(threshold * word_num) >= ceil(t1)) thres = ceil(threshold * word_num); else thres = ceil(t1);
	candid.clear();
	int* str_hash = new int[word_num];
	for (int i = 0; i < word_num; i++)
	{
		word_list[i] = word[i];
		if (word_list[i] != -1)
			candid_list.push_back(make_pair(jaccard_list[word_list[i]]->size(), word_list[i]));
	}
	int* candid_set = new int[context.size()];
	memset(candid_set, 0, sizeof(int) * context.size());
//	int f = thres - u1;
int f = (double)thres/(u1*log(M1) + 1.0);
	tmp = (word_num - (f));
	if (tmp <= 0 || thres == 0 || thres <= u)
	{
		for (int i = 0; i < word_num; i++)
		{
			if (word_list[i] != -1)
				for (int j = 0; j < jaccard_list[word_list[i]]->size(); j++)
				{
					candid_set[(*jaccard_list[word_list[i]])[j]] ++;
				}
		}
		for (int i = 0; i < context.size(); i++)
			if (candid_set[i] >= thres)
				candid.push_back(i);
	}
	else
	{
		sort(candid_list.begin(), candid_list.end(), cmp);
		candid_num.clear();
		for (int i = 0; i < tmp; i++)
			for (int j = 0; j < candid_list[i].first; j++)
			{
				if (candid_set[(*jaccard_list[candid_list[i].second])[j]] == 0)
					candid_num.push_back((*jaccard_list[candid_list[i].second])[j]);
				candid_set[(*jaccard_list[candid_list[i].second])[j]] ++;
			}
		int cnt;
		for (int i = 0; i < candid_num.size(); i++)
		{
			cnt = candid_set[candid_num[i]];
			if (cnt < (thres - f)) continue;
			for (int j = tmp; j < candid_list.size(); j++)
			{
				if (cnt >= thres)
					break;
				int h = 0, t = candid_list[j].first, mid;
				while (h < t - 1)
				{
					mid = (h + t) >> 1;
					if ((*jaccard_list[candid_list[j].second])[mid] <= candid_num[i]) h = mid; else t = mid;
				}
				if ((*jaccard_list[candid_list[j].second])[h] == candid_num[i]) cnt ++;
			}
			if (cnt >= thres)
				candid.push_back(candid_num[i]);
		}
	}
	sort(candid.begin(),candid.end());
	for (int i = 0; i < candid.size(); i++)
	{
		double jaccard = (verify_jaccard(word_list, word_num, context_hash[candid[i]], context_word_len[candid[i]], threshold));
		if (jaccard >= threshold)
			result.push_back(make_pair(candid[i], jaccard));
	}
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
	if (myabs(s2.length() - s1.length()) > threshold) return max_int;
	if (s2.length() < s1.length())
	{
		string tmp = s1;
		s1 = s2;
		s2 = tmp;
	}
	int** dp = new int*[2];
	int len1 = s1.length(), len2 = s2.length();
	dp[0] = new int[len2];
	dp[1] = new int[len2]; 
	for (int i = 0; i < len2; i++)
	{
		if (s1[0] == s2[i]) dp[0][i] = i; else dp[0][i] = i + 1;
		if (i > 0 && dp[0][i-1] + 1 < dp[0][i]) dp[0][i] = dp[0][i-1] + 1;
	}
	int now = 0;
//	cout << s1 << " " << s2  <<endl;
	if (s1.length() > 1)
		for (int i = 1; i < len1; i++)
		{
			now = now ^ 1;
			int l = 0,r = len2;
			if (i - threshold > l) l = i - threshold;
			if (i + threshold + 1 < r) r = i + threshold + 1;
			bool ok = true;
			for (int j = l; j < r; j++)
			{
				dp[now][j] = max_int;
				if (i > 0 && j > 0) //				if (okf(i-1, j-1, threshold)) 
				{
					if (s1[i] == s2[j])
					{
						if (dp[now ^ 1][j - 1] < dp [now][j]) dp[now][j] = dp[now ^ 1][j-1];
					}
					else 
					{
						if (dp[now ^ 1][j - 1] + 1 < dp[now][j]) dp[now][j] = dp[now ^ 1][j-1] + 1;
					}
				}
				if (i > 0 && myabs(j - i + 1) <= threshold && dp[now ^ 1][j] + 1 < dp[now][j]) //(okf(i-1, j, threshold))
					dp[now][j] = dp[now ^ 1][j] + 1;
				if (j > 0 && myabs(j - 1 - i) <= threshold && dp[now][j-1] + 1 < dp[now][j]) // (okf(i, j-1, threshold))
					dp[now][j] = dp[now][j-1] + 1;
				if (ok && dp[now][j] <= threshold) ok = false;
			}	
			if (ok) 
			{
				for (int j = 0; j < 2; j++)
					delete[] dp[j];
				delete[] dp;
				return max_int;
			}
		}
	int result = dp[now][len2 - 1];
	for (int j = 0; j < 2; j++)
		delete[] dp[j];
	delete[] dp;
	return result;
}


int SimSearcher::searchED(const char *query, unsigned threshold, vector<pair<unsigned, unsigned> > &result)
{
	result.clear(); 
	string str(query);
	int len = str.length();
	//calc threshold
	int thres = getmax(len - q_gram + 1 - threshold * q_gram, 0);
	candid.clear();
	candid_list.clear();
	if (thres == 0)
		for (int i = 0; i < context.size(); i++)
			candid.push_back(i);
	else
	{
		int word_num = len - q_gram + 1;
		int* word_list = new int[word_num];
		for (int i = 0; i < word_num; i++)
		{
			string word = str.substr(i, q_gram);
			word_list[i] = search(word, ed_hash);
			if (word_list[i] != -1)
				candid_list.push_back(make_pair(ed_list[word_list[i]]->size(), word_list[i]));
		}
		int* candid_set = new int[context.size()];
		memset(candid_set, 0, sizeof(int) * context.size());
		int f = (double)thres/(u*log(M) + 1.0);
		int tmp = (word_num - f);
		if (tmp <= 0 || thres == 0 || thres <= u)
		{
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
		}
		else
		{
			sort(candid_list.begin(), candid_list.end(), cmp);
			candid_num.clear();
			for (int i = 0; i < tmp; i++)
			{
				for (int j = 0; j < candid_list[i].first; j++)
				{
					if (candid_set[(*ed_list[candid_list[i].second])[j]] == 0)
						candid_num.push_back((*ed_list[candid_list[i].second])[j]);
					candid_set[(*ed_list[candid_list[i].second])[j]] ++;
				}
			}
			int cnt;
			for (int i = 0; i < candid_num.size(); i++)
			{
				cnt = candid_set[candid_num[i]];
				if (cnt < thres - f) continue;
				if (myabs(context[candid_num[i]].length() - str.length()) > threshold) continue;
				for (int j = tmp; j < candid_list.size(); j++)
				{
					if (cnt >= thres) break;
					int h = 0, t = candid_list[j].first, mid;
					while (h < t - 1)
					{
						mid = (h + t) >> 1;
						if ((*ed_list[candid_list[j].second])[mid] <= candid_num[i]) h = mid; else t = mid;
					}
					if ((*ed_list[candid_list[j].second])[h] == candid_num[i]) cnt ++;
				}	
				if (cnt >= thres)
					candid.push_back(candid_num[i]);
			}
		}
		delete[] word_list;
		delete[] candid_set;
	}
	sort(candid.begin(),candid.end());
	for (int i = 0; i < candid.size(); i++)
	{
		int num = verify_ed(str, context[candid[i]], threshold);
		if (num <= threshold)
			result.push_back(make_pair(candid[i], num));
	}
	return SUCCESS;
}

