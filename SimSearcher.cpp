#include "SimSearcher.h"
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstring>
const int mode = 1000003;
const int mode2 = 10000003;
const int hash_size = 1000005;
const int max_int = 0xffff;
const double u = 0.0085;
const double u1 = 0.0085;
double M = 0;
double M1 = 0;
using namespace std;

void k_big(pair<int,int>* arr,int low,int high, int k)
{
	int tmp = low + rand()%(high - low+1);
	pair <int, int> pivot = arr[tmp]; arr[tmp] = arr[low]; arr[low] = pivot;
	int high_tmp = high;
	int low_tmp = low;
	while(low < high){
		while (low < high && arr[high] >= pivot)
		{
			--high;
		}
		arr[low] = arr[high];
		while (low < high && arr[low] <= pivot)
		{
			++low;
		}
		arr[high] = arr[low];
	}
	arr[low] = pivot;

	if (low == k - 1)
	{
		return;
	}else if(low > k - 1)
	{
		k_big(arr,low_tmp,low-1,k);
	}else
	{
		k_big(arr,low+1,high_tmp,k);
	}
}

SimSearcher::SimSearcher()
{
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
void split(const char* str, vector<string>& ret)  
{  
	ret.clear();
	string s(str);
    size_t last = 0;  
	for (int i = 0; i < s.length(); i++)
		if (s[i] == ' ') 
		{
			ret.push_back(s.substr(last, (i - last)));
			last = i + 1;
		}
	ret.push_back(s.substr(last, (s.length() - last)));
}  

int myhash(const char* str, int* dataset, int len)
{
	int res = 0;
	int res2 = 0;
	for (int i = 0; i < len; i++)
	{
		res = (res * 29 + (int)str[i]);
		if (res >= mode) res = res % mode;
		res2 = (res2 * 39 + (int)str[i]);
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

int search(const char* str, int* dataset, int len)
{
	int res = 0;
	int res2 = 0;
	for (int i = 0; i < len; i++)
	{
		res = (res * 29 + (int)str[i]);
		if (res >= mode) res = res % mode;
		res2 = (res2 * 39 + (int)str[i]);
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

int searcha(const char* str, int* dataset, int s, int len)
{
	int res = 0;
	int res2 = 0;
	for (int i = s; i < s + len; i++)
	{
		res = (res * 29 + (int)str[i]);
		if (res >= mode) res = res % mode;
		res2 = (res2 * 39 + (int)str[i]);
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
void SimSearcher::convert(const char* str, vector<int>& word, int x = 0)
{
	w.clear();
	tmp.clear();
	split(str, w);
	for (int i = 0; i < w.size(); i++)
	{
		if (x == 0)
			tmp.push_back(myhash(w[i].c_str(), jaccard_hash, w[i].length()));
		else
			tmp.push_back(search(w[i].c_str(), jaccard_hash, w[i].length()));
	}
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
	context_pre.clear();
	context.clear();
	while (getline(fin,str))
	{
		context_pre.push_back(str);
		context.push_back(str.c_str());
	}
	candid = new int[context.size()];
	candid_num = new int[context.size()];
	candid_list = new pair<int,int> [context.size()];
	candid_set = new int[context.size()];
	candid_ys = new int[context.size()];
	for (int i = 0; i < context.size(); i++)
		candid_ys[i] = -1;
	//build jaccard index
	//build invert list of jaccard
	//clear
	clear();
	jaccard_hash = new int[hash_size];
	jaccard_list = new vector<int>*[hash_size];
	context_hash = new int*[context.size()];
	context_word_len = new int[context.size()];
	context_len = new int[context.size()];
	memset(jaccard_hash, 0, sizeof(int) * hash_size);
	q = 5;
	for (int i = 0; i < hash_size; i++)
	{
		jaccard_list[i] = NULL;
	}
	vector<int> word;
	for (int i = 0; i < context.size(); i++)
	{
		word.clear();
		convert(context[i], word); 
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
		context_len[i] = context_pre[i].length();
		int r = context_pre[i].length() - q + 1;
		for (int j = 0; j < r; j++)
		{
			const char* str = context_pre[i].substr(j, q).c_str();
			int mark = myhash(str, ed_hash, q);
			if (ed_list[mark] == NULL)
			{
				ed_list[mark] = new vector<int>;
				ed_list[mark]->clear();
			}
			ed_list[mark]->push_back(i);
			if (ed_list[mark]->size() > M) M = ed_list[mark]->size();
		}
	}		
	return SUCCESS;
}

inline int cmp(pair<int, int> a, pair<int, int> b)
{
	return a.first < b.first;
}

int SimSearcher::searchJaccard(const char *query, double threshold, vector<pair<unsigned, double> > &result)
{
	result.clear();
	//find candidate
	vector<int> word;
	word.clear();
	convert(query, word, 1);
	int word_num = word.size();
	int* word_list = new int[word_num];
	int tmp = 0;
	double t1 = (double) ((word_num + min_context_len) * threshold) / (1.0 + threshold);
	int thres;
	if (ceil(threshold * word_num) >= ceil(t1)) thres = ceil(threshold * word_num); else thres = ceil(t1);
	int l_c = 0, l_cn=0, l_cl = 0;
	int* str_hash = new int[word_num];
	for (int i = 0; i < word_num; i++)
	{
		word_list[i] = word[i];
		if (word_list[i] != -1)
		{
			candid_list[l_cl] = (make_pair(jaccard_list[word_list[i]]->size(), word_list[i]));
			l_cl ++;
		}
	}
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
			{
				candid[l_c] = (i);
				l_c ++;
			}
	}
	else
	{
		k_big(candid_list, 0, l_cl-1, tmp);
		for (int i = 0; i < tmp; i++)
			for (int j = 0; j < candid_list[i].first; j++)
			{
				if (candid_ys[(*jaccard_list[candid_list[i].second])[j]] == -1)
				{
					candid_num[l_cn] = ((*jaccard_list[candid_list[i].second])[j]);
					candid_set[l_cn] = 0;
					candid_ys[((*jaccard_list[candid_list[i].second])[j])] = l_cn;
					l_cn ++;
				}
				candid_set[candid_ys[(*jaccard_list[candid_list[i].second])[j]]] ++;
			}
		int cnt;
		for (int i = 0; i < l_cn; i++)
		{
			int tt = candid_num[i];
			candid_ys[tt] = -1;
			cnt = candid_set[i];
			if (cnt < (thres - f)) continue;
			for (int j = tmp; j < l_cl; j++)
			{
				if (cnt >= thres) break;
				int h = 0, t = candid_list[j].first, mid;
				while (h < t - 1)
				{
					mid = (h + t) >> 1;
					if ((*jaccard_list[candid_list[j].second])[mid] <= tt) h = mid; else t = mid;
				}
				if ((*jaccard_list[candid_list[j].second])[h] == tt) cnt ++;
			}
//			cout << cnt << " " << thres<< endl;
			if (cnt >= thres)
			{
				candid[l_c] = (tt);
				l_c ++;
			}
		}
	}
	sort(candid, candid + l_c);
	for (int i = 0; i < l_c; i++)
	{
		double jaccard = (verify_jaccard(word_list, word_num, context_hash[candid[i]], context_word_len[candid[i]], threshold));
		if (jaccard >= threshold)
			result.push_back(make_pair(candid[i], jaccard));
	}
	delete[] word_list;
	return SUCCESS;
}

inline bool okf(int x, int y, int threshold)
{
	if (x >=0 && y >= 0 && myabs(y - x) <= threshold)
		return true;
	else
		return false;
}

int verify_ed(const char* s1, int len1, const char* s2, int len2, int threshold)
{
	if (myabs(len2 - len1) > threshold) return max_int;
	if (len2 < len1)
	{ 
		const char* tmp = s1;s1 = s2;s2 = tmp;
		int t = len1; len1 = len2; len2 = t;
	}
	int** dp = new int*[2];
	dp[0] = new int[len2];
	dp[1] = new int[len2]; 
	int t;
	if (len2 <= threshold) t = len2; else t = threshold + 1;
	for (int i = 0; i < t; i++)
	{
		if (s1[0] == s2[i]) dp[0][i] = i; else dp[0][i] = i + 1;
		if (i > 0 && dp[0][i-1] + 1 < dp[0][i]) dp[0][i] = dp[0][i-1] + 1;
	}
	int now = 0;
//	cout << s1 << " " << s2  <<endl;
	if (len1 > 1)
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
				return max_int;
		}
	int result = dp[now][len2 - 1];
	return result;
}



int SimSearcher::searchED(const char *query, unsigned threshold, vector<pair<unsigned, unsigned> > &result)
{
	result.clear(); 
	int len = strlen(query);
	//calc threshold
	int thres = getmax(len - q_gram + 1 - threshold * q_gram, 0);
	int word_num = len - q_gram + 1;
	int l_c = 0, l_cn=0, l_cl = 0;
	if (thres == 0)
		for (int i = 0; i < context.size(); i++)
		{
			candid[l_c] = (i);
			l_c ++;
		}
	else
	{
		int* word_list = new int[word_num];
		for (int i = 0; i < word_num; i++)
		{
			word_list[i] = searcha(query, ed_hash, i, q_gram);
			if (word_list[i] != -1)
			{
				candid_list[l_cl] = (make_pair(ed_list[word_list[i]]->size(), word_list[i]));
				l_cl ++;
			}
		}
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
					candid[l_c] = (i);
					l_c ++;
				}
			}
		}
		else
		{
			k_big(candid_list, 0, l_cl-1, tmp);
			for (int i = 0; i < tmp; i++)
			{
				for (int j = 0; j < candid_list[i].first; j++)
				{
					if (candid_ys[(*ed_list[candid_list[i].second])[j]] == -1)
					{
						candid_num[l_cn] = ((*ed_list[candid_list[i].second])[j]);
						candid_set[l_cn] = 0;
						candid_ys[((*ed_list[candid_list[i].second])[j])] = l_cn;
						l_cn ++;
					}
					candid_set[candid_ys[(*ed_list[candid_list[i].second])[j]]] ++;
				}
			}
			int cnt;
			for (int i = 0; i < l_cn; i++)
			{
				int tt = candid_num[i];
				candid_ys[tt] = -1;
				cnt = candid_set[i];
				if (cnt < thres - f) continue;
				if (myabs(context_len[candid_num[i]] - len) > threshold) continue;
				for (int j = tmp; j < l_cl; j++)
				{
					if (cnt >= thres) break;
					int h = 0, t = candid_list[j].first, mid;
					while (h < t - 1)
					{
						mid = (h + t) >> 1;
						if ((*ed_list[candid_list[j].second])[mid] <= tt) h = mid; else t = mid;
					}
					if ((*ed_list[candid_list[j].second])[h] == tt) cnt ++;
				}	
				if (cnt >= thres)
				{
					candid[l_c] = tt;
					l_c ++;
				}
			}
		}
		delete[] word_list;
	}
	sort(candid, candid + l_c);
	for (int i = 0; i < l_c; i++)
	{
		int num = verify_ed(query, len, context[candid[i]], context_len[candid[i]], threshold);
		if (num <= threshold)
			result.push_back(make_pair(candid[i], num));
	}
	return SUCCESS;
}

