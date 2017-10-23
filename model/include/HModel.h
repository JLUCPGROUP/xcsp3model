/*
 * HModel.h
 *
 *  Created on: 2017年9月26日
 *      Author: leezear
 */

#ifndef HMODEL_H_
#define HMODEL_H_
#define _EXPORTING
#ifdef _EXPORTING
#define API_DECLSPEC __declspec(dllexport)
#else
#define API_DECLSPEC __declspec(dllimport)
#endif

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>

using namespace std;

namespace cudacp {

class API_DECLSPEC HVar {
public:
	int id;
	string name;
	vector<int> vals;
	unordered_map<int, int> val_map;
	vector<int> anti_map;
	const int std_min = 0;
	const int std_max;
	HVar(const int id, const string name, const int min_val, const int max_val);
	HVar(const int id, const string name, vector<int>& v);
	~HVar();
	void Show();
private:
};

class API_DECLSPEC HTab {
public:
	int id;
	string name;
	bool semantics;
	vector<HVar*> scope;
	vector<vector<int>> tuples;
	bool isSTD = false;
	HTab(const int id, const bool sem, vector<vector<int>>& ts,
		vector<HVar*>& scp);
	HTab(HTab* t, vector<HVar*>& scp);
	int GetAllSize() const;
	void GetSTDTuple(vector<int>& src_tuple, vector<int>& std_tuple);
	void GetORITuple(vector<int>& std_tuple, vector<int>& ori_tuple);
	bool SAT(vector<int>& t);
	bool SAT_STD(vector<int>& t);
	void Show();
	void GetTuple(int idx, vector<int>& t, vector<int>& t_idx);
private:
	//临时变量
	vector<int> tmp_t_;
	API_DECLSPEC friend ostream& operator<<(ostream &os, const vector<HVar*>& a);
};

class API_DECLSPEC HModel {
public:
	vector<HVar*> vars;
	vector<HTab*> tabs;
	unordered_map<string, HVar*> var_n_;
	//	unordered_map<int, HVar*> var_i_;
	HModel();
	virtual ~HModel();
	void AddVar(const int id, const string name, const int min_val,
		const int max_val);
	void AddVar(const int id, const string name, vector<int>& v);
	void AddTab(const int id, const bool sem, vector<vector<int>>& ts, vector<HVar*>& scp);
	void AddTab(const int id, const bool sem, vector<vector<int>>& ts, vector<string>& scp);
	void AddTabAsPrevious(HTab* t, vector<string>& scp);
	int max_domain_size() const { return mds_; }
	void Show();
private:
	size_t mds_ = 0;
	size_t mas_ = 0;

};

} /* namespace cudacp */

#endif /* HMODEL_H_ */
