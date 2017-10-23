/*
 * HModel.cpp
 *
 *  Created on: 2017年9月26日
 *      Author: leezear
 */

#include "HModel.h"

namespace cudacp {
////////////////////////////////////////////////////////////////////
HVar::HVar(const int id, const string name, const int min_val,
	const int max_val) :
	id(id), name(name), std_max(max_val - min_val) {
	int j = 0;
	const int size = max_val - min_val + 1;
	vals.resize(size);
	anti_map.resize(size);
	for (int i = min_val; i <= max_val; ++i) {
		val_map[i] = j;
		vals[j] = j;
		anti_map[j] = i;
		++j;
	}
}

HVar::HVar(const int id, const string name, vector<int>& v) :
	id(id), name(name), std_max(v.size() - 1) {
	vals.resize(vals.size());
	for (size_t i = 0; i < vals.size(); ++i) {
		val_map[v[i]] = i;
		vals[i] = i;
		anti_map[i] = v[i];
	}
}

void HVar::Show() {
	cout << id << "| " << name << ": ";
	for (size_t i = 0; i < vals.size(); ++i) {
		cout << vals[i] << "[" << val_map[i] << "] ";
	}
	cout << endl;
}

HVar::~HVar() {}
////////////////////////////////////////////////////////////////////
HTab::HTab(const int id, const bool sem, vector<vector<int>>& ts, vector<HVar*>& scp) :
	id(id), semantics(sem), scope(scp) {
	unsigned long all_size = 1;
	for (auto i : scp)
		all_size *= i->vals.size();
	unsigned long sup_size;

	if (!sem)
		sup_size = all_size - ts.size();
	else
		sup_size = ts.size();
	vector<int> ori_t_(scope.size());
	vector<int> std_t_(scope.size());
	tmp_t_.resize(scope.size());
	tuples.resize(sup_size, vector<int>(scope.size()));

	if (sem) {
		for (size_t i = 0; i < sup_size; i++) {
			GetSTDTuple(ts[i], std_t_);
			tuples[i] = std_t_;
		}
	}
	else {
		int j = 0;
		for (int i = 0; (i < all_size) && (j <= sup_size); ++i) {
			GetTuple(i, ori_t_, std_t_);
			if (find(ts.begin(), ts.end(), ori_t_) == ts.end())
				tuples[j++] = std_t_;
		}
	}

	semantics = true;
	isSTD = true;
	//Show();
	//	tuples = ts;
}

HTab::HTab(HTab * t, vector<HVar *>& scp) :
	id(t->id + 1),
	semantics(t->semantics),
	scope(scp),
	isSTD(true) {
	tuples = t->tuples;
}

//void HTab::GetSTDTuple(vector<int>& tuple) {
//	for (size_t i = 0; i < tuple.size(); ++i)
//		tuple[i] = scope[i]->val_map[tuple[i]];
//}
void HTab::GetSTDTuple(vector<int>& src_tuple, vector<int>& std_tuple) {
	for (size_t i = 0; i < src_tuple.size(); ++i)
		std_tuple[i] = scope[i]->val_map[src_tuple[i]];
}

void HTab::GetORITuple(vector<int>& std_tuple, vector<int>& ori_tuple) {
	for (size_t i = 0; i < std_tuple.size(); ++i)
		ori_tuple[i] = scope[i]->anti_map[std_tuple[i]];
}

int HTab::GetAllSize() const {
	int size = 1;
	for (auto v : scope)
		size *= v->vals.size();
	return size;
}

bool HTab::SAT(vector<int>& t) {
	return find(tuples.begin(), tuples.end(), t) == tuples.end();
}

bool HTab::SAT_STD(vector<int>& t) {
	GetORITuple(t, tmp_t_);
	return find(tuples.begin(), tuples.end(), tmp_t_) == tuples.end();
}

void HTab::GetTuple(int idx, vector<int>& src_t, vector<int>& std_t) {
	for (int i = (scope.size() - 1); i >= 0; --i) {
		HVar* v = scope[i];
		std_t[i] = idx % v->vals.size();
		src_t[i] = v->anti_map[std_t[i]];
		idx /= v->vals.size();
	}
}

//void HTab::GetTuple(int idx, vector<int>& t) {
//	for (int i = (scope.size() - 1); i >= 0; --i) {
//		HVar* v = scope[i];
//		t[i] = idx % v->vals.size();
//		idx /= v->vals.size();
//	}
//}

ostream & operator<<(ostream &os, const vector<HVar*>& a) {
	for (auto v : a)
		os << v->id << "[" << v->name << "] ";
	return os;
}

void HTab::Show() {
	const string sem = semantics ? "supports" : "conflicts";
	cout << "id: " << id << " semantics: " << sem << " size: " << tuples.size()
		<< " arity:" << scope.size() << " scope = {" << scope << endl;
	for (auto t : tuples) {
		cout << "( ";
		for (auto i : t)
			cout << i << " ";
		cout << ")";
	}
	cout << endl;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
HModel::HModel() {
	// TODO Auto-generated constructor stub

}

void HModel::AddVar(const int id, const string name, const int min_val,
	const int max_val) {
	HVar* var = new HVar(id, name, min_val, max_val);
	var_n_[name] = var;
	//	if ((size_t) id >= vars.size())
	//		vars.reserve(id + 1);
	//	vars[id] = var;
	vars.push_back(var);
	mds_ = max(mds_, var->vals.size());
}

void HModel::AddVar(const int id, const string name, vector<int>& v) {
	HVar* var = new HVar(id, name, v);
	var_n_[name] = var;
	//	if ((size_t) id >= vars.size())
	//		vars.reserve(id + 1);
	//	vars[id] = var;
	vars.push_back(var);
	mds_ = max(mds_, var->vals.size());
}

void HModel::AddTab(const int id, const bool sem, vector<vector<int>>& ts, vector<HVar*>& scp) {
	HTab* t = new HTab(id, sem, ts, scp);
	tabs.push_back(t);
	mas_ = max(mas_, t->scope.size());
}

void HModel::AddTab(const int id, const bool sem, vector<vector<int>>& ts, vector<string>& scp) {
	vector<HVar*> scope(scp.size());
	for (size_t i = 0; i < scp.size(); ++i)
		scope[i] = var_n_[scp[i]];
	AddTab(id, sem, ts, scope);
}

void HModel::AddTabAsPrevious(HTab* t, vector<string>& scp) {
	vector<HVar*> scope(scp.size());
	for (size_t i = 0; i < scp.size(); ++i)
		scope[i] = var_n_[scp[i]];
	HTab* nt = new HTab(t, scope);
	tabs.push_back(nt);
	mas_ = max(mas_, nt->scope.size());
}

void HModel::Show() {
	cout << "--------------Variables--------------" << endl;
	cout << "size: " << vars.size() << "\tmax domain size :" << mds_ << endl;
	for (auto v : vars)
		v->Show();
	cout << "-------------Constraints-------------" << endl;
	cout << "size: " << tabs.size() << "\tmax arity size :" << mas_ << endl;
	for (auto t : tabs)
		t->Show();
}

HModel::~HModel() {
	for (auto i : vars)
		delete i;
	vars.clear();

	for (auto i : tabs)
		delete i;
	tabs.clear();
}

} /* namespace cudacp */
