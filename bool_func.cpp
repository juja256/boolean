#include "core.h"
#include <math.h>
#include <stdlib.h>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif
BooleanFunction::BooleanFunction(): n(0), m(0), table(0), function_ptr(0), cached_FS(0), cached_WS(0) {

}

BooleanFunction::BooleanFunction(u32 n_, u32 m_, BooleanVector* table_): n(n_), m(m_), table(table_), function_ptr(0) {
	cached_FS = new double*[m_];
	cached_WS = new int*[m_];
	for (u32 i = 0;i < m_;i++) {
		cached_FS[i] = 0;
		cached_WS[i] = 0;
	}
	this->in_dim = 0x00000001 << n_;
}

BooleanFunction::BooleanFunction(u32 n_, u32 m_, BooleanVector (*function_ptr_)(const BooleanVector&)): n(n_), m(m_), function_ptr(function_ptr_) {
	cached_FS = new double*[m_];
	cached_WS = new int*[m_];
	for (u32 i = 0;i < m_;i++) {
		cached_FS[i] = 0;
		cached_WS[i] = 0;
	}
	this->in_dim = (u64)1 << n_;
	this->table = new BooleanVector[in_dim];
	
	for (u64 i = 0; i < in_dim; i++) {
		this->table[i] = this->function_ptr(BooleanVector((u32)i));
	}
}

BooleanFunction::BooleanFunction(const BooleanFunction& f) {
	this->in_dim = f.in_dim;
	this->n = f.n;
	this->m = f.m;
	this->cached_FS = f.cached_FS;
	this->cached_WS = f.cached_WS;
	this->table = f.table;
}

BooleanFunction::~BooleanFunction() {
	if (this->table)
		delete[] this->table;
	for (u32 i = 0; i < this->m; i++) {
		if (this->cached_FS[i])
			delete[] this->cached_FS[i];
		if (this->cached_WS[i])
			delete this->cached_WS[i];
	}
	delete[] this->cached_FS;
	delete[] this->cached_WS;
}

BooleanVector BooleanFunction::GetCordinateVector(u32 index) {
	BooleanVector res(this->in_dim);
	for (u64 i = 0; i < this->in_dim; i++) {
		res.SetBit(i, this->table[i][index]);
	}
	return res;
}

BooleanVector BooleanFunction::GetAlgebraicNormalForm(u32 index) {
	BooleanVector vec = this->GetCordinateVector(index);
	return ANF(vec, this->n);
}

u32 BooleanFunction::GetAlgebraicDegree(u32 index) {
	BooleanVector anf = this->GetAlgebraicNormalForm(index);
	u32 max = 0;
	u32 hw;
	for (int i = 0; i < anf.GetSize(); i++) {
		hw = HW(i);
		if (anf[i] && (hw > max))
			max = hw;
	}
	return max;
}

u32 BooleanFunction::GetAlgebraicDegree() {
	u32 max = 0, a;
	for (u32 i = 0;i < this->m; i++) {
		a = GetAlgebraicDegree(i);
		if (a > max)
			max = a;
	}
	return max;
}

int* BooleanFunction::GetWalshSpectrum(u32 i) {
	if (!this->cached_WS[i])
		this->cached_WS[i] = WAT(this->GetCordinateVector(i), this->n);
	return this->cached_WS[i];
}

double* BooleanFunction::GetFourierSpectrum(u32 i) {
	if (this->cached_FS[i])
		return this->cached_FS[i];
	int* v = FFT(this->GetCordinateVector(i), this->n);
	double* v2 = new double[this->in_dim];
	for (u64 j = 0; j < this->in_dim; j++)
		v2[j] = 1. * v[j] / this->in_dim;
	delete[] v;
	this->cached_FS[i] = v2;
	return v2;
}

u32 BooleanFunction::GetDissballance(u32 i) {
	if (!this->cached_WS[i]) {
		GetWalshSpectrum(i);
	}
	return abs(cached_WS[i][0]);
}

u64 BooleanFunction::GetUnlinearity(u32 i) {
	if (!this->cached_WS[i])
		GetWalshSpectrum(i);
	u64 max=0;
	for (u64 j = 0; j < this->in_dim; j++) {
		if (abs(cached_WS[i][j]) > max)
			max = abs(cached_WS[i][j]);
	}
	return ((this->in_dim) >> 1) - max / 2;
}

int BooleanFunction::GetCorrelationImmunityLevel(u32 i) {
	if (!this->cached_WS[i])
		GetWalshSpectrum(i);
	bool fl = true;
	for (int k = n; k >= 0; k--) {
		for (u32 j = 0;j < this->in_dim; j++) {
			if ((HW(j) <= k) && (this->cached_WS[i][j])) {
				fl = false;
				break;
			}
		}
		if (fl)
			return k;
		else
			fl = true;
	}
	return -1;
}

BooleanVector BooleanFunction::Eval(const BooleanVector& vec) {
	if (vec.GetBlocks() != 1) throw BooleanException(EVL_EXC);
	return this->table[vec.GetInt(0)];
}

u32 BooleanFunction::GetErrorExpandingCoefficient(u32 index, u32 v) {
	u32 s = 0;
	BooleanVector x;
	BooleanVector e;
	e.SetBit(v, 1);
	for (u32 i = 0; i < this->in_dim; i++) {
		s += (this->Eval(x)[index] ^ this->Eval(x ^ e)[index]);
		++x;
	}
	return s;
}

u32 BooleanFunction::GetErrorExpandingCoefficientAverage(u32 v) {
	u32 s = 0;
	BooleanVector x;
	BooleanVector e;
	e.SetBit(v, 1);
	for (u32 i = 0; i < this->in_dim; i++) {
		s += HW(this->Eval(x) ^ this->Eval(x ^ e));
		++x;
	}
	return s;
}

bool BooleanFunction::GetAvalancheEffectZeroLevel(u32 index) {
	for (u32 i = 0; i < this->n; i++) {
		u32 n_ = 1 << (this->n - 1);
		if (this->GetErrorExpandingCoefficient(index, i) != n_)
			return false;
	}
	return true;
}

bool BooleanFunction::GetAvalancheEffectZeroLevel() {
	for (u32 i = 0; i < this->m; i++) {
		if (!GetAvalancheEffectZeroLevel(i))
			return false;
	}
	return true;
}

bool BooleanFunction::GetAvalancheEffectAverage() {
	u32 s = 0;
	for (u32 i = 0; i < this->n; i++) {
		s = this->GetErrorExpandingCoefficientAverage(i);
		if (s != this->m * (this->in_dim >> 1))
			return false;
	}
	return true;
}

void BooleanFunction::Dump(const char* file) {
	std::ofstream of(file);
	of << "Truth table:\n";
	for (u32 i=0; i<this->in_dim; i++) {
		of << i << ": ";
		for (u32 j=0; j<this->m; j++) {
			of << this->table[i][j];
		}
		of << "\n";
	}
}

BooleanVector* BooleanFunction::Derivative(const BooleanVector& a) {
	BooleanVector* t = new BooleanVector[this->in_dim];
	BooleanVector x;
	for (u32 i = 0; i < this->in_dim; i++) {
		t[i] = this->Eval(x) ^ this->Eval(x ^ a);
		++x;
	}
	return t;
}

u32 GetMostFrequentOutputFrequency(BooleanVector* table, u32 in_dim, u32 out_dim) {
	u32* t = new u32[out_dim]();
	u32 k = 0;
	u32 max = 0;
	for (int i = 0; i < in_dim; i++) {
		k = ++t[table[i].GetInt(0)];
		if (k > max) max = k;
	}
	delete[] t;
	return max;
}

double BooleanFunction::GetMaximumDifferentialProbability() {
	u32 out_dim = 0x00000001 << this->m;
	BooleanVector a((u32)1);
	u32 max = 0;
	u32 t;
	BooleanVector* D;
	
	for (int i = 1; i < this->in_dim; i++) {
		D = this->Derivative(a);
		t = GetMostFrequentOutputFrequency(D, this->in_dim, out_dim);
		if (t > max)
			max = t;
		++a;
		delete[] D;
	}
	return (double)max / this->in_dim;
}

