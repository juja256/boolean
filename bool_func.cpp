#include "core.h"
#include <math.h>

BooleanFunction::BooleanFunction(u32 n_, u32 m_, BooleanVector* table_): n(n_), m(n_), table(table_) {
	cached_FS = new double*[m_];
	cached_WS = new int*[m_];
	for (u32 i = 0;i < m_;i++) {
		cached_FS[i] = 0;
		cached_WS[i] = 0;
	}
	this->in_dim = 0x00000001 << (n_ - 1);
}

BooleanFunction::BooleanFunction(u32 n_, u32 m_, BooleanVector(*function_ptr)(BooleanVector)): n(n_), m(m_) {
	cached_FS = new double*[m_];
	cached_WS = new int*[m_];
	for (u32 i = 0;i < m_;i++) {
		cached_FS[i] = 0;
		cached_WS[i] = 0;
	}
	this->in_dim = 0x00000001 << (n_ - 1);
	this->table = new BooleanVector[in_dim];
	
	for (u64 i = 0; i < in_dim; i++) {
		this->table[i] = (*function_ptr)(BooleanVector((u32)i));
	}
}
BooleanFunction::~BooleanFunction() {
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
	this->cached_WS[i] = WAT(this->GetCordinateVector(i), this->n);
	return this->cached_WS[i];
}

double* BooleanFunction::GetFourierSpectrum(u32 i) {
	int* v = FFT(this->GetCordinateVector(i), this->n);
	double* v2 = new double[this->in_dim];
	for (u64 i = 0; i < this->in_dim; i++)
		v2[i] = 1. * v[i] / this->in_dim;
	delete[] v;
	this->cached_FS[i] = v2;
	return v2;
}

u32 BooleanFunction::GetDissballance(u32 i) {
	if (!this->cached_WS[i]) {
		GetWalshSpectrum(i);
	}
	return cached_WS[i][0];
}

u32 BooleanFunction::GetUnlinearity(u32 i) {
	if (!this->cached_WS[i])
		GetWalshSpectrum(i);
	u32 max=0;
	for (u64 j = 0; j < this->in_dim; j++) {
		if (abs(cached_WS[i][j]) > max)
			max = abs(cached_WS[i][j]);
	}
	return ((this->in_dim) >> 1) - max / 2;
}