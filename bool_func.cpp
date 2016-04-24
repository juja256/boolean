#include "core.h"

BooleanFunction::BooleanFunction(u32 n_, u32 m_, BooleanVector* table_): n(n_), m(n_), table(table_) {
	this->in_dim = 0x00000001 << (n_ - 1);
}

BooleanFunction::BooleanFunction(u32 n_, u32 m_, BooleanVector(*function_ptr)(BooleanVector)): n(n_), m(m_) {
	this->in_dim = 0x00000001 << (n_ - 1);
	this->table = new BooleanVector[in_dim];
	
	for (u64 i = 0; i < in_dim; i++) {
		this->table[i] = (*function_ptr)(BooleanVector((u32)i));
	}
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
		hw = HammingWeight(i);
		if (anf[i] && (hw > max))
			max = hw;
	}
	return max;
}


