#include "core.h"

BooleanVector ANF(const BooleanVector& outs, u32 n) { // Fast Mebius Transform +++++
	BooleanVector res(outs);
	u64 s = outs.GetSize();
	u32 step, step2;
	for (u32 i = 1; i <= n; i++) {
		step = 1 << i;
		for (u32 j = 0; j < s; j += step) {
			step2 = (step >> 1);
			for (u32 k = j; k < j + step2; k += 1) {
				res.SetBit(k + step2, res[k + step2] ^ res[k]);
			}
		}
	}
	return res;
}

int* FFT(const BooleanVector& outs, u32 n) {
	u64 s = outs.GetSize();
	int* spectre = new int[s];
	for (u32 i = 0; i < s; i++)
		spectre[i] = outs[i];

	u32 step, step2;
	int u0, u1;
	for (u32 i = 1; i <= n; i++) {
		step = 1 << i;
		for (u32 j = 0; j < s; j += step) {
			step2 = (step >> 1);
			for (u32 k = j; k < j + step2; k += 1) {
				u0 = spectre[k];
				u1 = spectre[k + step2];
				spectre[k] = u0 + u1;
				spectre[k + step2] = u0 - u1;
			}
		}
	}
	return spectre;
}

int* WAT(const BooleanVector& outs, u32 n) {
	u64 s = outs.GetSize();
	int* spectre = new int[s];
	for (u32 i = 0; i < s; i++)
		spectre[i] = outs[i] ? -1 : 1;

	u32 step, step2;
	int u0, u1;
	for (u32 i = 1; i <= n; i++) {
		step = 1 << i;
		for (u32 j = 0; j < s; j += step) {
			step2 = (step >> 1);
			for (u32 k = j; k < j + step2; k += 1) {
				u0 = spectre[k];
				u1 = spectre[k + step2];
				spectre[k] = u0 + u1;
				spectre[k + step2] = u0 - u1;
			}
		}
	}
	return spectre;
}

u32 HW(u32 i) {
	i = i - ((i >> 1) & 0x55555555);
	i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

u32 HW(const BooleanVector& v) {
	u32 w=0;
	for (u32 i = 0; i < v.GetBlocks(); i++) {
		w += HW(v.GetInt(i));
	}
	return w;
}