#include "core.h"

BooleanVector ANF(BooleanVector& outs, u32 n) { // Fast Mebius Transform +++++
	BooleanVector res(outs);
	u64 s = outs.size;
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

u32 HammingWeight(u32 i) {
	i = i - ((i >> 1) & 0x55555555);
	i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}