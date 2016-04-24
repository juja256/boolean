#include "core.h"
#include <iostream>
#define PRINT_VECTOR(v) for (u64 i=0;i<(v).GetBlocks()*32;i++) std::cout<<(v)[i];std::cout<<"\n";

int main() {
	BooleanVector v7((u32)0x00000004);
	PRINT_VECTOR(v7)
	v7.SetBit(0, 1);
	PRINT_VECTOR(v7)
	u32 v_8[1] = {248};
	BooleanVector v8(v_8, 8);
	PRINT_VECTOR(ANF(v8, 3))
	int * fs = WAT(v8, 3);
	for (u32 i = 0;i < 8;i++)
		std::cout << fs[i] << " ";



	delete[] fs;
	system("PAUSE");
	return 0;
}
