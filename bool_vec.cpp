#include "core.h"

#define UN_INV 0x80000000
#define MAX_U32 0xffffffff

#ifdef _DEBUG
#include <iostream>
#define PRINT_VECTOR(v) for (u64 i=0;i<(v).GetBlocks()*32;i++) std::cout<<(v)[i];std::cout<<"\n";
#endif

BooleanException::BooleanException(): code(0) {}

BooleanException::BooleanException(u32 code_): code(code_) {}

BooleanVector::BooleanVector(): size(0), blocks(0), a(0), heap_alloc(false) {

}

BooleanVector::BooleanVector(u32 vec): size(32), blocks(1), heap_alloc(true) {
	this->a = new u32[1];
	this->a[0] = vec;
}

BooleanVector::BooleanVector(u32* c, u64 size_): size(size_), a(c), heap_alloc(false) {
	this->blocks = (size_ % 32 == 0) ? size_ / 32 : size_ / 32 + 1;
}

BooleanVector::BooleanVector(u64 size_): size(size_), heap_alloc(true) {
	if (size_ % 32 != 0) {
		this->blocks = size_/32 + 1;
	}
	else {
		this->blocks = size/32;
	}
	this->a = new u32[this->blocks]; 
	for (u32 i = 0; i<this->blocks; i++)
		this->a[i] = 0;
}

BooleanVector::BooleanVector(const BooleanVector& vec) { // proper copy
	
	this->size = vec.size;
	this->blocks = vec.blocks;
	this->a = new u32[vec.blocks];
	this->heap_alloc = true;
	for (u32 i=0; i<vec.blocks; i++) {
		this->a[i] = vec.a[i];
	}
}

BooleanVector& BooleanVector::operator=(const BooleanVector& vec) {
	if ((this->size != 0) && (this->heap_alloc))
		delete[] (this->a);
	this->size = vec.size;
	this->blocks = vec.blocks;
	this->a = new u32[vec.blocks];
	this->heap_alloc = true;
	for (u32 i=0; i<vec.blocks; i++)
		this->a[i] = vec.a[i];
	return *this;
}

BooleanVector::~BooleanVector() {
	if ((this->size != 0) && (this->heap_alloc))  
		delete[] (this->a);
}

BooleanVector BooleanVector::GetInverse() {
	BooleanVector res(*this);
	for (u32 i=0; i<this->blocks; i++) {
		this->a[i] = ~(this->a[i]);
	}
	return res;
}

u64 BooleanVector::GetSize() { return this->size; }
u32 BooleanVector::GetBlocks() { return this->blocks; }

bool BooleanVector::operator[](u32 i) {
	return (this->a[i/32] & ( 1 << (i%32) )) >> (i%32);	
}

BooleanVector BooleanVector::operator^(const BooleanVector& vec) const { 
	u32 b = (vec.blocks > this->blocks) ? this->blocks : vec.blocks;
	int b_ = (vec.blocks - this->blocks);

	u64 res_size = (vec.size > this->size) ? vec.size : this->size;
	BooleanVector res(res_size);
	for (u32 i=0; i<b; i++) {
		res.a[i] = (this->a[i]) ^ (vec.a[i]);
	}
	if (b_>0)
		for (u32 i=b; i<vec.blocks; i++) {
			res.a[i] = vec.a[i];
		}
	else if (b_<0)
		for (u32 i=b; i<this->blocks; i++) {
			res.a[i] = this->a[i];
		}
	return res;
}

BooleanVector operator<<(const BooleanVector& v, u32 p) {
	if (p == 0) {
		return BooleanVector(v);
	}
	u32 deg = v.Deg();
	int offset = deg - v.blocks * 32 + p;
	BooleanVector res;
	if (offset >= 0) {
		BooleanVector v2((u64)(v.blocks * 32 + offset));
		res = v ^ v2;	
	}
	else {
		res = v;
	}

	int buf = 0;
	for (int i=0; i<res.blocks; i++) {
		u32 cur = res.a[i];
		res.a[i] <<= p;
		res.a[i] += buf;
		buf = (cur & (MAX_U32 << 32 - p)) >> 32 -p;
	}
	return res;
}

u32 BooleanVector::Deg() const { // Danger !!1
	u32 t=1;
	for (int i=(this->blocks)-1; i>=0; i--) {
		if (this->a[i] == 0) 
			continue;
		else {
			u32 cur = this->a[i];
			u32 c = 32;
			while (c) {
				if ((UN_INV & cur) >> 31)
					return c + i*32 -1;
				cur <<= 1;
				--c;
			}
		}
	}
	return 0;
}

BooleanVector BooleanVector::operator*(const BooleanVector& vec) {
	//u64 res_size = vec.blocks + this->blocks;
	BooleanVector res((u64)0);
	for (u32 i = 0; i < this->blocks; i++) {
		for (u32 j = 0; j<32; j++) {
			if ((*this)[i*32 + j]) {
				res = res ^ (vec << (i*32 + j));
			}
		}
	}
	return res;	
}

BooleanVector operator%(const BooleanVector& v1, const BooleanVector& v2) {
	u32 k;
	u64 n = v2.Deg();
	BooleanVector res(v1);
	while (res.Deg() >= n) {
		k = res.Deg();
		res = res ^ (v2 << k-n);
		//PRINT_VECTOR(res)
	}
	return res;
}

BooleanVector BooleanVector::Pow(u32 p, const BooleanVector& v) {
	BooleanVector e(p);
	BooleanVector b((u32)0x00000001);
	BooleanVector c(*this);
	u32 s = e.Deg();
	for (u32 i = 0; i <= s; i++) {
		if (e[i]) {
			b = (b * c) % v;
		}
		c = (c * c) % v;
	}
	return b;
}

bool BooleanVector::SetBit(u64 index, u32 b) {
	u32 block = index / 32;
	u32 bit = index % 32;
	if (index >= this->size) throw BooleanException(FLW_EXC);
	this->a[block] ^= (-b ^ this->a[block]) & (1 << bit);
	return b;
}
