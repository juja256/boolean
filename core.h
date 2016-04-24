#ifndef CORE_H
#define CORE_H

#define XOR_EXC 0x00000001
#define RES_EXC 0x00000002
#define FLW_EXC 0x00000003

typedef unsigned long long u64;
typedef unsigned int u32;
typedef unsigned char u8;

class BooleanException {
	u32 code;
public:
	BooleanException();
	BooleanException(u32 code);
};

class BooleanVector {
	u64 size;
	u32 blocks;
	u32* a;
	bool heap_alloc;
public:
	explicit BooleanVector(u32 a);
	explicit BooleanVector(u64 a);
	BooleanVector(u32* c, u64 size);
	BooleanVector(const BooleanVector& vec);
	BooleanVector();
	~BooleanVector();
	
	u32 GetBlocks();
	u64 GetSize(); 
	u32 Deg() const; 
	BooleanVector GetInverse();
	bool operator[](u32 i); 
	BooleanVector operator^(const BooleanVector& v) const; 
	BooleanVector operator*(const BooleanVector& v);
	BooleanVector& operator=(const BooleanVector& v);
	
	BooleanVector Pow(u32 p, const BooleanVector& v);
	bool SetBit(u64 index, u32 b);

	friend BooleanVector operator%(const BooleanVector& v1, const BooleanVector& v2);
	friend BooleanVector operator<<(const BooleanVector& v, u32 p);
	friend BooleanVector ANF(BooleanVector& outs, u32 n);
	friend int* FourierSeries(BooleanVector& outs);
	friend int* WalshAdamarSeries(BooleanVector& outs);
};

BooleanVector operator<<(const BooleanVector& v, u32 p);
BooleanVector operator%(const BooleanVector& v1, const BooleanVector& v2);
BooleanVector ANF(BooleanVector& outs, u32 n);
u32 HammingWeight(u32 i);
int* FourierSeries(BooleanVector& outs);
int* WalshAdamarSeries(BooleanVector& outs);

class BooleanFunction {
	u32 n;
	u32 m;
	BooleanVector* table;
	u64 in_dim;
public:
	BooleanFunction(u32 n_, u32 m_, BooleanVector* table_);
	BooleanFunction(u32 n_, u32 m_, BooleanVector (*function_ptr)(BooleanVector));
	
	BooleanVector GetCordinateVector(u32 i);
	int* GetFourierSpectrum(u32 i);
	int* GetWalshSpectrum(u32 i);
	BooleanVector GetAlgebraicNormalForm(u32 i);
	u32 GetAlgebraicDegree(u32 i);
	int GetUnlinearity(u32 i);
	int GetDissballance(u32 i);
	int GetErrorExpandingLevel(u32 i);
	bool GetAvalancheEffectZeroLevel();
	bool GetAvalancheEffectAverage();
};

#endif