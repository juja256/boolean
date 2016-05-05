#ifndef CORE_H
#define CORE_H

#define XOR_EXC 0x00000001
#define RES_EXC 0x00000002
#define FLW_EXC 0x00000003
#define EVL_EXC 0x00000004

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
	
	u32 GetBlocks() const;
	u64 GetSize() const; 
	u32 Deg() const; 
	BooleanVector GetInverse();
	u32 GetInt(u32 index) const;
	bool operator[](u32 i) const; 
	BooleanVector operator^(const BooleanVector& v) const; 
	BooleanVector operator*(const BooleanVector& v);
	BooleanVector& operator=(const BooleanVector& v);
	void operator++();
	BooleanVector Pow(u32 p, const BooleanVector& v) const;
	bool SetBit(u64 index, u32 b);
	void Annulate();
	friend BooleanVector operator%(const BooleanVector& v1, const BooleanVector& v2);
	friend BooleanVector operator<<(const BooleanVector& v, u32 p);
};

BooleanVector operator<<(const BooleanVector& v, u32 p);
BooleanVector operator%(const BooleanVector& v1, const BooleanVector& v2);

class BooleanFunction {
	u32 n;
	u32 m;
	BooleanVector* table;
	u64 in_dim;
	double** cached_FS;
	int** cached_WS;
	BooleanVector(*function_ptr)(const BooleanVector&);
public:
	BooleanFunction();
	BooleanFunction(u32 n_, u32 m_, BooleanVector* table_);
	BooleanFunction(u32 n_, u32 m_, BooleanVector (*function_ptr)(const BooleanVector&));
	BooleanFunction(const BooleanFunction&);
	~BooleanFunction();
	BooleanVector GetCordinateVector(u32 i);
	BooleanVector Eval(const BooleanVector& vec);
	double* GetFourierSpectrum(u32 i);
	int* GetWalshSpectrum(u32 i);
	BooleanVector GetAlgebraicNormalForm(u32 i);
	u32 GetAlgebraicDegree();
	u32 GetAlgebraicDegree(u32 i);
	u64 GetUnlinearity(u32 i);
	u32 GetDissballance(u32 i);
	int GetCorrelationImmunityLevel(u32 i);
	u32 GetErrorExpandingCoefficient(u32 index, u32 v);
	u32 GetErrorExpandingCoefficientAverage(u32 v);
	bool GetAvalancheEffectZeroLevel(u32 i);
	bool GetAvalancheEffectZeroLevel();
	bool GetAvalancheEffectAverage();
	BooleanVector* Derivative(const BooleanVector& a);
	double GetMaximumDifferentialProbability();
	void Dump(const char* file);
};

/* algs */
BooleanVector ANF(const BooleanVector& outs, u32 n);
int* FFT(const BooleanVector& outs, u32 n);
int* WAT(const BooleanVector& outs, u32 n);
u32 HW(u32 i);
u32 HW(const BooleanVector& v);

#endif
