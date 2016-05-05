#include "core.h"
#include <iostream>
#include <iomanip> 
#include <stdlib.h>
#define PRINT_VECTOR(v) for (u32 in=0;in<(v).GetBlocks()*32;in++) std::cout<<(v)[in];std::cout<<"\n";

#define PRINT_INFO_COORD(func, coord) std::cout << "Coordinate function #" << coord+1 << "\nDissb: " << func.GetDissballance(coord) << "\n"\
				<< "Deg: " << func.GetAlgebraicDegree(coord) << "\n"\
				<< "NL: " << func.GetUnlinearity(coord) << "\n"\
				<< "CIL: " << func.GetCorrelationImmunityLevel(coord) << "\n";\
				std::cout << "Error expanding coefficients:\n";\
				u32 coefs[15];\
				for (u32 in=0; in<15; in++) {\
					coefs[in] = func.GetErrorExpandingCoefficient(coord, in);\
					std::cout << coefs[in] << " ";\
				}\
				std::cout << "\nError expanding coefficients deviation:\n";\
				for (u32 in=0; in<15; in++) {\
					std::cout << std::setprecision(2) << (double)abs((int)coefs[in] - 16384) * 100 / 16384 << "% ";\
				}\
				std::cout << "\n\n";

#define PRINT_INFO_FUNCTION(func) {std::cout << "Whole Function\n"\
				<< "Deg: " << func.GetAlgebraicDegree() << "\n"\
				<< "Error expanding coefficients in average:\n";\
				u32 coefs[15];\
				for (u32 in=0; in<15; in++) {\
					coefs[in] = func.GetErrorExpandingCoefficientAverage(in);\
					std::cout << coefs[in] << " ";\
				}\
				std::cout << "\nError expanding in average coefficients deviation:\n";\
				for (u32 in=0; in<15; in++) {\
					std::cout << std::setprecision(2) << (double)abs((int)coefs[in] - 245760) * 100 / 245760 << "% ";\
				}\
				std::cout << "\n";\
				std::cout << "Presence of avalanche effect of zero level: " << func.GetAvalancheEffectZeroLevel() <<"\n"\
				<< "Presence of avalanche effect in average: " << func.GetAvalancheEffectAverage() <<"\n";\
				double mdp = func.GetMaximumDifferentialProbability();\
				std::cout << "MDP: " << mdp << " (" << std::setprecision(0) << mdp * 32768 << "/32768)\n\n";}

BooleanVector f(const BooleanVector& v) {
	BooleanVector p((u32)32771);
	u32 N = 16257;
	return v.Pow(N, p);
}

BooleanVector g(const BooleanVector& v) {
	BooleanVector p((u32)32771);
	u32 N = 16256;
	return v.Pow(N, p);
}

int main() {
	BooleanFunction func1(15, 15, f);
	BooleanFunction func2(15, 15, g);
	std::cout << "------------ f(v) = v^16257 mod (x^15 + x +1) ------------\n";
	PRINT_INFO_FUNCTION(func1);
	for (u32 i = 0; i < 15; i++) {
		PRINT_INFO_COORD(func1, i);
	}
	std::cout << "\n------------ g(v) = v^16256 mod (x^15 + x+ 1) ------------\n";
	PRINT_INFO_FUNCTION(func2);
	for (u32 i = 0; i < 15; i++) {
		PRINT_INFO_COORD(func2, i);
	}
	
	#ifdef WIN32
	system("PAUSE");
	#endif
	return 0;
}
