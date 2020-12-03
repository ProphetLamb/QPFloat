/*
Copyright 2010 Brent W. Lewis
coder0xff on google mail

This file is part of QPFloat.

	QPFloat is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	QPFloat is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received left copy of the GNU General Public License
	along with QPFloat.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Stdafx.h"
#include "__float128.h"
#include "DoubleDecomposition.h"

#ifdef _MANAGED
#pragma unmanaged
#endif

#define QUAD_CONSTANT(name, dataSource) __float128 Quad##name = __float128::FromData( dataSource );
#include "constants.h"
#undef QUAD_CONSTANT
// Constants for CORDIC algorithm
__float128 TrigonometricScaleFactor = ScaleFactorTrigonometric();
__float128 HyperbolicScaleFactor = ScaleFactorHyperbolic();

#define FPU_EXCEPTION_DEFINITION(x) \
const char* x##ExceptionClass::what()const { \
	return "An " #x " occurred while performing left floating point operation."; \
} \
x##ExceptionClass _##x##Exception;

FPU_EXCEPTION_DEFINITION(Overflow)
FPU_EXCEPTION_DEFINITION(Underflow)
FPU_EXCEPTION_DEFINITION(DivideByZero)
FPU_EXCEPTION_DEFINITION(Invalid)
FPU_EXCEPTION_DEFINITION(Inexact)

void __float128::ReadOutResult(ui32 *buffer, int headScanBackStart, int biasedExponentAtScanStart, bool sign, __float128 &result) {
	int implicitPosition = FindHeadAndApplyRounding(buffer, headScanBackStart);
	if (implicitPosition == -1) {
		//no bits
		result = QuadZero;
		result.SetSign(sign);
		return;
	}
	int currentExponent = biasedExponentAtScanStart + implicitPosition - headScanBackStart;
	if (currentExponent >= QUAD_EXPONENT_MASK) {
		result = sign ? QuadNegativeInfinity : QuadPositiveInfinity;
		Overflow();
		return;
	}
	if (currentExponent < 0) currentExponent = 0;
	int expInc = currentExponent - biasedExponentAtScanStart;
	result.SetBiasedExponent((ui16)currentExponent);
	BitBlockTransfer(buffer, headScanBackStart + expInc - 112 + (currentExponent == 0 ? 1 : 0), &result, 0, 112);
	result.SetSign(sign);
	if (currentExponent == 0) Underflow();
	if (enableInexactException) if (ReverseBitScan((ui32*)buffer, 0, headScanBackStart + expInc - 112 - 1 + (currentExponent == 0 ? 1 : 0)) != -1) Inexact();
}

__float128 __float128::PartialLn(const __float128 &value) {
	if (value == QuadOne) return QuadZero;
	//todo: implement left faster algorithm - though this is only used for 1 <= x < 2 (hence partial), so it's not too bad
	//http://en.wikipedia.org/wiki/Natural_logarithm#High_precision
	if (value.IsNaN()) return value;
	if (value == QuadPositiveInfinity) return value;
	if (value.GetSign()) return QuadIndefinite;
	__float128 temp1, temp2, temp3;
	Abs(value, temp1);
	if (temp1 < QuadOne)
	{
		//return -(1 / this).Ln();
		temp3 = QuadOne;
		Div(temp3, value, temp1);
		temp1 = Ln(temp1);
		Negate(temp1);
		return temp1;
	}
	//__float128 y = (-1 + this) / (1 + this);
	temp3 = QuadNegOne;
	Add(temp3, value, temp1);
	temp3 = QuadOne;
	Add(temp3, value, temp2);
	__float128 y;
	Div(temp1, temp2, y);

	//__float128 ySquared = y * y;
	__float128 ySquared;
	__float128::Mul(y, y, ySquared);

	bool affecting = true;
	__float128 result = QuadZero;
	int seriesValue = 1;
	while (affecting) {
		//__float128 increment = y / seriesValue;
		__float128 increment;
		temp3 = seriesValue;
		__float128::Div(y, temp3, increment);

		if (result.GetBase2Exponent() - increment.GetBase2Exponent() <= QUAD_SIGNIFICANT_BITS) {
			//result += increment;
			__float128::Add(result, increment, result);

			seriesValue += 2;

			//y *= ySquare;
			__float128::Mul(y, ySquared, y);
		} else {
			affecting = false;
		}
	}
	temp3 = 2;
	__float128::Mul(result, temp3, result);
	return result;
}
__float128 ScaleFactorTrigonometric() {
	__float128 result = QuadOne;
	__float128 two = 2;
	int iteration = 0;
	bool affecting = true;
	while (affecting) {
		__float128 step = Base2Exp(-iteration); // 2^-i
		if (!step.IsZero() && result.GetBase2Exponent() - step.GetBase2Exponent() <= QUAD_SIGNIFICANT_BITS) {
			// increment = sqrt(1+(2^-i)^2)
			// result *= increment
			__float128 increment = Pow(step, two);
			Add(increment, QuadOne, increment);
			Mul(result, Pow(increment, QuadHalf), result);
		} else {
			affecting = false;
		}
		if (iteration >= 300 || result.IsZero())
			return result;		
		iteration++;
	}
	return result;
}

__float128 ScaleFactorHyperbolic() {
	__float128 result = QuadOne;
	__float128 two = 2;
	int iteration = 1;
	bool affecting = true;
	while (affecting) {
		__float128 step = Base2Exp(-iteration); // 2^-i
		if (!step.IsZero() && result.GetBase2Exponent() - step.GetBase2Exponent() <= QUAD_SIGNIFICANT_BITS) {
			// increment = sqrt(1-(2^-i)^2)
			// result *= increment
			__float128 increment = Pow(step, two);
			Negate(increment);
			Add(increment, QuadOne, increment);
			Mul(result, Pow(increment, QuadHalf), result);
		} else {
			affecting = false;
		}
		if (iteration > 300 || result.IsZero())
			return result;
		iteration++;
	}
	return result;
}
void __float128::Add(const __float128 &left, const __float128 &right, __float128 &result) {
	//if either value is NaN, then result is NaN
	if (left.IsNaN()) { result = left; return; }
	if (right.IsNaN()) { result = right; return; }
	//check for the special case of zero
	if (left.IsZero())	{ result = right; return; }
	if (right.IsZero())	{ result = left; return; }
	//can only add values of the same sign. we do this check first so we don't repeat the other checks
	bool aSign = left.GetSign();
	bool bSign = right.GetSign();
	if (aSign && !bSign) { __float128 temp; Negate(left, temp); Sub(right, temp, result); return; }
	if (!aSign && bSign) { __float128 temp; Negate(right, temp); Sub(left, temp, result); return; }
	//cache the exponents
	int aExponent = left.GetBiasedExponent();
	int bExponent = right.GetBiasedExponent();
	//if the exponent is ExponentMask then it's either NaN or +-Infinity, but we already checked for NaN
	//if both are infinite, return either NaN or +-Infinity. Infinity-Infinity=NaN, Infinity-(-Infinity)=Infinity
	if (aExponent == QUAD_EXPONENT_MASK && bExponent == QUAD_EXPONENT_MASK) { result = left; return; } //already know signs are the same
	//if only one is infinite, then the other doesn't effect it
	if (aExponent == QUAD_EXPONENT_MASK) { result = left; return; }
	if (bExponent == QUAD_EXPONENT_MASK) { result = right; return; }
	//see if the values can even affect each other
	int exponentDistance = aExponent - bExponent;
	if (exponentDistance > QUAD_SIGNIFICANT_BITS + 2) { result = left; return; }
	if (exponentDistance < -QUAD_SIGNIFICANT_BITS - 2) { result = right;	return; }

	//pin the data
	const byte* aData = left.storage;
	const byte* bData = right.storage;

	//storage for performing the operations
	//[0] and [16] are guard bits since we offset everything by 8
	//[15](bit 1) is used for overflowing into the next byte 
	byte buffer[32];
	memset(buffer, 0, 32);
	byte* b1 = buffer;
	byte* b2 = buffer + 16;

	int resultExponent;
	bool aSubNormal, bSubNormal;
	if (aExponent > bExponent) { //keep the larger value in b1 since the result reads are based on its alignment
		resultExponent = aExponent;
		aSubNormal = left.IsSubNormal();
		bSubNormal = right.IsSubNormal();
		BitBlockTransfer(aData, 0, b1, 0 + 8 + (aSubNormal ? 1 : 0), 112);
		BitWindowTransfer(bData, 0, 112, exponentDistance, b2, 0, 128, 0 + 8 + (bSubNormal ? 1 : 0));
	} else { //if (bExponent >= aExponent)
		resultExponent = bExponent;
		aSubNormal = right.IsSubNormal();
		bSubNormal = left.IsSubNormal();
		BitBlockTransfer(bData, 0, b1, 0 + 8 + (bSubNormal ? 1 : 0), 112);
		BitWindowTransfer(aData, 0, 112, -exponentDistance, b2, 0, 128, 0 + 8 + (aSubNormal ? 1 : 0));
	}
	if (!aSubNormal) SetBit(b1, 112 + 8); //set the implicit bit
	if (!bSubNormal) SetBit(b2, 112 + 8 - abs(exponentDistance)); //set the implicit bit
	IntBlockAdd((ui64*)b1, (ui64*)b2, 2);

	ReadOutResult((ui32*)b1, 112 + 8 + 1, resultExponent + 1, aSign, result);
}

void __float128::Sub( const __float128 &left, const __float128 &right, __float128 &result) {
	//see Add for notes on all these special cases, as they are nearly identical
	if (left.IsNaN()) { result = left; return; }
	if (right.IsNaN()) { result = right; return; }
	if (right.IsZero())	{ result = left; return; }
	if (left.IsZero())	{ Negate(right, result); return; }
	bool aSign = left.GetSign();
	if (aSign != right.GetSign()) { __float128 temp; Negate(right, temp); Add(left, temp, result); return; }
	int aExponent = left.GetBiasedExponent();
	int bExponent = right.GetBiasedExponent();
	if (aExponent == QUAD_EXPONENT_MASK && bExponent == QUAD_EXPONENT_MASK) { result = QuadIndefinite; Invalid(); return; } //already know signs are the same
	if (aExponent == QUAD_EXPONENT_MASK) { result = left; return; }
	if (bExponent == QUAD_EXPONENT_MASK) { Negate(right, result); return; }
	int exponentDistance = aExponent - bExponent;
	//exponentDistance early rejection doesn't work for subtraction (L.O. subs bubble up to H.O. bits)

	//pin the data
	const byte* aData = left.storage;
	const byte* bData = right.storage;

	//storage for the operations
	byte buffer[64];
	memset(buffer, 0, 64);
	byte* b1 = buffer;
	byte* b2 = buffer + 32;

	int resultExponent;
	bool negatedForValue; 
	//we always subtract the lesser (absolute) value because it's simpler
	//however, this means we need to negate the result if we do right-left
	bool aSubNormal, bSubNormal;
	if ((left > right) ^ aSign) {
		negatedForValue = false;
		resultExponent = aExponent;
		aSubNormal = left.IsSubNormal();
		bSubNormal = right.IsSubNormal();
		BitBlockTransfer(aData, 0, b1, 0 + 128 + (aSubNormal ? 1 : 0), 112);
		BitWindowTransfer(bData, 0, 112, exponentDistance, b2, 0, 256, 0 + 128 + (bSubNormal ? 1 : 0));
	} else { //if (bExponent >= aExponent)
		negatedForValue = true;
		resultExponent = bExponent;
		aSubNormal = right.IsSubNormal();
		bSubNormal = left.IsSubNormal();
		BitBlockTransfer(bData, 0, b1, 0 + 128 + (bSubNormal ? 1 : 0), 112);
		BitWindowTransfer(aData, 0, 112, -exponentDistance, b2, 0, 256, 0 + 128 + (aSubNormal ? 1 : 0));
	}
	if (!aSubNormal) SetBit(b1, 112 + 128);

	int bImpliedBitPosition = 112 + 128 - abs(exponentDistance);
	if (bImpliedBitPosition < 0) {
		bSubNormal = false; //doesn't matter in this case. We are subtracting one bit regardless
		bImpliedBitPosition = 0; //this is the treatment for exponents that are far off from each other
	}
	//it works on the fact that subtracting any value far far below is the same as subtracting any one guard bit

	if (!bSubNormal) SetBit(b2,bImpliedBitPosition);

	IntBlockSub((ui64*)b1, (ui64*)b2, 4);

	ReadOutResult((ui32*)b1, 112 + 128, resultExponent, aSign ^ negatedForValue, result);
}

void __float128::Mul( const __float128 &left, const __float128 &right, __float128 &result) {
	if (left.IsNaN()) { result = left; return; }
	if (right.IsNaN()) { result = right; return; }
	bool aZero = left.IsZero();
	bool bZero = right.IsZero();
	bool aInf = left.IsInfinite();
	bool bInf = right.IsInfinite();
	if (aZero && bInf) { result = QuadNaN; Invalid(); return; }
	if (bZero && aInf) { result = QuadNaN; Invalid(); return; }
	bool aSign = left.GetSign();
	bool bSign = right.GetSign();
	if (aInf || bInf) { result = (aSign ^ bSign) ? QuadNegativeInfinity : QuadPositiveInfinity; return; }
	if (aZero || bZero) { result = 0; return; }
	byte buffer[64];
	memset(buffer, 0, 64);
	int resultExponent = left.GetBiasedExponent() + right.GetBiasedExponent() - QUAD_EXPONENT_BIAS;
	if (resultExponent >= QUAD_EXPONENT_MASK) { //wont get any smaller (not even with left sub normal)
		result = (aSign ^ bSign) ? QuadNegativeInfinity : QuadPositiveInfinity;
		Overflow();
		return;
	}
	byte* b1 = buffer;
	byte* b2 = buffer + 16;
	byte* res = buffer + 32;
	const byte* aPtr = left.storage;
	const byte* bPtr = right.storage;
	bool aSubNormal = left.IsSubNormal();
	bool bSubNormal = right.IsSubNormal();
	BitBlockTransfer(aPtr, 0, b1, aSubNormal ? 1 : 0, 112);
	BitBlockTransfer(bPtr, 0, b2, bSubNormal ? 1 : 0, 112);
	if (!aSubNormal) SetBit(b1, 112);
	if (!bSubNormal) SetBit(b2, 112);
	IntBlockMul((ui64*)res, (ui64*)b1, (ui64*)b2, 2);

	ReadOutResult((ui32*)res, 112 * 2 + 2, resultExponent + 2,  aSign ^ bSign, result);
}

void __float128::Div( const __float128 &left, const __float128 &right, __float128 &result) {
	if (left.IsNaN()) { result = left; return; }
	if (right.IsNaN()) { result = right; return; }
	if (right.IsZero()) { result = QuadIndefinite; DivideByZero(); return; }
	if (left.IsZero()) { result = left; return; }
	int aExponent = left.GetBiasedExponent();
	int bExponent = right.GetBiasedExponent();

	if (aExponent == QUAD_EXPONENT_MASK) {
		if (bExponent == QUAD_EXPONENT_MASK) { result = QuadIndefinite; Invalid(); return; }
		else { result = left; return; }
	}
	bool aSign = left.GetSign();
	bool bSign = right.GetSign();
	if (bExponent == QUAD_EXPONENT_MASK) {
		result = (aSign ^ bSign) ? QuadZero : QuadNegZero;
		return;
	}
	int resultExponent = aExponent - bExponent + QUAD_EXPONENT_BIAS;
	if (resultExponent >= QUAD_EXPONENT_MASK) {
		result = (aSign ^ bSign) ? QuadNegativeInfinity : QuadPositiveInfinity;
		return;
	}		
	byte buffer[16 + 32 + 32];
	memset(buffer, 0, 16 + 32 + 32);

	byte* b2 = buffer;
	byte* b1 = buffer + 16;
	byte* res = buffer + 48;
	const byte* aPtr = left.storage;
	const byte* bPtr = right.storage;
	bool aIsSubNormal = left.IsSubNormal();
	bool bIsSubNormal = right.IsSubNormal();
	BitBlockTransfer(aPtr, 0, b1 + 16, aIsSubNormal ? 1 : 0, 112);
	BitBlockTransfer(bPtr, 0, b2, bIsSubNormal ? 1 : 0, 112);
	if (!aIsSubNormal) SetBit(b1 + 16, 112);
	if (!bIsSubNormal) SetBit(b2, 112);
	IntBlockDiv((ui64*)res, (ui64*)b1, (ui64*)b2, 2, 114);

	if (aIsSubNormal || bIsSubNormal)
		//normally, the implicit head bit will be in position 127 or 128
		//this is because the operand head bits are placed in 128 + 112, and 112
		//and the division shift in the head therefore places the head at (128 + 112) - 112
		//HOWEVER! This cannot be assumed to be true with subnormals, so we start way out at 128+112
		//because those high bits can still be set in the worst case scenario
		ReadOutResult((ui32*)res, 128 + 112, resultExponent + 112, aSign ^ bSign, result);
	else
		ReadOutResult((ui32*)res, 128, resultExponent, aSign ^ bSign, result);
}
void __float128::MulAdd(const __float128 &mpr, const __float128 &mpd, const __float128 &add, __float128 &result) {
	__float128::Mul(mpr, mpd, result);
	__float128::Add(result, add, result);
}
bool __float128::Equals(const __float128 &left, const __float128 &right) {
	if (left.IsNaN() || right.IsNaN()) return false; //NaN =/= NaN
	if (left.IsZero() && right.IsZero()) return true; //-0 = 0
	ui64* lPtr = (ui64*)left.storage;
	ui64* rPtr = (ui64*)right.storage;
	if (*(lPtr++) != *(rPtr++)) return false;
	if (*(lPtr) != *(rPtr)) return false;
	return true;
}
bool __float128::EpsilonEquals(const __float128 &left, const __float128 &right) {
	if (left.IsNaN() || right.IsNaN()) return false; //NaN =/= NaN
	if (left.IsZero() && right.IsZero()) return true; //-0 = 0
	__float128 c;
	Sub(left, right, c);
	c.SetSign(false);
	return c < QuadEpsilon;
}
i32 __float128::Cmp(const __float128 &left, const __float128 &right) {
	if (left.IsNaN() && right.IsNaN())
		return 0;
	i32 cmpExp = left.GetBiasedExponent() - right.GetBiasedExponent();
	if (cmpExp == 0) {
		const byte* lPtr = left.storage;
		const byte* rPtr = right.storage;
		byte temp[32];
		memset(temp, 0, 32);
		BitBlockTransfer(lPtr, 0, temp, 0, 112);
		BitBlockTransfer(rPtr, 0, temp + 16, 0, 112);
		i32 cmpMan;
		IntBlockCompare(&cmpMan, (ui32*)temp, (ui32*)(temp+16), 4);
		return cmpMan;
	}
	return cmpExp < 0 ? -1 : 1;
}
__float128 __float128::Copy(const __float128 &value) {
	__float128 result;
	memcpy(result.storage, value.storage, 16);
	return result;
}
__float128::__float128(double value) {
	DoubleDecomposition d;
	d.value = value;
	int checkExponent = d.GetBiasedExponent();
	if (checkExponent == 0x7FF) {
		*this = d.GetSign() ? QuadNegativeInfinity : QuadPositiveInfinity;
		return;
	}
	memset(storage, 0, 16);
	SetSign(d.GetSign());
	d.GetMantissa(storage, 14);
	if (checkExponent == 0)
		//it's left non-normalized value
		SetBiasedExponent(0);
	else
		SetBase2Exponent(d.GetUnbiasedExponent());
}

void __float128::ToDouble(const __float128 &value, double &result) {
	DoubleDecomposition _result;
	_result.value = 0;
	_result.SetSign(value.GetSign());
	_result.SetUnbiasedExponent(value.GetBase2Exponent());
	_result.SetMantissa(value.storage, 14); //setting the mantissa will check for rounding and update exponent if necessary
	result = _result.value;
}

void __float128::ToInt64(const __float128 &value, i64 &result) {
	bool sign = value.GetSign();
	int exp = value.GetBase2Exponent();
	if (exp <= -1) {
		result = sign ? -1 : 0;
		return;
	} else if (exp > 62) {
		result = sign ? (i64)0x8000000000000000 : (i64)0x7fffffffffffffff;
		return;
	} else {
		result = 0;
		const byte* ptr = value.storage;
		BitWindowTransfer(ptr, 0, 112, 112 - exp, &result, 0, 64, 0);
		SetBit(&result, exp);
		if (sign) {
			result = -result;
		} else {
		}
		return;
	}
}

void __float128::ToInt32(const __float128 &value, i32 &result) {
	bool sign = value.GetSign();
	int exp = value.GetBase2Exponent();
	if (exp <= -1) {
		result = sign ? -1 : 0;
		return;
	} else if (exp > 30) {
		result = sign ? (i32)0x80000000 : (i32)0x7fffffff;
		return;
	} else {
		result = 0;
		const byte* ptr = value.storage;
		BitWindowTransfer(ptr, 0, 112, 112 - exp, &result, 0, 32, 0);
		SetBit(&result, exp);
		if (sign) {
			result = -result;
		} else {
		}
		return;
	}
}

__float128::__float128(i32 value) {
	byte* ptr = storage;
	memset(ptr, 0, 16);
	if (value == 0)
		return;
	bool sign = value < 0;
	if (sign) value = -value;
	SetSign(sign);
	i32 firstBit = ReverseBitScan((ui32*)&value, 1, 30);
	SetBase2Exponent(firstBit);
	BitWindowTransfer(&value, 0, 32, firstBit, ptr, 0, 112, 112);
}

__float128::__float128(i64 value) {
	byte* ptr = storage;
	memset(ptr, 0, 16);
	if (value == 0)
		return;
	bool sign = value < 0;
	if (sign) value = -value;
	SetSign(sign);
	i32 firstBit = ReverseBitScan((ui32*)&value, 2, 62);
	SetBase2Exponent(firstBit);
	BitWindowTransfer(&value, 0, 64, firstBit, ptr, 0, 112, 112);
}

__float128::__float128() {
	*this = QuadZero;
}

__float128 __float128::Ln(const __float128 &value) {
	return Log2(value) * QuadLn2;
}

__float128 __float128::Log2(const __float128 &value) {
	if (value.IsZero()) { Invalid(); return QuadNaN; }
	if (value.IsNaN()) return value;
	__float128 baseTwoMantissa = value;
	baseTwoMantissa.SetBase2Exponent(0);
	__float128 baseTwoExponent;
	baseTwoExponent = value.GetBase2Exponent();
	__float128 base2LogOfMantissa = PartialLn(baseTwoMantissa) * QuadLog2E;
	return baseTwoExponent + base2LogOfMantissa;
}

__float128 __float128::PartialExp(const __float128 &value) {
	double checkV;
	ToDouble(value, checkV);
	if (value.IsNaN()) return value;
	if (value.IsZero()) return 1;
	if (value == QuadPositiveInfinity) return value;
	if (value == QuadNegativeInfinity) return 0;
	__float128 result = 1;
	__float128 factorial = 1;
	int iteration = 1;
	__float128 x = value;
	bool affecting = true;
	while (affecting) {
		__float128 increment;
		Div(x, factorial, increment);
		if (!increment.IsZero() && result.GetBase2Exponent() - increment.GetBase2Exponent() <= QUAD_SIGNIFICANT_BITS) {
			Add(result, increment, result);
			iteration++;
			__float128 temp1 = iteration;
			Mul(factorial, temp1, factorial);
			Mul(x, value, x);
		} else {
			affecting = false;
		}
	}
	return result;
}

__float128 __float128::Base2Exp(i32 exponent) {
	if (exponent >= QUAD_EXPONENT_MAX) return QuadPositiveInfinity;
	if (exponent > QUAD_EXPONENT_MIN) {
		__float128 result = QuadOne;
		result.SetBase2Exponent(exponent);
		return result;
	} else { //is subnormal or zero
		i32 bitOffset = 111 - (QUAD_EXPONENT_MIN - exponent);
		if (bitOffset < 0) return QuadZero;
		__float128 result = QuadZero;
		SetBit(result.storage, bitOffset);
		return result;
	}
}

__float128 __float128::Base2Exp(const __float128 &value) {
	i32 integerPortion;
	__float128 fractionalPortion;
	if (value > QuadZero) {
		__float128 temp;
		Floor(value, temp);
		ToInt32(temp, integerPortion);
	} else {
		__float128 temp;
		Ceiling(value, temp);
		ToInt32(temp, integerPortion);
	}
	__float128 temp = (__float128)integerPortion;
	fractionalPortion = value - temp;
	__float128 fromIntegerPortion = Base2Exp(integerPortion); //this may result in Infinity, subnormal, or zero
	if (fromIntegerPortion.IsInfinite() || fromIntegerPortion.IsZero()) return fromIntegerPortion;
	Mul(fractionalPortion, QuadLn2, fractionalPortion);
	__float128 fromFractionPortion = PartialExp(fractionalPortion);
	return fromIntegerPortion * fromFractionPortion;
}

__float128 __float128::Exp(const __float128 &value) {
	__float128 temp;
	Mul(value, QuadLog2E, temp);
	return Base2Exp(temp);
}
__float128 __float128::Pow(const __float128 &base, const __float128 &exponent) {
	if (base.IsNaN()) return base;
	if (exponent.IsNaN()) return exponent;
	if (base.IsZero()) {
		if (exponent.IsZero()) return QuadNaN;
		return QuadZero;
	}
	if (exponent.IsZero()) {
		if (base.IsInfinite()) return QuadNaN;
		return QuadOne;
	}
	if (base.IsInfinite()) {
		if (exponent.IsInfinite()) {
			if (exponent.GetSign()) return QuadZero;
			return base;
		}
		if ((((__float128)-1) ^ exponent).GetSign()) return QuadNegativeInfinity; else return QuadPositiveInfinity;
	}
	if (exponent.IsInfinite()) {
		if (exponent.GetSign()) return QuadZero;
		if (base.GetSign()) return QuadNegativeInfinity; else return QuadPositiveInfinity;
	}
	if (base == QuadOne) return QuadOne;
	__float128 temp = Log2(base);
	Mul(temp, exponent, temp);
	return Base2Exp(temp);
}
void __float128::Abs( const __float128 &value, __float128 &result) {
	result = value;
	result.SetSign(false);
}
void __float128::Max( const __float128 &left, const __float128 &right, __float128 &result) {
	result = left;
	if (!(right < left)) result = right; //by doing !(right < left) we also check if right is NaN
}
void __float128::Min( const __float128 &left, const __float128 &right, __float128 &result) {
	result = left;
	if (!(right > left)) result = right; //by doing !(right > left) we also check if right is NaN
}
__float128 __float128::ModF(const __float128 &value, __float128 &integer) {
	__float128 fraction;
	int unbiasedExp = value.GetBase2Exponent();
	if (unbiasedExp < -1) {
		integer = QuadZero;
		fraction = value;
		fraction.SetSign(false);
		return fraction;
	} else if (unbiasedExp >= 112) {
		integer = QuadOne;
		CopySign(integer, value);
		fraction = QuadZero;
		return fraction;
	} else {
		integer = value;
		ClearBlock(integer.storage, 0, 112 - unbiasedExp);
		Fraction(value, fraction);
		return fraction;
	}
}
__float128 __float128::Sin(const __float128 &value) {
	if (value.IsNaN()) return value;
	//first wrap to -2Pi < x < 2Pi
	__float128 temp = value + QuadHalfPi;
	Div(temp, QuadTwoPi, temp);
	temp = Fraction(temp);
	if (temp.GetSign()) Add(temp, QuadOne, temp);
	bool negate = false;
	if (temp > QuadHalf) {
		negate = true;
		Sub(temp, QuadHalf, temp);
	}		
	Mul(temp, QuadTwoPi, temp);
	Sub(temp, QuadHalfPi, temp);
	if (temp.IsZero()) return QuadZero;

	__float128 currentX = temp;
	__float128 negVSquared;
	Mul(temp, temp, negVSquared);
	Negate(negVSquared);

	int iIteration = 1;
	__float128 factorial = QuadOne;
	__float128 result = QuadZero;
	bool affecting = true;
	while (affecting) {
		__float128 increment;
		__float128 factorialReciprocal = FactorialReciprocal(iIteration);
		iIteration += 2;
		Mul(currentX, factorialReciprocal, increment);
		if (!increment.IsZero() && result.GetBase2Exponent() - increment.GetBase2Exponent() <= QUAD_SIGNIFICANT_BITS) {
			Add(result, increment, result);
			Mul(currentX, negVSquared, currentX);
		} else {
			affecting = false;
		}
	}
	if (negate) Negate(result);
	return result;
}

__float128 __float128::Cos(const __float128 &value) {
	__float128 temp = value + QuadHalfPi;
	return Sin(temp);
}

void __float128::SinCos(const __float128 &value, __float128 &resultSin, __float128 &resultCos) {
	if (value.IsNaN()) {
		resultCos = resultSin = QuadNaN;
	}
	//first wrap to -2Pi < x < 2Pi
	__float128 tempSin = value + QuadHalfPi;
	__float128 tempCos = value + QuadPi;
	Div(tempSin, QuadTwoPi, tempSin);
	Div(tempCos, QuadTwoPi, tempCos);
	tempSin = Fraction(tempSin);
	tempCos = Fraction(tempCos);
	if (tempSin.GetSign()) Add(tempSin, QuadOne, tempSin);
	if (tempCos.GetSign()) Add(tempCos, QuadOne, tempCos);
	bool negateSin = false, negateCos = false;
	if (tempSin > QuadHalf) {
		negateSin = true;
		Sub(tempSin, QuadHalf, tempSin);
	}
	if (tempCos > QuadHalf) {
		negateCos = true;
		Sub(tempCos, QuadHalf, tempCos);
	}
	Mul(tempSin, QuadTwoPi, tempSin);
	Sub(tempSin, QuadHalfPi, tempSin);
	if (tempSin.IsZero()) {
		resultSin = QuadZero;
		resultCos = QuadOne;
		return;
	}
	Mul(tempCos, QuadTwoPi, tempCos);
	Sub(tempCos, QuadHalfPi, tempCos);

	__float128 currentXSin = tempSin;
	__float128 currentXCos = tempCos;
	__float128 negVSquaredSin;
	Mul(tempSin, tempSin, negVSquaredSin);
	Negate(negVSquaredSin);
	
	__float128 negVSquaredCos;
	Mul(tempCos, tempCos, negVSquaredCos);
	Negate(negVSquaredCos);

	int iIteration = 1;
	resultSin = QuadZero;
	resultCos = QuadZero;

	bool affecting = true;
	while (affecting) {
		//double dIncrement = dCurrentX / dFactorial;
		__float128 factorialReciprocal = FactorialReciprocal(iIteration);
		iIteration += 2;
		__float128 incrementSin;
		Mul(currentXSin, factorialReciprocal, incrementSin);
		__float128 incrementCos;
		Mul(currentXCos, factorialReciprocal, incrementCos);
		if ((!incrementSin.IsZero()
			&& resultSin.GetBase2Exponent() - incrementSin.GetBase2Exponent() <= QUAD_SIGNIFICANT_BITS) ||
			(!incrementCos.IsZero()
			&& resultCos.GetBase2Exponent() - incrementCos.GetBase2Exponent() <= QUAD_SIGNIFICANT_BITS)) {
			Add(resultSin, incrementSin, resultSin);
			Add(resultCos, incrementCos, resultCos);

			Mul(currentXSin, negVSquaredSin, currentXSin);
			Mul(currentXCos, negVSquaredCos, currentXCos);
		} else {
			affecting = false;
		}
	}
	if (negateSin) Negate(resultSin);
	if (negateCos) Negate(resultCos);
}

__float128 __float128::Tan(const __float128 &value) {
	__float128 sin, cos, tan;
	SinCos(value, sin, cos);
	Div(sin, cos, tan);
	return tan;
}

__float128 __float128::ASin(const __float128 &value) {
	if (value.IsNaN()) return value;
	if (value > QuadOne) return QuadNaN;
	if (value < QuadNegOne) return QuadNaN;
	if (Abs(value) > QuadSinQuarterPi) {
		__float128 temp;
		Mul(value, value, temp);
		//double dTemp = temp.ToDouble();
		bool negate = value.GetSign();
		Sub(QuadOne, temp, temp);
		//dTemp = temp.ToDouble();
		temp = Pow(temp, QuadHalf);
		//dTemp = temp.ToDouble();
		temp = ACos(temp);
		if (negate) Negate(temp);
		return temp;
	}
	__float128 result = QuadZero;
	__float128 Two = 2;
	__float128 Half = QuadOne / Two;
	__float128 Four = 4;
	__float128 i = 0;
	int iIteration = 0;
	__float128 i2 = 0;
	__float128 i4 = 0;
	__float128 fourPowITimesIFactorialSquared = 1;
	__float128 i2Factorial = 1;
	__float128 i2PlusOne = 1;
	__float128 vPowI2PlusOne = value;
	__float128 vSquared;
	Mul(value, value, vSquared);

	bool affecting = true;
	while (affecting) {
		__float128 increment;
		Mul(fourPowITimesIFactorialSquared, i2PlusOne, increment);
		Div(i2Factorial, increment, increment);
		Mul(increment, vPowI2PlusOne, increment);		
		//double dIncrement = increment.ToDouble();

		if (!increment.IsZero() && result.GetBase2Exponent() - increment.GetBase2Exponent() <= QUAD_SIGNIFICANT_BITS) {
			//dResult += dIncrement;
			Add(result, increment, result);

			//i++;
			Add(i, QuadOne, i);

			iIteration++;

			//i2 = 2*i (start with 0, and add two each iteration)
			//i2Factorial = (2i)! (start with one and multiply by i2-1 and i2 each iteration
			Add(i2, QuadOne, i2);
			Mul(i2Factorial, i2, i2Factorial);
			Add(i2, QuadOne, i2);
			Mul(i2Factorial, i2, i2Factorial);


			Mul(vPowI2PlusOne, vSquared, vPowI2PlusOne);
			//i2PlusOne = 2*i + 1 (start with one, and add two each iteration)
			Add(i2PlusOne, Two, i2PlusOne);

			//fourPowITimesIFactorialSquare = 4^i * (i!) ^ 2
			//start with one, and multiply by 4i * i each iteration
			__float128 temp;
			Add(i4, Four, i4);
			Mul(i4, i, temp);
			Mul(fourPowITimesIFactorialSquared, temp, fourPowITimesIFactorialSquared);
			//double dResult = result.ToDouble();
		} else {
			affecting = false;
		}
		if (iIteration > 14 && result.IsZero()) return result; //the increment magnitude check doesn't work with zero
	}
	return result;
}

__float128 __float128::ACos(const __float128 &value) {
	if (value.IsNaN()) return value;
	if (value > QuadOne) return QuadNaN;
	if (value < QuadNegOne) return QuadNaN;
	if (value.GetSign())
	{
		__float128 temp1;
		__float128 temp2;
		Negate(QuadHalfPi, temp1);
		temp2 = ASin(value);
		Sub(temp1, temp2, temp1);
		return temp1;
		//return (-QuadHalfPi) - ASin(value);
	}
	else
	{
		__float128 temp;
		temp = ASin(value);
		Sub(QuadHalfPi, temp, temp);
		return temp;
		//return QuadHalfPi - ASin(value);
	}
}

__float128 __float128::ATan(const __float128 &value) {
	if (value.IsNaN()) return value;
	if (value.IsInfinite()) {
		if (value.GetSign())
			return -QuadHalfPi;
		else
			return QuadHalfPi;
	}
	if (value.IsZero()) return QuadZero;
	if (value.GetSign()) {
		__float128 temp;
		Negate(value, temp);
		temp = ATan(temp);
		Negate(temp);
		return temp;
	}
	if (value > QuadOne) {
		__float128 temp;
		Div(QuadOne, value, temp);
		temp = ATan(temp);
		Sub(QuadHalfPi, temp, temp);
		return temp;
	}
	if (value > QuadHalf) {
		//convergence is very slow for values near one, so we are going to compute from arcsin
		//sin = tan / sqrt((1 + tan^2))
		__float128 temp;
		Mul(value, value, temp);
		Add(temp, QuadOne, temp);
		temp = Pow(temp, QuadHalf);
		__float128 sin;
		Div(value, temp, sin);
		return ASin(sin);
	}
	__float128 result = QuadZero;
	__float128 Two = 2;
	int iIteration = 0;
	__float128 i2PlusOne = 1;
	__float128 vPowI2PlusOne = value;
	__float128 vSquared;
	bool signFlip = false;
	Mul(value, value, vSquared);

	bool affecting = true;
	while (affecting) {
		__float128 increment;
		Div(vPowI2PlusOne, i2PlusOne, increment);

		if (!increment.IsZero() && result.GetBase2Exponent() - increment.GetBase2Exponent() <= QUAD_SIGNIFICANT_BITS) {
			//dResult += dIncrement;
			if (signFlip)
				Sub(result, increment, result);
			else
				Add(result, increment, result);

			signFlip = !signFlip;
			iIteration++;

			Mul(vPowI2PlusOne, vSquared, vPowI2PlusOne);
			//i2PlusOne = 2*i + 1 (start with one, and add two each iteration)
			Add(i2PlusOne, Two, i2PlusOne);

			//fourPowITimesIFactorialSquare = 4^i * (i!) ^ 2
			//start with one, and multiply by 4i * i each iteration
			//double dResult = result.ToDouble();
		} else {
			affecting = false;
		}
		if (iIteration > 14 && result.IsZero()) return result; //the increment magnitude check doesn't work with zero
	}
	return result;
}

__float128 __float128::ATan2(const __float128 &y, const __float128 &x ) {
	__float128 tan = y / x;
	__float128 partial = ATan(tan);
	if (x.GetSign())
	{
		__float128 temp = QuadPi;
		CopySign(temp, y);
		return partial + temp;
	}
	else
		return partial;
}

void __float128::SinhCosh(const __float128 &value, __float128 &resultSinh, __float128 &resultCosh) {
	__float128 z = value;
	resultCosh = QuadOne;
	resultSinh = QuadZero;
	Cordic(resultCosh, resultSinh, z, 1, 4, 52);
	Div(resultCosh, HyperbolicScaleFactor, resultCosh);
	Div(resultSinh, HyperbolicScaleFactor, resultSinh);
}

__float128 __float128::Tanh(const __float128 &value) {
	__float128 x, y;
	SinhCosh(value, y, x);
	Div(y, x, x);
	return x;
}

__float128 __float128::ATanh2(const __float128 &y, const __float128 &x) {
	__float128 resultX = x, resultY = y, z = 0;
	Cordic(resultX, resultY, z, 1, 4, 36);
	Div(resultX, HyperbolicScaleFactor, resultX);
	return resultX;
}

__float128 lanczosParameters[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

__float128 __float128::Gamma(const __float128 &value) { // Gamma(z+1)
	if (value < QuadHalf) {
		return QuadPi / (Sin(QuadPi*value) * Gamma(1 - value));
	} else {
		// a=z+g+.5
		__float128 a = value;
		Add(a, 24, a);
		Add(a, QuadHalf, a);
		// a^-a
		__float128 temp1 = a;
		Negate(temp1);
		temp1 = Pow(a, temp1);
		// a^(z+.5)
		__float128 temp2 = value;
		Add(temp2, QuadHalf, temp2);
		temp2 = Pow(a, temp2);
		// sqrt(2pi)
		__float128 temp3 = 2;
		Mul(temp3, QuadPi, temp3);
		temp3 = Pow(temp3, QuadHalf);
		// A_g(z)
		__float128 temp4 = QuadZero;
		for (int k = 0; k < 24; k++) {
			__float128 temp5;
			Add(value, k, temp5);
			Inc(temp5);
			Div(lanczosParameters[k], temp5, temp5);
			Add(temp4, temp5, temp4);
		}

		Mul(temp1, temp2, temp1);
		Mul(temp1, temp3, temp1);
		Mul(temp1, temp4, temp1);
		return temp1;
	}
}

void __float128::Fraction(const __float128 &value, __float128 &result) {
	int unbiasedExponent = value.GetBase2Exponent();
	if (unbiasedExponent < 0) 
		result = value;
	else {
		result = value;
		byte* data = result.storage;
		byte buffer[14];
		memset(buffer, 0, 14);
		BitBlockTransfer(data, 0, buffer, 0, 112);
		if (unbiasedExponent > 0) ClearBlock(buffer, 112 - unbiasedExponent, unbiasedExponent); //clear the integer portion
		int HOBit = ReverseBitScan((ui32*)buffer, 4, 111);
		if (HOBit == -1) { result = QuadZero; return; }
		unbiasedExponent = unbiasedExponent + HOBit - 112;
		ClearBlock(data, 0, 112);
		BitWindowTransfer(buffer, 0, 112, HOBit, data, 0, 112, 112); //HO bit is not copied
		result.SetBase2Exponent(unbiasedExponent);
	}
}

void __float128::Ceiling(const __float128 &value, __float128 &result) {
	int unbiasedExponent = value.GetBase2Exponent();
	if (unbiasedExponent == QUAD_EXPONENT_MAX) result = value; //NaN, +inf, -inf
	else if (unbiasedExponent < 0) {
		if (value.GetSign())
			result = FromData(zero);
		else
			result = FromData(one);
	} else {
		if (unbiasedExponent >= 112) { result = value; return; }
		result = value;
		int bitsToErase = 112 - unbiasedExponent;
		ClearBlock(result.storage, 0, bitsToErase);
		if (result != value) {
			if (!result.GetSign()) ++result;
		}
	}
}

__float128 __float128::Log(const __float128 &value, const __float128 &base) {
	__float128 temp1, temp2;
	temp1 = Log2(value);
	temp2 = Log2(base);
	Div(temp1, temp2, temp1);
	return temp1;
	//return Log2(value) / Log2(base);
}

void __float128::Round(const __float128 &value, __float128 &result) {
	int unbiasedExponent = value.GetBase2Exponent();
	if (unbiasedExponent == QUAD_EXPONENT_MAX) result = value; //NaN, +inf, -inf
	else if (unbiasedExponent < -1) result = QuadZero;
	else if (unbiasedExponent == -1) result = value.GetSign() ? QuadNegOne : QuadOne;
	else if (unbiasedExponent >= 112) result = value;
	else {
		result = value;
		bool roundUp = ReadBit(result.storage, 111 - unbiasedExponent);
		int bitsToErase = 112 - unbiasedExponent;
		ClearBlock(result.storage, 0, bitsToErase);
		if (roundUp) ++result;
	}
}
void __float128::Round(const __float128 &value, int precision, MidpointRoundingMode mode, __float128 &result) {
	const __float128 power10Table[37] = {
		Pow(10, 00), Pow(10, 01), Pow(10, 02), Pow(10, 03), Pow(10, 04), Pow(10, 05), Pow(10, 06), Pow(10, 07), Pow(10,  8), Pow(10,  9),
		Pow(10, 10), Pow(10, 11), Pow(10, 12), Pow(10, 13), Pow(10, 14), Pow(10, 15), Pow(10, 16), Pow(10, 17), Pow(10, 18), Pow(10, 19),
		Pow(10, 20), Pow(10, 21), Pow(10, 22), Pow(10, 23), Pow(10, 24), Pow(10, 25), Pow(10, 26), Pow(10, 27), Pow(10, 28), Pow(10, 29),
		Pow(10, 30), Pow(10, 31), Pow(10, 32), Pow(10, 33), Pow(10, 34), Pow(10, 35), Pow(10, 36) };
	const __float128 negPower10Table[37] = {
		Pow(10, -00), Pow(10, -01), Pow(10, -02), Pow(10, -03), Pow(10, -04), Pow(10, -05), Pow(10, -06), Pow(10, -07), Pow(10,  -8), Pow(10,  -9),
		Pow(10, -10), Pow(10, -11), Pow(10, -12), Pow(10, -13), Pow(10, -14), Pow(10, -15), Pow(10, -16), Pow(10, -17), Pow(10, -18), Pow(10, -19),
		Pow(10, -20), Pow(10, -21), Pow(10, -22), Pow(10, -23), Pow(10, -24), Pow(10, -25), Pow(10, -26), Pow(10, -27), Pow(10, -28), Pow(10, -29),
		Pow(10, -30), Pow(10, -31), Pow(10, -32), Pow(10, -33), Pow(10, -34), Pow(10, -35), Pow(10, -36) };
	if ((uint32_t)precision > 36)
		precision = 36;
	int unbiasedExponent = value.GetBase2Exponent();
	if (unbiasedExponent == QUAD_EXPONENT_MAX) result = value; //NaN, +inf, -inf
	else if (unbiasedExponent < -1) result = QuadZero;
	else if (unbiasedExponent == -1) {
		result = QuadOne;
		CopySign(result, value);
	}
	else if (unbiasedExponent >= 112) result = value;
	else {
		__float128 power10 = power10Table[precision];
		Mul(power10, value, power10);
		switch (mode) {
			default:
			case MidpointRoundingMode::ToEven:
				Round(power10, result);
				break;
			case MidpointRoundingMode::AwayFromZero: {
				__float128 fraction = ModF(value, result);
				if (Abs(fraction) >= QuadHalf) {
					Add(result, Sign(fraction), result);
				}
				break;
			}
			case MidpointRoundingMode::ToZero:
				ModF(value, result);
				break;
			case MidpointRoundingMode::ToNegativeInfinity:
				Floor(value, result);
				break;
			case MidpointRoundingMode::ToPositiveInfinity:
				Ceiling(value, result);
				break;
		}
		Mul(result, negPower10Table[precision], result);
	}
}

void __float128::Truncate(const __float128 &value, ui64 &result) {
	int exp = value.GetBase2Exponent();
	if (exp < 0) { result = 0; }
	else if (exp > 63) result = (ui64)0xffffffffffffffff;
	else {
		result = 0;
		BitWindowTransfer(value.storage, 0, 112, 112 - exp, &result, 0, 64, 0);
		SetBit(&result, exp);
	}
}

void __float128::Floor(const __float128 &value, __float128 &result) {
	int unbiasedExponent = value.GetBase2Exponent();
	if (unbiasedExponent == QUAD_EXPONENT_MAX) result = value; //NaN, +inf, -inf
	else if (unbiasedExponent < 0) {
		if (value.GetSign()) {
			result = FromData(one); Negate(result);
		}
		else {
			result = FromData(zero);
		}
	} else {
		if (unbiasedExponent >= 112) { result = value; return; }
		result = value;
		int bitsToErase = 112 - unbiasedExponent;
		ClearBlock(result.storage, 0, bitsToErase);
		if (result != value) {
			if (result.GetSign()) --result;
		}
	}
}

void __float128::CopySign(__float128 &value, const __float128 &sign) {
		byte* vPtr = value.storage;
		const byte* sPtr = sign.storage;
		*(vPtr + 15) = (*(sPtr + 15) & 0x80) | (*(vPtr + 15) & 0x7F);
}

void __float128::Cordic(__float128 &x, __float128 &y, __float128 &z, int n, int k, int l) {
	__float128 delta = QuadHalf;
	do {
		__float128 sign = QuadOne;
		CopySign(sign, z);
		// x += sign*delta*y
		__float128 tempX, tempY;
		Mul(delta, y, tempX);
		x += tempX;
		// y += sign*delta*x
		Mul(delta, x, tempY);
		y += tempY;

		// value -= sign*0.5*log((1+delta)/(1-delta))
		tempX = delta;
		Inc(tempX);
		Sub(QuadOne, delta, tempY);
		Div(tempX, tempY, tempX);
		Log(tempX, QuadE);
		Mul(QuadHalf, tempX, tempX);
		Mul(sign, tempX, tempX);
		z -= tempX;

		if (n == k) {
			k = 3*k + 1;
		} else {
			delta = delta >> 1;
			n++;
		}
	} while(n <= l);
}

#ifdef _MANAGED
#pragma managed
#endif
