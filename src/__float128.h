/*
Copyright 2010 Brent W. Lewis
coder0xff on hotmail

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

#include "Helpers.h"
#include <exception>

#pragma warning(disable: 4949)
#pragma unmanaged
#pragma warning(default: 4949)

#pragma region Exceptions
#define FPU_EXCEPTION_DECLARATION(x) \
class x##ExceptionClass : public std::exception \
{ \
 virtual const char* what() const throw(); \
} extern _##x##Exception; \
inline void x##() { if (enable##x##Exception) throw _##x##Exception; }

FPU_EXCEPTION_DECLARATION(Overflow)
FPU_EXCEPTION_DECLARATION(Underflow)
FPU_EXCEPTION_DECLARATION(DivideByZero)
FPU_EXCEPTION_DECLARATION(Invalid)
FPU_EXCEPTION_DECLARATION(Inexact)
#undef FPU_EXCEPTION_DECLARATION
#pragma endregion

struct __float128;
#define QUAD_CONSTANT(name, dataSource) extern __float128 Quad##name ;
#include "constants.h"
#undef QUAD_CONSTANT

struct __float128
{
private:
	byte storage[16];
	inline void SetBiasedExponent(ui16 value)
	{
		value &= QUAD_EXPONENT_MASK;
		ui16* ptr = (ui16*)this + 7;
		*ptr &= ~QUAD_EXPONENT_MASK;
		*ptr |= value;
	}
	inline ui16 GetBiasedExponent() const
	{
		return *((ui16*)this + 7) & QUAD_EXPONENT_MASK;
	}
	static void ReadOutResult(ui32 *buffer, int headScanBackStart, int biasedExponentAtScanStart, bool sign, __float128 &result);
	//the partial functions implement the Taylor series for the function
	//they are slow on large values, because they have no optimizations
	//but they are needed by the optimized functions - Ln, and Exp
	static __float128 PartialLn(const __float128 &value);
	static __float128 PartialExp(const __float128 &value);
	static __float128 Factorial(i32 value)
	{
		if (value < 0) return QuadNaN;
		if (value > MAX_FACTORIAL) return QuadPositiveInfinity;
		return factorials[value];
	}
	static __float128 FactorialReciprocal(i32 value)
	{
		if (value < 0) return QuadNaN;
		if (value > MAX_FACTORIAL) return QuadPositiveInfinity;
		return factorialReciprocals[value];
	}
	static void CordicTan(__float128 &x, __float128 &y, __float128 &z, __float128 delta, int n, int k);
	static void CordicArctan(__float128 &x, __float128 &y, __float128 &z, __float128 delta, int n, int k);
public:
	__float128();
	static inline __float128& FromData(const byte *data)
	{
		return *((__float128*)data);
	}
	inline bool GetSign() const
	{
		const byte* ptr = storage;
		return (*(ptr + 15) & 0x80) != 0;
	}
	inline void SetSign(bool value)
	{
		byte* ptr = storage;
		if (value)
			*(ptr + 15) |= 0x80;
		else
			*(ptr + 15) &= 0x7f;
	}
	inline i32 GetBase2Exponent() const {
		return (i32)GetBiasedExponent() - QUAD_EXPONENT_BIAS;
	}
	inline void SetBase2Exponent(i32 value) {
		if (value >= QUAD_EXPONENT_MAX) {
			*this = GetSign() ? QuadNegativeInfinity : QuadPositiveInfinity;
			return;
		}
		if (value < QUAD_EXPONENT_MIN) value = QUAD_EXPONENT_MIN;
		SetBiasedExponent((ui16)(value + QUAD_EXPONENT_BIAS));
	}
	inline bool IsZero() const {
		ui32* ptr = (ui32*)storage;
		if (*(ptr + 3) != 0)
			if (*(ptr + 3) != 0x80000000) return false;
		if (*ptr != 0) return false;
		ptr++;
		if (*ptr != 0) return false;
		ptr++;
		if (*ptr != 0) return false;
		return true;
	}
	inline bool IsNaN() const {
		if (GetBiasedExponent() != 0x7FFF) return false;
		i32* ptr = (i32*)storage;
		if (*ptr != 0) return true;
		ptr++;
		if (*ptr != 0) return true;
		ptr++;
		if (*ptr != 0) return true;
		ptr++;
		if (*(ui16*)ptr != 0) return true;
		return false;
	}
	inline bool IsInfinite() const {
		ui32* ptr = (ui32*)storage;
		if (*ptr != 0) return false;
		ptr++;
		if (*ptr != 0) return false;
		ptr++;
		if (*ptr != 0) return false;
		ptr++;
		if ((*ptr & 0x7FFFFFFF) != 0x7FFF0000) return false;
		return true;
	}
	inline bool IsSubNormal() const {
		return GetBiasedExponent() == 0 && !IsZero();
	}
	inline bool IsNormal() const {
		ui16 bExp = GetBiasedExponent();
		return bExp != 0x7FFF && bExp != 0;
	}
	inline static void Negate( __float128 &value)
	{
		*(value.storage + 15) ^= 0x80;
	}
	inline static void Negate(const __float128 &value, __float128 &result )
	{
		result = value;
		Negate(result);
	}
	__float128 operator-() const
	{
		__float128 copy = *this;
		Negate(copy);
		return copy;
	}
	inline static int Sign(const __float128 &value) {
		return value.IsZero() ? 0 : value.GetSign() ? -1 : 1;
	}
	static void Quantize(__float128 &value, const __float128 &exponent) {
		value.SetBase2Exponent(exponent.GetBase2Exponent());
	}
	static void Add(const __float128 &left, const __float128 &right, __float128 &result);
	static void Sub(const __float128 &left, const __float128 &right, __float128 &result);
	static void Mul(const __float128 &left, const __float128 &right, __float128 &result);
	static void Div(const __float128 &left, const __float128 &right, __float128 &result);
	static void FMulAdd(const __float128 &mpr, const __float128 &mpd, const __float128 &add, __float128 &result);
	inline static void Inc(__float128 &value) {
		Add(value, QuadOne, value);
	}
	inline static void Dec(__float128 &value) {
		__float128 temp = QuadOne;
		Sub(value, temp, value);
	}
	static bool Eql(const __float128 &left, const __float128 &right);
	inline bool Eql(const __float128 &other) const {
		return Eql(*this, other);
	}
	static bool EpsEql(const __float128 &left, const __float128 &right);
	inline bool EpsEql(const __float128 &other) const {
		return EpsEql(*this, other);
	}
	static i32 Cmp(const __float128 &left, const __float128 &right);
	inline i32 Cmp(const __float128 &other) const {
		return Cmp(*this, other);
	}
	__float128 operator+(const __float128 &other) const {
		__float128 result;
		Add(*this, other, result);
		return result;
	}
	inline __float128& operator+=(const __float128 &other) {
		Add(*this, other, *this);
		return *this;
	}
	inline __float128 operator-(const __float128 &other) const {
		__float128 result;
		Sub(*this, other, result);
		return result;
	}
	inline  __float128& operator-=(const __float128 &other) {
		Add(*this, other, *this);
		return *this;
	}
	inline __float128 operator*(const __float128 &other) const {
		__float128 result;
		Mul(*this, other, result);
		return result;
	}
	inline  __float128& operator*=(const __float128 &other) {
		Add(*this, other, *this);
		return *this;
	}
	inline __float128 operator/(const __float128 &other) const {
		__float128 result;
		Div(*this, other, result);
		return result;
	}
	inline  __float128& operator/=(const __float128 &other) {
		Div(*this, other, *this);
		return *this;
	}
	inline bool operator==(const __float128 &other) const {
		return Eql(*this, other);
	}
	inline bool operator!=(const __float128 &other) const {
		return !Eql(*this, other);
	}
	inline bool operator >(const __float128 &other) const {
		if ((*this).IsNaN() && other.IsNaN())
			return false;
		return Cmp(*this, other) > 0;
	}
	inline bool operator <(const __float128 &other) const {
		if ((*this).IsNaN() && other.IsNaN())
			return false;
		return Cmp(*this, other) < 0;
	}
	inline bool operator>=(const __float128 &other) const {
		if (Eql(*this, other)) return true; //ensure that equals test happens first because it's much faster
		else return *this > other;
	}
	inline bool operator<=(const __float128 &other) const {
		if (Eql(*this, other)) return true; //ensure that equals test happens first because it's much faster
		else return *this < other;
	}
	inline __float128& operator++() {
		__float128 temp = QuadOne;
		Add(temp, *this, *this);
		return *this;
	}
	inline __float128& operator--() {
		__float128 temp = QuadNegOne;
		Add(temp, *this, *this);
		return *this;
	}
	inline __float128 operator<<(i32 shift) const {
		__float128 temp = *this;
		temp.SetBase2Exponent(temp.GetBase2Exponent() + shift);
		return temp;
	}
	inline __float128 operator>>(i32 shift) const {
		__float128 temp = *this;
		temp.SetBase2Exponent(temp.GetBase2Exponent() - shift);
		return temp;
	}
	static inline __float128 Copy(const __float128 &value);
	__float128(double value);
	static void ToDouble(const __float128 &value, double &result);
	__float128(i64 value);
	static void ToInt64(const __float128 &value, i64 &result);
	__float128(i32 value);
	static void ToInt32(const __float128 &value, i32 &result);
	static void Abs(const __float128 &value, __float128 &result);
	static inline __float128 Abs(const __float128 &value) {
		__float128 result;
		Abs(value, result);
		return result;
	}
	static void Max(const __float128 &left, const __float128 &right, __float128 &result);
	static inline __float128 Max(const __float128 &left, const __float128 &right) {
		__float128 result;
		Max(left, right, result);
		return result;
	}
	static void Min(const __float128 &left, const __float128 &right, __float128 &result);
	static inline __float128 Min(const __float128 &left, const __float128 &right) {
		__float128 result;
		Min(left, right, result);
		return result;
	}
	static __float128 Ln(const __float128 &value);
	static __float128 Exp(const __float128 &value);
	static __float128 Base2Exp(i32 value);
	static __float128 Base2Exp(const __float128 &value);
	static __float128 Pow(const __float128 &base, const __float128 &exponent);
	inline __float128 operator^(const __float128 &right) const { return Pow(*this, right); }
	static __float128 Log(const __float128 &value, const __float128 &base);
	static __float128 Log2(const __float128 &value);
	static void Ceiling(const __float128 &value, __float128 &result);
	static inline __float128 Ceiling(const __float128 &value)
	{
		__float128 result;
		Ceiling(value, result);
		return result;
	}
	static void Floor(const __float128 &value, __float128 &result);
	static inline __float128 Floor(const __float128 &value)
	{
		__float128 result;
		Floor(value, result);
		return result;
	}
	static void Round(const __float128 &value, __float128 &result);
	static inline __float128 Round(const __float128 &value)
	{
		__float128 result;
		Round(value, result);
		return result;
	}
	enum MidpointRoundingMode : int {
        ToEven = 0, // Round
        AwayFromZero = 1,
        ToZero = 2, // Trunc
        ToNegativeInfinity = 3, // Floor
        ToPositiveInfinity = 4 // Ceiling
	};
	static void Round(const __float128 &value, int precision, MidpointRoundingMode mode, __float128 &result);
	static inline __float128 Round(const __float128 &value, int precision, MidpointRoundingMode mode) {
		__float128 result;
		Round(value, precision, mode, result);
		return result;
	}
	//result = Floor(Abs(value))
	static void Truncate(const __float128 &value, ui64 &result);
	static inline ui64 Truncate(const __float128 &value)
	{
		ui64 result;
		Truncate(value, result);
		return result;
	}
	//result = Abs(value) - Floor(Abs(value))
	static void Fraction(const __float128 &value, __float128 &result);
	static inline __float128 Fraction(const __float128 &value)
	{
		__float128 result;
		Fraction(value, result);
		return result;
	} 
	static __float128 ModF(const __float128 &value, __float128 &integer);
	static __float128 Sin(const __float128 &value);
	static __float128 Cos(const __float128 &value);
	static void SinCos(const __float128 &value, __float128 &resultSin, __float128 &resultCos);
	static __float128 Tan(const __float128 &value);
	static __float128 ASin(const __float128 &value);
	static __float128 ACos(const __float128 &value);
	static __float128 ATan(const __float128 &value);
	static __float128 ATan2(const __float128 &y, const __float128 &x);
	static void SinhCosh(const __float128 &value, __float128 &resultSinh, __float128 &resultCosh);
	static __float128 Tanh(const __float128 &value);
	static __float128 ATanh2(const __float128 &y, const __float128 &x);
	static inline void CopySign(__float128 &value, const __float128 &sign);
	static __float128 Gamma(const __float128 &value);
	static __float128 Factorial(const __float128 &value) { return Gamma(value + QuadOne); }
};

#pragma warning(disable: 4949)
#pragma managed
#pragma warning(default: 4949)