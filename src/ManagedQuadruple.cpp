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

	You should have received a copy of the GNU General Public License
	along with QPFloat.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "stdafx.h"

#ifdef _MANAGED

#include <locale>
#include "ManagedQuadruple.h"
#include "DoubleDecomposition.h"
#include "__float128.h"

using namespace System::Collections::Generic;

namespace System
{
	array<byte>^ Quadruple::mantissa::get()
	{
		array<byte>^ results = gcnew array<byte>(QUAD_MANTISSA_BYTES);
		pin_ptr<byte> buf = &results[0];
		pin_ptr<byte> ptr = storage;
		BitBlockTransfer(ptr, 0, buf, QUAD_MANTISSA_BYTES * 8 - QUAD_SIGNIFICANT_BITS, QUAD_SIGNIFICANT_BITS); //place in H.O. bits
		return results;
	}
	void Quadruple::mantissa::set(array<byte>^ value)
	{
		int bitCount = value->Length * 8;
		pin_ptr<byte> buf = &value[0];
		pin_ptr<byte> ptr = storage;
		if (bitCount < QUAD_SIGNIFICANT_BITS) ClearBlock(ptr, 0, QUAD_SIGNIFICANT_BITS);
		BitWindowTransfer(buf, 0, bitCount, bitCount, ptr, 0, QUAD_SIGNIFICANT_BITS, QUAD_SIGNIFICANT_BITS); //read H.O. bits to H.O. bits
	}

	//reads bits from the buffer and stores the result in result, adjusting the exponent, performing round-to-even, and handling subnormals
	void Quadruple::ReadOutResult(ui32* buffer, int headScanBackStart, int biasedExponentAtScanStart, bool sign, Quadruple %result )
	{
		pin_ptr<byte> rPtr = result.storage;
		int implicitPosition = FindHeadAndApplyRounding(buffer, headScanBackStart);
		if (implicitPosition == -1)
		{
			//no bits
			result = 0;
			result.IsSigned = sign;
			return;
		}
		int currentExponent = biasedExponentAtScanStart + implicitPosition - headScanBackStart;
		if (currentExponent >= QUAD_EXPONENT_MASK)
		{
			if (sign)
				result = NegativeInfinity;
			else
				result = PositiveInfinity;
			Overflow();
			return;
		}
		if (currentExponent < 0) currentExponent = 0;
		int expInc = currentExponent - biasedExponentAtScanStart;
		result.biasedExponent = (ui16)currentExponent;
		BitBlockTransfer(buffer, headScanBackStart + expInc - 112, rPtr, 0, 112);
		result.IsSigned = sign;
		if (currentExponent == 0) Underflow();
		if (EnableInexactException) if (ReverseBitScan((ui32*)buffer, 0, headScanBackStart + expInc - 112 - 1) != -1) Inexact();
	}

	void Quadruple::WriteGroupedString(System::Text::StringBuilder^ builder, Quadruple q, String^ decimalSeparator, array<int>^ groupSizes, String^ groupSeparator, int decimalDigits)
	{
		//if we are just appending zeros after the decimal, then put them in a temporary buffer
		System::Text::StringBuilder^ nonZeroWaitCache = gcnew System::Text::StringBuilder();

		Quadruple currentValue = q;
		if (currentValue.IsSubNormal)
			currentValue *= Pow(10, 4931);
		Quadruple::Abs(currentValue, currentValue);
		Quadruple ten = 10;
		Quadruple temp = Log(currentValue, ten);
		int decimalDigitPlace = (int)Floor(temp);
		if (currentValue / (ten ^ decimalDigitPlace) >= ten)
		{
			//an error occurred calculating the logarithm - specifically numeric rounding error
			decimalDigitPlace++;
		}
		if (decimalDigitPlace < 0)
			decimalDigitPlace = 0; //start from the ones place at minimum
		int currentDecimalDigitPlace = decimalDigitPlace;
		currentValue = currentValue / (ten ^ currentDecimalDigitPlace);
		// Numeric groups
		int numericGroup = 0;
		int numericGroupIndex = 0;
		//it's not excessively large or small, so we'll display it without scientific notation
		while (currentValue != Quadruple::Zero || currentDecimalDigitPlace >= 0)
		{
			int digitValue = (int)currentValue;
			if (currentDecimalDigitPlace == -1) nonZeroWaitCache->Append(decimalSeparator);
			nonZeroWaitCache->Append(digitValue);
			if (currentDecimalDigitPlace >= 0 || digitValue != 0)
			{
				builder->Append(nonZeroWaitCache);
				nonZeroWaitCache->Clear();
			}
			if (decimalDigitPlace < 0)
			{
				if (numericGroupIndex == groupSizes[numericGroup]) // Move to next number group
				{
					nonZeroWaitCache->Append(groupSeparator);
					numericGroupIndex = 0;
					if (numericGroup == groupSizes->Length - 1) // Loop to first number group
						numericGroup = 0;
					else
						numericGroup++;
				}
				numericGroupIndex++;
			}
			currentValue = Quadruple::Fraction(currentValue);
			Mul(currentValue, ten, currentValue);
			currentDecimalDigitPlace--;
		}
		if (decimalDigitPlace < decimalDigits)
			builder->Append(gcnew String('0', decimalDigits - decimalDigitPlace));
	}
	void Quadruple::WriteSignedCurrencyString(System::Text::StringBuilder^ builder, Quadruple q, System::Globalization::NumberFormatInfo^ format, int precision)
	{
		if (q.IsSigned)
		{
			switch(format->CurrencyNegativePattern)
			{
				case 0:
					builder->Append("(")->Append(format->CurrencySymbol);
					break;
				case 1:
					builder->Append(format->NegativeSign)->Append(format->CurrencySymbol);
					break;
				case 2:
					builder->Append(format->CurrencySymbol)->Append(format->NegativeSign);
					break;
				case 3:
					builder->Append(format->CurrencySymbol);
					break;
				case 4:
					builder->Append("(");
					break;
				case 5:
					builder->Append(format->NegativeSign);
					break;
				case 8:
					builder->Append(format->NegativeSign);
					break;
				case 9:
					builder->Append(format->NegativeSign)->Append(format->CurrencySymbol)->Append(" ");
					break;
				case 11:
					builder->Append(format->CurrencySymbol)->Append(" ");
					break;
				case 12:
					builder->Append(format->CurrencySymbol)->Append(" ")->Append(format->NegativeSign);
					break;
				case 14:
					builder->Append("(")->Append(format->CurrencySymbol)->Append(" ");
					break;
				case 15:
					builder->Append("(");
					break;
			}
		}
		else
		{
			switch(format->CurrencyPositivePattern)
			{
				case 0:
					builder->Append(format->CurrencySymbol);
					break;
				case 2:
					builder->Append(format->CurrencySymbol)->Append(" ");
					break;
			}
		}
		WriteGroupedString(builder, q, format->CurrencyDecimalSeparator, format->CurrencyGroupSizes, format->CurrencyGroupSeparator, -1);
		if (q.IsSigned)
		{
			switch(format->CurrencyNegativePattern)
			{
				case 0:
					builder->Append(")");
					break;
				case 3:
					builder->Append(format->NegativeSign);
					break;
				case 4:
					builder->Append(format->CurrencySymbol)->Append(")");
					break;
				case 5:
					builder->Append(format->CurrencySymbol);
					break;
				case 6:
					builder->Append(format->NegativeSign)->Append(format->CurrencySymbol);
					break;
				case 7:
					builder->Append(format->CurrencySymbol)->Append(format->NegativeSign);
				case 8:
					builder->Append(" ")->Append(format->CurrencySymbol);
					break;
				case 10:
					builder->Append(" ")->Append(format->CurrencySymbol)->Append(format->NegativeSign);
					break;
				case 11:
					builder->Append(format->NegativeSign);
					break;
				case 13:
					builder->Append(format->NegativeSign)->Append(" ")->Append(format->CurrencySymbol);
					break;
				case 14:
					builder->Append(")");
					break;
				case 15:
					builder->Append(" ")->Append(format->CurrencySymbol)->Append(")");
					break;
			}
		}
		else
		{
			switch(format->CurrencyPositivePattern)
			{
				case 1:
					builder->Append(format->CurrencySymbol);
					break;
				case 3:
					builder->Append(" ")->Append(format->CurrencySymbol);
					break;
			}
		}
	}
	void Quadruple::WriteSingedDecimalString(System::Text::StringBuilder^ builder, Quadruple q, System::Globalization::NumberFormatInfo^ format, int precision)
	{
		//if we are just appending zeros after the decimal, then put them in a temporary buffer
		System::Text::StringBuilder^ nonZeroWaitCache = gcnew System::Text::StringBuilder();
		if (q.IsSigned)
			builder->Append(format->NegativeSign);
		
		Quadruple currentValue = q;
		if (currentValue.IsSubNormal)
			currentValue *= Pow(10, 4931);
		Quadruple::Abs(currentValue, currentValue);
		Quadruple ten = 10;
		Quadruple temp = Log(currentValue, ten);
		int decimalDigitPlace = (int)Floor(temp);
		if (currentValue / (ten ^ decimalDigitPlace) >= ten)
		{
			//an error occurred calculating the logarithm - specifically numeric rounding error
			decimalDigitPlace++;
		}
		if (decimalDigitPlace < 0)
			decimalDigitPlace = 0; //start from the ones place at minimum
		int currentDecimalDigitPlace = decimalDigitPlace;
		currentValue = currentValue / (ten ^ currentDecimalDigitPlace);
		//it's not excessively large or small, so we'll display it without scientific notation
		while (currentValue != Quadruple::Zero || currentDecimalDigitPlace >= 0)
		{
			int digitValue = (int)currentValue;
			if (currentDecimalDigitPlace == -1) nonZeroWaitCache->Append(format->NumberDecimalSeparator);
			nonZeroWaitCache->Append(digitValue);
			if (currentDecimalDigitPlace >= 0 || digitValue != 0)
			{
				builder->Append(nonZeroWaitCache);
				nonZeroWaitCache->Clear();
			}
			currentValue = Quadruple::Fraction(currentValue);
			Mul(currentValue, ten, currentValue);
			currentDecimalDigitPlace--;
		}
		int integerDigits = builder->Length - decimalDigitPlace != 0 ? 1 + decimalDigitPlace : 0;
		if (integerDigits < precision)
			builder->Insert(0, gcnew String(format->NativeDigits[0][0], precision - integerDigits));
	}
	void Quadruple::WriteExponentialString(System::Text::StringBuilder^ builder, Quadruple q, String^ exponentSign, System::Globalization::NumberFormatInfo^ format, int precision)
	{
		Quadruple currentValue = q;
		int scientificExponent = 0;
		if (currentValue.IsSubNormal)
		{
			currentValue *= Pow(10, 4931);
			scientificExponent = -4931;
		}
		Quadruple::Abs(currentValue, currentValue);
		Quadruple ten = 10;
		Quadruple temp = Log(currentValue, ten);
		int decimalDigitPlace = (int)Floor(temp);
		scientificExponent += decimalDigitPlace;
		Quadruple digitMultiplier = ten ^ decimalDigitPlace; //use a round to make sure it's accurate
		if (currentValue / digitMultiplier >= ten)
		{
			//an error occurred calculating the logarithm - specifically numeric rounding error
			decimalDigitPlace++;
			scientificExponent++;
		}
		int displayedDigitPlaces = 0;
		digitMultiplier = ten ^ decimalDigitPlace;
		currentValue = currentValue / digitMultiplier;
		WriteFixedString(builder, currentValue, format, precision); // Append significant digits
		builder->Append(exponentSign);
		if (scientificExponent > 0) builder->Append(format->PositiveSign); // Append exponent
		builder->Append(scientificExponent);
	}
	void Quadruple::WriteFixedString(System::Text::StringBuilder^ builder, Quadruple q, System::Globalization::NumberFormatInfo^ format, int precision)
	{
		//if we are just appending zeros after the decimal, then put them in a temporary buffer
		System::Text::StringBuilder^ nonZeroWaitCache = gcnew System::Text::StringBuilder();

		Quadruple currentValue = q;
		if (currentValue.IsSubNormal)
			currentValue *= Pow(10, 4931);
		Quadruple::Abs(currentValue, currentValue);
		Quadruple ten = 10;
		Quadruple temp = Log(currentValue, ten);
		int decimalDigitPlace = (int)Floor(temp);
		if (currentValue / (ten ^ decimalDigitPlace) >= ten)
		{
			//an error occurred calculating the logarithm - specifically numeric rounding error
			decimalDigitPlace++;
		}
		if (decimalDigitPlace < 0)
			decimalDigitPlace = 0; //start from the ones place at minimum
		int currentDecimalDigitPlace = decimalDigitPlace;
		currentValue = currentValue / (ten ^ currentDecimalDigitPlace);
		//it's not excessively large or small, so we'll display it without scientific notation
		while (currentValue != Quadruple::Zero || currentDecimalDigitPlace >= precision)
		{
			int digitValue = (int)currentValue;
			if (currentDecimalDigitPlace == -1) nonZeroWaitCache->Append(format->NumberDecimalSeparator);
			nonZeroWaitCache->Append(digitValue);
			if (currentDecimalDigitPlace >= 0 || digitValue != 0)
			{
				builder->Append(nonZeroWaitCache);
				nonZeroWaitCache->Clear();
			}
			currentValue = Quadruple::Fraction(currentValue);
			Mul(currentValue, ten, currentValue);
			currentDecimalDigitPlace--;
		}
		if (decimalDigitPlace < precision)
			builder->Append(gcnew String(format->NativeDigits[0][0], precision - decimalDigitPlace));
	}
	void Quadruple::WriteSignedGeneralString(System::Text::StringBuilder^ builder, Quadruple q, System::Globalization::NumberFormatInfo^ format, int precision)
	{
		//if we are just appending zeros after the decimal, then put them in a temporary buffer
		System::Text::StringBuilder^ nonZeroWaitCache = gcnew System::Text::StringBuilder();
		if (q.IsSigned)
			builder->Append(format->NegativeSign);

		Quadruple currentValue = q;
		int scientificExponent = 0;
		if (currentValue.IsSubNormal)
		{
			currentValue *= Pow(10, 4931);
			scientificExponent = -4931;
		}
		Quadruple::Abs(currentValue, currentValue);
		Quadruple ten = 10;
		Quadruple temp = Log(currentValue, ten);
		int decimalDigitPlace = (int)Floor(temp);
		scientificExponent += decimalDigitPlace;
		Quadruple digitMultiplier = ten ^ decimalDigitPlace; //use a round to make sure it's accurate
		if (currentValue / digitMultiplier >= ten)
		{
			//an error occurred calculating the logarithm - specifically numeric rounding error
			scientificExponent++;
		}
		int currentDecimalDigitPlace = decimalDigitPlace;
		if (scientificExponent < 40 && scientificExponent > -5)
		{
			//it's not excessively large or small, so we'll display it without scientific notation
			while (currentValue != Quadruple::Zero || currentDecimalDigitPlace >= 0)
			{
				int digitValue = (int)currentValue;
				if (currentDecimalDigitPlace == -1) nonZeroWaitCache->Append(format->NumberDecimalSeparator);
				nonZeroWaitCache->Append(digitValue);
				if (currentDecimalDigitPlace >= 0 || digitValue != 0)
				{
					builder->Append(nonZeroWaitCache);
					nonZeroWaitCache->Clear();
				}
				currentValue = Quadruple::Fraction(currentValue);
				Mul(currentValue, ten, currentValue);
				currentDecimalDigitPlace--;
			}
			if (decimalDigitPlace < precision)
				builder->Append(gcnew String(format->NativeDigits[0][0], precision - decimalDigitPlace));
		}
		else
		{
			WriteFixedString(builder, currentValue, format, 0); // Append significant digits
			builder->Append("e");
			if (scientificExponent > 0) builder->Append(format->PositiveSign); // Append exponent
			builder->Append(scientificExponent);
		}
	}
	void Quadruple::WriteSignedNumericString(System::Text::StringBuilder^ builder, Quadruple q, System::Globalization::NumberFormatInfo^ format, int precision)
	{
		if (q.IsSigned)
			builder->Append(format->NegativeSign);
		WriteGroupedString(builder, q, format->NumberDecimalSeparator, format->NumberGroupSizes, format->NumberGroupSeparator, -1);
	}
	void Quadruple::WriteSignedPrecentString(System::Text::StringBuilder^ builder, Quadruple q, System::Globalization::NumberFormatInfo^ format, int precision)
	{
		WriteGroupedString(builder, q*100, format->PercentDecimalSeparator, format->PercentGroupSizes, format->PercentGroupSeparator, -1);
	}
	void Quadruple::WriteSignedRoundTripString(System::Text::StringBuilder^ builder, Quadruple q, System::Globalization::NumberFormatInfo^ format)
	{
		if (q.IsSigned)
			builder->Append(format->NegativeSign);
		WriteExponentialString(builder, q, "e", format, 37);
	}

	void Quadruple::Add( Quadruple %a, Quadruple %b, Quadruple %result )
	{
		pin_ptr<byte> aPtr = a.storage;
		pin_ptr<byte> bPtr = b.storage;
		pin_ptr<byte> rPtr = result.storage;
		__float128::Add(*(__float128*)aPtr, *(__float128*)bPtr, *(__float128*)rPtr);
	}

	void Quadruple::Sub( Quadruple %a, Quadruple %b, Quadruple %result )
	{
		pin_ptr<byte> aPtr = a.storage;
		pin_ptr<byte> bPtr = b.storage;
		pin_ptr<byte> rPtr = result.storage;
		__float128::Sub(*(__float128*)aPtr, *(__float128*)bPtr, *(__float128*)rPtr);
	}

	void Quadruple::Mul( Quadruple %a, Quadruple %b, Quadruple %result )
	{
		pin_ptr<byte> aPtr = a.storage;
		pin_ptr<byte> bPtr = b.storage;
		pin_ptr<byte> rPtr = result.storage;
		__float128::Mul(*(__float128*)aPtr, *(__float128*)bPtr, *(__float128*)rPtr);
	}

	void Quadruple::Div( Quadruple %a, Quadruple %b, Quadruple %result )
	{
		//minimize managed to native transitions by just doing one
		pin_ptr<byte> aPtr = a.storage;
		pin_ptr<byte> bPtr = b.storage;
		pin_ptr<byte> rPtr = result.storage;
		__float128::Div(*(__float128*)aPtr, *(__float128*)bPtr, *(__float128*)rPtr);
	}

	bool Quadruple::operator==( Quadruple a, Quadruple b )
	{
		pin_ptr<byte> aPtr = a.storage;
		pin_ptr<byte> bPtr = b.storage;
		return (*(__float128*)aPtr) == (*(__float128*)bPtr);
	}

	bool Quadruple::operator!=( Quadruple a, Quadruple b )
	{
		pin_ptr<byte> aPtr = a.storage;
		pin_ptr<byte> bPtr = b.storage;
		return (*(__float128*)aPtr) != (*(__float128*)bPtr);
	}

	bool Quadruple::operator>( Quadruple a, Quadruple b )
	{
		pin_ptr<byte> aPtr = a.storage;
		pin_ptr<byte> bPtr = b.storage;
		return (*(__float128*)aPtr) > (*(__float128*)bPtr);
	}

	bool Quadruple::operator<( Quadruple a, Quadruple b )
	{
		pin_ptr<byte> aPtr = a.storage;
		pin_ptr<byte> bPtr = b.storage;
		return (*(__float128*)aPtr) < (*(__float128*)bPtr);
	}

	Quadruple::operator Quadruple( Double v )
	{
		__float128 temp = v;
		return *(Quadruple*)&temp;
	}

	Quadruple::operator Double(Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		double result;
		__float128::ToDouble(*(__float128*)vPtr, result);
		return result;
	}

	Quadruple::operator float(Quadruple v)
	{
		return (float)((Double)v);
	}

	Quadruple::operator Quadruple( i64 v )
	{
		__float128 temp = v;
		return *(Quadruple*)&temp;
	}

	Quadruple::operator i64(Quadruple v)
	{
		pin_ptr<byte> vPtr = v.storage;
		i64 result;
		__float128::ToInt64(*(__float128*)vPtr, result);
		return result;
	}

	Quadruple::operator Quadruple( i32 v )
	{
		__float128 temp = v;
		return *(Quadruple*)&temp;
	}

	Quadruple::operator i32(Quadruple v)
	{
		pin_ptr<byte> vPtr = v.storage;
		i32 result;
		__float128::ToInt32(*(__float128*)vPtr, result);
		return result;
	}

	Quadruple::operator Quadruple(__float128 v)
	{
		Quadruple result;
		pin_ptr<byte> resultPtr = result.storage;
		*(__float128*)resultPtr = v;
		return result;
	}

	Quadruple::operator __float128(Quadruple v)
	{
		__float128 result;
		pin_ptr<byte> vPtr = v.storage;
		result = *(__float128*)vPtr;
		return result;
	}

	void Quadruple::Abs( Quadruple %v, Quadruple %result )
	{
		result = v;
		result.IsSigned = false;
	}

	void Quadruple::Max( Quadruple %a, Quadruple %b, Quadruple% result)
	{
		result = a;
		if (!(b < a)) result = b; //by doing !(b < a) we also check if b is NaN
	}

	void Quadruple::Min( Quadruple %a, Quadruple %b, Quadruple% result)
	{
		result = a;
		if (!(b > a)) result = b; //by doing !(b > a) we also check if b is NaN
	}

	Quadruple Quadruple::Ln( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Ln(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::Base2Exp( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Base2Exp(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::Exp( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Exp(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::Pow( Quadruple base, Quadruple exponent )
	{
		pin_ptr<byte> basePtr = base.storage;
		pin_ptr<byte> expPtr = exponent.storage;
		__float128 result = __float128::Pow(*(__float128*)basePtr, *(__float128*)expPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::Log( Quadruple v, Quadruple base )
	{
		pin_ptr<byte> vPtr = v.storage;
		pin_ptr<byte> basePtr = base.storage;
		__float128 result = __float128::Log(*(__float128*)vPtr, *(__float128*)basePtr);
		return *(Quadruple*)&result;	
	}

	Quadruple Quadruple::Log2( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Log2(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	void Quadruple::Ceiling( Quadruple %v, Quadruple %result )
	{
		pin_ptr<byte> vPtr = v.storage;
		pin_ptr<byte> resultPtr = result.storage;
		__float128::Ceiling(*(__float128*)vPtr, *(__float128*)resultPtr);
	}

	void Quadruple::Floor( Quadruple %v, Quadruple %result )
	{
		pin_ptr<byte> vPtr = v.storage;
		pin_ptr<byte> resultPtr = result.storage;
		__float128::Floor(*(__float128*)vPtr, *(__float128*)resultPtr);
	}

	void Quadruple::Truncate( Quadruple %v, ui64 %result )
	{
		pin_ptr<byte> vPtr = v.storage;
		result = __float128::Truncate(*(__float128*)vPtr);
	}

	void Quadruple::Round( Quadruple %v, Quadruple %result )
	{
		pin_ptr<byte> vPtr = v.storage;
		pin_ptr<byte> resultPtr = result.storage;
		__float128::Round(*(__float128*)vPtr, *(__float128*)resultPtr);
	}

	void Quadruple::Fraction( Quadruple %v, Quadruple %result )
	{
		pin_ptr<byte> vPtr = v.storage;
		pin_ptr<byte> rPtr = result.storage;
		__float128::Fraction(*(__float128*)vPtr, *(__float128*)rPtr);
	}

	System::Quadruple Quadruple::Sin( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Sin(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	System::Quadruple Quadruple::Cos( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Cos(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	System::Quadruple Quadruple::Tan( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Tan(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	System::Quadruple Quadruple::ASin( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::ASin(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	System::Quadruple Quadruple::ACos( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::ACos(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	System::Quadruple Quadruple::ATan( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::ATan(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	System::Quadruple Quadruple::ATan2( Quadruple y, Quadruple x )
	{
		pin_ptr<byte> xPtr = x.storage;
		pin_ptr<byte> yPtr = y.storage;
		__float128 result = __float128::ATan2(*(__float128*)yPtr, *(__float128*)xPtr);
		return *(Quadruple*)&result;
	}

	String^ Quadruple::ToString()
	{
		return ToString("g");
	}

	String^ Quadruple::ToString( String^ format, IFormatProvider^ provider )
	{
		System::Text::StringBuilder^ builder = gcnew System::Text::StringBuilder();
		System::Globalization::NumberFormatInfo^ formatProvider = FormatProviderToNumberFormat(provider);
		int precision = format->Length > 1 ? Int32::Parse(format->Remove(0,1)) : -1;
		switch (tolower(format[0]))
		{
			case 'c':
				WriteSignedCurrencyString(builder, *this, formatProvider, precision);
				break;
			case 'd':
				WriteSingedDecimalString(builder, *this, formatProvider, precision);
				break;
			case 'e':
				if (this->IsSigned)
					builder->Append(formatProvider->NegativeSign);
				WriteExponentialString(builder, *this, gcnew String(format[0], 1), formatProvider, precision);
				break;
			case 'f':
				if (this->IsSigned)
					builder->Append(formatProvider->NegativeSign);
				WriteFixedString(builder, *this, formatProvider, precision);
				break;
			case 'g':
				WriteSignedGeneralString(builder, *this, formatProvider, precision);
				break;
			case 'n':
				WriteSignedNumericString(builder, *this, formatProvider, precision);
				break;
			case 'p':
				WriteSignedPrecentString(builder, *this, formatProvider, precision);
				break;
			case 'r':
				WriteSignedRoundTripString(builder, *this, formatProvider);
				break;
			default:
				throw gcnew System::FormatException();
		}
		return builder->ToString();
	}

	String^ Quadruple::ToString(String^ format)
	{
		return ToString(format, System::Globalization::CultureInfo::CurrentCulture->NumberFormat);
	}
	String^ Quadruple::ToString(IFormatProvider^ provider)
	{
		return ToString("g", provider);
	}

	Quadruple Quadruple::FromString(String^ str, IFormatProvider^ format)
	{
		System::Globalization::NumberFormatInfo^ provider = FormatProviderToNumberFormat(format);

		str = str->Trim();
		System::Text::StringBuilder^ builder = gcnew System::Text::StringBuilder(str);
		Quadruple result = 0;
		bool negative = false;
		if (builder->default[0] == provider->NegativeSign[0])
		{
			negative = true;
			builder->Remove(0, 1);
		}
		Quadruple ten = 10;
		int postDecimalDigits = 0;
		bool postDecimal = false;
		while (str->Length > 0)
		{
			Char digit = builder->default[0];
			builder->Remove(0, 1);
			if (digit >= '0' && digit <= '9')
			{
				int digitValue = (int)Char::GetNumericValue(digit);
				result *= ten;
				result += digitValue;
				if (postDecimal) postDecimalDigits++;
			}
			else if (digit == provider->NumberDecimalSeparator[0])
			{
				postDecimal = true;
			}
			else if (digit == 'e' || digit == 'E')
			{
				if (builder->default[0] == provider->PositiveSign[0]) builder->Remove(0, 1);
				postDecimalDigits -= System::Convert::ToInt32(builder->ToString());
				break;
			}
			else
			{
				throw gcnew FormatException();
			}
		}
		Quadruple temp = Quadruple::Pow(ten, -postDecimalDigits);
		Quadruple::Mul(result, temp, result);
		result.IsSigned = negative;
		return result;
	}

	int Quadruple::CompareTo(Quadruple other)	//	AK
	{
		pin_ptr<byte> aPtr = storage;
		pin_ptr<byte> bPtr = other.storage;
		return (*(__float128*)aPtr) < (*(__float128*)bPtr) ? -1
			: (*(__float128*)aPtr) > (*(__float128*)bPtr) ? +1
			: 0;
	}

	bool Quadruple::Equals(Quadruple other)
	{
		pin_ptr<byte> s = this->storage;
		pin_ptr<byte> o = other.storage;
		ui64* sPtr = (ui64*)s;
		ui64* oPtr = (ui64*)o;
		if (*(sPtr++) != *(oPtr++)) return false;
		if (*(sPtr) != *(oPtr)) return false;
		return true;
	}

}

#endif