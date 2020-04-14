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

#include <ctype.h>
#include "ManagedQuadruple.h"
#include "ManagedQuadrupleStringFactory.h"
#include "DoubleDecomposition.h"
#include "__float128.h"

using namespace System::Collections::Generic;
using namespace System::Globalization;

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

	Quadruple Quadruple::Sin( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Sin(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::Cos( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Cos(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::Tan( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Tan(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::ASin( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::ASin(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::ACos( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::ACos(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::ATan( Quadruple v )
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::ATan(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::ATan2( Quadruple y, Quadruple x )
	{
		pin_ptr<byte> xPtr = x.storage;
		pin_ptr<byte> yPtr = y.storage;
		__float128 result = __float128::ATan2(*(__float128*)yPtr, *(__float128*)xPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::Sinh(Quadruple v)
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Sinh(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::Cosh(Quadruple v)
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Cosh(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::Tanh(Quadruple v)
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Tanh(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::Coth(Quadruple v)
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::Coth(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::ASinh(Quadruple v)
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::ASinh(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::ACosh(Quadruple v)
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::ACosh(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::ATanh(Quadruple v)
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::ATanh(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}

	Quadruple Quadruple::ACoth(Quadruple v)
	{
		pin_ptr<byte> vPtr = v.storage;
		__float128 result = __float128::ACoth(*(__float128*)vPtr);
		return *(Quadruple*)&result;
	}
	
	String^ Quadruple::ToString()
	{
		NumberFormatInfo^ nfi = CultureInfo::CurrentCulture->NumberFormat;
		return ToString("G", nfi);
	}

	String^ Quadruple::ToString( String^ format, IFormatProvider^ provider )
	{
		//verfiy valid provider
		NumberFormatInfo^ nfi;
		if (provider->GetType()->IsAssignableFrom(CultureInfo::typeid)) nfi = ((CultureInfo^)provider)->NumberFormat;
		else if (provider->GetType()->IsAssignableFrom(NumberFormatInfo::typeid)) nfi = ((NumberFormatInfo^)provider);
		else throw gcnew ArgumentException("Invalid type of parameter 'provider'. Globalization::CultureInfo or Globalization::NumberFormatInfo must be assignable from the type of 'provider'.");

		wchar_t symbol;
		int precision = -1;
		if (format == nullptr || (format = format->Trim())->Length == 0) throw gcnew ArgumentNullException("format", "The parameter 'format' may not be null, and has to contain a format specifier symbol.");
		else if (format->Length == 1) symbol = format[0];
		else
		{
			symbol = format[0];
			precision = Convert::ToInt32(format->Remove(0, 1));
			if (precision < 0) throw gcnew FormatException("The precision specifier of the format string may not be negativ.");
		}
		symbol = toupper(symbol);

		String^ result;
		switch (symbol)
		{
		case 'B': //bytes in decimal
			int i;
			for (i = 0; i < 16; i++)
			{
				result += storage[i] + ",";
			}
			while (i < precision) //left pad bytes
			{
				result->Insert(0, "0,");
			}
			return result->Remove(result->Length - 1, 1); //remove tailing comma

		case 'C': //decimal digits provided by provider ; fixed notation
			result = QuadrupleStringFactory(Abs(*this), nfi->CurrencyDecimalSeparator, String::Empty, String::Empty, nullptr
				, nfi->PositiveInfinitySymbol, nfi->NegativeInfinitySymbol, nfi->NaNSymbol, nfi->NativeDigits[0][0]
				, nfi->CurrencyGroupSeparator, nfi->CurrencyGroupSizes[0], precision == -1 ? nfi->CurrencyDecimalDigits : precision, -1, false, true);
			if (IsSigned)
			{
				String^ currencyNegativFormat;
				switch (nfi->CurrencyNegativePattern) //{0}:n ; {1}:$
				{
				case 0: currencyNegativFormat = "({1}{0})"; break;
				case 4: currencyNegativFormat = "({0}{1})"; break;
				case 14: currencyNegativFormat = "({1} {0})"; break;
				case 15: currencyNegativFormat = "({0} {1})"; break;
				}
				if (currencyNegativFormat != nullptr) return String::Format(currencyNegativFormat, result, nfi->CurrencySymbol);
				switch (nfi->CurrencyNegativePattern) //{0}:n ; {1}:$ ; {2}:-
				{
				case 1: currencyNegativFormat = "{2}{1}{0}"; break;
				case 2: currencyNegativFormat = "{1}{2}{0}"; break;
				case 3: currencyNegativFormat = "{1}{0}{2}"; break;
				case 5: currencyNegativFormat = "{2}{0}{1}"; break;
				case 6: currencyNegativFormat = "{0}{2}{1}"; break;
				case 7: currencyNegativFormat = "{0}{1}{2}"; break;
				case 8: currencyNegativFormat = "{2}{0} {1}"; break;
				case 9: currencyNegativFormat = "{2}{1} {0}"; break;
				case 10: currencyNegativFormat = "{0} {1}{2}"; break;
				case 11: currencyNegativFormat = "{1} {0}{2}"; break;
				case 12: currencyNegativFormat = "{1} {2}{0}"; break;
				case 13: currencyNegativFormat = "{0}{2} {1}"; break;
				}
				if (currencyNegativFormat != nullptr) return String::Format(currencyNegativFormat, result, nfi->CurrencySymbol, nfi->NegativeSign);
				else throw gcnew ArgumentException("Invalid value '" + nfi->CurrencyNegativePattern + "' for CurrencyNegativePattern provided by 'provider'. Valid are values from 0 to 15.");
			}
			else
			{
				String^ currencyPositivFormat;
				switch (nfi->CurrencyPositivePattern) //{0}:n ; {1}:$
				{
				case 0: currencyPositivFormat = "{1}{0}"; break;
				case 1: currencyPositivFormat = "{0}{1}"; break;
				case 2: currencyPositivFormat = "{1} {0}"; break;
				case 3: currencyPositivFormat = "{0} {1}"; break;
				default: throw gcnew ArgumentException("Invalid value '" + nfi->CurrencyPositivePattern + "' for CurrencyPositivePattern provided by 'provider'. Valid are values from 0 to 3.");
				}
				return String::Format(currencyPositivFormat, result, nfi->CurrencySymbol);
			}
			break;

		case 'E': //no groups ; scientific notation
			return QuadrupleStringFactory(*this, nfi->NumberDecimalSeparator, nfi->NegativeSign, nfi->PositiveSign, format->Substring(0, 1)
				, nfi->PositiveInfinitySymbol, nfi->NegativeInfinitySymbol, nfi->NaNSymbol, nfi->NativeDigits[0][0]
				, nullptr, -1, precision, -1, false, false);

		case 'F': //no groups ; default decimal digits or digits in format string ; fixed notation
			return QuadrupleStringFactory(*this, nfi->NumberDecimalSeparator, nfi->NegativeSign, nfi->PositiveSign, nullptr
				, nfi->PositiveInfinitySymbol, nfi->NegativeInfinitySymbol, nfi->NaNSymbol, nfi->NativeDigits[0][0]
				, nullptr, -1, precision == -1 ? nfi->NumberDecimalDigits : precision, -1, false, true);

		case 'G': //lenght provided by format ; no groups ; auto notation
			return QuadrupleStringFactory(*this, nfi->NumberDecimalSeparator, nfi->NegativeSign, nfi->PositiveSign, "E"
				, nfi->PositiveInfinitySymbol, nfi->NegativeInfinitySymbol, nfi->NaNSymbol, nfi->NativeDigits[0][0]
				, nullptr, -1, -1, precision, true, false);

		case 'N': //no decimal places ; groups ; fixed notation
			result = QuadrupleStringFactory(Abs(*this), nfi->NumberDecimalSeparator, nfi->NegativeSign, nfi->PositiveSign, nullptr
				, nfi->PositiveInfinitySymbol, nfi->NegativeInfinitySymbol, nfi->NaNSymbol, nfi->NativeDigits[0][0]
				, nfi->NumberGroupSeparator, nfi->NumberGroupSizes[0], precision == -1 ? nfi->NumberDecimalDigits : precision, -1, false, true);
			if (IsSigned)
			{
				String^ numberNegativFormat;
				switch (nfi->NumberNegativePattern) //{0}:n
				{
				case 0: numberNegativFormat = "({0})"; break;
				case 1: numberNegativFormat = "-{0}"; break;
				case 2: numberNegativFormat = "- {0}"; break;
				case 3: numberNegativFormat = "{0}-"; break;
				case 4: numberNegativFormat = "{0} -"; break;
				default: throw gcnew ArgumentException("Invalid value '" + nfi->NumberNegativePattern + "' for CurrencyNegativePattern provided by 'provider'. Valid are values from 0 to 4.");
				}
				return String::Format(numberNegativFormat, result);
			}
			else return result;

		case 'P': //2 decimal digits ; groups ; 2 decimal places ; fixed notation
			result = QuadrupleStringFactory(Abs(*this), nfi->PercentDecimalSeparator, String::Empty, String::Empty, nullptr
				, nfi->PositiveInfinitySymbol, nfi->NegativeInfinitySymbol, nfi->NaNSymbol, nfi->NativeDigits[0][0]
				, nfi->PercentGroupSeparator, nfi->PercentGroupSizes[0], precision == -1 ? nfi->PercentDecimalDigits : precision, -1, false, true);
			if (IsSigned)
			{
				String^ percentNegativFormat;
				switch (nfi->PercentNegativePattern) //{0}:n ; {1}:%
				{
				case 0: percentNegativFormat = "-{0} {1}"; break;
				case 1: percentNegativFormat = "-{0}{1}"; break;
				case 2: percentNegativFormat = "-{1}{0}"; break;
				case 3: percentNegativFormat = "{1}-{0}"; break;
				case 4: percentNegativFormat = "{1}{0}-"; break;
				case 5: percentNegativFormat = "{0}-{1}"; break;
				case 6: percentNegativFormat = "{0}{1}-"; break;
				case 7: percentNegativFormat = "-{1} {0}"; break;
				case 8: percentNegativFormat = "{0} {1}-"; break;
				case 9: percentNegativFormat = "{1} {0}-"; break;
				case 10: percentNegativFormat = "{1} -{0}"; break;
				case 11: percentNegativFormat = "{0}- {1}"; break;
				default: throw gcnew ArgumentException("Invalid value '" + nfi->PercentNegativePattern + "' for PercentNegativePattern provided by 'provider'. Valid are values from 0 to 11.");
				}
				return String::Format(percentNegativFormat, result, nfi->PercentSymbol);
			}
			else
			{
				String^ percentPositivFormat;
				switch (nfi->PercentPositivePattern) //{0}:n ; {1}:%
				{
				case 0: percentPositivFormat = "{0} {1}"; break;
				case 1: percentPositivFormat = "{0}{1}"; break;
				case 2: percentPositivFormat = "{1}{0}"; break;
				case 3: percentPositivFormat = "{1} {0}"; break;
				default: throw gcnew ArgumentException("Invalid value '" + nfi->PercentPositivePattern + "' for PercentNegativePattern provided by 'provider'. Valid are values from 0 to 3.");
				}
				return String::Format(percentPositivFormat, result, nfi->PercentSymbol);
			}
			return result;
		case 'R':
			//serialize with default precision
			result = QuadrupleStringFactory(*this, nfi->NumberDecimalSeparator, nfi->NegativeSign, nfi->PositiveSign, "E"
				, nfi->PositiveInfinitySymbol, nfi->NegativeInfinitySymbol, nfi->NaNSymbol, nfi->NativeDigits[0][0]
				, nullptr, -1, 34, -1, true, false);
			if (IsInfinite || IsNaN)
				return result;
			if (Quadruple::FromString(result) != *this) //if parsing didnot reconstructed accurately
			{
				//increased # of decimal digits
				result = QuadrupleStringFactory(*this, nfi->NumberDecimalSeparator, nfi->NegativeSign, nfi->PositiveSign, "E"
					, nfi->PositiveInfinitySymbol, nfi->NegativeInfinitySymbol, nfi->NaNSymbol, nfi->NativeDigits[0][0]
					, nullptr, -1, 36, -1, true, false);
			}
			return result;
		default: throw gcnew FormatException("Unknown format specifier '" + (wchar_t)format[0] + "'. Supported are 'B', 'C', 'E', 'F', 'G', 'N', 'P', 'R'.");
		}
	}

	String^ Quadruple::ToString( String^ format)
	{
		return ToString(format, CultureInfo::CurrentCulture->NumberFormat);
	}

	String^ Quadruple::ToString(IFormatProvider^ provider )
	{
		return ToString("G", provider);
	}

	Quadruple Quadruple::FromString( String^ str )
	{
		str = str->Trim();
		if (str->StartsWith("0b", true, CultureInfo::InvariantCulture))
		{
			if (str->Length != 130)
				throw gcnew ArgumentException("The binary string 'str' must provide exactly 128bits.");
		}
		else
		{
			System::Text::StringBuilder^ s = gcnew System::Text::StringBuilder(str);
			Quadruple result = 0;
			bool negative = false;
			if (s->default[0] == '-')
			{
				negative = true;
				s->Remove(0, 1);
			}
			Quadruple ten = 10;
			int postDecimalDigits = 0;
			bool postDecimal = false;
			while (s->Length > 0)
			{
				Char digit = s->default[0];
				s->Remove(0, 1);
				if (digit >= '0' && digit <= '9')
				{
					int digitValue = (int)Char::GetNumericValue(digit);
					result *= ten;
					result += digitValue;
					if (postDecimal) postDecimalDigits++;
				}
				else if (digit == '.')
				{
					postDecimal = true;
				}
				else if (digit == 'e' || digit == 'E')
				{
					if (s->default[0] == '+') s->Remove(0, 1);
					postDecimalDigits -= Convert::ToInt32(s->ToString());
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
	}

	int Quadruple::CompareTo(Quadruple x)	//	AK
	{
		pin_ptr<byte> aPtr = storage;
		pin_ptr<byte> bPtr = x.storage;
		return (*(__float128*)aPtr) < (*(__float128*)bPtr) ? -1
			: (*(__float128*)aPtr) > (*(__float128*)bPtr) ? +1
			: 0;
	}
}

#endif