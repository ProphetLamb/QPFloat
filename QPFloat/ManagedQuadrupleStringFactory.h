#pragma once
#include "ManagedQuadruple.h"

using namespace System;
using namespace System::Text;

static String^ QuadrupleStringFactory(Quadruple value, String^ decimalSeperatorSign, String^ negativeSign, String^ positiveSign, String^ exponentSign
	, String^ infinity, String^ negativInfinity, String^ nan, wchar_t zero
	, String^ groupSeperator, int groupSize, int fixedDecimalDigits, int fixedLength
	, bool autoFixedNotation, bool forceFixedNotation)
{
	if (value.IsInfinite)
	{
		if (value.IsSigned) return negativInfinity;
		else return infinity;
	}
	if (value.IsNaN) return nan;

	StringBuilder^ result = gcnew StringBuilder();
	//if we are just appending zeros after the decimal, then put them in a temporary buffer
	StringBuilder^ nonZeroWaitCache = gcnew StringBuilder();
	bool sign = value.IsSigned;

	if (value.IsZero) //neg zero must be handled seperatly
	{
		result->Append(zero);
		if (fixedLength != -1)
		{
			result->Append(decimalSeperatorSign);
			for (int i = 0; i < fixedLength; i++) result->Append(zero);
		}
		else if (fixedDecimalDigits != -1)
		{
			result->Append(decimalSeperatorSign);
			for (int i = 0; i < fixedDecimalDigits; i++) result->Append(zero);
		}
		if (!forceFixedNotation && !autoFixedNotation)
		{
			result->Append(exponentSign);
			result->Append(positiveSign);
			result->Append(zero);
		}
		if (sign) result->Insert(0, negativeSign);
		return result->ToString();
	}

	int scientificExponent = 0;
	if (value.IsSubNormal)
	{
		value *= Quadruple::Pow(10, 4931);
		scientificExponent = -4931;
	}
	Quadruple::Abs(value, value);
	Quadruple ten = 10;
	Quadruple temp = Quadruple::Log(value, ten);
	int decimalDigitPlace = (int)Quadruple::Floor(temp);
	scientificExponent += decimalDigitPlace;
	Quadruple digitMultiplier = ten ^ decimalDigitPlace; //use a round to make sure it's accurate
	if (value / digitMultiplier >= ten)
	{
		//an error occurred calculating the logarithm - specifically numeric rounding error
		decimalDigitPlace++;
		scientificExponent++;
	}

	if (forceFixedNotation || autoFixedNotation && scientificExponent < 34 && scientificExponent > -5 && fixedLength > 0 && decimalDigitPlace < fixedLength) //to ensure fixed length use scientific notation if the decimal part is too long.
	{
		int maximumDigits; //enforce maximum length
		if (fixedLength > 0) maximumDigits = fixedLength; 
		else if (fixedDecimalDigits > -1) maximumDigits = decimalDigitPlace + fixedDecimalDigits + 1;
		else maximumDigits = 34;
		int displayedDigitPlaces = 0;
		bool decimalSeperatorQueued = false;
		int displayedDecimalDigitPlaces = 0;

		if (decimalDigitPlace < 0) decimalDigitPlace = 0; //start from the ones place at minimum
		digitMultiplier = ten ^ decimalDigitPlace;
		value = value / digitMultiplier;
		while ((value != Quadruple::Zero || decimalDigitPlace >= 0) && displayedDigitPlaces < maximumDigits)
		{
			int digitValue = (int)value;
			if (decimalDigitPlace == -1)
			{
				nonZeroWaitCache->Append(decimalSeperatorSign);
				decimalSeperatorQueued = true;
			}
			nonZeroWaitCache->Append(digitValue);
			if (decimalDigitPlace >= 0 || digitValue != 0)
			{
				if (decimalSeperatorQueued) displayedDecimalDigitPlaces+=nonZeroWaitCache->Length;
				result->Append(nonZeroWaitCache);
				nonZeroWaitCache->Clear();
			}
			value = Quadruple::Fraction(value);
			Quadruple::Mul(value, ten, value);
			decimalDigitPlace--;
			displayedDigitPlaces++;
		}
		//apply group seperators
		int digits = result->Length - displayedDecimalDigitPlaces;
		if (groupSize > 0 && digits > groupSize)
		{
			while ((digits -= groupSize) >= 1)
			{
				result->Insert(digits, groupSeperator);
			}
		}
		//increase length if nessessary
		int paddingZeros = 0;
		if (fixedLength > 0 && result->Length - displayedDecimalDigitPlaces <= fixedLength)
		{
			paddingZeros = fixedLength - result->Length - displayedDecimalDigitPlaces;
		}
		if (fixedDecimalDigits > 0 && displayedDecimalDigitPlaces < fixedDecimalDigits)
		{
			paddingZeros = fixedDecimalDigits - displayedDecimalDigitPlaces;
		}
		if (displayedDecimalDigitPlaces == 0 && paddingZeros > 0)
		{
			result->Append(decimalSeperatorSign);
		}
		while (paddingZeros-- > 0)
		{
			result->Append(zero);
		}
	}
	else
	{
		digitMultiplier = ten ^ decimalDigitPlace;
		value = value / digitMultiplier;
		//recursive call to fixed notation
		result->Append(QuadrupleStringFactory(value, decimalSeperatorSign, negativeSign, positiveSign, nullptr
			, nullptr, nullptr, nullptr, zero, groupSeperator, groupSize, fixedDecimalDigits, fixedLength, false, true));
		result->Append(exponentSign);
		if (scientificExponent >= 0)
		{
			result->Append(positiveSign);
			result->Append(scientificExponent);
		}
		else
		{
			result->Append(negativeSign);
			result->Append(-scientificExponent);
		}
	}
	if (sign) result->Insert(0, negativeSign);
	return result->ToString();
}
