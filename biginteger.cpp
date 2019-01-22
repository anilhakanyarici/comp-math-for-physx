#include "biginteger.h"
#include <stdio.h>
#include "stdexcept"
#include <memory>
#include <algorithm>
#include <locale>

#include "drbg.h"

using namespace mp;

struct InternalBasicOperators
{
    static uint *add(uint *left, int leftLength, uint *right, int rightLength, int &resultLength)
    {
        uint *operands[2] { left, right };
        int lengths[2] { leftLength, rightLength };
        int lgtr = leftLength < rightLength;

        int max = lengths[lgtr];
        int min = lengths[!lgtr];
        uint *maxOperand = operands[lgtr];
        uint *minOperand = operands[!lgtr];

        uint *result = new uint[max + 1];
        ulong carry = 0;
        for (int i = 0; i < min; ++i){
            ulong sum = (ulong)minOperand[i] + (ulong)maxOperand[i] + carry;
            carry = sum >> 32;
            result[i] = (uint)(sum);
        }

        if (carry > 0){
            for (int i = min; i < max; ++i){
                ulong sum = (ulong)maxOperand[i] + carry;
                carry = sum >> 32;
                result[i] = (uint)(sum);
            }
            result[max] = (uint)(carry);
            resultLength = max + (carry > 0);
        }
        else{
            for (int i = min; i < max; ++i)
                result[i] = maxOperand[i];
            resultLength = max;
        }
        return result;
    }
    static uint *addSingle(uint *left, int leftLength, uint right, int &resultLength)
    {
        ulong sum = (ulong)left[0] + right;
        uint carry = sum >> 32;
        uint* result = new uint[leftLength + 1];
        result[0] = (uint)sum;
        for(int i = 1; i < leftLength; ++i){
            sum = (ulong)carry + left[i];
            result[i] = (uint)sum;
            carry = sum >> 32;
        }
        result[leftLength] = carry;
        resultLength = leftLength + (result[leftLength] != 0);
        return result;
    }
    static uint *divrem(uint *numerator, int numeratorLength, uint *denominator, int denominatorLength, int &quotientLength, int &remainderLength)
    {
        if(denominatorLength == 1){
            remainderLength = 1;
            uint divisor = denominator[0];
            uint *q = new uint[numeratorLength];
            q[numeratorLength - 1] = 0;

            ulong dividend = numerator[numeratorLength - 1];
            int qPos = numeratorLength - 1;
            int rPos = qPos;
            if (dividend >= divisor)
            {
                ulong quot = dividend / divisor;
                q[qPos--] = (uint)quot;
                numerator[rPos] = (uint)(dividend % divisor);
            }
            else
                qPos--;
            rPos--;
            while (rPos > -1)
            {
                int rPosPlusOne = rPos + 1;
                dividend = ((ulong)numerator[rPosPlusOne] << 32) | numerator[rPos];
                ulong quot = dividend / divisor;
                q[qPos--] = (uint)quot;
                numerator[rPosPlusOne] = 0;
                numerator[rPos--] = (uint)(dividend % divisor);
            }
            if (q[numeratorLength - 1] == 0)
                quotientLength = numeratorLength - 1;
            else
                quotientLength = numeratorLength;

            return q;
        } else {
            int numLastU = numeratorLength - 1;
            int opLDiff = numLastU - (denominatorLength - 1);
            quotientLength = opLDiff;
            for (int iu = numLastU; ; iu--)
            {
                if (iu < opLDiff)
                {
                    quotientLength = 1 + quotientLength;
                    break;
                }
                if (denominator[iu - opLDiff] != numerator[iu])
                {
                    if (denominator[iu - opLDiff] < numerator[iu])
                        quotientLength = 1 + quotientLength;
                    break;
                }
            }

            uint *quotient = new uint[quotientLength];

            uint denFirst = denominator[denominatorLength - 1];
            uint denSecond = denominator[denominatorLength - 2];
            int leftShiftBit = InternalBasicOperators::countOfZeroBitStart(denFirst);
            int rightShiftBit = 32 - leftShiftBit;
            if (leftShiftBit > 0)
            {
                denFirst = (denFirst << leftShiftBit) | (denSecond >> rightShiftBit);
                denSecond <<= leftShiftBit;
                if (denominatorLength > 2)
                    denSecond |= denominator[denominatorLength - 3] >> rightShiftBit;
            }

            for (int uInd = quotientLength; --uInd >= 0;)
            {
                uint hiNumDig = (uInd + denominatorLength <= numLastU) ? numerator[uInd + denominatorLength] : 0;

                ulong currNum = ((ulong)hiNumDig << 32) | numerator[uInd + denominatorLength - 1];
                uint nextNum = numerator[uInd + denominatorLength - 2];
                if (leftShiftBit > 0)
                {
                    currNum = (currNum << leftShiftBit) | (nextNum >> rightShiftBit);
                    nextNum <<= leftShiftBit;
                    if (uInd + denominatorLength >= 3)
                        nextNum |= numerator[uInd + denominatorLength - 3] >> rightShiftBit;
                }
                ulong rQuot = currNum / denFirst;
                ulong rRem = (uint)(currNum % denFirst);
                if (rQuot > 0xFFFFFFFF)
                {
                    rRem += denFirst * (rQuot - 0xFFFFFFFF);
                    rQuot = 0xFFFFFFFF;
                }
                while (rRem <= 0xFFFFFFFF && rQuot * denSecond > (((ulong)((uint)rRem) << 32) | nextNum))
                {
                    rQuot--;
                    rRem += denFirst;
                }
                if (rQuot > 0)
                {
                    ulong borrow = 0;
                    for (int u = 0; u < denominatorLength; ++u)
                    {
                        borrow += denominator[u] * rQuot;
                        uint uSub = (uint)borrow;
                        borrow >>= 32;
                        if (numerator[uInd + u] < uSub)
                            ++borrow;
                        numerator[uInd + u] -= uSub;
                    }
                    if (hiNumDig < borrow)
                    {
                        uint uCarry = 0;
                        for (int iu2 = 0; iu2 < denominatorLength; ++iu2)
                        {
                            uCarry = InternalBasicOperators::addCarry(&numerator[uInd + iu2], denominator[iu2], uCarry);
                        }
                        --rQuot;
                    }
                    numLastU = uInd + denominatorLength - 1;
                }
                quotient[uInd] = (uint)rQuot;
            }

            remainderLength = denominatorLength;

            while (numerator[remainderLength - 1] == 0 && remainderLength > 1)
                remainderLength = remainderLength - 1;
            while (quotient[quotientLength - 1] == 0 && quotientLength > 1)
                quotientLength = quotientLength - 1;
            return quotient;
        }
    }
    static void rem(uint *numerator, int numeratorLength, uint *denominator, int denominatorLength, int &remainderLength)
    {
        if(denominatorLength == 1){
            remainderLength = 1;
            uint divisor = denominator[0];

            ulong dividend = numerator[numeratorLength - 1];
            int rPos = numeratorLength - 1;
            if (dividend >= divisor)
                numerator[rPos] = (uint)(dividend % divisor);
            rPos--;
            while (rPos > -1)
            {
                int rPosPlusOne = rPos + 1;
                dividend = ((ulong)numerator[rPosPlusOne] << 32) | numerator[rPos];
                numerator[rPosPlusOne] = 0;
                numerator[rPos--] = (uint)(dividend % divisor);
            }
        } else {
            int numLastU = numeratorLength - 1;
            int opLDiff = numLastU - (denominatorLength - 1);
            int quotientLength = opLDiff;
            for (int iu = numLastU; ; iu--)
            {
                if (iu < opLDiff)
                {
                    quotientLength = 1 + quotientLength;
                    break;
                }
                if (denominator[iu - opLDiff] != numerator[iu])
                {
                    if (denominator[iu - opLDiff] < numerator[iu])
                        quotientLength = 1 + quotientLength;
                    break;
                }
            }

            uint denFirst = denominator[denominatorLength - 1];
            uint denSecond = denominator[denominatorLength - 2];
            int leftShiftBit = InternalBasicOperators::countOfZeroBitStart(denFirst);
            int rightShiftBit = 32 - leftShiftBit;
            if (leftShiftBit > 0)
            {
                denFirst = (denFirst << leftShiftBit) | (denSecond >> rightShiftBit);
                denSecond <<= leftShiftBit;
                if (denominatorLength > 2)
                    denSecond |= denominator[denominatorLength - 3] >> rightShiftBit;
            }

            for (int uInd = quotientLength; --uInd >= 0;)
            {
                uint hiNumDig = (uInd + denominatorLength <= numLastU) ? numerator[uInd + denominatorLength] : 0;

                ulong currNum = ((ulong)hiNumDig << 32) | numerator[uInd + denominatorLength - 1];
                uint nextNum = numerator[uInd + denominatorLength - 2];
                if (leftShiftBit > 0)
                {
                    currNum = (currNum << leftShiftBit) | (nextNum >> rightShiftBit);
                    nextNum <<= leftShiftBit;
                    if (uInd + denominatorLength >= 3)
                        nextNum |= numerator[uInd + denominatorLength - 3] >> rightShiftBit;
                }
                ulong rQuot = currNum / denFirst;
                ulong rRem = (uint)(currNum % denFirst);
                if (rQuot > 0xFFFFFFFF)
                {
                    rRem += denFirst * (rQuot - 0xFFFFFFFF);
                    rQuot = 0xFFFFFFFF;
                }
                while (rRem <= 0xFFFFFFFF && rQuot * denSecond > (((ulong)((uint)rRem) << 32) | nextNum))
                {
                    rQuot--;
                    rRem += denFirst;
                }
                if (rQuot > 0)
                {
                    ulong borrow = 0;
                    for (int u = 0; u < denominatorLength; ++u)
                    {
                        borrow += denominator[u] * rQuot;
                        uint uSub = (uint)borrow;
                        borrow >>= 32;
                        if (numerator[uInd + u] < uSub)
                            ++borrow;
                        numerator[uInd + u] -= uSub;
                    }
                    if (hiNumDig < borrow)
                    {
                        uint uCarry = 0;
                        for (int iu2 = 0; iu2 < denominatorLength; ++iu2)
                        {
                            uCarry = InternalBasicOperators::addCarry(&numerator[uInd + iu2], denominator[iu2], uCarry);
                        }
                        rQuot--;
                    }
                    numLastU = uInd + denominatorLength - 1;
                }
            }
            remainderLength = denominatorLength;
            while (numerator[remainderLength - 1] == 0 && remainderLength > 1)
                remainderLength = remainderLength - 1;
        }
    }
    static uint *mul(uint *left, int leftLength, uint *right, int rightLength, int &resultLength)
    {
        if (leftLength > rightLength)
        {
            int tmp = leftLength;
            leftLength = rightLength;
            rightLength = tmp;
            uint *tmpb = left;
            left = right;
            right = tmpb;
        }

        resultLength = leftLength + rightLength;
        uint *result = new uint[resultLength];

        for(int i = 0; i < resultLength; ++i)
            result[i] = 0;

        for (int i = 0; i < leftLength; ++i)
        {
            if (left[i] == 0) continue;

            ulong carry = 0;
            for (int j = 0, k = i; j < rightLength; ++j, ++k)
            {
                ulong val = ((ulong)left[i] * right[j]) + result[k] + carry;
                result[k] = (uint)val;
                carry = (val >> 32);
            }
            result[i + rightLength] = (uint)carry;
        }
        while (result[resultLength - 1] == 0 && resultLength > 1)
            resultLength = resultLength - 1;
        return result;
    }
    static uint *mulSingle(uint *left, int leftLength, uint right, int &resultLength)
    {
        resultLength = leftLength + 1;
        uint *result = new uint[resultLength];

        result[0] = 0;
        ulong carry;
        for (int i = 0; i < leftLength; ++i)
        {
            carry = 0;
            ulong val = ((ulong)left[i] * (ulong)right) + (ulong)result[i];
            result[i] = (uint)val;
            carry = (val >> 32);
            result[i + 1] = (uint)carry;
        }

        while (result[resultLength - 1] == 0 && resultLength > 1)
            resultLength = resultLength - 1;
        return result;
    }
    static uint *sub(uint *left, int leftLength, uint *right, int rightLength, int &resultLength)
    {
        uint *result = new uint[leftLength];
        for(int i = leftLength - rightLength; i < leftLength; ++i)
            result[i] = 0;

        int carry = 0;
        int i = 0;
        long diff;
        for ( ; i < rightLength; ++i)
        {
            diff = (long)left[i] - (long)right[i] + carry;
            result[i] = (uint)(diff);
            carry = (diff >> 63);
        }
        for ( ; carry && i < leftLength; ++i)
        {
            diff = (long)left[i] + carry;
            result[i] = (uint)(diff);
            carry = (diff >> 63);
        }
        for ( ; i < leftLength; ++i)
            result[i] = left[i];

        resultLength = leftLength;
        while (result[resultLength - 1] == 0 && resultLength > 1)
            resultLength = resultLength - 1;
        return result;
    }

    static uint *shiftRight(uint *digits, int digitLength, int shift, int &resultLength)
    {
        if (shift == 0) {
            uint *clone = new uint[digitLength];
            for(int i = 0; i < digitLength; ++i)
                clone[i] = digits[i];
            return clone;
        }

        int fullShift = shift >> 5;
        int remShift = shift & 31;

        int predictedLen = (digitLength) - fullShift;
        uint *result = new uint[predictedLen];
        for(int i = 0; i < predictedLen; ++i)
            result[i] = digits[fullShift + i];
        if(remShift > 0){
            result[0] = result[0] >> remShift;
            for(int i = 1; i < predictedLen; ++i){
                result[i - 1] |= result[i] << (32 - remShift);
                result[i] = result[i] >> remShift;
            }
        }
        if(result[predictedLen - 1] == 0 && predictedLen > 1)
            predictedLen--;
        resultLength = predictedLen;
        return result;
    }
    static uint *shiftLeft(uint *digits, int digitLength, int shift, int &resultLength)
    {
        if (shift == 0) {
            uint *clone = new uint[digitLength];
            for(int i = 0; i < digitLength; ++i)
                clone[i] = digits[i];
            return clone;
        }

        int fullShift = shift >> 5;
        int remShift = shift & 31;
        int needRemShift = remShift > 0;
        resultLength = digitLength + fullShift + needRemShift;
        uint *result = new uint[resultLength];
        result[resultLength - 1] = 0;

        for(int i = 0; i < digitLength; ++i)
            result[i + fullShift] = digits[i];
        if(remShift > 0){
            for(int i = resultLength - 1; i > fullShift; --i){
                uint dig = result[i] << remShift;
                dig |= result[i - 1] >> (32 - remShift);
                result[i] = dig;
            }
            result[fullShift] = result[fullShift] << remShift;
        }
        for(int i = 0; i < fullShift; ++i)
            result[i] = 0;
        if(result[resultLength - 1] == 0 && resultLength > 1)
            resultLength--;

        return result;
    }
    static uint *And(uint *left, int leftLength, uint *right, int rightLength, int &resultLength)
    {
        int min = leftLength > rightLength ? rightLength : leftLength;

        uint *result = new uint[min];
        for (int i = 0; i < min; ++i)
            result[i] = left[i] & right[i];

        resultLength = min;
        while (result[resultLength - 1] == 0 && resultLength > 1)
            resultLength = resultLength - 1;
        return result;
    }
    static uint *Or(uint *left, int leftLength, uint *right, int rightLength, int &resultLength)
    {
        uint *operands[2] { left, right };
        int lengths[2] { leftLength, rightLength };
        int lgtr = leftLength < rightLength;

        int max = lengths[lgtr];
        int min = lengths[!lgtr];
        uint *maxOperand = operands[lgtr];

        uint *result = new uint[max];
        int i = 0;
        for ( ; i < min; ++i)
            result[i] = left[i] | right[i];
        for ( ; i < max; ++i)
            result[i] = maxOperand[i];
        resultLength = max;
        return result;
    }
    static uint *Xor(uint *left, int leftLength, uint *right, int rightLength, int &resultLength)
    {
        uint *operands[2] { left, right };
        int lengths[2] { leftLength, rightLength };
        int lgtr = leftLength < rightLength;

        int max = lengths[lgtr];
        int min = lengths[!lgtr];
        uint *maxOperand = operands[lgtr];

        uint *result = new uint[max];
        int i = 0;
        for ( ; i < min; ++i)
            result[i] = left[i] ^ right[i];
        for ( ; i < max; ++i)
            result[i] = maxOperand[i] ^ 0;
        resultLength = max;
        while (result[resultLength - 1] == 0 && resultLength > 1)
            resultLength = resultLength - 1;
        return result;
    }
    static uint *Not(uint *digits, int digitLength, int &resultLength)
    {
        uint *result = new uint[digitLength];
        for(int i = 0; i < digitLength; ++i)
            result[i] = ~digits[i];
        resultLength = digitLength;
        while (result[resultLength - 1] == 0 && resultLength > 1)
            resultLength = resultLength - 1;
        return result;
    }

    static int compare(uint *left, int leftLength, uint *right, int rightLength)
    {
        int c = (int)(leftLength > rightLength);
        c = -(int)(leftLength < rightLength) * ((c + 1) & 1) + c;
        for(int i = leftLength - 1; !c && (i >= 0); --i)
            c = (int)(left[i] > right[i]) - (int)(left[i] < right[i]);
        return c;
    }
    static std::string toString(uint *digits, int digitLength, int sign)
    {
        if (sign == 0 || digitLength == 0) {
            return "0";
        }
        else if (digitLength == 1 && sign == 1){
            char *str = new char[20];
            sprintf(str, "%d", digits[0]);
            std::string u16str = str;
            delete[] str;
            return u16str;
        }

        const uint kuBase = 1000000000; // 10^9
        int cuMax = digitLength * 10 / 9 + 2;
        uint *rguDst = new uint[cuMax];
        int cuDst = 0;
        for(int i = 0; i < cuMax; ++i)
            rguDst[i] = 0;

        for (int iuSrc = digitLength; --iuSrc >= 0;){
            uint uCarry = digits[iuSrc];
            for (int iuDst = 0; iuDst < cuDst; ++iuDst){
                ulong uuRes = ((ulong)rguDst[iuDst] << 32) | uCarry;
                rguDst[iuDst] = (uint)(uuRes % kuBase);
                uCarry = (uint)(uuRes / kuBase);
            }
            if (uCarry != 0){
                rguDst[cuDst++] = uCarry % kuBase;
                uCarry /= kuBase;
                if (uCarry != 0)
                    rguDst[cuDst++] = uCarry;
            }
        }
        int cchMax = cuDst * 9;
        int rgchBufSize = cchMax + 2 + (sign == -1);
        int ichDst = cchMax;

        std::shared_ptr<char[]> rgchRef(new char[rgchBufSize]);
        char *rgch = rgchRef.get();

        for (int iuDst = 0; iuDst < cuDst - 1; ++iuDst){
            uint uDig = rguDst[iuDst];
            for (int cch = 9; --cch >= 0; ){
                int ascii = (48 + uDig % 10);
                rgch[--ichDst] = (char)ascii;
                uDig /= 10;
            }
        }
        for (uint uDig = rguDst[cuDst - 1]; uDig != 0; ){
            int ascii = (48 + uDig % 10);
            rgch[--ichDst] = (char)ascii;
            uDig /= 10;
        }
        delete[] rguDst;
        if (sign == -1){
            rgch[--ichDst] = '-';
        }
        rgch += ichDst;
        rgch[cchMax - ichDst] = '\0';
        return std::string(rgch);
    }
    static QString toQString(uint *digits, int digitLength, int sign)
    {
        if (sign == 0 || digitLength == 0) {
            return "0";
        }
        else if (digitLength == 1 && sign == 1){
            char *str = new char[20];
            sprintf(str, "%d", digits[0]);
            QString u16str = str;
            delete[] str;
            return u16str;
        }

        const uint kuBase = 1000000000; // 10^9
        int cuMax = digitLength * 10 / 9 + 2;
        uint *rguDst = new uint[cuMax];
        int cuDst = 0;
        for(int i = 0; i < cuMax; ++i)
            rguDst[i] = 0;

        for (int iuSrc = digitLength; --iuSrc >= 0;){
            uint uCarry = digits[iuSrc];
            for (int iuDst = 0; iuDst < cuDst; ++iuDst){
                ulong uuRes = ((ulong)rguDst[iuDst] << 32) | uCarry;
                rguDst[iuDst] = (uint)(uuRes % kuBase);
                uCarry = (uint)(uuRes / kuBase);
            }
            if (uCarry != 0){
                rguDst[cuDst++] = uCarry % kuBase;
                uCarry /= kuBase;
                if (uCarry != 0)
                    rguDst[cuDst++] = uCarry;
            }
        }
        int cchMax = cuDst * 9;
        int rgchBufSize = cchMax + 2 + (sign == -1);
        int ichDst = cchMax;

        std::shared_ptr<char[]> rgchRef(new char[rgchBufSize]);
        char *rgch = rgchRef.get();

        for (int iuDst = 0; iuDst < cuDst - 1; ++iuDst){
            uint uDig = rguDst[iuDst];
            for (int cch = 9; --cch >= 0; ){
                int ascii = (48 + uDig % 10);
                rgch[--ichDst] = (char)ascii;
                uDig /= 10;
            }
        }
        for (uint uDig = rguDst[cuDst - 1]; uDig != 0; ){
            int ascii = (48 + uDig % 10);
            rgch[--ichDst] = (char)ascii;
            uDig /= 10;
        }
        delete[] rguDst;
        if (sign == -1){
            rgch[--ichDst] = '-';
        }
        rgch += ichDst;
        rgch[cchMax - ichDst] = '\0';
        return QString(rgch);
    }

    static uint *unsignedParse(char *value, int charLen, int &digitLength)
    {
        int offset = charLen & 7;
        int base108Len = (charLen >> 3);
        uint lastDigit = 0;

        ++base108Len;
        for(int i = 0; i < offset; ++i){
            lastDigit *= 10;
            lastDigit += value[i] - 48;
        }
        uint *base108Digits = new uint[base108Len];
        int uiLast = base108Len - 1;
        base108Digits[uiLast] = lastDigit;
        for(int i = uiLast - 1; i >= 0; --i){
            uint curDigit = 0;
            offset += 8;
            for(int j = offset - 8; j < offset; ++j){
                curDigit *= 10;
                curDigit += value[j] - 48;
            }
            base108Digits[i] = curDigit;
        }

        const uint cBase = 100000000u;
        uint *digits = new uint[1];
        digits[0] = 0;
        int resultLength = 1;
        uint *tmp = digits;
        digits = InternalBasicOperators::addSingle(digits, resultLength, base108Digits[base108Len - 1], resultLength);
        delete[] tmp;
        for(int i = base108Len - 2; i >= 0; --i){
            tmp = digits;
            digits = InternalBasicOperators::mulSingle(digits, resultLength, cBase, resultLength);
            delete[] tmp;
            tmp = digits;
            uint base108dig = base108Digits[i];
            digits = InternalBasicOperators::addSingle(digits, resultLength, base108dig, resultLength);
            delete[] tmp;
        }
        delete[] base108Digits;
        digitLength = resultLength;
        return digits;
    }
    static long unsignedBitsLength(uint *digits, int digitLength)
    {
        uint uiLast = digits[digitLength - 1];
        return 32 * (long)digitLength - InternalBasicOperators::countOfZeroBitStart(uiLast);
    }

    static uint addCarry(uint *u1, uint u2, uint uCarry)
    {
        ulong uu = (ulong)*u1 + u2 + uCarry;
        *u1 = (uint)uu;
        return (uint)(uu >> 32);
    }
    static int countOfZeroBitStart(uint u)
    {
        int cbit = (u == 0);
        int f = (u & 0xFFFF0000) == 0;
        cbit += f << 4;
        f = ((u << cbit) & 0xFF000000) == 0;
        cbit += f << 3;
        f = ((u << cbit) & 0xF0000000) == 0;
        cbit += f << 2;
        f = ((u << cbit) & 0xC0000000) == 0;
        cbit += f << 1;
        f = ((u << cbit) & 0x80000000) == 0;
        cbit += f;
        return cbit;
    }
    static int getBit(uint *digits, int digitLength, long bitPosition) //Tested OK.
    {
        int digitPos = (int)(bitPosition / 32);
        if (digitLength <= digitPos)
            return 0;

        int smallBitPos = (int)(bitPosition & 31);
        return (int)((digits[digitPos] >> smallBitPos) & 1);
    }

    static long signedBitsLength(uint *digits, int digitLength, int sign) //Tested OK.
    {
        if (sign == 0)
            return 1;
        if (digitLength == 1 && digits[0] == 0)
            return 1;

        uint lastDigit = digits[digitLength - 1];
        unsigned char lastByte = 0;
        int bitsLength = digitLength * 32;

        if ((lastByte = (unsigned char)(lastDigit >> 24)) != 0) { }
        else if ((lastByte = (unsigned char)(lastDigit >> 16)) != 0) { bitsLength -= 8; }
        else if ((lastByte = (unsigned char)(lastDigit >> 8)) != 0) { bitsLength -= 16; }
        else if ((lastByte = (unsigned char)(lastDigit)) != 0) { bitsLength -= 24; }

        if ((lastByte >> 7) == 1 && sign == -1)
            bitsLength += 8;
        return bitsLength;
    }

    static uint *fromUnsignedBytes(unsigned char *data, int data_length, bool bigEndian, int &digitLength) //Tested OK.
    {
        digitLength = data_length/ 4;
        if ((data_length & 3) > 0)
            ++digitLength;

        uint *digits = new uint[digitLength];

        if (bigEndian)
        {
            int digitPos = digitLength - 1;
            int dataPos = 0;

            int nullDataLength = data_length & 3;
            if (nullDataLength == 1)
            {
                digits[digitPos--] = data[dataPos++];
            }
            else if (nullDataLength == 2)
            {
                uint digit = 0;
                digit |= (uint)(data[dataPos++] << 8);
                digit |= (uint)(data[dataPos++]);
                digits[digitPos--] = digit;
            }
            else if (nullDataLength == 3)
            {
                uint digit = 0;
                digit |= (uint)(data[dataPos++] << 16);
                digit |= (uint)(data[dataPos++] << 8);
                digit |= (uint)(data[dataPos++]);
                digits[digitPos--] = digit;
            }

            while (digitPos > -1)
            {
                uint current = 0;
                current |= (uint)(data[dataPos++] << 24);
                current |= (uint)(data[dataPos++] << 16);
                current |= (uint)(data[dataPos++] << 8);
                current |= (uint)(data[dataPos++]);
                digits[digitPos--] = current;
            }
        }
        else
        {
            int digitPos = 0;
            int dataPos = 0;
            int lastDigitPos = digitLength - 1;
            while (digitPos < lastDigitPos)
            {
                uint current = 0;
                current |= (uint)(data[dataPos++]);
                current |= (uint)(data[dataPos++] << 8);
                current |= (uint)(data[dataPos++] << 16);
                current |= (uint)(data[dataPos++] << 24);
                digits[digitPos++] = current;
            }

            int nullDataLength = data_length & 3;

            if (nullDataLength == 1)
            {
                digits[lastDigitPos] = data[dataPos++];
            }
            else if (nullDataLength == 2)
            {
                uint digit = 0;
                digit |= (uint)(data[dataPos++]);
                digit |= (uint)(data[dataPos++] << 8);
                digits[lastDigitPos] = digit;
            }
            else if (nullDataLength == 3)
            {
                uint digit = 0;
                digit |= (uint)(data[dataPos++]);
                digit |= (uint)(data[dataPos++] << 8);
                digit |= (uint)(data[dataPos++] << 16);
                digits[lastDigitPos] = digit;
            }
            else if (nullDataLength == 0)
            {
                uint digit = 0;
                digit |= (uint)(data[dataPos++]);
                digit |= (uint)(data[dataPos++] << 8);
                digit |= (uint)(data[dataPos++] << 16);
                digit |= (uint)(data[dataPos++] << 24);
                digits[lastDigitPos] = digit;
            }
        }
        InternalBasicOperators::trim(digits, digitLength);
        return digits;
    }

    static void trim(uint *digits, int &digitLength)
    {
        while (digits[digitLength - 1] == 0 && digitLength > 1)
            digitLength--;
    }
};

struct BigInteger::pimpl
{
    int _digit_length;
    uint *_digits;
    int _sign;
    bool _weak = false;

    ~pimpl()
    {
        if(this->_weak)
            return;
        delete[] this->_digits;
    }
};

BigInteger::BigInteger()
{
    this->_pimpl = std::shared_ptr<pimpl>(new pimpl);
    this->_pimpl->_sign = 0;
    this->_pimpl->_digits = new uint[1];
    this->_pimpl->_digits[0] = 0;
    this->_pimpl->_digit_length = 1;
}
BigInteger::BigInteger(int value)
{
    this->_pimpl = std::shared_ptr<pimpl>(new pimpl);
    this->_pimpl->_sign = (value >> 31);
    this->_pimpl->_sign = this->_pimpl->_sign + (int)(value > 0);
    value *= this->_pimpl->_sign;
    this->_pimpl->_digit_length = 1;
    this->_pimpl->_digits = new uint[1];
    this->_pimpl->_digits[0] = (uint)value;
}
BigInteger::BigInteger(uint value)
{
    this->_pimpl = std::shared_ptr<pimpl>(new pimpl);
    this->_pimpl->_sign = value > 0;
    this->_pimpl->_digit_length = 1;
    this->_pimpl->_digits = new uint[1];
    this->_pimpl->_digits[0] = (uint)value;
}
BigInteger::BigInteger(long value)
{
    this->_pimpl = std::shared_ptr<pimpl>(new pimpl);
    this->_pimpl->_sign = (value >> 63);
    this->_pimpl->_sign = this->_pimpl->_sign + (int)(value > 0);
    value *= this->_pimpl->_sign;
    int vgtuim = (int)(value > 0xffffffff);
    this->_pimpl->_digit_length = 1 + vgtuim;
    this->_pimpl->_digits = new uint[this->_pimpl->_digit_length];
    this->_pimpl->_digits[vgtuim] = (uint)(value >> 32);
    this->_pimpl->_digits[0] = (uint)value;
}
BigInteger::BigInteger(ulong value)
{
    this->_pimpl = std::shared_ptr<pimpl>(new pimpl);
    this->_pimpl->_sign = value > 0;
    int vgtuim = (int)(value > 0xffffffff);
    this->_pimpl->_digit_length = 1 + vgtuim;
    this->_pimpl->_digits = new uint[this->_pimpl->_digit_length];
    this->_pimpl->_digits[vgtuim] = (uint)(value >> 32);
    this->_pimpl->_digits[0] = (uint)value;
}

BigInteger::BigInteger(char *bytes, long bytes_length, int sign, bool big_endian)
{
    this->_pimpl = std::shared_ptr<pimpl>(new pimpl);
    this->_pimpl->_digits = InternalBasicOperators::fromUnsignedBytes(reinterpret_cast<unsigned char*>(bytes), bytes_length, big_endian, this->_pimpl->_digit_length);
    this->_pimpl->_sign = sign;
}

BigInteger::BigInteger(const QByteArray &bytes, int sign, bool big_endian)
{
    this->_pimpl = std::shared_ptr<pimpl>(new pimpl);
    this->_pimpl->_digits = InternalBasicOperators::fromUnsignedBytes(reinterpret_cast<unsigned char*>(const_cast<char*>(bytes.data())), bytes.size(), big_endian, this->_pimpl->_digit_length);
    this->_pimpl->_sign = sign;
}

//BigInteger::BigInteger(const BigInteger &value)
//{
//    this->_pimpl = std::shared_ptr<pimpl>(new pimpl);
//    this->_pimpl->_digit_length = value._pimpl->_digit_length;
//    this->_pimpl->_sign = value._pimpl->_sign;
//    this->_pimpl->_digits = new uint[this->_pimpl->_digit_length];
//    for(int i = 0; i < this->_pimpl->_digit_length; ++i){
//        this->_pimpl->_digits[i] = value._pimpl->_digits[i];
//    }
//}

//BigInteger &BigInteger::operator =(const BigInteger &value)
//{
//    this->_pimpl = std::shared_ptr<pimpl>(new pimpl);
//    this->_pimpl->_digit_length = value._pimpl->_digit_length;
//    this->_pimpl->_sign = value._pimpl->_sign;
//    this->_pimpl->_digits = new uint[this->_pimpl->_digit_length];
//    for(int i = 0; i < this->_pimpl->_digit_length; ++i){
//        this->_pimpl->_digits[i] = value._pimpl->_digits[i];
//    }
//    return *this;
//}

BigInteger BigInteger::twoPowerOf(const long &exponent)
{
    long digitlen = (exponent >> 5) + 1;
    uint *digits = new uint[digitlen];
    for(int i = 0; i < digitlen; ++i)
        digits[i] = 0;
    digits[digitlen - 1] = 1 << (exponent & 31);
    return BigInteger(digits, digitlen, 1, true);
}

bool BigInteger::isZero() const { return this->_pimpl->_sign == 0; }
bool BigInteger::isOne() const { return this->_pimpl->_digit_length == 1 && this->_pimpl->_digits[0] == 1 && this->_pimpl->_sign == 1; }
bool BigInteger::isMinusOne() const { return this->_pimpl->_digit_length == 1 && this->_pimpl->_digits[0] == 1 && this->_pimpl->_sign == -1; }
bool BigInteger::isNegative() const { return this->_pimpl->_sign == -1; }
bool BigInteger::isPositive() const { return this->_pimpl->_sign == 1; }
bool BigInteger::isPowerOfTwo() const
{
    int uiLast = this->_pimpl->_digit_length - 1;
    uint uLast = this->_pimpl->_digits[uiLast];

    if ((uLast & (uLast - 1)) != 0)
        return false;

    for(int i = 0; i < uiLast; ++i)
        if(this->_pimpl->_digits[i] != 0)
            return false;

    return true;
}
bool BigInteger::isEven() const { return (this->_pimpl->_digits[0] & 1) == 0; }
bool BigInteger::isOdd() const { return (this->_pimpl->_digits[0] & 1); }
int BigInteger::sign() const { return this->_pimpl->_sign; }
int BigInteger::digitLength() const { return this->_pimpl->_digit_length; }

long BigInteger::bytesLength() const
{
    long bitlen = this->bitsLength();
    return (bitlen >> 3) + (int)((bitlen & 7) > 0);
}
long BigInteger::bitsLength() const { return InternalBasicOperators::signedBitsLength(this->_pimpl->_digits, this->_pimpl->_digit_length, this->_pimpl->_sign); }
uint BigInteger::firstDigit() const { return this->_pimpl->_digits[0]; }
uint *BigInteger::data() const { return this->_pimpl->_digits; }

char *BigInteger::bytesData() const
{
    return reinterpret_cast<char *>(this->_pimpl->_digits);
}

BigInteger &BigInteger::add(const BigInteger &value)
{
    if(this->_pimpl->_sign == value._pimpl->_sign){
        if(this->_pimpl->_sign == 0){
            delete[] this->_pimpl->_digits;
            this->_pimpl->_digits = new uint[1] { 0 };
            return *this;
        }
        uint *digits = InternalBasicOperators::add(this->_pimpl->_digits, this->_pimpl->_digit_length, value._pimpl->_digits, value._pimpl->_digit_length, this->_pimpl->_digit_length);
        delete[] this->_pimpl->_digits;
        this->_pimpl->_digits = digits;
        return *this;
    } else {
        int c = InternalBasicOperators::compare(this->_pimpl->_digits, this->_pimpl->_digit_length, value._pimpl->_digits, value._pimpl->_digit_length);
        uint *digits;
        if(c == 1){
            digits = InternalBasicOperators::sub(this->_pimpl->_digits, this->_pimpl->_digit_length, value._pimpl->_digits, value._pimpl->_digit_length, this->_pimpl->_digit_length);
            delete[] this->_pimpl->_digits;
            this->_pimpl->_digits = digits;
        } else if(c == -1){
            digits = InternalBasicOperators::sub(value._pimpl->_digits, value._pimpl->_digit_length, this->_pimpl->_digits, this->_pimpl->_digit_length, this->_pimpl->_digit_length);
            delete[] this->_pimpl->_digits;
            this->_pimpl->_digits = digits;
            this->_pimpl->_sign = value._pimpl->_sign;
        } else {
            delete[] this->_pimpl->_digits;
            this->_pimpl->_digits = new uint[1] { 0 };
            this->_pimpl->_digit_length = 1;
            this->_pimpl->_sign = 0;
        }
        return *this;
    }
}

BigInteger &BigInteger::sub(const BigInteger &value)
{
    if(this->_pimpl->_sign == -value._pimpl->_sign){
        if(this->_pimpl->_sign == 0){
            delete[] this->_pimpl->_digits;
            this->_pimpl->_digits = new uint[1] { 0 };
            return *this;
        }
        uint *digits = InternalBasicOperators::add(this->_pimpl->_digits, this->_pimpl->_digit_length, value._pimpl->_digits, value._pimpl->_digit_length, this->_pimpl->_digit_length);
        delete[] this->_pimpl->_digits;
        this->_pimpl->_digits = digits;
        return *this;
    } else {
        int c = InternalBasicOperators::compare(this->_pimpl->_digits, this->_pimpl->_digit_length, value._pimpl->_digits, value._pimpl->_digit_length);
        uint *digits;
        if(c == 1){
            digits = InternalBasicOperators::sub(this->_pimpl->_digits, this->_pimpl->_digit_length, value._pimpl->_digits, value._pimpl->_digit_length, this->_pimpl->_digit_length);
            delete[] this->_pimpl->_digits;
            this->_pimpl->_digits = digits;
        } else if(c == -1){
            digits = InternalBasicOperators::sub(value._pimpl->_digits, value._pimpl->_digit_length, this->_pimpl->_digits, this->_pimpl->_digit_length, this->_pimpl->_digit_length);
            delete[] this->_pimpl->_digits;
            this->_pimpl->_digits = digits;
            this->_pimpl->_sign = -value._pimpl->_sign;
        } else {
            delete[] this->_pimpl->_digits;
            this->_pimpl->_digits = new uint[1] { 0 };
            this->_pimpl->_digit_length = 1;
            this->_pimpl->_sign = 0;
        }
        return *this;
    }
}

BigInteger &BigInteger::mul(const BigInteger &value)
{
    this->_pimpl->_sign = this->_pimpl->_sign * value._pimpl->_sign;
    if(this->_pimpl->_sign == 0){
        delete[] this->_pimpl->_digits;
        this->_pimpl->_digits = new uint[1] { 0 };
        this->_pimpl->_digit_length = 1;
        return *this;
    }
    uint *digits = InternalBasicOperators::mul(this->_pimpl->_digits, this->_pimpl->_digit_length, value._pimpl->_digits, value._pimpl->_digit_length, this->_pimpl->_digit_length);
    delete[] this->_pimpl->_digits;
    this->_pimpl->_digits = digits;
    return *this;
}

BigInteger &BigInteger::div(const long &divisor)
{
    long rem;
    return this->divrem(divisor, rem);
}
BigInteger &BigInteger::div(const ulong &divisor)
{
    ulong rem;
    return this->divrem(divisor, rem);
}
BigInteger &BigInteger::div(const BigInteger &divisor)
{
    if(divisor.isZero())
        throw std::runtime_error("Cannot divide a number by zero.");
    BigInteger rem;
    return this->divrem(divisor, rem);
}
BigInteger &BigInteger::rem(long modulus)
{
    if(modulus == 0)
        throw std::runtime_error("Cannot divide a number by zero.");
    if(this->_pimpl->_sign == 0)
        return *this;
    int modulus_sign = (modulus >> 63);
    modulus_sign = modulus_sign + (int)(modulus > 0);
    modulus *= modulus_sign;
    int vgtuim = (int)(modulus > 0xffffffff);
    int modulus_digitLength = 1 + vgtuim;
    uint *modulus_digits = new uint[modulus_digitLength];
    modulus_digits[vgtuim] = (uint)(modulus >> 32);
    modulus_digits[0] = (uint)modulus;
    this->_pimpl->_sign = this->_pimpl->_sign * modulus_sign;
    if(this->_pimpl->_digit_length >= modulus_digitLength)
        InternalBasicOperators::rem(this->_pimpl->_digits, this->_pimpl->_digit_length, modulus_digits, modulus_digitLength, this->_pimpl->_digit_length);
    return *this;
}
BigInteger &BigInteger::rem(const ulong &modulus)
{
    if(modulus == 0)
        throw std::runtime_error("Cannot divide a number by zero.");
    if(this->_pimpl->_sign == 0)
        return *this;
    int vgtuim = (int)(modulus > 0xffffffff);
    int modulus_digitLength = 1 + vgtuim;
    uint *modulus_digits = new uint[modulus_digitLength];
    modulus_digits[vgtuim] = (uint)(modulus >> 32);
    modulus_digits[0] = (uint)modulus;
    if(this->_pimpl->_digit_length >= modulus_digitLength)
        InternalBasicOperators::rem(this->_pimpl->_digits, this->_pimpl->_digit_length, modulus_digits, modulus_digitLength, this->_pimpl->_digit_length);
    return *this;
}
BigInteger &BigInteger::rem(const BigInteger &modulus)
{
    if(modulus.isZero())
        throw std::runtime_error("Cannot divide a number by zero.");
    if(this->_pimpl->_sign == 0)
        return *this;
    this->_pimpl->_sign = this->_pimpl->_sign * modulus._pimpl->_sign;
    if(this->_pimpl->_digit_length >= modulus._pimpl->_digit_length)
        InternalBasicOperators::rem(this->_pimpl->_digits, this->_pimpl->_digit_length, modulus._pimpl->_digits, modulus._pimpl->_digit_length, this->_pimpl->_digit_length);
    return *this;
}
BigInteger &BigInteger::divrem(const long &divisor, long &remainder)
{
    BigInteger r, div = divisor;
    this->divrem(div, r);
    remainder = (long)r;
    return *this;
}
BigInteger &BigInteger::divrem(const ulong &divisor, ulong &remainder)
{
    BigInteger r, div = divisor;
    this->divrem(div, r);
    remainder = (ulong)r;
    remainder += divisor;
    remainder %= divisor;
    return *this;
}
BigInteger &BigInteger::divrem(const BigInteger &divisor, BigInteger &remainder)
{
    if(divisor._pimpl->_sign == 0)
        throw std::runtime_error("Cannot divide a number by zero.");
    if(this->_pimpl->_sign == 0){
        if(this->_pimpl->_digits != remainder._pimpl->_digits){
            delete[] remainder._pimpl->_digits;
            remainder._pimpl->_digits = new uint[1] { 0 };
            remainder._pimpl->_digit_length = 1;
            remainder._pimpl->_sign = 0;
        }
        return *this;
    } else {
        this->_pimpl->_sign = this->_pimpl->_sign * divisor._pimpl->_sign;
        remainder._pimpl->_sign = this->_pimpl->_sign;
        int c = InternalBasicOperators::compare(this->_pimpl->_digits, this->_pimpl->_digit_length, divisor._pimpl->_digits, divisor._pimpl->_digit_length);
        if(c == -1){
            if(this->_pimpl->_digits != remainder._pimpl->_digits){
                delete[] remainder._pimpl->_digits;
                remainder._pimpl->_digits = this->_pimpl->_digits;
                remainder._pimpl->_digit_length = this->_pimpl->_digit_length;
                this->_pimpl->_digits = new uint[1] { 0 };
                this->_pimpl->_digit_length = 1;
                this->_pimpl->_sign = 0;
            }
            return *this;
        } else if (c == 0){
            if(this->_pimpl->_digits == remainder._pimpl->_digits){
                delete[] this->_pimpl->_digits;
                this->_pimpl->_digits = new uint[1] { 0 };
                this->_pimpl->_digit_length = 1;
            } else {
                delete[] this->_pimpl->_digits;
                this->_pimpl->_digits = new uint[1] { 1 };
                this->_pimpl->_digit_length = 1;
                delete[] remainder._pimpl->_digits;
                remainder._pimpl->_digits = new uint[1] { 0 };
                remainder._pimpl->_digit_length = 1;
                remainder._pimpl->_sign = 0;
            }
            return *this;
        }else {
            int quotientLength, remainderLength;
            uint *quot = InternalBasicOperators::divrem(this->_pimpl->_digits, this->_pimpl->_digit_length, divisor._pimpl->_digits, divisor._pimpl->_digit_length, quotientLength, remainderLength);
            if(this->_pimpl->_digits == remainder._pimpl->_digits){
                delete[] quot;
                this->_pimpl->_digit_length = remainderLength;
                return *this;
            } else {
                delete[] remainder._pimpl->_digits;
                remainder._pimpl->_digits = this->_pimpl->_digits;
                remainder._pimpl->_digit_length = remainderLength;
                this->_pimpl->_digits = quot;
                this->_pimpl->_digit_length = quotientLength;
                remainder._pimpl->_sign = (!(remainder._pimpl->_digit_length == 1 && !remainder._pimpl->_digits[0])) * remainder._pimpl->_sign;
            }
            return *this;
        }
    }
}

BigInteger BigInteger::add(const BigInteger &left, const BigInteger &right)
{
    if(left._pimpl->_sign == right._pimpl->_sign){
        if(left._pimpl->_sign == 0){
            uint *digits = new uint[1] { 0 };
            return BigInteger(digits, 1, 0, true);
        }
        int resultLength;
        uint *digits = InternalBasicOperators::add(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length, resultLength);
        return BigInteger(digits, resultLength, left._pimpl->_sign, resultLength);
    } else {
        int c = InternalBasicOperators::compare(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length);
        if(c == 1){
            int resultLength;
            uint *digits = InternalBasicOperators::sub(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length, resultLength);
            return BigInteger(digits, resultLength, left._pimpl->_sign, true);
        } else if(c == -1){
            int resultLength;
            uint *digits = InternalBasicOperators::sub(right._pimpl->_digits, right._pimpl->_digit_length, left._pimpl->_digits, left._pimpl->_digit_length, resultLength);
            return BigInteger(digits, resultLength, right._pimpl->_sign, true);
        } else {
            uint *digits = new uint[1] { 0 };
            return BigInteger(digits, 1, 0, true);
        }
    }
}
BigInteger BigInteger::sub(const BigInteger &left, const BigInteger &right)
{
    if(left._pimpl->_sign == -right._pimpl->_sign){
        if(left._pimpl->_sign == 0){
            return BigInteger(0u);
        }
        int resultLength;
        uint *digits = InternalBasicOperators::add(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length, resultLength);
        return BigInteger(digits, resultLength, left._pimpl->_sign, true);
    } else {
        int c = InternalBasicOperators::compare(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length);
        if(c == 1){
            int resultLength;
            uint *digits = InternalBasicOperators::sub(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length, resultLength);
            return BigInteger(digits, resultLength, left._pimpl->_sign, true);
        } else if(c == -1){
            int resultLength;
            uint *digits = InternalBasicOperators::sub(right._pimpl->_digits, right._pimpl->_digit_length, left._pimpl->_digits, left._pimpl->_digit_length, resultLength);
            return BigInteger(digits, resultLength, -right._pimpl->_sign, true);
        } else {
            return BigInteger(0u);
        }
    }
}
BigInteger BigInteger::mul(const BigInteger &left, const BigInteger &right)
{
    int result_sign = left._pimpl->_sign * right._pimpl->_sign;
    if(result_sign == 0)
        return BigInteger(0u);

    int resultLength;
    uint *digits = InternalBasicOperators::mul(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length, resultLength);
    return BigInteger(digits, resultLength, result_sign, true);
}
BigInteger BigInteger::div(const BigInteger &dividend, const BigInteger &divisor)
{
    if(divisor._pimpl->_sign == 0)
        throw std::runtime_error("Cannot divide a number by zero.");
    if(dividend._pimpl->_sign == 0){
        return BigInteger(0u);
    } else {
        int c = InternalBasicOperators::compare(dividend._pimpl->_digits, dividend._pimpl->_digit_length, divisor._pimpl->_digits, divisor._pimpl->_digit_length);
        int result_sign = dividend._pimpl->_sign * divisor._pimpl->_sign;
        if(c == -1)
            return BigInteger(0u);
        else if (c == 0)
            return BigInteger(result_sign);
        else {
            int quotientLength, remainderLength;
            uint *remainder = new uint[dividend._pimpl->_digit_length];
            for(int i = 0; i < dividend._pimpl->_digit_length; ++i)
                remainder[i] = dividend._pimpl->_digits[i];
            uint *quot = InternalBasicOperators::divrem(remainder, dividend._pimpl->_digit_length, divisor._pimpl->_digits, divisor._pimpl->_digit_length, quotientLength, remainderLength);
            delete[] remainder;
            return BigInteger(quot, quotientLength, result_sign, true);
        }
    }
}
BigInteger BigInteger::rem(const BigInteger &dividend, const BigInteger &divisor)
{
    if(divisor.isZero())
        throw std::runtime_error("Cannot divide a number by zero.");
    if(dividend._pimpl->_sign == 0)
        return BigInteger(0u);
    if(dividend._pimpl->_digit_length >= divisor._pimpl->_digit_length){
        uint *remainder = new uint[dividend._pimpl->_digit_length];
        for(int i = 0; i < dividend._pimpl->_digit_length; ++i)
            remainder[i] = dividend._pimpl->_digits[i];
        int remainderLength;
        InternalBasicOperators::rem(remainder, dividend._pimpl->_digit_length, divisor._pimpl->_digits, divisor._pimpl->_digit_length, remainderLength);
        return BigInteger(remainder, remainderLength, dividend._pimpl->_sign * divisor._pimpl->_sign, true);
    } else {
        BigInteger rem = dividend.clone();
        rem._pimpl->_sign *= divisor._pimpl->_sign;
        return rem;
    }
}
BigInteger BigInteger::divrem(const BigInteger &dividend, const BigInteger &divisor, BigInteger &remainder)
{
    if(dividend._pimpl->_digits == remainder._pimpl->_digits)
        return dividend.clone().divrem(divisor, remainder);
    if(divisor._pimpl->_sign == 0)
        throw std::runtime_error("Cannot divide a number by zero.");
    if(dividend._pimpl->_sign == 0){
        return BigInteger(0u);
    } else {
        int result_sign = dividend._pimpl->_sign * divisor._pimpl->_sign;
        remainder._pimpl->_sign = result_sign;
        delete[] remainder._pimpl->_digits;
        remainder._pimpl->_digits = new uint[dividend._pimpl->_digit_length];
        for(int i = 0; i < dividend._pimpl->_digit_length; ++i)
            remainder._pimpl->_digits[i] = dividend._pimpl->_digits[i];
        if(dividend._pimpl->_digit_length >= divisor._pimpl->_digit_length) {
            int quotientLength;
            uint *quot = InternalBasicOperators::divrem(remainder._pimpl->_digits, dividend._pimpl->_digit_length, divisor._pimpl->_digits, divisor._pimpl->_digit_length, quotientLength, remainder._pimpl->_digit_length);
            return BigInteger(quot, quotientLength, result_sign, true);
        } else {
            remainder._pimpl->_digit_length = dividend._pimpl->_digit_length;
            return BigInteger(0u);
        }
    }
}
BigInteger BigInteger::bitOr(const BigInteger &left, const BigInteger &right)
{
    int resultLength;
    uint *digits = InternalBasicOperators::Or(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length, resultLength);
    return BigInteger(digits, resultLength, (int)(resultLength != 1 || digits[0] != 0), true);
}
BigInteger BigInteger::bitAnd(const BigInteger &left, const BigInteger &right)
{
    int resultLength;
    uint *digits = InternalBasicOperators::And(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length, resultLength);
    return BigInteger(digits, resultLength, (int)(resultLength != 1 || digits[0] != 0), true);
}
BigInteger BigInteger::bitXor(const BigInteger &left, const BigInteger &right)
{
    int resultLength;
    uint *digits = InternalBasicOperators::Xor(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length, resultLength);
    return BigInteger(digits, resultLength, (int)(resultLength != 1 || digits[0] != 0), true);
}

BigInteger &BigInteger::bitshift(const int &shift)
{
    uint *digits;
    if(shift > 0)
        digits = InternalBasicOperators::shiftLeft(this->_pimpl->_digits, this->_pimpl->_digit_length, shift, this->_pimpl->_digit_length);
    else
        digits = InternalBasicOperators::shiftRight(this->_pimpl->_digits, this->_pimpl->_digit_length, -shift, this->_pimpl->_digit_length);
    delete[] this->_pimpl->_digits;
    this->_pimpl->_digits = digits;
    return *this;
}

BigInteger &BigInteger::rightShift(const int &shift)
{
    return this->bitshift(-shift);
}

BigInteger &BigInteger::leftShift(const int &shift)
{
    return this->bitshift(shift);
}
BigInteger &BigInteger::bitAnd(const BigInteger &value)
{
    uint *digits = InternalBasicOperators::And(this->_pimpl->_digits, this->_pimpl->_digit_length, value._pimpl->_digits, value._pimpl->_digit_length, this->_pimpl->_digit_length);
    delete[] this->_pimpl->_digits;
    this->_pimpl->_digits = digits;
    return *this;
}
BigInteger &BigInteger::bitOr(const BigInteger &value)
{
    uint *digits = InternalBasicOperators::Or(this->_pimpl->_digits, this->_pimpl->_digit_length, value._pimpl->_digits, value._pimpl->_digit_length, this->_pimpl->_digit_length);
    delete[] this->_pimpl->_digits;
    this->_pimpl->_digits = digits;
    return *this;
}
BigInteger &BigInteger::bitXor(const BigInteger &value)
{
    uint *digits = InternalBasicOperators::Xor(this->_pimpl->_digits, this->_pimpl->_digit_length, value._pimpl->_digits, value._pimpl->_digit_length, this->_pimpl->_digit_length);
    delete[] this->_pimpl->_digits;
    this->_pimpl->_digits = digits;
    return *this;
}
BigInteger &BigInteger::complement()
{
    uint *digits = InternalBasicOperators::Not(this->_pimpl->_digits, this->_pimpl->_digit_length, this->_pimpl->_digit_length);
    delete[] this->_pimpl->_digits;
    this->_pimpl->_digits = digits;
    return *this;
}
BigInteger &BigInteger::decrease()
{
    return this->sub(1l);
}
BigInteger &BigInteger::increase()
{
    return this->add(1l);
}
BigInteger &BigInteger::negate()
{
    this->_pimpl->_sign *= -1;
    return *this;
}
BigInteger &BigInteger::square()
{
    return this->mul(*this);
}

BigInteger &BigInteger::cube()
{
    BigInteger sq = this->clone().mul(*this);
    return this->mul(sq);
}
BigInteger &BigInteger::abs()
{
    this->_pimpl->_sign *= this->_pimpl->_sign;
    return *this;
}
BigInteger &BigInteger::pow(uint exponent)
{
    if (exponent == 1)
        return *this;
    if (exponent == 0){
        delete[] this->_pimpl->_digits;
        this->_pimpl->_digits = new uint[1] { 1 };
        this->_pimpl->_digit_length = 1;
        this->_pimpl->_sign = 1;
    }
    BigInteger value = this->clone();
    if (!(exponent & 1)){
        delete[] this->_pimpl->_digits;
        this->_pimpl->_digits = new uint[1] { 1 };
        this->_pimpl->_digit_length = 1;
        this->_pimpl->_sign = 1;
    }

    exponent >>= 1;
    while (exponent != 0)
    {
        value.square();
        if ((exponent & 1) == 1)
            this->mul(value);
        exponent >>= 1;
    }
    return *this;
}

void BigInteger::assign(const uint &value)
{
    this->_pimpl->_digit_length = 1;
    this->_pimpl->_digits[0] = value;
    this->_pimpl->_sign = value > 0;
}

int BigInteger::bit(const long &pos) const
{
    if (this->_pimpl->_sign == 0)
        return 0;
    if (pos > this->bitsLength())
        return 0;
    int bit = InternalBasicOperators::getBit(this->_pimpl->_digits, this->_pimpl->_digit_length, pos);
    if (this->_pimpl->_sign == 1)
        return bit;
    else
        return ~bit;
}

std::string BigInteger::toString() const
{
    if(this->_pimpl->_digit_length == 1){
        char buffer[50];
        int len = sprintf(&buffer[1], "%u", this->_pimpl->_digits[0]);
        buffer[len + 1] = 0;
        buffer[0] = '-';
        return std::string((const char *)&buffer[(this->_pimpl->_sign + 2) >> 1]);
    } else {
        return InternalBasicOperators::toString(this->_pimpl->_digits, this->_pimpl->_digit_length, this->_pimpl->_sign);
    }
}

QString BigInteger::toQString() const
{
    if(this->_pimpl->_digit_length == 1){
        char buffer[50];
        int len = sprintf(&buffer[1], "%u", this->_pimpl->_digits[0]);
        buffer[len + 1] = 0;
        buffer[0] = '-';
        return QString((const char *)&buffer[(this->_pimpl->_sign + 2) >> 1]);
    } else {
        return InternalBasicOperators::toQString(this->_pimpl->_digits, this->_pimpl->_digit_length, this->_pimpl->_sign);
    }
}

BigInteger BigInteger::clone() const
{
    uint *digits = new uint[this->_pimpl->_digit_length];
    for(int i = 0; i < this->_pimpl->_digit_length; ++i)
        digits[i] = this->_pimpl->_digits[i];
    return BigInteger(digits, this->_pimpl->_digit_length, this->_pimpl->_sign, true);
}

BigInteger BigInteger::parse(const char *chars)
{
    char *value = (char *)chars;
    int sign = 1;
    if(chars[0] == '-'){
        sign = -1;
        value = (char *)&chars[1];
    }
    int count = 0;
    char c = value[0];
    while(c != 0){
        if(c < 48 || 57 < c)
            throw std::runtime_error("Used invalid numeric literal for parsing string.");
        ++count;
        c = value[count];
    }

    int zeroTrim = 0;
    c = value[0];
    while(c == 48){
        ++zeroTrim;
        c = value[zeroTrim];
    }
    char *start = &value[zeroTrim];
    count -= zeroTrim;
    int digitLength;
    uint *digits = InternalBasicOperators::unsignedParse(start, count, digitLength);
    if(digitLength == 1 && digits[0] == 0)
        sign = 0;
    return BigInteger(digits, digitLength, sign, true);
}
BigInteger BigInteger::parse(const std::string &value)
{
    return BigInteger::parse(value.c_str());
}

BigInteger BigInteger::parseFromHex(std::string hex)
{
    bool isNeg = hex[0] == ('-');
    if (isNeg)
        hex = hex.substr(1, hex.size() - 1);

    std::transform(hex.begin(), hex.end(), hex.begin(), ::tolower);

    if (hex.rfind("0x") == 0)
        hex = hex.substr(2, hex.size() - 2);

    int byteLength = (hex.size() + 1) / 2;
    unsigned char *bytes = new unsigned char[byteLength];
    for(int i = 0; i < byteLength; ++i)
        bytes[i] = 0;

    std::string zero = "0";
    if ((hex.size() & 1) == 1)
        hex = zero.append(hex);
    for (int i = 0; i < byteLength; ++i)
    {
        char c1 = hex[2 * i];
        char c2 = hex[2 * i + 1];
        if (c1 > 96)
            bytes[i] |= (unsigned char)((c1 - 87) << 4);
        else
            bytes[i] |= (unsigned char)((c1 - 48) << 4);

        if (c2 > 96)
            bytes[i] |= (unsigned char)((c2 - 87));
        else
            bytes[i] |= (unsigned char)((c2 - 48));
    }
    int digitLength;
    uint *digits = InternalBasicOperators::fromUnsignedBytes(bytes, byteLength, true, digitLength);
    delete[] bytes;
    if (digitLength == 1 && digits[0] == 0)
        return 0;
    return BigInteger(digits, digitLength, isNeg ? -1 : 1, true);
}

BigInteger BigInteger::random(const long &bytes_length)
{
    QByteArray bytes = DRBG::fromRandomSeed().generate(bytes_length);
    return BigInteger(bytes.data(), bytes_length, 1);
}

void BigInteger::swap(BigInteger &l, BigInteger &r)
{
    BigInteger t = r;
    r = l.clone();
    l = t.clone();
}

std::vector<int> BigInteger::nonAdjacentForm(const int &window) const
{
    BigInteger d = this->clone();
    int modulus = 1 << window;

    std::vector<int> naf;
    naf.resize(d.bitsLength() + 1, 0);

    int modMinOne = modulus - 1;
    int halfOfModulus = modulus >> 1;
    for (int i = 0; !d.isZero(); ++i)
    {
        if (d.isOdd())
        {
            int mod = (int)d.firstDigit() & modMinOne; //d mod 2 ^ w

            if (mod >= halfOfModulus)
            {
                int inc = modulus - mod;
                naf[i] = -inc;
                d += inc;
            }
            else
            {
                naf[i] = mod;
                d -= (uint)mod;
            }
        }
        d = d >> 1;
    }
    return naf;
}

int BigInteger::compare(const BigInteger &left, long right)
{
    int right_sign = (right >> 63);
    right_sign = right_sign + (int)(right > 0);
    right *= right_sign;
    int vgtuim = (int)(right > 0xffffffff);
    int right_digitLength = 1 + vgtuim;
    uint *right_digits = new uint[right_digitLength];
    right_digits[vgtuim] = (uint)(right >> 32);
    right_digits[0] = (uint)right;

    if(left._pimpl->_sign < right_sign)
        return -1;
    if(left._pimpl->_sign > right_sign)
        return 1;
    return InternalBasicOperators::compare(left._pimpl->_digits, left._pimpl->_digit_length, right_digits, right_digitLength);

}
int BigInteger::compare(const BigInteger &left, ulong right)
{
    int right_sign = right > 0;
    int vgtuim = (int)(right > 0xffffffff);
    int right_digitLength = 1 + vgtuim;
    uint *right_digits = new uint[right_digitLength];
    right_digits[vgtuim] = (uint)(right >> 32);
    right_digits[0] = (uint)right;

    if(left._pimpl->_sign < right_sign)
        return -1;
    if(left._pimpl->_sign > right_sign)
        return 1;
    return InternalBasicOperators::compare(left._pimpl->_digits, left._pimpl->_digit_length, right_digits, right_digitLength);
}
int BigInteger::compare(const BigInteger &left, const BigInteger &right)
{
    if(left._pimpl->_sign > right._pimpl->_sign)
        return 1;
    if(left._pimpl->_sign < right._pimpl->_sign)
        return -1;
    return InternalBasicOperators::compare(left._pimpl->_digits, left._pimpl->_digit_length, right._pimpl->_digits, right._pimpl->_digit_length);
}

BigInteger BigInteger::pow(const BigInteger &value, uint exponent)
{
    if (exponent == 1)
        return value;
    if (exponent == 0)
        return BigInteger(1u);

    BigInteger res = BigInteger(1u);
    if (exponent & 1)
        res = value;

    BigInteger vsq = value.clone();
    exponent >>= 1;
    while (exponent != 0)
    {
        vsq.square();
        if ((exponent & 1) == 1)
            res.mul(vsq);
        exponent >>= 1;
    }
    return res;
}
BigInteger BigInteger::factorial(uint value)
{
    if (value == 0 || value == 1)
        return BigInteger(1u);

    uint cbit = (value == 0);
    uint f = (value & 0xFFFF0000) == 0;
    cbit += f << 4;
    f = ((value << cbit) & 0xFF000000) == 0;
    cbit += f << 3;
    f = ((value << cbit) & 0xF0000000) == 0;
    cbit += f << 2;
    f = ((value << cbit) & 0xC0000000) == 0;
    cbit += f << 1;
    f = ((value << cbit) & 0x80000000) == 0;
    cbit += f;

    uint t = (32 - cbit);
    uint *a = new uint[t];
    uint shift = 0;

    for (uint i = 0; i < t; ++i){
        if ((value & 1) == 1){
            value = a[i] = value >> 1;
            shift += value;
        }
        else{
            a[i] = (value = value >> 1) - 1;
            shift += value;
        }
    }
    BigInteger inner = BigInteger(1u), T = BigInteger(1u);
    uint abovePi, belowPi;
    for (uint i = 1; i < t; ++i){
        abovePi = a[i - 1];
        belowPi = a[i] + 1;
        for (uint j = belowPi; j <= abovePi; ++j)
            inner.mul((j << 1) | 1);
        inner.pow(i);
        T.mul(inner);
        inner.assign(1u);
    }
    delete[] a;
    T.bitshift((int)shift);
    return T;
}
BigInteger BigInteger::factorize(BigInteger &value)
{
    if (value.isZero())
        return BigInteger(0u);
    if (value.isOne())
        return BigInteger(1u);

    BigInteger a = 5u;
    a.divrem(value, a);
    BigInteger b = 26u;
    b.divrem(value, b);
    while (true)
    {
        BigInteger c = a < b ? BigInteger::sub(b, a) : BigInteger::sub(a, b);
        BigInteger d = BigInteger::gcd(c, value);
        if (d.isOne())
        {
            a.square().increase().rem(value);
            b.square().increase().rem(value);
            b.square().increase().rem(value);
            continue;
        }
        else{
            value.div(d);
            return d;
        }
    }
}
BigInteger BigInteger::gcd(const BigInteger &left, const BigInteger &right)
{
    if (right.isZero())
        return left.clone();
    else{
        BigInteger rem = left % right; //BigInteger::rem(left, right);
        return BigInteger::gcd(right, rem);
    }
}
BigInteger BigInteger::lcm(const BigInteger &left, const BigInteger &right)
{
    BigInteger gcd = BigInteger::gcd(left, right);
    BigInteger leftFac = left / gcd;
    BigInteger rightFac = right / gcd;
    gcd.mul(leftFac).mul(rightFac);
    return gcd;
}
BigInteger BigInteger::modpow(const BigInteger &value, const BigInteger &exponent, const BigInteger &modulus)
{
    if(modulus.isNegative())
        throw std::runtime_error("modulus cannot be negative");
    if(value.isZero())
        return BigInteger(0u);
    if(value.isNegative() && exponent.isNegative())
        throw std::runtime_error("value and exponent cannot be negative currently");
    BigInteger temp_value = value.clone();
    if(exponent.isNegative())
        temp_value = BigInteger::modInverse(temp_value, modulus);

    BigInteger result = BigInteger(1u);
    long bitLength = InternalBasicOperators::unsignedBitsLength(exponent._pimpl->_digits, exponent._pimpl->_digit_length);
    for (long i = 0; i < bitLength; ++i)
    {
        int bit = 0;
        int digitPos = (int)(i >> 5);
        if(digitPos < exponent._pimpl->_digit_length)
            bit = (exponent._pimpl->_digits[digitPos] >> (i & 31)) & 1;
        if (temp_value._pimpl->_digit_length > modulus._pimpl->_digit_length)
            temp_value.divrem(modulus, temp_value);
        if (bit == 1)
        {
            result.mul(temp_value);
            if (result._pimpl->_digit_length > modulus._pimpl->_digit_length)
                result.divrem(modulus, result);
        }
        temp_value.square();
    }
    result.divrem(modulus, result);
    return result;
}
BigInteger BigInteger::squareroot(const BigInteger &value)
{
    if(value.isNegative())
        throw std::runtime_error("value cannot be a negative");
    BigInteger xf = BigInteger(1u);
    BigInteger xl = BigInteger(1u);
    do
    {
        xf = xl;
        xl = value;
        xl.div(xf).add(xf).bitshift(-1);
    } while (!(xf == xl || xf.sub(xl).isMinusOne()));
    if(xf.isMinusOne())
        xl.decrease();
    return xl;
}
BigInteger BigInteger::modInverse(const BigInteger &value, const BigInteger &modulus)
{
    if (value.isOne())
        return BigInteger(1u);

    BigInteger x1 = BigInteger::zero(), x2 = modulus.clone(), y1 = BigInteger::one(), y2 = value.clone();

    BigInteger t1, t2, q = BigInteger::divrem(x2, y2, t2);
    q.negate();
    t1 = q.clone();

    while (!y2.isOne())
    {
        if (t2.isZero())
            return BigInteger::zero();

        x1 = y1.clone(); x2 = y2.clone();
        y1 = t1.clone(); y2 = t2.clone();
        q = BigInteger::divrem(x2, y2, t2);

        t1 = x1 - q * y1;
    }
    if (y1.sign() == -1)
        return y1 + modulus;
    else
        return y1;
}

BigInteger BigInteger::square(const BigInteger &value)
{
    return value.clone().square();
}

BigInteger::operator int() const
{
    uint f_digit = this->_pimpl->_digits[0];
    f_digit -= (f_digit > 0x7FFFFFFF) * 0x7FFFFFFF;
    return (int)f_digit * this->_pimpl->_sign;
}
BigInteger::operator uint() const
{
    uint f_digit = this->_pimpl->_digits[0];
    return (uint)(f_digit * this->_pimpl->_sign);
}
BigInteger::operator long() const
{
    if(this->_pimpl->_digit_length > 1){
        ulong dig = ((ulong *)this->_pimpl->_digits)[0];
        dig -= (dig > 0x7FFFFFFFFFFFFFFF) * 0x7FFFFFFFFFFFFFFF;
        return (long)dig * this->_pimpl->_sign;
    } else {
        long f_digit = this->_pimpl->_digits[0];
        return f_digit * this->_pimpl->_sign;
    }
}
BigInteger::operator ulong() const
{
    if(this->_pimpl->_digit_length > 1){
        ulong dig = ((ulong *)this->_pimpl->_digits)[0];
        return (ulong)dig * this->_pimpl->_sign;
    } else {
        ulong f_digit = this->_pimpl->_digits[0];
        return f_digit * this->_pimpl->_sign;
    }
}

BigInteger::BigInteger(uint *digits, int length, int sign, bool internal)
{
    this->_pimpl = std::shared_ptr<pimpl>(new pimpl);

    if(internal){
        if(length == 0){
            delete[] digits;
            this->_pimpl->_digit_length = 1;
            this->_pimpl->_digits = new uint[1] { 0 };
            this->_pimpl->_sign = 0;
        } else {
            if(length == 1 && digits[0] == 0)
                this->_pimpl->_sign = 0;
            else
                this->_pimpl->_sign = sign;
            if(sign == 0){
                this->_pimpl->_digit_length = 1;
                this->_pimpl->_digits = new uint[1] { 0 };
                this->_pimpl->_sign = 0;
            } else{
                this->_pimpl->_digits = digits;
                this->_pimpl->_digit_length = length;
            }
        }
    } else {
        if(length == 0){
            this->_pimpl->_sign = 0;
            this->_pimpl->_digit_length = length;
            this->_pimpl->_digits = new uint[1] { 0 };
        } else {
            while(digits[length - 1] == 0 && length > 1)
                length--;
            if(length == 1 && digits[0] == 0)
                this->_pimpl->_sign = 0;
            else
                this->_pimpl->_sign = sign;
            if(sign == 0){
                this->_pimpl->_digit_length = 1;
                this->_pimpl->_digits = new uint[1] { 0 };
                this->_pimpl->_sign = 0;
            } else {
                this->_pimpl->_digit_length = length;
                this->_pimpl->_digits = new uint[length];
                for(int i = 0; i < length; ++i)
                    this->_pimpl->_digits[i] = digits[i];
            }
        }
    }
}
