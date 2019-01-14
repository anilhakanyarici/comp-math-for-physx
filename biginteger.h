#ifndef BIGINTEGER_H
#define BIGINTEGER_H

#include <string>

namespace mp
{

class BigInteger final
{
public:
    BigInteger();
    BigInteger(int value);
    BigInteger(uint value);
    BigInteger(long value);
    BigInteger(ulong value);

    BigInteger(const BigInteger &big);

    ~BigInteger() { delete[] this->_digits; }

    bool isZero() const { return this->_sign == 0; }
    bool isOne() const { return this->_digitLength == 1 && this->_digits[0] == 1 && this->_sign == 1; }
    bool isMinusOne() const { return this->_digitLength == 1 && this->_digits[0] == 1 && this->_sign == -1; }
    bool isNegative() const { return this->_sign == -1; }
    bool isPositive() const { return this->_sign == 1; }
    bool isPowerOfTwo() const;
    bool isEven() const { return (this->_digits[0] & 1) == 0; }
    bool isOdd() const { return (this->_digits[0] & 1); }
    int getSign() const { return this->_sign; }
    int getDigitLength() const { return this->_digitLength; }
    uint *c_data() const{ return this->_digits; }

    inline BigInteger &add(int value) { return this->add((long)value); }
    inline BigInteger &sub(int value) { return this->sub((long)value); }
    inline BigInteger &mul(int value) { return this->mul((long)value); }
    inline BigInteger &div(int divisor) { return this->div((long)divisor); }
    inline BigInteger &rem(int modulus) { return this->rem((long)modulus); }
    BigInteger &divrem(int divisor, int &remainder) { long rem; this->divrem((long)divisor, rem); remainder = (int)rem; return *this; }
    inline BigInteger &add(uint value) { return this->add((ulong)value); }
    inline BigInteger &sub(uint value) { return this->sub((ulong)value); }
    inline BigInteger &mul(uint value) { return this->mul((ulong)value); }
    inline BigInteger &div(uint divisor) { return this->div((ulong)divisor); }
    inline BigInteger &rem(uint modulus) { return this->rem((ulong)modulus); }
    BigInteger &divrem(uint divisor, uint &remainder) { ulong rem; this->divrem((ulong)divisor, rem); remainder = (uint)rem; return *this; }
    BigInteger &add(long value);
    BigInteger &sub(long value);
    BigInteger &mul(long value);
    BigInteger &div(long divisor);
    BigInteger &rem(long modulus);
    BigInteger &divrem(long divisor, long &remainder);
    BigInteger &add(ulong value);
    BigInteger &sub(ulong value);
    BigInteger &mul(ulong value);
    BigInteger &div(ulong divisor);
    BigInteger &rem(ulong modulus);
    BigInteger &divrem(ulong divisor, ulong &remainder);
    BigInteger &add(const BigInteger &value);
    BigInteger &sub(const BigInteger &value);
    BigInteger &mul(const BigInteger &value);
    BigInteger &div(const BigInteger &divisor);
    BigInteger &rem(const BigInteger &modulus);
    BigInteger &divrem(const BigInteger &divisor, BigInteger &remainder);

    static BigInteger add(const BigInteger &left, const BigInteger &right);
    static BigInteger sub(const BigInteger &left, const BigInteger &right);
    static BigInteger mul(const BigInteger &left, const BigInteger &right);
    static BigInteger div(const BigInteger &dividend, const BigInteger &divisor);
    static BigInteger rem(const BigInteger &dividend, const BigInteger &divisor);
    static BigInteger divrem(const BigInteger &dividend, const BigInteger &divisor, BigInteger &remainder);
    static BigInteger bitOr(const BigInteger &left, const BigInteger &right);
    static BigInteger bitAnd(const BigInteger &left, const BigInteger &right);
    static BigInteger bitXor(const BigInteger &left, const BigInteger &right);

    BigInteger &bitshift(int shift);
    BigInteger &bitAnd(const BigInteger &value);
    BigInteger &bitOr(const BigInteger &value);
    BigInteger &bitXor(const BigInteger &value);
    BigInteger &complement();
    BigInteger &decrease();
    BigInteger &increase();
    BigInteger &negate();
    BigInteger &square();
    BigInteger &abs();
    BigInteger &pow(uint exponent);

    void assign(uint value);
    BigInteger clone() const;
    static BigInteger parse(const char *chars);
    static BigInteger parse(const std::string &value);

    static inline int compare(int left, const BigInteger &right) { return -BigInteger::compare(right, (long)left); }
    static inline int compare(uint left, const BigInteger &right) { return -BigInteger::compare(right, (ulong)left); }
    static inline int compare(long left, const BigInteger &right) { return -BigInteger::compare(right, left); }
    static inline int compare(ulong left, const BigInteger &right) { return -BigInteger::compare(right, left); }
    static inline int compare(const BigInteger &left, int right) { return BigInteger::compare(left, (long)right); }
    static inline int compare(const BigInteger &left, uint right) { return BigInteger::compare(left, (ulong)right); }
    static int compare(const BigInteger &left, long right);
    static int compare(const BigInteger &left, ulong right);
    static int compare(const BigInteger &left, const BigInteger &right);

    BigInteger &operator =(int value);
    BigInteger &operator =(uint value);
    BigInteger &operator =(long value);
    BigInteger &operator =(ulong value);

    BigInteger &operator =(const BigInteger &big);

    friend BigInteger operator +(const BigInteger &left, int right) { BigInteger r = right; return BigInteger::add(left, r); }
    friend BigInteger operator -(const BigInteger &left, int right) { BigInteger r = right; return BigInteger::sub(left, r); }
    friend BigInteger operator *(const BigInteger &left, int right) { BigInteger r = right; return BigInteger::mul(left, r); }
    friend BigInteger operator /(const BigInteger &left, int right) { BigInteger r = right; return BigInteger::div(left, r); }
    friend BigInteger operator %(const BigInteger &left, int right) { BigInteger r = right; return BigInteger::rem(left, r); }
    friend BigInteger operator +(const BigInteger &left, uint right) { BigInteger r = right; return BigInteger::add(left, r); }
    friend BigInteger operator -(const BigInteger &left, uint right) { BigInteger r = right; return BigInteger::sub(left, r); }
    friend BigInteger operator *(const BigInteger &left, uint right) { BigInteger r = right; return BigInteger::mul(left, r); }
    friend BigInteger operator /(const BigInteger &left, uint right) { BigInteger r = right; return BigInteger::div(left, r); }
    friend BigInteger operator %(const BigInteger &left, uint right) { BigInteger r = right; return BigInteger::rem(left, r); }
    friend BigInteger operator +(const BigInteger &left, long right) { BigInteger r = right; return BigInteger::add(left, r); }
    friend BigInteger operator -(const BigInteger &left, long right) { BigInteger r = right; return BigInteger::sub(left, r); }
    friend BigInteger operator *(const BigInteger &left, long right) { BigInteger r = right; return BigInteger::mul(left, r); }
    friend BigInteger operator /(const BigInteger &left, long right) { BigInteger r = right; return BigInteger::div(left, r); }
    friend BigInteger operator %(const BigInteger &left, long right) { BigInteger r = right; return BigInteger::rem(left, r); }
    friend BigInteger operator +(const BigInteger &left, ulong right) { BigInteger r = right; return BigInteger::add(left, r); }
    friend BigInteger operator -(const BigInteger &left, ulong right) { BigInteger r = right; return BigInteger::sub(left, r); }
    friend BigInteger operator *(const BigInteger &left, ulong right) { BigInteger r = right; return BigInteger::mul(left, r); }
    friend BigInteger operator /(const BigInteger &left, ulong right) { BigInteger r = right; return BigInteger::div(left, r); }
    friend BigInteger operator %(const BigInteger &left, ulong right) { BigInteger r = right; return BigInteger::rem(left, r); }
    inline friend BigInteger operator +(const BigInteger &left, const BigInteger &right) { return BigInteger::add(left, right); }
    inline friend BigInteger operator -(const BigInteger &left, const BigInteger &right) { return BigInteger::sub(left, right); }
    inline friend BigInteger operator *(const BigInteger &left, const BigInteger &right) { return BigInteger::mul(left, right); }
    inline friend BigInteger operator /(const BigInteger &left, const BigInteger &right) { return BigInteger::div(left, right); }
    inline friend BigInteger operator %(const BigInteger &left, const BigInteger &right) { return BigInteger::rem(left, right); }

    inline BigInteger &operator +=(int value){ return this->add(value); }
    inline BigInteger &operator +=(uint value){ return this->add(value); }
    inline BigInteger &operator +=(long value){ return this->add(value); }
    inline BigInteger &operator +=(ulong value){ return this->add(value); }
    inline BigInteger &operator +=(const BigInteger &value){ return this->add(value); }
    inline BigInteger &operator -=(int value){ return this->sub(value); }
    inline BigInteger &operator -=(uint value){ return this->sub(value); }
    inline BigInteger &operator -=(long value){ return this->sub(value); }
    inline BigInteger &operator -=(ulong value){ return this->sub(value); }
    inline BigInteger &operator -=(const BigInteger &value){ return this->sub(value); }
    inline BigInteger &operator *=(int value){ return this->mul(value); }
    inline BigInteger &operator *=(uint value){ return this->mul(value); }
    inline BigInteger &operator *=(long value){ return this->mul(value); }
    inline BigInteger &operator *=(ulong value){ return this->mul(value); }
    inline BigInteger &operator *=(const BigInteger &value){ return this->mul(value); }
    inline BigInteger &operator /=(int value){ return this->div(value); }
    inline BigInteger &operator /=(uint value){ return this->div(value); }
    inline BigInteger &operator /=(long value){ return this->div(value); }
    inline BigInteger &operator /=(ulong value){ return this->div(value); }
    inline BigInteger &operator /=(const BigInteger &value){ return this->div(value); }
    inline BigInteger &operator %=(int value){ return this->rem(value); }
    inline BigInteger &operator %=(uint value){ return this->rem(value); }
    inline BigInteger &operator %=(long value){ return this->rem(value); }
    inline BigInteger &operator %=(ulong value){ return this->rem(value); }
    inline BigInteger &operator %=(const BigInteger &value){ return this->rem(value); }

    inline friend BigInteger operator ^(const BigInteger &left, const BigInteger &right) { return BigInteger::bitXor(left, right); }
    inline friend BigInteger operator |(const BigInteger &left, const BigInteger &right) { return BigInteger::bitOr(left, right); }
    inline friend BigInteger operator &(const BigInteger &left, const BigInteger &right) { return BigInteger::bitAnd(left, right); }
    inline friend BigInteger operator >>(BigInteger value, int shift) { return value.bitshift(-shift); }
    inline friend BigInteger operator <<(BigInteger value, int shift) { return value.bitshift(shift); }
    inline friend BigInteger operator ~(BigInteger value) { return value.complement(); }

    inline BigInteger &operator ^=(const BigInteger &value) { return this->bitXor(value); }
    inline BigInteger &operator &=(const BigInteger &value) { return this->bitAnd(value); }
    inline BigInteger &operator |=(const BigInteger &value) { return this->bitOr(value); }
    inline BigInteger &operator >>=(int shift) { return this->bitshift(-shift); }
    inline BigInteger &operator <<=(int shift) { return this->bitshift(shift); }

    inline friend BigInteger operator -(BigInteger value) { return value.negate(); }
    inline friend BigInteger operator +(BigInteger value) { return value.abs(); }
    inline BigInteger operator --() { return this->decrease(); }
    inline BigInteger operator ++() { return this->increase(); }
    inline BigInteger operator --(int) { return operator --(); }
    inline BigInteger operator ++(int) { return operator ++(); }

    inline friend bool operator <(int left, const BigInteger &right) { return BigInteger::compare(left, right) == -1; }
    inline friend bool operator >(int left, const BigInteger &right) { return BigInteger::compare(left, right) == 1; }
    inline friend bool operator ==(int left, const BigInteger &right) { return BigInteger::compare(left, right) != 0; }
    inline friend bool operator !=(int left, const BigInteger &right) { return BigInteger::compare(left, right) == 0; }
    inline friend bool operator <=(int left, const BigInteger &right) { return BigInteger::compare(left, right) < 1; }
    inline friend bool operator >=(int left, const BigInteger &right) { return BigInteger::compare(left, right) > -1; }
    inline friend bool operator <(uint left, const BigInteger &right) { return BigInteger::compare(left, right) == -1; }
    inline friend bool operator >(uint left, const BigInteger &right) { return BigInteger::compare(left, right) == 1; }
    inline friend bool operator ==(uint left, const BigInteger &right) { return BigInteger::compare(left, right) != 0; }
    inline friend bool operator !=(uint left, const BigInteger &right) { return BigInteger::compare(left, right) == 0; }
    inline friend bool operator <=(uint left, const BigInteger &right) { return BigInteger::compare(left, right) < 1; }
    inline friend bool operator >=(uint left, const BigInteger &right) { return BigInteger::compare(left, right) > -1; }
    inline friend bool operator <(long left, const BigInteger &right) { return BigInteger::compare(left, right) == -1; }
    inline friend bool operator >(long left, const BigInteger &right) { return BigInteger::compare(left, right) == 1; }
    inline friend bool operator ==(long left, const BigInteger &right) { return BigInteger::compare(left, right) != 0; }
    inline friend bool operator !=(long left, const BigInteger &right) { return BigInteger::compare(left, right) == 0; }
    inline friend bool operator <=(long left, const BigInteger &right) { return BigInteger::compare(left, right) < 1; }
    inline friend bool operator >=(long left, const BigInteger &right) { return BigInteger::compare(left, right) > -1; }
    inline friend bool operator <(ulong left, const BigInteger &right) { return BigInteger::compare(left, right) == -1; }
    inline friend bool operator >(ulong left, const BigInteger &right) { return BigInteger::compare(left, right) == 1; }
    inline friend bool operator ==(ulong left, const BigInteger &right) { return BigInteger::compare(left, right) != 0; }
    inline friend bool operator !=(ulong left, const BigInteger &right) { return BigInteger::compare(left, right) == 0; }
    inline friend bool operator <=(ulong left, const BigInteger &right) { return BigInteger::compare(left, right) < 1; }
    inline friend bool operator >=(ulong left, const BigInteger &right) { return BigInteger::compare(left, right) > -1; }
    inline bool operator <(int value) const { return BigInteger::compare(*this, value) == -1; }
    inline bool operator >(int value) const { return BigInteger::compare(*this, value) == 1; }
    inline bool operator !=(int value) const { return BigInteger::compare(*this, value) != 0; }
    inline bool operator ==(int value) const { return BigInteger::compare(*this, value) == 0; }
    inline bool operator >=(int value) const { return BigInteger::compare(*this, value) < 1; }
    inline bool operator <=(int value) const { return BigInteger::compare(*this, value) > -1; }
    inline bool operator <(uint value) const { return BigInteger::compare(*this, value) == -1; }
    inline bool operator >(uint value) const { return BigInteger::compare(*this, value) == 1; }
    inline bool operator !=(uint value) const { return BigInteger::compare(*this, value) != 0; }
    inline bool operator ==(uint value) const { return BigInteger::compare(*this, value) == 0; }
    inline bool operator >=(uint value) const { return BigInteger::compare(*this, value) < 1; }
    inline bool operator <=(uint value) const { return BigInteger::compare(*this, value) > -1; }
    inline bool operator <(long value) const { return BigInteger::compare(*this, value) == -1; }
    inline bool operator >(long value) const { return BigInteger::compare(*this, value) == 1; }
    inline bool operator !=(long value) const { return BigInteger::compare(*this, value) != 0; }
    inline bool operator ==(long value) const { return BigInteger::compare(*this, value) == 0; }
    inline bool operator >=(long value) const { return BigInteger::compare(*this, value) < 1; }
    inline bool operator <=(long value) const { return BigInteger::compare(*this, value) > -1; }
    inline bool operator <(ulong value) const { return BigInteger::compare(*this, value) == -1; }
    inline bool operator >(ulong value) const { return BigInteger::compare(*this, value) == 1; }
    inline bool operator !=(ulong value) const { return BigInteger::compare(*this, value) != 0; }
    inline bool operator ==(ulong value) const { return BigInteger::compare(*this, value) == 0; }
    inline bool operator >=(ulong value) const { return BigInteger::compare(*this, value) < 1; }
    inline bool operator <=(ulong value) const { return BigInteger::compare(*this, value) > -1; }
    inline bool operator <(const BigInteger &value) const { return BigInteger::compare(*this, value) == -1; }
    inline bool operator >(const BigInteger &value) const { return BigInteger::compare(*this, value) == 1; }
    inline bool operator !=(const BigInteger &value) const { return BigInteger::compare(*this, value) != 0; }
    inline bool operator ==(const BigInteger &value) const { return BigInteger::compare(*this, value) == 0; }
    inline bool operator <=(const BigInteger &value) const { return BigInteger::compare(*this, value) < 1; }
    inline bool operator >=(const BigInteger &value) const { return BigInteger::compare(*this, value) > -1; }

    operator int() const;
    operator uint() const;
    operator long() const;
    operator ulong() const;
    std::string toString() const;

    static BigInteger pow(BigInteger value, uint exponent);
    static BigInteger factorial(uint value);
    static BigInteger factorize(BigInteger &value);
    static BigInteger gcd(const BigInteger &left, const BigInteger &right);
    static BigInteger lcm(const BigInteger &left, const BigInteger &right);
    static BigInteger modpow(BigInteger value, const BigInteger &exponent, const BigInteger &modulus);
    static BigInteger squareroot(const BigInteger &value);
    static BigInteger modInverse(const BigInteger &value, const BigInteger &modulus);

private:
    int _digitLength;
    uint *_digits;
    int _sign;

    BigInteger(uint* digits, int length, int sign, bool internal);
};
}



#endif // BIGINTEGER_H
