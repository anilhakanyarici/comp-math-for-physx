#ifndef BIGINTEGER_H
#define BIGINTEGER_H

#include <QString>

#include <vector>
#include <memory>

namespace mp
{

class BigInteger final
{
    struct pimpl;
    std::shared_ptr<pimpl> _pimpl;

public:
    BigInteger();
    BigInteger(int value);
    BigInteger(uint value);
    BigInteger(long value);
    BigInteger(ulong value);
    BigInteger(char *bytes, long bytes_length, int sign, bool big_endian = false);
    BigInteger(const QByteArray &bytes, int sign,  bool big_endian = false);

//    BigInteger(const BigInteger &value);
//    BigInteger &operator =(const BigInteger &value);

    static inline BigInteger zero() { return BigInteger(); }
    static inline BigInteger one() { return BigInteger(1u); }
    static BigInteger twoPowerOf(const long &exponent);

    bool isZero() const;
    bool isOne() const;
    bool isMinusOne() const;
    bool isNegative() const;
    bool isPositive() const;
    bool isPowerOfTwo() const;
    bool isEven() const;
    bool isOdd() const;
    int sign() const;
    int digitLength() const;
    long bytesLength() const;
    long bitsLength() const;
    uint firstDigit() const;
    uint *data() const;
    char *bytesData() const;

    inline BigInteger &add(const int &value) { return this->add((long)value); }
    inline BigInteger &sub(const int &value) { return this->sub((long)value); }
    inline BigInteger &mul(const int &value) { return this->mul((long)value); }
    inline BigInteger &div(const int &divisor) { return this->div((long)divisor); }
    inline BigInteger &rem(const int &modulus) { return this->rem((long)modulus); }
    BigInteger &divrem(const int &divisor, int &remainder) { long rem; this->divrem((long)divisor, rem); remainder = (int)rem; return *this; }
    inline BigInteger &add(const uint &value) { return this->add((ulong)value); }
    inline BigInteger &sub(const uint &value) { return this->sub((ulong)value); }
    inline BigInteger &mul(const uint &value) { return this->mul((ulong)value); }
    inline BigInteger &div(const uint &divisor) { return this->div((ulong)divisor); }
    inline BigInteger &rem(const uint &modulus) { return this->rem((ulong)modulus); }
    BigInteger &divrem(const uint &divisor, uint &remainder) { ulong rem; this->divrem((ulong)divisor, rem); remainder = (uint)rem; return *this; }
    inline BigInteger &add(const long &value) { return this->add(BigInteger(value)); }
    inline BigInteger &sub(const long &value) { return this->sub(BigInteger(value)); }
    inline BigInteger &mul(const long &value) { return this->mul(BigInteger(value)); }
    BigInteger &div(const long &divisor);
    BigInteger &rem(long modulus);
    BigInteger &divrem(const long &divisor, long &remainder);
    inline BigInteger &add(const ulong &value) { return this->add(BigInteger(value)); }
    inline BigInteger &sub(const ulong &value) { return this->sub(BigInteger(value)); }
    inline BigInteger &mul(const ulong &value) { return this->mul(BigInteger(value)); }
    BigInteger &div(const ulong &divisor);
    BigInteger &rem(const ulong &modulus);
    BigInteger &divrem(const ulong &divisor, ulong &remainder);
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

    BigInteger &bitshift(const int &shift);
    BigInteger &rightShift(const int &shift);
    BigInteger &leftShift(const int &shift);
    BigInteger &bitAnd(const BigInteger &value);
    BigInteger &bitOr(const BigInteger &value);
    BigInteger &bitXor(const BigInteger &value);
    BigInteger &complement();
    BigInteger &decrease();
    BigInteger &increase();
    BigInteger &negate();
    BigInteger &square();
    BigInteger &cube();
    BigInteger &abs();
    BigInteger &pow(uint exponent);

    void assign(const uint &value);
    int bit(const long &pos) const;
    BigInteger clone() const;
    std::vector<int> nonAdjacentForm(const int &window) const;
    static BigInteger parse(const char *chars);
    static BigInteger parse(const std::string &value);
    static BigInteger parseFromHex(std::string hex);
    static BigInteger random(const long &bytes_length);
    static void swap(BigInteger &l, BigInteger &r);

    static inline int compare(int left, const BigInteger &right) { return -BigInteger::compare(right, (long)left); }
    static inline int compare(uint left, const BigInteger &right) { return -BigInteger::compare(right, (ulong)left); }
    static inline int compare(long left, const BigInteger &right) { return -BigInteger::compare(right, left); }
    static inline int compare(ulong left, const BigInteger &right) { return -BigInteger::compare(right, left); }
    static inline int compare(const BigInteger &left, int right) { return BigInteger::compare(left, (long)right); }
    static inline int compare(const BigInteger &left, uint right) { return BigInteger::compare(left, (ulong)right); }
    static int compare(const BigInteger &left, long right);
    static int compare(const BigInteger &left, ulong right);
    static int compare(const BigInteger &left, const BigInteger &right);

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
    inline friend BigInteger operator >>(const BigInteger &value, int shift) { return value.clone().bitshift(-shift); }
    inline friend BigInteger operator <<(const BigInteger &value, int shift) { return value.clone().bitshift(shift); }
    inline friend BigInteger operator ~(const BigInteger &value) { return value.clone().complement(); }

    inline BigInteger &operator ^=(const BigInteger &value) { return this->bitXor(value); }
    inline BigInteger &operator &=(const BigInteger &value) { return this->bitAnd(value); }
    inline BigInteger &operator |=(const BigInteger &value) { return this->bitOr(value); }
    inline BigInteger &operator >>=(int shift) { return this->bitshift(-shift); }
    inline BigInteger &operator <<=(int shift) { return this->bitshift(shift); }

    inline friend BigInteger operator -(BigInteger value) { return value.clone().negate(); }
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

    inline friend bool operator <(const BigInteger &left, int value) { return BigInteger::compare(left, value) == -1; }
    inline friend bool operator >(const BigInteger &left, int value) { return BigInteger::compare(left, value) == 1; }
    inline friend bool operator !=(const BigInteger &left, int value) { return BigInteger::compare(left, value) != 0; }
    inline friend bool operator ==(const BigInteger &left, int value) { return BigInteger::compare(left, value) == 0; }
    inline friend bool operator >=(const BigInteger &left, int value) { return BigInteger::compare(left, value) < 1; }
    inline friend bool operator <=(const BigInteger &left, int value) { return BigInteger::compare(left, value) > -1; }
    inline friend bool operator <(const BigInteger &left, uint value) { return BigInteger::compare(left, value) == -1; }
    inline friend bool operator >(const BigInteger &left, uint value) { return BigInteger::compare(left, value) == 1; }
    inline friend bool operator !=(const BigInteger &left, uint value) { return BigInteger::compare(left, value) != 0; }
    inline friend bool operator ==(const BigInteger &left, uint value) { return BigInteger::compare(left, value) == 0; }
    inline friend bool operator >=(const BigInteger &left, uint value) { return BigInteger::compare(left, value) < 1; }
    inline friend bool operator <=(const BigInteger &left, uint value) { return BigInteger::compare(left, value) > -1; }
    inline friend bool operator <(const BigInteger &left, long value) { return BigInteger::compare(left, value) == -1; }
    inline friend bool operator >(const BigInteger &left, long value) { return BigInteger::compare(left, value) == 1; }
    inline friend bool operator !=(const BigInteger &left, long value) { return BigInteger::compare(left, value) != 0; }
    inline friend bool operator ==(const BigInteger &left, long value) { return BigInteger::compare(left, value) == 0; }
    inline friend bool operator >=(const BigInteger &left, long value) { return BigInteger::compare(left, value) < 1; }
    inline friend bool operator <=(const BigInteger &left, long value) { return BigInteger::compare(left, value) > -1; }
    inline friend bool operator <(const BigInteger &left, ulong value) { return BigInteger::compare(left, value) == -1; }
    inline friend bool operator >(const BigInteger &left, ulong value) { return BigInteger::compare(left, value) == 1; }
    inline friend bool operator !=(const BigInteger &left, ulong value) { return BigInteger::compare(left, value) != 0; }
    inline friend bool operator ==(const BigInteger &left, ulong value) { return BigInteger::compare(left, value) == 0; }
    inline friend bool operator >=(const BigInteger &left, ulong value) { return BigInteger::compare(left, value) < 1; }
    inline friend bool operator <=(const BigInteger &left, ulong value) { return BigInteger::compare(left, value) > -1; }
    inline friend bool operator <(const BigInteger &left, const BigInteger &value) { return BigInteger::compare(left, value) == -1; }
    inline friend bool operator >(const BigInteger &left, const BigInteger &value) { return BigInteger::compare(left, value) == 1; }
    inline friend bool operator !=(const BigInteger &left, const BigInteger &value) { return BigInteger::compare(left, value) != 0; }
    inline friend bool operator ==(const BigInteger &left, const BigInteger &value) { return BigInteger::compare(left, value) == 0; }
    inline friend bool operator <=(const BigInteger &left, const BigInteger &value) { return BigInteger::compare(left, value) < 1; }
    inline friend bool operator >=(const BigInteger &left, const BigInteger &value) { return BigInteger::compare(left, value) > -1; }

    operator int() const;
    operator uint() const;
    operator long() const;
    operator ulong() const;
    std::string toString() const;
    QString toQString() const;

    static BigInteger pow(const BigInteger &value, uint exponent);
    static BigInteger factorial(uint value);
    static BigInteger factorize(BigInteger &value);
    static BigInteger gcd(const BigInteger &left, const BigInteger &right);
    static BigInteger lcm(const BigInteger &left, const BigInteger &right);
    static BigInteger modpow(const BigInteger &value, const BigInteger &exponent, const BigInteger &modulus);
    static BigInteger squareroot(const BigInteger &value);
    static BigInteger modInverse(const BigInteger &value, const BigInteger &modulus);
    static BigInteger square(const BigInteger &value);

private:
    BigInteger(uint* digits, int length, int sign, bool internal);
};
}

#endif // BIGINTEGER_H
