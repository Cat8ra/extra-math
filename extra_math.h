#pragma once

#include <cmath>
#include <string>
#include <ctime>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

/*namespace Random {
    void initialize() { srand((unsigned int)time(NULL)); }
    char nextChar() { return rand(); }
    int nextInt() {
        return rand() + rand() * RAND_MAX + rand() * RAND_MAX * RAND_MAX;
    }
    unsigned int nextUint() {
        return (unsigned int)nextInt();
    }
    unsigned long long nextLong() {
        return nextInt() * 4294967296 + nextInt();
    }
    unsigned long long next_ulong() {
        return (unsigned long long)nextLong();
    }
};*/

struct Random {
    Random(long long seed = 0) {
        now = seed;
        first = true;
    }
    bool nextBool() {
        step();
        return now >> 19 % 2 != 0;
    }
    char nextChar() {
        step();
        return now >> 19;
    }
    short nextSInt() {
        step();
        return now >> 19;
    }
    unsigned short nextUSInt() {
        step();
        return now >> 19;
    }
    int nextInt() {
        step();
        return now >> 19;
    }
    unsigned int nextUInt() {
        step();
        return now >> 19;
    }
    long long nextLong() {
        step();
        return now;
    }
    unsigned long long nextULong() {
        step();
        return now;
    }
    double nextDouble() {
        step();
        return *(double*)(&now);
    }

private:
    void step() {
        now = now * a[((now % 3) + 3) % 3] + b[((now >> 7 % 3) + 4) % 3];
    }
    static const long long a[3];
    static const long long b[3];
    long long now;
    bool first;
};


const long long Random::b[3] = { 14294630557836824467, 14294630557836824469, 9585689890975426951 };
const long long Random::a[3] = { 15404481456978043987, 12899375936557105357, 7731693896872981757 };

/// <summary>
/// An implementation of complex number (using long double).
/// </summary>
struct Complex {
    Complex() {}
    Complex(long double n) { real = n; }
    Complex(long double a, long double b) { real = a; imaginary = b; }

    bool operator ==(Complex right) {
        return this->real == right.real && this->imaginary == right.imaginary;
    }
    bool operator !=(Complex right) {
        return this->real != right.real || this->imaginary != right.imaginary;
    }

    Complex operator =(Complex right) {
        if (this == &right)
            return *this;
        else {
            real = right.real;
            imaginary = right.imaginary;

            return *this;
        }
    }

    Complex operator +() {
        return *this;
    }
    Complex operator -() {
        Complex Res = Complex();
        Res.real = -real;
        Res.imaginary = -imaginary;
        return Res;
    }

    Complex operator +(Complex b) {
        Complex Res = Complex();

        Res.real = real + b.real;
        Res.imaginary = imaginary + b.imaginary;

        return Res;
    }
    Complex operator -(Complex b) {
        Complex Res = Complex();

        Res.real = real - b.real;
        Res.imaginary = imaginary - b.imaginary;

        return Res;
    }

    Complex operator *(Complex b) {
        Complex Res = Complex();

        Res.real = real * b.real - imaginary * b.imaginary;
        Res.imaginary = b.real * imaginary + real * b.imaginary;

        return Res;
    }
    Complex operator /(Complex b) {
        Complex Res = Complex();

        Res.real = (real * b.real + imaginary * b.imaginary) / (b.real * b.real + b.imaginary * b.imaginary);
        Res.imaginary = (imaginary * b.real - real * b.imaginary) / (b.real * b.real + b.imaginary * b.imaginary);

        return Res;
    }

    Complex operator +=(Complex b) {
        real += b.real;
        imaginary += b.imaginary;

        return *this;
    }
    Complex operator -=(Complex b) {
        real -= b.real;
        imaginary -= b.imaginary;

        return *this;
    }

    Complex operator *=(unsigned int b) {
        *this = *this * Complex(b);
        return *this;
    }
    Complex operator /=(unsigned int b) {
        *this = *this / Complex(b);
        return *this;
    }

    static long double abs(Complex z){
        return sqrtl(z.real * z.real + z.imaginary * z.imaginary);
    }
    std::string toString() {
        if (imaginary < 0)
            return std::to_string(real) + " - " + std::to_string(abs(imaginary)) + 'i';
        else
            return std::to_string(real) + " + " + std::to_string(imaginary) + 'i';
    }

    static Complex conjugated(Complex a) {
        return Complex(a.real, -a.imaginary);
    }

    long double real = 0;
    long double imaginary = 0;

};
//Unsigned ultra long realization
/// <summary>
/// An implementation of unsigned long arithmetic.
/// </summary>
class UltraLong {

    public:
        /// <summary>
        /// UltraLong constructor.
        /// </summary>
        /// <returns>
        /// Zero.
        /// </returns>
        UltraLong() { }
        /// <summary>
        /// UltraLong constructor.
        /// </summary>
        /// <param name = "n">
        /// : The value that initialises UltraLong.
        /// </param>
        UltraLong(unsigned int n) { value[0] = n; }
        UltraLong(SignedUltraLong x) {
            for (unsigned int i = 0; i < LENGTH; i++)
                this->value[i] = x.value[i];
        }
        /*TODO*/
        /// <summary>
        /// UltraLong constructor.
        /// </summary>
        /// <param name = "n">
        /// : The value that initialises UltraLong. If it is less than zero then
        /// UltraLong will be equal ~n + 1, so n + (-n) = 0.
        /// </param>
        UltraLong(int n) {
            if (n < 0) {
                value[0] = (unsigned int)(-n);
                *this = -*this;
            }
            else
                value[0] = (unsigned int)n;
        }
        UltraLong(long long n) {
            if (n < 0) {
                value[0] = (-n) % UINT_RANGE;
                value[1] = (unsigned int)((-n) / UINT_RANGE);
                *this = -*this;
            }
            else
                value[0] = n % UINT_RANGE;
                value[1] = (unsigned int)(n / UINT_RANGE);
        }
        /*TODO*/
        /// <summary>
        /// UltraLong constructor.
        /// </summary>
        /// <param name = "n">
        /// : The value that initialises UltraLong.
        /// </param>
        UltraLong(unsigned long long n) {
            value[0] = n % UINT_RANGE;
            value[1] = (unsigned int)(n / UINT_RANGE);
        }

        bool operator ==(UltraLong right) {
            if (this == &right)
                return true;
            else {
                for (unsigned int i = 0; i < LENGTH; i++)
                    if (this->value[i] != right.value[i]) return false;

                return true;
            }
        }
        bool operator !=(UltraLong right) {
            if (this == &right)
                return false;
            else {
                for (unsigned int i = 0; i < LENGTH; i++)
                    if (this->value[i] != right.value[i]) return true;

                return false;
            }
        }
        bool operator <(UltraLong right) {
            unsigned int i = LENGTH;
            while (i != 0) {
                i--;
                if (this->value[i] < right.value[i]) return true;
                if (this->value[i] > right.value[i]) return false;
            }

            return false;
        }
        bool operator >(UltraLong right){
            unsigned int i = LENGTH;
            while (i != 0) {
                i--;
                if (this->value[i] > right.value[i]) return true;
                if (this->value[i] < right.value[i]) return false;
            }

            return false;
        }
        bool operator >=(UltraLong right) {
            unsigned int i = LENGTH;
            while (i != 0) {
                i--;
                if (this->value[i] > right.value[i]) return true;
                if (this->value[i] < right.value[i]) return false;
            }

            return true;
        }
        bool operator <=(UltraLong right) {
            unsigned int i = LENGTH;
            while (i != 0) {
                i--;
                if (this->value[i] < right.value[i]) return true;
                if (this->value[i] > right.value[i]) return false;
            }

            return true;
        }

        UltraLong operator =(UltraLong right) {
            if (this == &right)
                return *this;
            else {
                for (unsigned int i = 0; i < LENGTH; i++)
                    this->value[i] = right.value[i];

                return *this;
            }
        }

        UltraLong operator -() {

            UltraLong Res = *this;

            for (unsigned int i = 0; i < LENGTH; i++)
                Res.value[i] = ~Res.value[i];

            Res++;

            return Res;

        }
        UltraLong operator +() {
            return *this;
        }
        UltraLong operator ++(int) {

            this->value[0]++;

            for (unsigned int i = 1; i < LENGTH; i++)
            {
                if (this->value[i - 1] != 0)
                    break;

                this->value[i]++;
            }

            return *this;

        }

        UltraLong operator +(UltraLong b) {

            UltraLong Res = UltraLong();

            unsigned int Overflow = 0;

            for (unsigned int i = 0; i < LENGTH; i++)
            {
                Res.value[i] = this->value[i] + b.value[i] + Overflow;

                if (Res.value[i] < this->value[i] || (Overflow == 1 && Res.value[i] == this->value[i]))
                    Overflow = 1;
                else
                    Overflow = 0;
            }

            return Res;

        }
        UltraLong operator -(UltraLong b) {

            UltraLong res = UltraLong();

            unsigned int overflow = 0;

            for (unsigned int i = 0; i < LENGTH; i++)
            {
                res.value[i] = this->value[i] + ~b.value[i] + overflow;

                if (res.value[i] < this->value[i] || (overflow == 1 && res.value[i] == this->value[i]))
                    overflow = 1;
                else
                    overflow = 0;
            }

            res++;

            return res;

        }

        UltraLong operator *(unsigned int b) {

            UltraLong res = UltraLong();

            unsigned long long overflow = 0,
                               overflow1 = 0;

            for (unsigned int i = 0; i < LENGTH; i++)
            {
                overflow1 = (unsigned long long)this->value[i] * b + overflow;

                if (overflow1 < overflow)
                    overflow = overflow1 + 1;
                else
                    overflow = overflow1;

                res.value[i] = overflow % UINT_RANGE;
                overflow /= UINT_RANGE;
            }

            return res;

        }
        UltraLong operator /(unsigned int b) {
            
            UltraLong res = UltraLong();

            unsigned int i = LENGTH;
            unsigned long long prev = 0;
            
            while (i != 0) {
                i--;
            
                prev += this->value[i];
                
                res.value[i] = (unsigned int)(prev / b);
                prev = (prev % b) * UINT_RANGE;
            }

            return res;

        }
        /*UltraLong operator *(UltraLong b) {
            UltraLong Res = UltraLong();
            for (unsigned int i = 0; i < LENGTH; i++)
            {
                for (unsigned int j = 0; j < LENGTH - i; j++)
                {
                    Res += (UltraLong(this->Value[i]) * b.Value[j]).SuperLeftShift(i + j);
                }
            }
            return Res;
        }*/
        /*UltraLong operator *(UltraLong b) {
            unsigned int ThisNumberLength = KaratsubaIterationLength;
            if (ThisNumberLength < 4) {
                UltraLong Res = UltraLong();
                for (unsigned int i = 0; i < ThisNumberLength; i++)
                {
                    for (unsigned int j = 0; j < ThisNumberLength - i; j++)
                    {
                        Res += (UltraLong(this->Value[i]) * b.Value[j]).SuperLeftShift(i + j);
                    }
                }
                return Res;
            }
            unsigned int RightNumberLength = ThisNumberLength / 2;
            unsigned int LeftNumberLength = ThisNumberLength - RightNumberLength;
            UltraLong LeftNumberA = UltraLong();
            UltraLong LeftNumberB = UltraLong();
            for (unsigned int i = 0; i < LeftNumberLength; i++)
            {
                LeftNumberA.Value[i] = this->Value[RightNumberLength + i];
                LeftNumberB.Value[i] = b.Value[RightNumberLength + i];
            }
            UltraLong RightNumberA = UltraLong();
            UltraLong RightNumberB = UltraLong();
            for (unsigned int i = 0; i < RightNumberLength; i++)
            {
                RightNumberA.Value[i] = this->Value[i];
                RightNumberB.Value[i] = b.Value[i];
            }
            KaratsubaIterationLength = LeftNumberLength;
            UltraLong ac = LeftNumberA * LeftNumberB;
            //KaratsubaIterationLength = ThisNumberLength;
            KaratsubaIterationLength = RightNumberLength;
            UltraLong bd = RightNumberA * RightNumberB;
            //KaratsubaIterationLength = ThisNumberLength;
            KaratsubaIterationLength = LeftNumberLength + 1;
            UltraLong abcd = (LeftNumberA + RightNumberA) * (RightNumberB + LeftNumberB);
            KaratsubaIterationLength = ThisNumberLength;
            UltraLong adbc = abcd - ac - bd;
            UltraLong Res = ac.SuperLeftShift(2 * RightNumberLength) + adbc.SuperLeftShift(RightNumberLength) + bd;
            return Res;
        }*/
        UltraLong operator *(UltraLong b) {

            /*if (*this == UltraLong::Zero || b == UltraLong::Zero) {
                return UltraLong();
            }*/


            lastMultOverflow = false;

            UltraLong res = UltraLong();

            unsigned long long temp1[4 * LENGTH] = { 0 },
                               temp2[4 * LENGTH] = { 0 },
                               temp0[4 * LENGTH] = { 0 };

            for (unsigned int i = 0; i < LENGTH; i++)
            {
                temp1[2 * i] = this->value[i] % SUINT_RANGE;
                temp1[2 * i + 1] = this->value[i] / SUINT_RANGE;
                temp2[2 * i] = b.value[i] % SUINT_RANGE;
                temp2[2 * i + 1] = b.value[i] / SUINT_RANGE;
            }

            multiply(temp1, temp2, temp0);

            for (unsigned int i = 1; i < 4 * LENGTH; i++)
            {
                temp0[i] += temp0[i - 1] / SUINT_RANGE;
                temp0[i - 1] %= SUINT_RANGE;
            }

            temp0[2 * LENGTH - 1] %= SUINT_RANGE;


            for (size_t i = 2 * LENGTH; i < 4 * LENGTH; i++)
                if (temp0[i] != 0) {
                    lastMultOverflow = true;
                    break;
                }

            for (unsigned int i = 0; i < LENGTH; i++)
            {
                res.value[i] = (unsigned int)(temp0[2 * i] + temp0[2 * i + 1] * SUINT_RANGE);
            }

            return res;

        }
        UltraLong operator /(UltraLong b) {

            if (b == Zero) { throw "UltraLong division by zero."; }
            if (b == One) { return *this; }

            if (*this < b) { return Zero; }

            std::pair<unsigned int, unsigned int> lead_b = lead(b);

            UltraLong l = (*this / (lead_b.first + 1)).superRightShift(lead_b.second);
            UltraLong r = *this;

            while (r - l > 1)
            {
                UltraLong m = UltraLong::middle(l, r);

                if (m * b > *this || lastMultOverflow) {
                    r = m;
                }
                else {
                    l = m;
                }
            }

            return l;

        }

        unsigned int operator %(unsigned int b) {
            
            unsigned int i = LENGTH;
            unsigned long long prev = 0;

            while (i != 0) {
                i--;
                prev *= UINT_RANGE;
                prev += this->value[i];
                prev = prev % b;
            }

            return (unsigned int)prev;

        }
        UltraLong operator %(UltraLong b) {
            return *this - *this / b * b;
        }

        UltraLong operator ^(UltraLong b) {
            UltraLong res = UltraLong();

            for (unsigned int i = 0; i < LENGTH; i++)
                res.value[i] = this->value[i] ^ b.value[i];

            return res;
        }
        UltraLong operator &(UltraLong b) {
            UltraLong res = UltraLong();

            for (unsigned int i = 0; i < LENGTH; i++)
                res.value[i] = this->value[i] & b.value[i];

            return res;
        }
        UltraLong operator |(UltraLong b) {
            UltraLong res = UltraLong();

            for (unsigned int i = 0; i < LENGTH; i++)
                res.value[i] = this->value[i] | b.value[i];

            return res;
        }
        UltraLong operator ~() {
            UltraLong res = UltraLong();

            for (unsigned int i = 0; i < LENGTH; i++)
                res.value[i] = ~this->value[i];

            return res;
        }

        UltraLong operator +=(UltraLong right) {
            *this = *this + right;
            return *this;
        }
        UltraLong operator -=(UltraLong right) {
            *this = *this - right;
            return *this;
        }

        UltraLong operator *=(UltraLong right) {
            *this = *this * right;
            return *this;
        }
        UltraLong operator *=(unsigned int right) {
            *this = *this * right;
            return *this;
        }
        UltraLong operator /=(UltraLong right) {
            *this = *this / right;
            return *this;
        }

        UltraLong operator %=(UltraLong right) {
            *this = *this % right;
            return *this;
        }

        UltraLong operator ^=(UltraLong right) {
            *this = *this ^ right;
            return *this;
        }
        UltraLong operator &=(UltraLong right) {
            *this = *this & right;
            return *this;
        }
        UltraLong operator |=(UltraLong right) {
            *this = *this | right;
            return *this;
        } 

        UltraLong operator >>(unsigned long long right) {

            unsigned long long n = right;

            UltraLong res = *this;
            res = res.superRightShift((unsigned int)(n / BITS_IN_UINT));

            n %= BITS_IN_UINT;
            unsigned int mod = 1u << n;
            res.value[0] >>= n;

            for (unsigned int i = 1; i < LENGTH; i++)
            {
                res.value[i - 1] += (res.value[i] % mod) << (BITS_IN_UINT - n);
                res.value[i] >>= n;
            }

            return res;

        }
        UltraLong operator <<(unsigned long long right) {

            unsigned long long n = right;

            UltraLong res = *this;
            res = res.superLeftShift((unsigned int)(n / BITS_IN_UINT));

            n %= BITS_IN_UINT;
            unsigned int mod = 1u << (BITS_IN_UINT - n);

            for (unsigned int i = LENGTH - 1; i > 0; i--)
            {
                res.value[i] <<= n;
                res.value[i] += res.value[i - 1] / mod;
            }

            return res;

        }

        bool isEven() {
            return this->value[0] % 2 == 0;
        }
        bool isOdd() {
            return this->value[0] == 1;
        }

        /// <summary>
        /// Modular exponentiation.
        /// </summary>
        /// <returns>The remainder after dividing a to b-th power by mod</returns>
        static UltraLong modPow(UltraLong a, UltraLong b, UltraLong mod) {
            if (mod == 0) { throw "UltraLong modulo by zero."; }
            if (mod == One) {
                return Zero;
            }
            if (a == Zero || a == One) {
                return a;
            }
            if (b == 0) {
                return 1;
            }
            if (b.isEven()) {
                return modPow(a * a % mod, b.rightShift(1), mod);
            }
            else {
                return modPow(a * a % mod, b.rightShift(1), mod) * a % mod;
            }

        }

        /// <summary>
        /// Cast UltraLong to unsigned int.
        /// </summary>
        /// <returns>The remainder after dividing UltraLong by 2**32.</returns>
        unsigned int toUint() {
            return this->value[0];
        }
        /// <summary>
        /// Cast UltraLong to unsigned long long.
        /// </summary>
        /// <returns>The remainder after dividing UltraLong by 2**64.</returns>
        unsigned long long toUlong() {
            return LENGTH > 1 ? this->value[1] * UINT_RANGE + this->value[0] : this->value[0];
        }

        long double toLongDouble() {
            long double ans = 0;
            unsigned int i = LENGTH;
            while (i != 0)
            {
                i--;
                ans += this->value[i] * pow(UINT_RANGE, i);
            }
            return ans;
        }

        /// <summary>
        /// Cast UltraLong to string.
        /// </summary>
        /// <returns>UltraLong in decimal botation.</returns>
        std::string toString() {

            std::string s = "";

            UltraLong n = *this;

            while (n > 0)
            {
                unsigned int digit = n % 10;
                s += std::to_string(digit);
                n = n / 10;
            }

            if (s == "") s = "0";

            std::string s1 = "";
            unsigned long long i = s.size();
            while (i != 0)
            {
                i--;
                s1 += s[i];
            }

            return s1;

        }

        static long double divide(UltraLong a, UltraLong b) {
            return a.toLongDouble() / b.toLongDouble();
        }

        static long double log(UltraLong n) {
            if (n == 0) {
                throw "Logarithm of zero is minus infinity.";
            }
            long double x = divide(n - 1, n + 1);
            long double x_2 = x * x;
            long double ans = 0;
            unsigned int s = 1;
            for (unsigned int i = 0; i < LOG_PRECISE; i++)
            {
                ans += x / s;
                x *= x_2;
                s += 2;
            }
            return ans * 2;
        }

        static UltraLong sqrt(UltraLong n) {

            if (n < 2) return n;

            UltraLong l = 0, r = n;

            while (r - l > 1) {

                UltraLong m = (l + r) / 2;

                if (m * m > n)
                    r = m;
                else
                    l = m;

            }

            return l;

        }

        static UltraLong parse(std::string s) {
            UltraLong ans = UltraLong();
            unsigned long long size = s.size();
            for (unsigned long long i = 0; i < size; i++)
            {
                ans *= 10;
                ans += s[i] - '0';
            }
            return ans;
        }

        friend class SignedUltraLong;

    protected:

        static const unsigned int LENGTH = 2;
        static const unsigned int UPPER_BOUND_LENGTH = 2;

        static const unsigned long long UINT_RANGE = 4294967296;
        static const unsigned long long SUINT_RANGE = 65536;
        static const unsigned int BITS_IN_UINT = 32;
        static const unsigned int LOG_PRECISE = 10;

        static long double PI;
        static bool precalc;

        //static inline unsigned int KaratsubaIterationLength = LENGTH;

        unsigned int value[LENGTH] = { 0 };

        static UltraLong One;
        static UltraLong Zero;
        static UltraLong MinusOne;

        static unsigned int rev[4 * UPPER_BOUND_LENGTH];
        static Complex wlen_pw[4 * UPPER_BOUND_LENGTH];

        static bool lastMultOverflow;

        static UltraLong middle(UltraLong a, UltraLong b) {
            return a + (b - a) / 2;
        }

        static std::pair<unsigned int, unsigned int> lead(UltraLong a) {
            unsigned int i = LENGTH;
            while (i != 0) {
                i--;
                if (a.value[i] != 0)
                    return { a.value[i], i };
            }
            return { 0, 0 };
        }

        static void swap(Complex* a, Complex* b) {

            if (a == b) { return; }

            Complex c = Complex();

            c = *a;
            *a = *b;
            *b = c;

        }

        static void calcRev(unsigned int n, unsigned int log_n) {

            for (unsigned int i = 0; i < n; ++i) {

                rev[i] = 0;

                for (unsigned int j = 0; j < log_n; ++j)
                    if (i & (1 << j))
                        rev[i] |= 1 << (log_n - 1 - j);

            }

        }

        static void multiply(const unsigned long long a[], const unsigned long long b[], unsigned long long res[]) {

            Complex fa[4 * UPPER_BOUND_LENGTH] = { Complex() },
                    fb[4 * UPPER_BOUND_LENGTH] = { Complex() }, 
                    fc[4 * UPPER_BOUND_LENGTH] = { Complex() }, 
                    fd[4 * UPPER_BOUND_LENGTH] = { Complex() };

            for (unsigned int i = 0; i < 2 * LENGTH; i++)
            {
                fa[i].real = (long double)a[i];
                fb[i].real = (long double)b[i];
                fc[i].real = fa[i].real;
                fc[i].imaginary = fb[i].real;
            }

            unsigned int n = 4 * UPPER_BOUND_LENGTH;

            fastFourierTransformation(fa, n, false),
            fastFourierTransformation(fb, n, false);

            fastFourierTransformation(fc, n, false);

            for (size_t i = 0; i < n; ++i) {
                Complex fcj = fc[i];
                Complex fcnj = Complex::conjugated(fc[4 * UPPER_BOUND_LENGTH - i - 1]);
                fd[i] = (fcj + fcnj);
            }

            for (size_t i = 0; i < n; ++i)
                fa[i] = fa[i] * fb[i];

            fastFourierTransformation(fa, n, true);

            for (size_t i = 0; i < 4 * LENGTH; ++i)
                res[i] = unsigned long long(fa[i].real + 0.5);
            
        }

        static void fastFourierTransformation(Complex* a, unsigned int n, bool reverse = false) {

            if (!precalc) {
                precalc = true;
                unsigned int lg_n = 0;
                unsigned int t_lg = 1;
                while (t_lg < n) { ++lg_n; t_lg <<= 1; }
                calcRev(n, lg_n);
            }

                for (unsigned int i = 0; i < n; ++i)
                    if (i < rev[i])
                        swap(&a[i], &a[rev[i]]);


            for (unsigned int len = 2; len <= n; len <<= 1) {

                long double ang = 2 * PI / len * (reverse ? -1 : +1);
                unsigned int len2 = len >> 1;

                Complex wlen(cos(ang), sin(ang));
                wlen_pw[0] = Complex(1, 0);
                for (unsigned int i = 1; i < len2; ++i)
                    wlen_pw[i] = wlen_pw[i - 1] * wlen;

                for (unsigned int i = 0; i < n; i += len) {
                    Complex t,
                           *pu = a + i,
                           *pv = a + i + len2,
                           *pu_end = a + i + len2,
                           *pw = wlen_pw;
                    for (; pu != pu_end; ++pu, ++pv, ++pw) {
                        t = *pv * *pw;
                        *pv = *pu - t;
                        *pu += t;
                    }
                }

            }

            if (reverse)
                for (unsigned int i = 0; i < n; ++i)
                    a[i] /= n;
        }

        UltraLong superLeftShift(unsigned int n) {
            if (n > LENGTH - 1)
                return 0;

            UltraLong res = UltraLong();

            unsigned int i = LENGTH - n;
            while (i != 0) {
                i--;
                res.value[i + n] = this->value[i];
            } 

            return res;

        }
        UltraLong superRightShift(unsigned int n) {

            UltraLong res = UltraLong();

            for (unsigned int i = 0; i < LENGTH - n; i++)
                res.value[i] = this->value[i + n];

            for (unsigned int i = LENGTH - n; i < LENGTH; i++)
                res.value[i] = 0;

            return res;

        }
        UltraLong rightShift(unsigned long long n) {

            UltraLong res = *this;
            res.superRightShift((unsigned int)(n / BITS_IN_UINT));

            n %= BITS_IN_UINT;
            unsigned int mod = 1u << n;
            res.value[0] >>= n;

            for (unsigned int i = 1; i < LENGTH; i++)
            {
                res.value[i - 1] += (res.value[i] % mod) << (BITS_IN_UINT - n);
                res.value[i] >>= n;
            }

            return res;

        }
};

UltraLong UltraLong::One = UltraLong(1);
UltraLong UltraLong::Zero = UltraLong(0);
UltraLong UltraLong::MinusOne = UltraLong(-1);
long double UltraLong::PI = 3.1415926535897932389l;
bool UltraLong::precalc = false;
unsigned int UltraLong::rev[4 * UltraLong::UPPER_BOUND_LENGTH];
Complex UltraLong::wlen_pw[4 * UltraLong::UPPER_BOUND_LENGTH];
bool UltraLong::lastMultOverflow = false;

class SignedUltraLong : UltraLong {
    SignedUltraLong(UltraLong x) {
        for (unsigned int i = 0; i < LENGTH; i++)
            this->value[i] = x.value[i];
    }
    bool operator <(SignedUltraLong right) {
        if (this->isNegative()) {
            if (right.isNegative())
                return (-*this) > (-right);
            return true;
        }
        if (right.isNegative()) {
            return false;
        }
        unsigned int i = LENGTH;
        while (i != 0) {
            i--;
            if (this->value[i] < right.value[i]) return true;
            if (this->value[i] > right.value[i]) return false;
        }

        return false;
    }
    bool operator >(SignedUltraLong right) {
        if (this->isNegative()) {
            if (right.isNegative())
                return (-*this) < (-right);
            return false;
        }
        if (right.isNegative()) {
            return true;
        }
        unsigned int i = LENGTH;
        while (i != 0) {
            i--;
            if (this->value[i] > right.value[i]) return true;
            if (this->value[i] < right.value[i]) return false;
        }

        return false;
    }
    bool operator >=(SignedUltraLong right) {
        if (this->isNegative()) {
            if (right.isNegative())
                return (-*this) <= (-right);
            return false;
        }
        if (right.isNegative()) {
            return true;
        }
        unsigned int i = LENGTH;
        while (i != 0) {
            i--;
            if (this->value[i] > right.value[i]) return true;
            if (this->value[i] < right.value[i]) return false;
        }

        return true;
    }
    bool operator <=(SignedUltraLong right) {
        if (this->isNegative()) {
            if (right.isNegative())
                return (-*this) >= (-right);
            return true;
        }
        if (right.isNegative()) {
            return false;
        }
        unsigned int i = LENGTH;
        while (i != 0) {
            i--;
            if (this->value[i] < right.value[i]) return true;
            if (this->value[i] > right.value[i]) return false;
        }

        return true;
    }
    SignedUltraLong operator /(unsigned int b) {

        if (this->isNegative()) {
            return -(-*this / b);
        }

        UltraLong res = UltraLong();

        unsigned int i = LENGTH;
        unsigned long long prev = 0;

        while (i != 0) {
            i--;

            prev += this->value[i];

            res.value[i] = (unsigned int)(prev / b);
            prev = (prev % b) * UINT_RANGE;
        }

        return res;

    }

    static SignedUltraLong parse(std::string s) {
        if (s == "")
            return SignedUltraLong(0);
        if (s[0] == '-')
            return -UltraLong::parse(s.substr(1));
        return UltraLong::parse(s);
    }

    std::string toString() {
        if (this->isNegative())
            return "-" + UltraLong(-*this).toString();
    }
protected:
    bool isNonnegative() {
        return this->value[LENGTH - 1] < (UINT_RANGE / 2);
    }
    bool isNegative() {
        return this->value[LENGTH - 1] >= (UINT_RANGE / 2);
    }
};

namespace Math {

    enum class PRIME_TESTS_OPTIMISE_LEVELS { O1, O2, O3 };
    const PRIME_TESTS_OPTIMISE_LEVELS PRIME_TESTS_OPTIMISE_LEVEL = PRIME_TESTS_OPTIMISE_LEVELS::O3;

    template <class T>
    T sqrt(T n) {
        return (T)sqrtl((long double)n);
    }

    template <>
    UltraLong sqrt(UltraLong n) {
        return UltraLong::sqrt(n);
    }

    template <class T>
    T MillerUpperBound_O1(T n) {
        return 2 * powl(logl(n), 2);
    }

    template <>
    UltraLong MillerUpperBound_O1<UltraLong>(UltraLong n) {
        return (unsigned long long)floorl(2 * powl(UltraLong::log(n), 2));
    }

    template <class T>
    T MillerUpperBound_O2(T n) {
        long double ans = logl(n);
        return  ans * logl(ans) / logl(2);
    }

    template <>
    UltraLong MillerUpperBound_O2<UltraLong>(UltraLong n) {
        long double ans = UltraLong::log(n);
        return (unsigned long long)floorl(ans * logl(ans) / logl(2));
    }

    template <class T>
    bool IsPrimePower(T n) {
        return false;
    }

    template <class T>
    T ModPow(T a, T b, T mod) {
        if (mod == 0) { throw "ModPow: mod is zero"; }
        if (a == 0 || mod == 1) { return 0; }
        if (b == 0) { return 1; }
        if (b % 2 == 1) { return ModPow<T>(a * a % mod, b / 2, mod) * a % mod; }
        return ModPow<T>(a * a % mod, b / 2, mod);
    }

    template <>
    UltraLong ModPow<UltraLong>(UltraLong a, UltraLong b, UltraLong mod) {
        return UltraLong::modPow(a, b, mod);
    }

    template <class T>
    bool IsStrongPseudoprime(T n, T base) {

        T exp = n - 1;
        T ost_1 = exp;
        T res = ModPow<T>(base, exp, n);

        if (res != 1)
            return false;

        while (true)
        {
            exp = exp / 2;
            res = ModPow<T>(base, exp, n);

            if (res == ost_1)
                return true;

            if (exp % 2 == 1)
            {
                res = ModPow<T>(base, exp, n);
                if (res == 1)
                    return true;

                break;
            }
        }

        return false;

    }

    template <class T>
    bool SmallIsPrime(T base) {

        if (base < 2) return false;

        T sq_b = sqrt<T>(base) + 1;

        for (T i = 2; i < sq_b; i++)
            if (base % i == 0)
                return false;

        return true;

    }

    template <class T>
    T SmallNextPrime(T base) {

        T x = base + 1;

        while (!SmallIsPrime<T>(x))
        {
            x++;
        }

        return x;

    }

    template <class T>
    bool MillerTest(T n, PRIME_TESTS_OPTIMISE_LEVELS optimise_level = PRIME_TESTS_OPTIMISE_LEVEL) {

        if (IsPrimePower<T>(n)) {
            return false;
        }

        T upper_bound;

        switch (optimise_level)
        {
        case PRIME_TESTS_OPTIMISE_LEVELS::O1: {
            upper_bound = MillerUpperBound_O1(n);
            break;
        }
        default:
            upper_bound = MillerUpperBound_O2(n);
            break;
        }

        T base = 2;

        while (base <= n) {

            if (!IsStrongPseudoprime<T>(n, upper_bound)) {
                return false;
            }

            base = SmallNextPrime<T>(base);

        }

        return true;

    }

    //TODO
    template <class T>
    T nextPrime(T n, PRIME_TESTS_OPTIMISE_LEVELS optimize_level = PRIME_TESTS_OPTIMISE_LEVEL) {

        if (n == (T)2) return (T)3;
        if (n < (T)2) return (T)2;

        T x = n + 2;

        while (!MillerTest(x, optimize_level))
            x += 2;

        return x;

    }

}


template <class T>
class Sequence {
    public:
        T next();
        Sequence<T> renew();
    private:
        T start;
        T now;
        unsigned long long iter;
};

template <class T>
class Arithmetic : public Sequence<T> {
    public:
        Arithmetic<T>(T start, T difference) {
            this->start = start;
            this->now = start;
            this->difference = difference;
            iter = 0;
        }
        T next() {
            iter++;
            T now1 = now;
            now += difference;
            return now1;
        }
        Arithmetic<T> renew() {
            return Arithmetic<T>(start, difference);
        }
    private:
        T start;
        T now;
        T difference;
        unsigned long long iter;
};

template <class T>
class Geometric : public Sequence<T> {
    public:
        Geometric(T start, T denominator) {
            this->start = start;
            this->now = start;
            this->denominator = denominator;
            iter = 0;
        }
        T next() {
            iter++;
            T now1 = now;
            now *= denominator;
            return now1;
        }
        Geometric<T> renew() {
            return Geometric<T>(start, denominator);
        }
    private:
        T start;
        T now;
        T denominator;
        unsigned long long iter;
};

template <class T>
class Fibonacci : public Sequence<T> {
    public:
        Fibonacci(T start = 0, T start1 = 1) {
            this->start = start;
            this->start1 = start1;
            iter = 0;
        }
        T next() {
            iter++;

            switch (iter) {
                case 1: { /*now = start;*/ return start; }
                case 2: { now = start1; prev = start; return start1; }
            }

            T now1 = now + prev;
            prev = now;
            now = now1;

            return now;
        }
        Fibonacci<T> renew() {
            return Fibonacci<T>(start, start1);
        }
    private:
        T start;
        T start1;
        T prev;
        T now;
        unsigned long long iter;
};

template <class T>
class Primes : public Sequence<T> {
public:
    Primes() {
        this->start = T(2);
        this->now = this->start;
        iter = 0;
    }
    T next() {
        iter++;

        T now1 = now;
        now = Math::nextPrime<T>(now);
        return now1;
    }
    Primes<T> renew() {
        return Primes<T>();
    }
private:
    T start;
    T prev;
    T now;
    unsigned long long iter;
};
//Arithmetic<int> nat = Arithmetic<int>(1, 1);

struct Permutation {
    //Gives identity permutation with length n
    Permutation(unsigned int length) {
        len = length;
        x.resize(len);
        for (unsigned int i = 0; i < len; i++)
            x[i] = i;
    }
    /*Permutation(int t[]) {
        int n = sizeof(t) / sizeof(int);
        for (int i = 0; i < n; i++){
            x[i] = t[i];
        }
    }*/

    Permutation operator *(Permutation b) {
        if (this->len != b.len)
            throw "Can't multiply permutations with different lengths.";

        Permutation ans = Permutation(this->len, 0);

        for (unsigned int i = 0; i < this->len; i++)
            ans.x[i] = b.x[this->x[i]];

        return ans;
    }

    Permutation inverse() {
        Permutation ans = Permutation(this->len, 0);
        for (unsigned int i = 0; i < this->len; i++)
        {
            ans.x[this->x[i]] = i;
        }
    }
private:
    unsigned int len;
    std::vector<unsigned int> x;
    Permutation(unsigned int length, unsigned int value) {
        len = length;
        x.assign(length, value);
    }
};

template <class T>
struct Matrix {

    Matrix<T>(const std::vector<std::vector<T>> &arr) {

        mx = arr;
        rows = mx.size();
        columns = rows == 0 ? 0 : mx[0].size();

        for (unsigned int i = 1; i < rows; i++)
            if (mx[i].size() < columns)
                columns = mx[i].size();
    
    }

    Matrix<T>(unsigned int r, unsigned int c, T &value = 0) {

        mx.resize(r);
        
        for (unsigned int i = 0; i < r; i++)
            mx[i].assign(c, value);
        
        columns = c;
        rows = r;

    }

    Matrix<T> operator =(Matrix<T> right) {
        if (this == &right)
            return *this;
        else {
            this->rows = right.rows();
            this->columns = right.columns();
            this->mx = right.mx();

            return *this;
        }
    }

    Matrix<T> operator +(Matrix<T> b) {

        if (this->columns != b.columns || this->rows != b.rows)
            throw "Matrixes sizes don't match (operator +).";

        Matrix<T> ans = Matrix<T>(b.rows, b.columns);

        for (unsigned int i = 0; i < b.rows; i++)
            for (unsigned int j = 0; j < b.columns; j++)
                ans.mx[i][j] = this->mx[i][j] + b.mx[i][j];

        return ans;

    }

    Matrix<T> operator -(Matrix<T> b) {

        if (this->columns != b.columns || this->rows != b.rows)
            throw "Matrixes sizes don't match (operator -).";

        Matrix<T> ans = Matrix<T>(b.rows, b.columns);

        for (unsigned int i = 0; i < b.rows; i++)
            for (unsigned int j = 0; j < b.columns; j++)
                ans.mx[i][j] = this->mx[i][j] - b.mx[i][j];

        return ans;

    }

    Matrix<T> operator *(Matrix<T> b) {

        if (this->columns != b.rows)
            throw "Matrixes sizes don't match (operator *).";

        Matrix<T> ans = Matrix<T>(this->rows, b.columns);

        for (unsigned int i = 0; i < this->rows; i++)
            for (unsigned int j = 0; j < b.columns; j++) 
                for (unsigned int k = 0; k < b.rows; k++)
                    ans.mx[i][j] += this->mx[i][k] * b.mx[k][j];

        return ans;

    }

    Matrix<T> operator *(T b) {

        Matrix<T> ans = Matrix<T>(rows, columns);

        for (unsigned int i = 0; i < rows; i++)
            for (unsigned int j = 0; j < columns; j++)
                ans.mx[i][j] = this->mx[i][j] * b;

        return ans;

    }

    static Matrix<T> randomMatrix(unsigned int columns, unsigned int rows, T(&generator)(void)) {

        Matrix<T> ans = Matrix<T>(rows, columns);

        for (unsigned int i = 0; i < rows; i++)
            for (unsigned int j = 0; j < columns; j++)
                ans.mx[i][j] = generator();

        return ans;

    }

    static Matrix<T> transpose(Matrix<T> a) {

        Matrix<T> ans = Matrix(a.columns, a.rows);

        for (unsigned int i = 0; i < a.rows; i++)
            for (unsigned int j = 0; j < a.columns; j++)    
                ans.mx[j][i] = a.mx[i][j];

        return ans;

    }

    static T determinant(Matrix<T> a) {

        if (a.columns != a.rows)
            throw "You can't found the determinant of non-square matrix.";

        T ans = (T)0;
        
        for (unsigned int i = 0; i < a.columns; i++)
        {
            T pre_ans = 1;
            for (unsigned int j = 0; j < a.rows; j++)
                pre_ans *= a.mx[j][(i + j) % a.columns];
            ans += pre_ans;
        }

        for (unsigned int i = 0; i < a.columns; i++)
        {
            T pre_ans = -1;
            for (unsigned int j = 0; j < a.rows; j++)
                pre_ans *= a.mx[a.rows - 1 - j][(i + j) % a.columns];
            ans += pre_ans;
        }

        return ans;

    }

    void print() {
        for (unsigned int i = 0; i < rows; i++){
            for (unsigned int j = 0; j < columns; j++)
                std::cout << mx[i][j] << ' ';
            std::cout << '\n';
        }
    }

    Matrix<T> complementaryMinor() {}

private:
    unsigned int rows, columns;
    std::vector<std::vector<T>> mx;

};

struct Polynomial {

public:

    Polynomial() { degree = 0; }

    Polynomial(std::vector<long double> coefficients) {

        unsigned int index = (unsigned int)coefficients.size(), 
                     power = 0;

        degree = index - 1;

        while (index != 0) {
            index--;

            if (coefficients[index] != 0)
                cfs[power] = coefficients[index];

            power++;
        }

    }

    Polynomial operator +(Polynomial b) {
        
        Polynomial ans = Polynomial();
        ans.cfs = this->cfs;
        
        for (std::pair<unsigned int, long double> each : b.cfs)
            if (ans.cfs.find(each.first) != ans.cfs.end())
                ans.cfs[each.first] += each.second;
            else
                ans.cfs[each.first] = each.second;
        
        
        ans.degree = this->degree > b.degree ? this->degree : b.degree;
        
        return ans;

    }

    Polynomial operator -(Polynomial b) {

        Polynomial ans = Polynomial();
        ans.cfs = this->cfs;
        
        for (std::pair<unsigned int, long double> each : b.cfs)
            if (ans.cfs.find(each.first) != ans.cfs.end())
                ans.cfs[each.first] -= each.second;
            else
                ans.cfs[each.first] = -each.second;
        
        ans.ClearZeroes();

        return ans;
        //TODO : degree
    }

    long double Value(long double x) {
        long double ans = 0, px = 1;
        for (unsigned int i = 0; i < degree + 1; i++){
            ans += px * cfs[i];
            px *= x;
        }
        return ans;
    }

    static Polynomial Derivative(Polynomial a) {
        
        Polynomial ans = Polynomial();
        
        for (std::pair<unsigned int, long double> each : a.cfs)
            if (each.first != 0)
                ans.cfs[each.first - 1] = each.second * each.first;
        
        ans.degree = std::max(0u, a.degree - 1);
        
        return ans;

    }

    struct Solutions
    {
    public:
        Solutions() {
            any_number = false;
        }
        bool IsAnyNumber() {
            return any_number;
        }
        void Add(long double sol) {
            ans.push_back(sol);
        }
        void SetAnyNumber(bool an) {
            any_number = an;
        }
        std::vector<long double> GetSolutions() {
            return ans;
        }
    protected:
        std::vector<long double> ans;
        bool any_number;
        
    };
    
    static Solutions Solve(Polynomial a) {
        switch (a.degree) {

        case 0: {
            Solutions ans = Solutions();
            ans.SetAnyNumber(true);
            return ans;
        }

        case 1:
            return SolveLinear(a);
        case 2:
            return SolveQuadratic(a);
        
        default: {
            Solutions ans = Solutions();
            std::vector<long double> der_zeroes = Solve(Polynomial::Derivative(a)).GetSolutions();
            if (der_zeroes.size() == 0){
                long double x = a.FindCriticalZero(-MAX_ROOT, MAX_ROOT);
                if (a.IsNearRoot(x))
                    ans.Add(x);
            }
            else{
                long double x = a.FindCriticalZero(-MAX_ROOT, der_zeroes[0]);
                if (a.IsNearRoot(x))
                    ans.Add(x);

                for (unsigned int i = 0; i < der_zeroes.size() - 1; i++){
                    x = a.FindZero(der_zeroes[i], der_zeroes[i + 1]);
                    if (a.IsNearRoot(x))
                        ans.Add(x);
                }

                x = a.FindCriticalZero(der_zeroes[der_zeroes.size() - 1], MAX_ROOT);
                if (a.IsNearRoot(x))
                    ans.Add(x);
            }
            return ans;
        }
        
        }    
    }

    void Print(std::string var) {

        bool first = true;
        std::string s = "";

        for (auto iter = this->cfs.rbegin(); iter != this->cfs.rend(); ++iter) {

            if (!first && iter->second > 0)
                s += '+';

            s += std::to_string(iter->second) + var + '^' + std::to_string(iter->first);
            first = false;

        }

        std::cout << s << '\n';

    }

private:

    std::map<unsigned int, long double> cfs;
    unsigned int degree;

    void ClearZeroes() {
        for (const std::pair<unsigned int, long double> &each : this->cfs)
            if (each.second == 0)
                this->cfs.erase(this->cfs.find(each.first));
    }

    static Solutions SolveLinear(Polynomial a) {
        Solutions ans = Solutions();
        ans.Add(-a.cfs[1]/a.cfs[0]);
        return ans;
    }

    static Solutions SolveQuadratic(Polynomial a) {
        Solutions ans = Solutions();
        long double discriminant = a.cfs[1] * a.cfs[1] - 4 * a.cfs[0] * a.cfs[2];
        
        if (discriminant > 0) {
            if (a.cfs[2] > 0) {
                ans.Add((-a.cfs[1] - sqrtl(discriminant)) / (2 * a.cfs[2]));
                ans.Add((-a.cfs[1] + sqrtl(discriminant)) / (2 * a.cfs[2]));
            }
            else {
                ans.Add((-a.cfs[1] + sqrtl(discriminant)) / (2 * a.cfs[2]));
                ans.Add((-a.cfs[1] - sqrtl(discriminant)) / (2 * a.cfs[2]));
            }
        }
        else if (discriminant == 0) {
            ans.Add(-a.cfs[1] / (2 * a.cfs[2]));
        }

        return ans;
    }

    long double FindZero(long double L, long double R) {
        if (Value(R) > Value(L))
            while (R - L > SOLVE_PRECISE){
                long double M = (L + R) / 2;
                if (Value(M) > 0)
                    R = M;
                else
                    L = M;
            }
        else
            while (R - L > SOLVE_PRECISE) {
                long double M = (L + R) / 2;
                if (Value(M) < 0)
                    R = M;
                else
                    L = M;
            }
        return L;
    }
    long double FindCriticalZero(long double L, long double R) {
        if (cfs[degree] * powl(L, degree) < cfs[degree] * powl(R, degree))
            while (R - L > SOLVE_PRECISE) {
                long double M = (L + R) / 2;
                if (Value(M) > 0 || isnan(Value(M)))
                    R = M;
                else
                    L = M;
            }
        else
            while (R - L > SOLVE_PRECISE) {
                long double M = (L + R) / 2;
                if (Value(M) < 0 || isnan(Value(M)))
                    R = M;
                else
                    L = M;
            }
        return L;
    }

    bool IsNearRoot(long double x) {
        return Value(x - SOLVE_PRECISE) * Value(x + SOLVE_PRECISE) <= 0;
    }

    static long double SOLVE_PRECISE;
    static long double MAX_ROOT;
};

long double Polynomial::SOLVE_PRECISE = 1e-13;
long double Polynomial::MAX_ROOT = 1e300;