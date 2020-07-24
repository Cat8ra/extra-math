#pragma once

#include <cmath>
#include <string>
#include <ctime>

/// <summary>
/// An implementation of complex number (using long double).
/// </summary>
struct Random {
    Random() { srand(time(NULL)); }
    Random(unsigned int n) { srand(n); }
    int next_int() {
        return rand() + rand() * RAND_MAX + rand() * RAND_MAX * RAND_MAX;
    }
    unsigned int next_uint() {
        return (unsigned int)next_int();
    }
    unsigned long long next_long() {
        return next_int() * 4294967296 + next_int();
    }
    unsigned long long next_ulong() {
        return (unsigned long long)next_long();
    }
};

struct Complex {
    Complex() {}
    Complex(long double n) { Real = n; }
    Complex(long double a, long double b) { Real = a; Imaginary = b; }
    
    bool operator ==(Complex right) {
        return this->Real == right.Real && this->Imaginary == right.Imaginary;
    }
    bool operator !=(Complex right) {
        return this->Real != right.Real || this->Imaginary != right.Imaginary;
    }

    Complex operator =(Complex right) {
        if (this == &right)
            return *this;
        else {
            Real = right.Real;
            Imaginary = right.Imaginary;

            return *this;
        }
    }

    Complex operator +() {
        return *this;
    }
    Complex operator -() {
        Complex Res = Complex();
        Res.Real = -Real;
        Res.Imaginary = -Imaginary;
        return Res;
    }
    
    Complex operator +(Complex b) {
        Complex Res = Complex();

        Res.Real = Real + b.Real;
        Res.Imaginary = Imaginary + b.Imaginary;

        return Res;
    }
    Complex operator -(Complex b) {
        Complex Res = Complex();

        Res.Real = Real - b.Real;
        Res.Imaginary = Imaginary - b.Imaginary;
        
        return Res;
    }

    Complex operator *(Complex b) {
        Complex Res = Complex();

        Res.Real = Real * b.Real - Imaginary * b.Imaginary;
        Res.Imaginary = b.Real * Imaginary + Real * b.Imaginary;

        return Res;
    }
    Complex operator /(Complex b) {
        Complex Res = Complex();

        Res.Real = (Real * b.Real + Imaginary * b.Imaginary) / (b.Real * b.Real + b.Imaginary * b.Imaginary);
        Res.Imaginary = (Imaginary * b.Real - Real * b.Imaginary) / (b.Real * b.Real + b.Imaginary * b.Imaginary);

        return Res;
    }

    Complex operator +=(Complex b) {
        Real += b.Real;
        Imaginary += b.Imaginary;

        return *this;
    }
    Complex operator -=(Complex b) {
        Real -= b.Real;
        Imaginary -= b.Imaginary;

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

    static Complex Abs(Complex z){
        return sqrt(z.Real * z.Real + z.Imaginary * z.Imaginary);
    }
    std::string ToString() {
        if (Imaginary < 0)
            return std::to_string(Real) + " - " + std::to_string(abs(Imaginary)) + 'i';
        else
            return std::to_string(Real) + " + " + std::to_string(Imaginary) + 'i';
    }

    long double Real = 0;
    long double Imaginary = 0;

};

//Unsigned ultra long realization
/// <summary>
/// An implementation of unsigned long arithmetic.
/// </summary>
struct UltraLong {

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
        UltraLong(unsigned int n) { Value[0] = n; }
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
                Value[0] = (unsigned int)(-n); 
                *this = -*this; 
            }
            else 
                Value[0] = (unsigned int)n;
        } 
        /*TODO*/
        /// <summary>
        /// UltraLong constructor.
        /// </summary>
        /// <param name = "n">
        /// : The value that initialises UltraLong.
        /// </param>
        UltraLong(unsigned long long n) {
            Value[0] = n % UINT_RANGE;
            Value[1] = (unsigned int)(n / UINT_RANGE);
        }

        bool operator ==(UltraLong right) {
            if (this == &right)
                return true;
            else {
                for (unsigned int i = 0; i < LENGTH; i++)
                    if (this->Value[i] != right.Value[i]) return false;

                return true;
            }
        }
        bool operator !=(UltraLong right) {
            if (this == &right)
                return false;
            else {
                for (unsigned int i = 0; i < LENGTH; i++)
                    if (this->Value[i] != right.Value[i]) return true;

                return false;
            }
        }
        bool operator <(UltraLong right) {
            for (unsigned int i = LENGTH - 1; i < LENGTH; i--){
                if (this->Value[i] < right.Value[i]) return true;
                if (this->Value[i] > right.Value[i]) return false;
            }

            return false;
        }
        bool operator >(UltraLong right) {
            for (unsigned int i = LENGTH - 1; i < LENGTH; i--){
                if (this->Value[i] > right.Value[i]) return true;
                if (this->Value[i] < right.Value[i]) return false;
            }

            return false;
        }
        bool operator >=(UltraLong right) {
            for (unsigned int i = LENGTH - 1; i < LENGTH; i--){
                if (this->Value[i] > right.Value[i]) return true;
                if (this->Value[i] < right.Value[i]) return false;
            }

            return true;
        }
        bool operator <=(UltraLong right) {
            for (unsigned int i = LENGTH - 1; i < LENGTH; i--) {
                if (this->Value[i] < right.Value[i]) return true;
                if (this->Value[i] > right.Value[i]) return false;
            }

            return true;
        }

        UltraLong operator =(UltraLong right) {
            if (this == &right)
                return *this;
            else {
                for (unsigned int i = 0; i < LENGTH; i++)
                    this->Value[i] = right.Value[i];

                return *this;
            }
        }
        
        UltraLong operator -() {

            UltraLong Res = *this;

            for (unsigned int i = 0; i < LENGTH; i++)
                Res.Value[i] = ~Res.Value[i];

            Res++;

            return Res;

        }
        UltraLong operator +() {
            return *this;
        }
        UltraLong operator ++(int) {

            this->Value[0]++;

            for (unsigned int i = 1; i < LENGTH; i++)
            {
                if (this->Value[i - 1] != 0)
                    break;

                this->Value[i]++;
            }

            return *this;

        }
        
        UltraLong operator +(UltraLong b) {

            UltraLong Res = UltraLong();

            unsigned int Overflow = 0;

            for (unsigned int i = 0; i < LENGTH; i++)
            {
                Res.Value[i] = this->Value[i] + b.Value[i] + Overflow;

                if (Res.Value[i] < this->Value[i] || (Overflow == 1 && Res.Value[i] == this->Value[i]))
                    Overflow = 1;
                else
                    Overflow = 0;
            }

            return Res;

        }
        UltraLong operator -(UltraLong b) {

            UltraLong Res = UltraLong();

            unsigned int Overflow = 0;

            for (unsigned int i = 0; i < LENGTH; i++)
            {
                Res.Value[i] = this->Value[i] + ~b.Value[i] + Overflow;

                if (Res.Value[i] < this->Value[i] || (Overflow == 1 && Res.Value[i] == this->Value[i]))
                    Overflow = 1;
                else
                    Overflow = 0;
            }

            Res++;

            return Res;

        }
        
        UltraLong operator *(unsigned int b) {

            UltraLong Res = UltraLong();

            unsigned long long Overflow = 0, 
                               Overflow1 = 0;

            for (unsigned int i = 0; i < LENGTH; i++)
            {
                Overflow1 = (unsigned long long)this->Value[i] * b + Overflow;

                if (Overflow1 < Overflow) 
                    Overflow = Overflow1 + 1;
                else 
                    Overflow = Overflow1;

                Res.Value[i] = Overflow % UINT_RANGE;
                Overflow /= UINT_RANGE;
            }

            return Res;

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

            UltraLong Res = UltraLong();

            unsigned long long temp1[2 * LENGTH] = { 0 },
                temp2[2 * LENGTH] = { 0 },
                temp0[2 * LENGTH] = { 0 };

            for (unsigned int i = 0; i < LENGTH; i++)
            {
                temp1[2 * i] = this->Value[i] % SUINT_RANGE;
                temp1[2 * i + 1] = this->Value[i] / SUINT_RANGE;
                temp2[2 * i] = b.Value[i] % SUINT_RANGE;
                temp2[2 * i + 1] = b.Value[i] / SUINT_RANGE;
            }

            multiply(temp1, temp2, temp0);

            for (unsigned int i = 1; i < 2 * LENGTH; i++)
            {
                temp0[i] += temp0[i - 1] / SUINT_RANGE;
                temp0[i - 1] %= SUINT_RANGE;
            }

            temp0[2 * LENGTH - 1] %= SUINT_RANGE;

            for (unsigned int i = 0; i < LENGTH; i++)
            {
                Res.Value[i] = (unsigned int)(temp0[2 * i] + temp0[2 * i + 1] * SUINT_RANGE);
            }

            return Res;

        }
        UltraLong operator /(UltraLong b) {

            if (b == Zero) { throw "UltraLong division by zero."; }
            if (b == One) { return *this; }

            if (*this < b) { return Zero; }

            UltraLong L = 0;
            UltraLong R = *this;

            while (R - L > 1)
            {
                UltraLong M = (L + R).RightShift(1);

                if (M * b > * this) {
                    R = M;
                }
                else {
                    L = M;
                }
            }

            return L;

        }
        
        UltraLong operator %(UltraLong b) {
            return *this - *this / b * b;
        }

        UltraLong operator ^(UltraLong b) {
            UltraLong Res = UltraLong();
            
            for (unsigned int i = 0; i < LENGTH; i++)
                Res.Value[i] = this->Value[i] ^ b.Value[i];
            
            return Res;
        }
        UltraLong operator &(UltraLong b) {
            UltraLong Res = UltraLong();

            for (unsigned int i = 0; i < LENGTH; i++)
                Res.Value[i] = this->Value[i] & b.Value[i];

            return Res;
        }
        UltraLong operator |(UltraLong b) {
            UltraLong Res = UltraLong();

            for (unsigned int i = 0; i < LENGTH; i++)
                Res.Value[i] = this->Value[i] | b.Value[i];

            return Res;
        }
        UltraLong operator ~() {
            UltraLong Res = UltraLong();

            for (unsigned int i = 0; i < LENGTH; i++)
                Res.Value[i] = ~this->Value[i];

            return Res;
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

        UltraLong operator >>(UltraLong right) {
        
            unsigned long long n = right.ToULONG();

            UltraLong Res = *this;
            Res = Res.SuperRightShift((unsigned int)(n / BITS_IN_UINT));

            n %= BITS_IN_UINT;
            unsigned int mod = 1u << n;
            Res.Value[0] >>= n;

            for (unsigned int i = 1; i < LENGTH; i++)
            {
                Res.Value[i - 1] += (Res.Value[i] % mod) << (BITS_IN_UINT - n);
                Res.Value[i] >>= n;
            }

            return Res;
        
        }
        UltraLong operator <<(UltraLong right) {

            unsigned long long n = right.ToULONG();

            UltraLong Res = *this;
            Res = Res.SuperLeftShift((unsigned int)(n / BITS_IN_UINT));

            n %= BITS_IN_UINT;
            unsigned int mod = 1u << (BITS_IN_UINT - n);

            for (unsigned int i = LENGTH - 1; i > 0; i--)
            {
                Res.Value[i] <<= n;
                Res.Value[i] += Res.Value[i - 1] / mod;
            }

            return Res;

        }

        bool isEven() {
            return this->Value[0] == 0;
        }
        bool isOdd() {
            return this->Value[0] == 1;
        }

        /// <summary>
        /// Modular exponentiation.
        /// </summary>
        /// <returns>The remainder after dividing a to b-th power by mod</returns>
        static UltraLong ModPow(UltraLong a, UltraLong b, UltraLong mod) {
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
                return ModPow(a * a % mod, b.RightShift(1), mod);
            }
            else {
                return ModPow(a * a % mod, b.RightShift(1), mod) * a % mod;
            }

        }
        
        /// <summary>
        /// Cast UltraLong to unsigned int.
        /// </summary>
        /// <returns>The remainder after dividing UltraLong by 2**32.</returns>
        unsigned int ToUINT() {
            return this->Value[0];
        }
        /// <summary>
        /// Cast UltraLong to unsigned long long.
        /// </summary>
        /// <returns>The remainder after dividing UltraLong by 2**64.</returns>
        unsigned long long ToULONG() {
            return LENGTH > 1 ? this->Value[1] * UINT_RANGE + this->Value[0] : this->Value[0];
        }

        long double ToLongDouble() {
            long double ans = 0;
            for (unsigned i = LENGTH - 1; i < LENGTH; i--)
            {
                ans += this->Value[i] * pow(UINT_RANGE, i);
            }
            return ans;
        }

        /// <summary>
        /// Cast UltraLong to string.
        /// </summary>
        /// <returns>UltraLong in decimal botation.</returns>
        std::string ToString() {

            std::string s = "";

            while (*this > 0)
            {
                unsigned int digit = (*this % 10).ToUINT();
                s += std::to_string(digit);
                *this = *this / 10;
            }

            if (s == "") s = "0";

            std::string s1 = "";
            for (unsigned long long i = s.size(); i <= s.size(); i--) {
                s1 += s[i];
            }

            return s1;

        }

        static long double Divide(UltraLong a, UltraLong b) {
            static const unsigned int BASE = 1000000000;
            UltraLong div = a / b;
            long double ans = div.ToLongDouble();
            UltraLong mod = a - div * b;
            UltraLong div2 = a * BASE / b - div * BASE;
            ans += div2.ToLongDouble() / BASE;
            UltraLong div3 = a * BASE * BASE / b - div * BASE * BASE - div2 * BASE;
            ans += div3.ToLongDouble() / BASE / BASE;
            return ans;
        }

        static long double Log(UltraLong n) {
            if (n == 0) {
                throw "Logarithm of zero is minus infinity.";
            }
            long double x = ((n - 1) / (n + 1)).ToLongDouble();
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

        static UltraLong Sqrt(UltraLong n) {

            if (n < 2) return n;
        
            UltraLong L = 0, R = n;
            
            while (R - L > 1) {
                
                UltraLong M = (L + R).RightShift(1);
                
                if (M * M > n)
                    R = M;
                else 
                    L = M;
                
            }

            return L;

        }

    private:
        
        static const unsigned int LENGTH = 32;
        static const unsigned int UPPER_BOUND_LENGTH = 32;

        static const unsigned long long UINT_RANGE = 4294967296;
        static const unsigned long long SUINT_RANGE = 65536;
        static const unsigned int BITS_IN_UINT = 32;
        static const unsigned int LOG_PRECISE = 10;

        static inline long double PI = 3.1415926535897932389l;
        static inline bool precalc = false;

        //static inline unsigned int KaratsubaIterationLength = LENGTH; 
        
        unsigned int Value[LENGTH] = { 0 };
        
        static UltraLong One;
        static UltraLong Zero;
        static UltraLong MinusOne;

        static inline unsigned int rev[2 * UPPER_BOUND_LENGTH] = { 0 };
        static inline Complex wlen_pw[2 * UPPER_BOUND_LENGTH] = { Complex() };
        
        static void swap(Complex* a, Complex* b) {

            if (a == b) { return; }

            Complex c = Complex();

            c = *a;
            *a = *b;
            *b = c;

        }

        static void calc_rev(unsigned int n, unsigned int log_n) {

            for (unsigned int i = 0; i < n; ++i) {

                rev[i] = 0;

                for (unsigned int j = 0; j < log_n; ++j)
                    if (i & (1 << j))
                        rev[i] |= 1 << (log_n - 1 - j);

            }

        }

        static void multiply(const unsigned long long a[], const unsigned long long b[], unsigned long long res[]) {

            Complex fa[2 * UPPER_BOUND_LENGTH] = { Complex() }, 
                    fb[2 * UPPER_BOUND_LENGTH] = { Complex() };

            for (unsigned int i = 0; i < 2 * LENGTH; i++)
            {
                fa[i].Real = (long double)a[i];
                fb[i].Real = (long double)b[i];
            }

            unsigned int n = 2 * UPPER_BOUND_LENGTH;

            FastFourierTransformation(fa, n, false), 
            FastFourierTransformation(fb, n, false);

            for (size_t i = 0; i < n; ++i)
                fa[i] = fa[i] * fb[i];

            FastFourierTransformation(fa, n, true);


            for (size_t i = 0; i < 2 * LENGTH; ++i)
                res[i] = unsigned long long(fa[i].Real + 0.5);

        }
        
        static void FastFourierTransformation(Complex* a, unsigned int n, bool reverse = false) {
            
            if (!precalc) {
                precalc = true;
                unsigned int lg_n = 0;
                unsigned int t_lg = 1;
                while (t_lg < n) { ++lg_n; }
                calc_rev(n, lg_n);
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
        
        UltraLong SuperLeftShift(unsigned int n) {

            UltraLong Res = UltraLong();

            for (unsigned int i = LENGTH - 1 - n; i < LENGTH; i--)
                Res.Value[i + n] = this->Value[i];

            return Res;

        }
        UltraLong SuperRightShift(unsigned int n) {

            UltraLong Res = UltraLong();

            for (unsigned int i = 0; i < LENGTH - n; i++)
                Res.Value[i] = this->Value[i + n];

            for (unsigned int i = LENGTH - n; i < LENGTH; i++)
                Res.Value[i] = 0;

            return Res;

        }
        UltraLong RightShift(unsigned long long n) {

            UltraLong Res = *this;
            Res.SuperRightShift((unsigned int)(n / BITS_IN_UINT));

            n %= BITS_IN_UINT;
            unsigned int mod = 1u << n;
            Res.Value[0] >>= n;

            for (unsigned int i = 1; i < LENGTH; i++)
            {
                Res.Value[i - 1] += (Res.Value[i] % mod) << (BITS_IN_UINT - n);
                Res.Value[i] >>= n;
            }

            return Res;

        }
};

UltraLong UltraLong::One = UltraLong(1);
UltraLong UltraLong::Zero = UltraLong(0);
UltraLong UltraLong::MinusOne = UltraLong(-1);

enum PRIME_TESTS_OPTIMISE_LEVELS {O1, O2, O3};
const inline PRIME_TESTS_OPTIMISE_LEVELS PRIME_TESTS_OPTIMISE_LEVEL = O3;

template <class T>
T Sqrt(T n) {
    return (T)sqrtl((long double)n);
}

template <>
UltraLong Sqrt(UltraLong n) {
    return UltraLong::Sqrt(n);
}

template <class T>
T MillerUpperBound_O1(T n) {
    return 2 * powl(logl(n), 2);
}

template <>
UltraLong MillerUpperBound_O1<UltraLong>(UltraLong n) {
    return (unsigned long long)floorl(2 * powl(UltraLong::Log(n), 2));
}

template <class T>
T MillerUpperBound_O2(T n) {
    long double ans = logl(n);
    return  ans * logl(ans) / logl(2);
}

template <>
UltraLong MillerUpperBound_O2<UltraLong>(UltraLong n) {
    long double ans = UltraLong::Log(n);
    return (unsigned long long)floorl(ans * logl(ans) / logl(2));
}

template <class T>
bool IsPrimePower(T n) {
    return false;
}

template <class T>
T ModPow(T a, T b, T mod){
    if (mod == 0) { throw "ModPow: mod is zero"; }
    if (a == 0 || mod == 1) { return 0; }
    if (b == 0) { return 1; }
    if (b % 2 == 1) { return ModPow<T>(a * a % mod, b / 2, mod) * a % mod; }
    return ModPow<T>(a * a % mod, b / 2, mod);
}

template <>
UltraLong ModPow<UltraLong>(UltraLong a, UltraLong b, UltraLong mod) {
    return UltraLong::ModPow(a, b, mod);
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
            res = UltraLong.ModPow(base, exp, n);
            if (res == 1)
                return true;

            break;
        }
    }

    return false;

}

template <class T>
T SmallNextPrime(T base) {

    T x = base + 1;

    while (!SmallIsPrime(x))
    {
        x++;
    }

    return x;

}

template <class T>
bool MillerTest(T n, PRIME_TESTS_OPTIMISE_LEVELS optimise_level = O2) {

    if (IsPrimePower<T>(n)) {
        return false;
    }

    T upper_bound;
    
    switch (optimise_level)
    {
    case O1: {
        upper_bound = MillerUpperBound_O1(n); 
        break;
    }
    default:
        upper_bound = MillerUpperBound_O2(n);
        break;
    }

    T base = 2;

    while (base <= n) {
    
        if (!IsStrongPseudoprime<T>(n, base)) {
            return false;
        }

        base = SmallNextPrime<T>(base);

    }

    return true;

}

template <class T>
bool SmallIsPrime(T base) {
    
    if (base < 2) return false;

    T sq_b = Sqrt<T>(base) + 1;
    
    for (T i = 2; i < sq_b; i++)
        if (base % i == 0)
            return false;

    return true;

}

template <class T>
T NextPrime(PRIME_TESTS_OPTIMISE_LEVELS optimize_level = O3) {
    if (*this == (T)2) { return (T)3; }
    if (*this < (T)2) { return (T)2; }
    switch (optimize_level)
    {
    case O1: {
        
    }
    default:
        break;
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
        now = NextPrime<T>(now);
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