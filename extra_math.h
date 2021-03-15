/// @section License
/// Copyright (c) 2020-2021 Launer Maksim
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in all
/// copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
/// SOFTWARE.

#pragma once

#include <cmath>
#include <string>
#include <ctime>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>
#include <unordered_set>
#include <unordered_map>

struct Random {
public:
    Random(long long seed);
    bool nextBool();
    char nextChar();
    short nextSInt();
    unsigned short nextUSInt();
    int nextInt();
    unsigned int nextUInt();
    long long nextLong();
    unsigned long long nextULong();
    double nextDouble();

private:
    void step();
    static const long long a[3];
    static const long long b[3];
    long long now;
};

/// <summary>
/// An implementation of complex number (using long double).
/// </summary>
struct Complex {
    Complex();
    Complex(long double n);
    Complex(long double a, long double b);

    bool operator ==(Complex right);
    bool operator !=(Complex right);

    Complex operator =(Complex right);

    Complex operator +();
    Complex operator -();

    Complex operator +(Complex b);
    Complex operator -(Complex b);

    Complex operator *(Complex b);
    Complex operator /(Complex b);

    Complex operator +=(Complex b);
    Complex operator -=(Complex b);

    Complex operator *=(unsigned int b);
    Complex operator /=(unsigned int b);

    static long double abs(Complex z);
    std::string toString();

    static Complex conjugated(Complex a);

    long double real = 0;
    long double imaginary = 0;
};

//Unsigned ultra long realization
/// <summary>
/// An implementation of unsigned long arithmetic.
/// </summary>
class UnsignedUltraLong {

public:
    /// <summary>
    /// UltraLong constructor.
    /// </summary>
    /// <returns>
    /// Zero.
    /// </returns>
    UnsignedUltraLong();
    /// <summary>
    /// UltraLong constructor.
    /// </summary>
    /// <param name = "n">
    /// : The value that initialises UltraLong.
    /// </param>
    UnsignedUltraLong(unsigned int n);
    /*UltraLong(SignedUltraLong x) {
        for (unsigned int i = 0; i < LENGTH; i++)
            this->value[i] = x.value[i];
    }*/
    /*TODO*/
    /// <summary>
    /// UltraLong constructor.
    /// </summary>
    /// <param name = "n">
    /// : The value that initialises UltraLong. If it is less than zero then
    /// UltraLong will be equal ~n + 1, so n + (-n) = 0.
    /// </param>
    UnsignedUltraLong(int n);
    UnsignedUltraLong(long long n);
    /*TODO*/
    /// <summary>
    /// UltraLong constructor.
    /// </summary>
    /// <param name = "n">
    /// : The value that initialises UltraLong.
    /// </param>
    UnsignedUltraLong(unsigned long long n);

    bool operator ==(UnsignedUltraLong right);
    bool operator !=(UnsignedUltraLong right);
    bool operator <(UnsignedUltraLong right);
    bool operator >(UnsignedUltraLong right);
    bool operator >=(UnsignedUltraLong right);
    bool operator <=(UnsignedUltraLong right);

    UnsignedUltraLong operator =(UnsignedUltraLong right);

    UnsignedUltraLong operator -();
    UnsignedUltraLong operator +();
    UnsignedUltraLong operator ++(int);

    UnsignedUltraLong operator +(UnsignedUltraLong b);
    UnsignedUltraLong operator -(UnsignedUltraLong b);

    UnsignedUltraLong operator *(unsigned int b);
    UnsignedUltraLong operator /(unsigned int b);
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
    UnsignedUltraLong operator *(UnsignedUltraLong b);
    UnsignedUltraLong operator /(UnsignedUltraLong b);

    unsigned int operator %(unsigned int b);
    UnsignedUltraLong operator %(UnsignedUltraLong b);

    UnsignedUltraLong operator ^(UnsignedUltraLong b);
    UnsignedUltraLong operator &(UnsignedUltraLong b);
    UnsignedUltraLong operator |(UnsignedUltraLong b);
    UnsignedUltraLong operator ~();

    UnsignedUltraLong operator +=(UnsignedUltraLong right);
    UnsignedUltraLong operator -=(UnsignedUltraLong right);
    UnsignedUltraLong operator *=(UnsignedUltraLong right);
    UnsignedUltraLong operator *=(unsigned int right);
    UnsignedUltraLong operator /=(UnsignedUltraLong right);

    UnsignedUltraLong operator %=(UnsignedUltraLong right);

    UnsignedUltraLong operator ^=(UnsignedUltraLong right);
    UnsignedUltraLong operator &=(UnsignedUltraLong right);
    UnsignedUltraLong operator |=(UnsignedUltraLong right);

    UnsignedUltraLong operator >>(unsigned long long right);
    UnsignedUltraLong operator <<(unsigned long long right);

    bool isEven();
    bool isOdd();

    /// <summary>
    /// Modular exponentiation.
    /// </summary>
    /// <returns>The remainder after dividing a to b-th power by mod</returns>
    static UnsignedUltraLong modPow(UnsignedUltraLong a, UnsignedUltraLong b, UnsignedUltraLong mod);

    /// <summary>
    /// Cast UltraLong to unsigned int.
    /// </summary>
    /// <returns>The remainder after dividing UltraLong by 2**32.</returns>
    unsigned int toUint();
    /// <summary>
    /// Cast UltraLong to unsigned long long.
    /// </summary>
    /// <returns>The remainder after dividing UltraLong by 2**64.</returns>
    unsigned long long toUlong();

    long double toLongDouble();

    /// <summary>
    /// Cast UltraLong to string.
    /// </summary>
    /// <returns>UltraLong in decimal botation.</returns>
    std::string toString();

    static long double divide(UnsignedUltraLong a, UnsignedUltraLong b);

    static long double log(UnsignedUltraLong n);

    static UnsignedUltraLong sqrt(UnsignedUltraLong n);

    static UnsignedUltraLong parse(std::string s);

    friend class UltraLong;

protected:

    static UnsignedUltraLong middle(UnsignedUltraLong a, UnsignedUltraLong b);

    static std::pair<unsigned int, unsigned int> lead(UnsignedUltraLong a);

    static void swap(Complex* a, Complex* b);

    static void calcRev(unsigned int n, unsigned int log_n);

    static void multiply(const unsigned long long a[], const unsigned long long b[], unsigned long long res[]);

    static void fastFourierTransformation(Complex* a, unsigned int n, bool reverse);

    UnsignedUltraLong superLeftShift(unsigned int n);
    UnsignedUltraLong superRightShift(unsigned int n);
    UnsignedUltraLong rightShift(unsigned long long n);
    static const unsigned int LENGTH = 2;
    static const unsigned int UPPER_BOUND_LENGTH = 2;

    static const unsigned long long UINT_RANGE = 4294967296;
    static const unsigned long long SUINT_RANGE = 65536;
    static const unsigned int BITS_IN_UINT = 32;
    static const unsigned int LOG_PRECISE = 10;

    static long double PI;
    static bool precalc;
    static UnsignedUltraLong One;
    static UnsignedUltraLong Zero;
    static UnsignedUltraLong MinusOne;
    static unsigned int rev[4 * UPPER_BOUND_LENGTH];
    static Complex wlen_pw[4 * UPPER_BOUND_LENGTH];
    static bool lastMultOverflow;
    unsigned int value[LENGTH];
};

class UltraLong : public UnsignedUltraLong {
public:
    operator UnsignedUltraLong();
    UltraLong();
    UltraLong(unsigned int n);
    UltraLong(UnsignedUltraLong x);
    UltraLong(int n);
    UltraLong(long long n);
    UltraLong(unsigned long long n);
    bool operator <(UltraLong right);
    bool operator >(UltraLong right);
    bool operator >=(UltraLong right);
    bool operator <=(UltraLong right);
    UltraLong operator /(unsigned int b);
    static UltraLong parse(std::string s);

    std::string toString();
    static UnsignedUltraLong abs(UltraLong x);
protected:
    bool isNonnegative();
    bool isNegative();
};

namespace Math {

    enum class PRIME_TESTS_OPTIMISE_LEVELS { O1, O2, O3 };
    const PRIME_TESTS_OPTIMISE_LEVELS PRIME_TESTS_OPTIMISE_LEVE = Math::PRIME_TESTS_OPTIMISE_LEVELS::O3;

    template <class T>
    T sqrt(T n);

    template <>
    UnsignedUltraLong sqrt(UnsignedUltraLong n);

    template <class T>
    T MillerUpperBound_O1(T n);

    template <>
    UnsignedUltraLong MillerUpperBound_O1<UnsignedUltraLong>(UnsignedUltraLong n);

    template <class T>
    T MillerUpperBound_O2(T n);

    template <>
    UnsignedUltraLong MillerUpperBound_O2<UnsignedUltraLong>(UnsignedUltraLong n);

    template <class T>
    bool IsPrimePower(T n);

    template <class T>
    T ModPow(T a, T b, T mod);

    template <>
    UnsignedUltraLong ModPow<UnsignedUltraLong>(UnsignedUltraLong a, UnsignedUltraLong b, UnsignedUltraLong mod);

    template <class T>
    bool IsStrongPseudoprime(T n, T base);

    template <class T>
    bool SmallIsPrime(T base);

    template <class T>
    T SmallNextPrime(T base);

    template <class T>
    bool MillerTest(T n, PRIME_TESTS_OPTIMISE_LEVELS optimise_level = PRIME_TESTS_OPTIMISE_LEVEL);

    /*template <class T>
    T PrimeTest() {}*/

    //TODO
    template <class T>
    T nextPrime(T n, PRIME_TESTS_OPTIMISE_LEVELS optimize_level = PRIME_TESTS_OPTIMISE_LEVEL);

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
        Arithmetic<T>(T start, T difference);
        T next();
        Arithmetic<T> renew();
    private:
        T start;
        T now;
        T difference;
        unsigned long long iter;
};

template <class T>
class Geometric : public Sequence<T> {
    public:
        Geometric(T start, T denominator);
        T next();
        Geometric<T> renew();
    private:
        T start;
        T now;
        T denominator;
        unsigned long long iter;
};

template <class T>
class Fibonacci : public Sequence<T> {
    public:
        Fibonacci<T>(T start = 0, T start1 = 1);
        T next();
        Fibonacci<T> renew();
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
    Primes();
    T next();
    Primes<T> renew();
private:
    T start;
    T prev;
    T now;
    unsigned long long iter;
};
//Arithmetic<int> nat = Arithmetic<int>(1, 1);

struct Permutation {
    //Gives identity permutation with length n
    Permutation(unsigned int length);
    /*Permutation(int t[]) {
        int n = sizeof(t) / sizeof(int);
        for (int i = 0; i < n; i++){
            x[i] = t[i];
        }
    }*/

    Permutation operator *(Permutation b);

    Permutation inverse();
private:
    unsigned int len;
    std::vector<unsigned int> x;
    Permutation(unsigned int length, unsigned int value);
};

template <class T>
struct Matrix {

    Matrix<T>(const std::vector<std::vector<T>> &arr);

    Matrix<T>(unsigned int r, unsigned int c, T &value);

    Matrix<T>(unsigned int r, unsigned int c);

    Matrix<T> operator =(Matrix<T> right);

    Matrix<T> operator +(Matrix<T> b);

    Matrix<T> operator -(Matrix<T> b);

    Matrix<T> operator *(Matrix<T> b);

    Matrix<T> operator *(T b);

    static Matrix<T> randomMatrix(unsigned int columns, unsigned int rows, T generator(void));

    static Matrix<T> transpose(Matrix<T> a);

    /*static T determinant(Matrix<T> a) {
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
    }*/

    void print();

    //Matrix<T> complementaryMinor() {}

private:
    unsigned int rows, columns;
    std::vector<std::vector<T>> mx;

};

struct Polynomial {

public:

    Polynomial();

    Polynomial(std::vector<long double> coefficients);

    Polynomial operator +(Polynomial b);

    Polynomial operator -(Polynomial b);

    long double Value(long double x);

    static Polynomial Derivative(Polynomial a);

    struct Solutions
    {
    public:
        Solutions();
        bool IsAnyNumber();
        void Add(long double sol);
        void SetAnyNumber(bool an);
        std::vector<long double> GetSolutions();
    protected:
        std::vector<long double> ans;
        bool any_number; 
    };
    
    static Solutions Solve(Polynomial a);

    /*void Print(std::string var) {
        bool first = true;
        std::string s = "";
        for (auto iter = this->cfs.rbegin(); iter != this->cfs.rend(); ++iter) {
            if (!first && iter->second > 0)
                s += '+';
            s += std::to_string(iter->second) + var + '^' + std::to_string(iter->first);
            first = false;
        }
        std::cout << s << '\n';
    }*/

private:

    std::unordered_map<unsigned int, long double> cfs;
    unsigned int degree;

    void ClearZeroes();

    static Solutions SolveLinear(Polynomial a);

    static Solutions SolveQuadratic(Polynomial a);

    long double FindZero(long double L, long double R);
    long double FindCriticalZero(long double L, long double R);

    bool IsNearRoot(long double x);

    static long double SOLVE_PRECISE;
    static long double MAX_ROOT;
};

class Graph {
public:
    Graph();
    Graph(unsigned int vertex);
    inline void addEdge(unsigned int first, unsigned int second);
    std::unordered_set<unsigned int> neighbours(unsigned int index);
    Graph(unsigned int vertex, std::vector<std::pair<unsigned int, unsigned int>> &edges);
    bool isTree();
protected:
    unsigned int vertex;
    std::vector<std::unordered_set<unsigned int>> edges;
    bool isTreeDFS(unsigned int v, unsigned int p);
    std::vector<unsigned int> color;
};

struct Vector3D {
public:
    Vector3D();
    Vector3D(long double x, long double y, long double z);
    long double abs();
    Vector3D operator +(Vector3D b);
    Vector3D operator -(Vector3D b);
    long double dot(Vector3D b);
    Vector3D operator * (long double b);
    Vector3D operator * (Vector3D b);
private:
    long double x, y, z;
};

struct Parser {
public:
    Parser();
    long double parse(std::string s, unsigned int shift = 0);
    std::string getError();
private:
    static int code;
    static unsigned int error_index;
    static std::string error;
    static void generate_error_message();
};

/*class MultiGraph : public Graph {
protected:
    unsigned int vertex;
    std::vector<std::unordered_multiset<unsigned int>> edges;
};*/
