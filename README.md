# extra-math
Some math that is not in standard library, for example, long arithmetic, a bit faster complex numbers and random numbers generator.

# Documentation
## UltraLong
### Description
`UltraLong` is an implementation of long arithmetic. There is `UnsignedUltraLong`, too. Before use you should change settings (else it'll be your fault). `UltraLong` assumes that size of `unsigned int` is 32 bits and size of `unsigned long long` is 64 bits. If it's not true then `UltraLong` may make errors in calculations. If so, you can always open an issue with size of `unsigned int` and `unsigned long long`.
### Settings
There are 5 private constants in `UltraLong` that are called "Settings":
|Name|Description|Assumed to be|
|:-: |:-         |:-:|
|`LENGTH`|Any `UltraLong` will be `LENGTH * BITS_IN_UINT` bits long|0 < `LENGTH` < 2<sup>30</sup> + 1|
|`UPPER_BOUND_LENGTH`|The least power of 2 that is equal or greater than `LENGTH`|1 < `UPPER_BOUND_LENGTH` < 2<sup>30</sup> + 1
|`UINT_RANGE`|`2^BITS_IN_UINT`, range of `unsigned int`| 2<sup>32</sup> |
|`SUINT_RANGE`|`2^(BITS_IN_UINT/2)`, range of `unsigned short int`| 2<sup>16</sup> |
|`BITS_IN_UINT`|Size of `unsigned int` in bits | 32|
### Functions
|Function|Is static|Description|
|:-:     |:-:      |:-         |
|`modPow(a, b, mod)`  | ✓ |Reminder after dividing `a` to `b`-th power by `mod`|
|`toUINT()`| ☓ |Reminder after dividing `UltraLong` by 2<sup>32</sup>|
|`toULONG()`| ☓ |Reminder after dividing `UltraLong` by 2<sup>64</sup>|
|`toString()`| ☓ |Returns `std::string` with `UltraLong` in decimal notation|
### Computational complexity
Let `n = LENGTH`
|Operation|Complexity|
|:-:      |:-:       |
|`+`|*O(n)*|
|`-`|*O(n)*|
|`*`|*O(n×*log *n×*log log *n)*|
|`/`|*O(n<sup>2</sup>×*log *n×*log log *n)*|
|`%`|*O(n<sup>2</sup>×*log *n×*log log *n)*|
|`ModPow`|*O(n<sup>3</sup>×*log *n×*log log *n)*|
|`ToUINT`|*O(1)*|
|`ToULONG`|*O(1)*|
|`ToString`|*O(n<sup>3</sup>×*log *n×*log log *n)*|

## Complex
### Description
The type that implements a complex number based on `long double`. Why not `complex <long double>`? Some compilers do not optimize the templates well.

## Random
### Description
Random numbers generator. Create an object with `Random rnd = Random()` or `Random rnd = Random(seed)` (`seed` is `long long` type) and use the functions below.
### Functions
|Function      |Return type           |
|:-:           |:-                    |
|`nextBool()`  | `bool`               |
|`nextChar()`  | `char`               |
|`nextSInt()`  | `short int`          |
|`nextUSInt()` | `unsigned short int` |
|`nextInt()`   | `int`                |
|`nextLong()`  | `long long`          |
|`nextUInt()`  | `unsigned int`       |
|`nextULong()` | `unsigned long long` |
|`nextDouble()`| `double`             |

DIEHARD test result (using `nextInt()`) is in the eponymous file.
## Sequence
### Description
An infinite sequence of numbers. Template class. Actually, it is some kind of abstract class or interface with several implementations that you can use out of the box.
### Functions
Let `T` be a type of number.
|Function |Return type  |Description|
|:-:      |:-:          |:-         |
|`next()` |`T`          |Returns the next number of the sequence.|
|`renew()`|`Sequence<T>`|`seq.renew()` returns `Sequence` like `seq`, but it starts from the first number.|
### Implementations
#### Arithmetic
An arithmetic progression. To create: `Arithmetic<T>(T a, T d)`, you will get the sequence: *a*, *a + d*, *a + 2d*, ...
#### Geometric
A geometric progression. To create: `Geometric<T>(T a, T d)`, you will get the sequence: *a*, *ad*, *ad<sup>2</sup>*, ...
#### Fibonacci
A Fibonacci sequence (each number is a sum of the two previous ones). To create: `Fibonacci<T>(T a, T b)`, you will get the sequence: *a*, *b*, *a + b*, *a + 2b*, *2a + 3b*, ...
By default, `a = 0` and `b = 1`.
#### Primes
A prime number sequence. **NOT TESTED.** To create: `Primes()`, you will get the sequence: *2*, *3*, *5*, *7*, ...
## Math
### Description
A namespace that contains some useful functions.
### Functions
|Function      |Description|Return type|
|:-:           |:-         |:-         |
|`sqrt<T>(n)`  |Returns the square root of a number (type `T`).|`T`|
|`ModPow<T>(n)`|Returns the square root of a number (type `T`).|`T`|

# Donate
If you decided to donate, think again. But if you haven't changed your mind yet, then I say thank you for supporting for a developing project! You can contact me on cat8ra@gmail.com