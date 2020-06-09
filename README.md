# extra-math
Some math that is not in standard library, for example, long arithmetic, a bit faster complex numbers and random numbers generator.

# Documentation
## UltraLong
### Description
`UltraLong` is a big unsigned number. Before use you should change settings (else it'll be your fault). `UltraLong` assumes that size of `unsigned int` is 32 bits and size of `unsigned long long` is 64 bits. If it's not true then `UltraLong` may make errors in calculations. If so, you can always open an issue with size of `unsigned int` and `unsigned long long` sizes.
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
|`ModPow(a, b, mod)`  | ✓ |Reminder after dividing `a` to `b`-th power by `mod`|
|`ToUINT()`| ☓ |Reminder after dividing `UltraLong` by 2<sup>32</sup>|
|`ToULONG()`| ☓ |Reminder after dividing `UltraLong` by 2<sup>64</sup>|
|`ToString()`| ☓ |Returns `std::string` with `UltraLong` in decimal notation|
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
Actually, the `Complex` type is very restricted now because it is used only for `UltraLong` realisation. Some of basic operators are missing.

## Random
### Description
Random numbers generator. Create an object with `Random rnd = Random()` and use the functions below.
### Functions
|Function|Return type|
|:-:     |:-  |
|`next_int()`| `int` |
|`next_long()`| `long long` |
|`next_uint()`| `unsigned int` |
|`next_ulong()`| `unsigned long long` |

# Donate
If you decided to donate, think again. But if you haven't changed your mind yet, then I say thank you for supporting for a developing project! You can contact me on cat8ra@gmail.com