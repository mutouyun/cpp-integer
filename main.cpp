#include "include/integer.hpp"

#include <iostream>

#define TEST_CASE_()                                           \
    std::cout << std::endl << __func__ << " ->:" << std::endl

void test_shift(void)
{
    TEST_CASE_();

    integer i = (uint32_t)-1;
    std::cout << i << std::endl;
    i <<= 66;
    std::cout << i << std::endl;
    i >>= 33;
    std::cout << i << std::endl;
}

void test_sub(void)
{
    TEST_CASE_();

    integer i = 123;
    std::cout << (i -= 321) << std::endl;
}

void test_mul(void)
{
    TEST_CASE_();

    integer i = 123;
    std::cout << (i *= 98) << std::endl;
    std::cout << (i *= (uint32_t)-1) << std::endl;
}

void test_div(void)
{
    TEST_CASE_();

    integer i = 123, j = 3;
    i *= (uint64_t)-1; // 123 * 18446744073709551615 = 2268949521066274848645
    j *= (uint32_t)-1; //   3 *           4294967295 =            12884901885
    std::cout << i << ", " << j << std::endl;
    std::cout << (i /= j) << std::endl;         // 176093659177
    std::cout << std::hex << i << std::endl;    // 2900000029
    std::cout << std::dec;
}

#include <cmath>

int main(void)
{
/*
    test_shift();
    test_sub();
    test_mul();
    test_div();
    std::cout << std::endl;

    std::cout << "Fibonacci Sequence:" << std::endl;
    integer a, b = 1;
    for (int i = 0; i < 1000; ++i)
    {
        std::cout << a << std::endl;
        integer c = a + b;
        a = std::move(b);
        b = std::move(c);
    }
*/
    auto is_prime = [](const auto& n)
    {
        using n_t = decltype(n);
        if (n == n_t(2)) return true;
        if (n & n_t(1))
        {
            auto s = static_cast<size_t>(std::sqrt(static_cast<size_t>(n)));
            for (size_t i = 3; i <= s; ++i)
            {
                if (i & 1)
                {
                    if (n % n_t(i)) continue;
                    return false;
                }
            }
            return true;
        }
        return false;
    };

    for (integer k = 2, c = 1; k < integer(100); ++k)
    {
        if (is_prime(k))
        {
            std::cout << c << " " << k << std::endl;
            c += 1;
        }
    }

    return 0;
}
