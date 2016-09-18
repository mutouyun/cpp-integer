/*
    cpp-integer - Code covered by the MIT License
    Author: mutouyun (http://darkc.at)
*/

#pragma once

#include <utility>   // std::swap, std::move
#include <deque>     // std::deque
#include <vector>    // std::vector
#include <tuple>     // std::tuple
#include <limits>    // std::numeric_limits
#include <algorithm> // std::min
#include <iostream>  // std::ostream, std::istream
#include <iomanip>   // std::setw
#include <cmath>     // std::power
#include <cassert>   // assert
#include <cstdint>   // uint32_t, ...

class integer
{
    using arr_t = std::deque<uint32_t>;

    enum { VALUE_DIGITS = std::numeric_limits<arr_t::value_type>::digits };

    bool  minus_ = false;
    arr_t values_;

public:
    integer(void)           = default;
    integer(const integer&) = default;
    integer(integer&&)      = default;
    integer(uint32_t n)
        : integer()
    {
        if (n == 0) return;
        values_.push_back(n);
    }
    integer(int32_t n)
        : integer(static_cast<uint32_t>(n < 0 ? (minus_ = true, -n) : n))
    {}
    integer(uint64_t n)
        : integer()
    {
        if (n == 0) return;
        values_.push_back(static_cast<arr_t::value_type>(n & 0xffffffffull));
        arr_t::value_type carry = static_cast<arr_t::value_type>(n >> 32);
        if (carry != 0) values_.push_back(carry);
    }
    integer(int64_t n)
        : integer(static_cast<uint64_t>(n < 0 ? (minus_ = true, -n) : n))
    {}

    operator uint32_t(void) const
    {
        if (values_.empty()) return 0;
        return static_cast<uint32_t>(values_[0]);
    }

    operator uint64_t(void) const
    {
        if (values_.empty()) return 0;
        if (values_.size() == 1) return static_cast<uint64_t>(values_[0]);
        uint64_t ret = values_[1];
        return ((ret <<= 32) |= values_[0]);
    }

    size_t size(void) const
    {
        return this->values_.size();
    }

    void swap(integer& y)
    {
        std::swap(this->minus_, y.minus_);
        this->values_.swap(y.values_);
    }

    integer abs(void) const
    {
        integer r = (*this);
        r.minus_ = false;
        return std::move(r);
    }

private:
    void clean_up(void)
    {
        while (this->size() > 0 && this->values_.back() == 0)
            this->values_.pop_back();
    }

    static int unsign_compare(const integer& x, const integer& y)
    {
        int r = 0;
        if (x.size() == y.size())
        {
            if (x.size() == 0) return 0;
            for (size_t i = x.size() - 1; i > 0; --i)
            {
                if (x.values_[i] != y.values_[i])
                {
                    r = (x.values_[i] < y.values_[i]) ? -1 : 1;
                    break;
                }
            }
            if (r == 0)
            {
                if (x.values_.front() != y.values_.front())
                    r = (x.values_.front() < y.values_.front()) ? -1 : 1;
            }
        }
        else r = (x.size() < y.size()) ? -1 : 1;
        return r;
    }

    integer& unsign_add(const integer& y)
    {
        if (this->size() < y.size()) this->values_.resize(y.size());
        arr_t::value_type carry = 0;
        for (size_t i = 0; i < this->size(); ++i)
        {
            uint64_t sum = this->values_[i];
            if (i < y.size()) sum += y.values_[i];
            else if (carry == 0) break;
            sum += carry;
            this->values_[i] = static_cast<arr_t::value_type>(sum & 0xffffffffull);
            carry            = static_cast<arr_t::value_type>(sum >> 32);
        }
        if (carry != 0) this->values_.push_back(carry);
        return (*this);
    }

    integer& unsign_sub(const integer& y)
    {
        arr_t::value_type borrow = 0;
        for (size_t i = 0; i < this->size(); ++i)
        {
            uint64_t sub = (i < y.size()) ? y.values_[i] : 0;
            sub += borrow;
            if (this->values_[i] < sub)
            {
                uint64_t num = 0x100000000ull + this->values_[i];
                this->values_[i] = static_cast<arr_t::value_type>(num - sub);
                borrow = 1;
            }
            else
            {
                this->values_[i] -= static_cast<arr_t::value_type>(sub);
                borrow = 0;
            }
        }
        clean_up();
        return (*this);
    }

    integer& unsign_mul(const integer& y)
    {
        integer ret;
        for (size_t j = 0; j < y.size(); ++j)
        {
            integer tmp;
            for (size_t i = 0; i < this->size(); ++i)
            {
                uint64_t mul = this->values_[i]; mul *= y.values_[j];
                tmp.unsign_add(integer{mul} <<= (VALUE_DIGITS * i));
            }
            ret.unsign_add(tmp <<= (VALUE_DIGITS * j));
        }
        return ((*this) = std::move(ret));
    }

    static std::tuple<uint64_t, integer> special_div(integer&& dd, const integer& ds)
    {
        static const auto u64_cast = [](const integer& x, size_t max)->uint64_t
        {
            size_t diff = max - x.size();
            if (diff >= 2) return 0;
            uint64_t ret = 0;
            if (!x.values_.empty())
            {
                ret = x.values_.back();
                if (diff == 0 && x.size() > 1)
                {
                    ((ret <<= VALUE_DIGITS) |= x.values_[x.values_.size() - 2]);
                }
            }
            return ret;
        };
        int cmp = unsign_compare(dd, ds);
        if (cmp == 0)
        {
            return std::tuple<uint64_t, integer>(1, 0);
        }
        else if (cmp < 0)
        {
            return std::tuple<uint64_t, integer>(0, std::move(dd));
        }
        else // dd.size() == ds.size()
        {
            size_t max = std::max(dd.size(), ds.size());
            uint64_t dd_64 = u64_cast(dd, max);
            uint64_t ds_64 = u64_cast(ds, max);
            assert(dd_64 && ds_64);
            uint64_t qu_64 = dd_64 / ds_64;
            integer re;
            integer ds_mul = ds;
            ds_mul.unsign_mul(qu_64);
            int cmp = 0, last_cmp = 0;
            const auto calc_rem = [&re](auto& dd, auto& ds_mul)
            {
                re = std::move(dd);
                re.unsign_sub(ds_mul);
            };
            while (true)
            {
                cmp = unsign_compare(dd, ds_mul);
                if (cmp < 0)
                {
                    --qu_64;
                    ds_mul.unsign_sub(ds);
                    if (last_cmp > 0)
                    {
                        calc_rem(dd, ds_mul);
                        break;
                    }
                }
                else if (cmp > 0)
                {
                    if (last_cmp < 0)
                    {
                        calc_rem(dd, ds_mul);
                        break;
                    }
                    ++qu_64;
                    ds_mul.unsign_add(ds);
                }
                else break;
                last_cmp = cmp;
            }
            return std::tuple<uint64_t, integer>(qu_64, std::move(re));
        }
    }

    integer& unsign_div(const integer& y, integer* rem)
    {
        static const auto cut = [](const integer& x, size_t n) -> std::tuple<integer, size_t>
        {
            if (n >= x.size()) return std::make_tuple(x, x.size());
            integer ret;
            size_t start_i = x.size() - 1, end = std::min(n, x.size()), k = 0;
            for (size_t i = start_i; k < end; --i, ++k)
                ret.values_.push_front(x.values_[i]);
            return std::make_tuple(std::move(ret), k);
        };
        auto cut_r = cut(*this, y.size());
        integer qu, re = std::move(std::get<0>(cut_r));
        size_t start = std::get<1>(cut_r);
        while (true)
        {
            auto div_r = special_div(std::move(re), y);
            qu.values_.push_front(static_cast<uint32_t>(std::get<0>(div_r)));
            re = std::move(std::get<1>(div_r));
            if (++start > this->size()) break;
            re.values_.push_front(this->values_[this->size() - start]); // carry
            re.clean_up();
        }
        this->swap(qu);
        clean_up();
        if (rem != nullptr) rem->swap(re);
        return (*this);
    }

public:
    static int compare(const integer& x, const integer& y)
    {
        if (x.minus_ != y.minus_) return x.minus_ ? -1 : 1;
        int r = unsign_compare(x, y);
        if (x.minus_) r = -r;
        return r;
    }

    integer& operator=(integer y)
    {
        this->swap(y);
        return (*this);
    }

    friend bool operator==(const integer& x, const integer& y)
    {
        return (integer::compare(x, y) == 0);
    }

    friend bool operator<(const integer& x, const integer& y)
    {
        return (integer::compare(x, y) < 0);
    }

    integer& operator>>=(size_t n)
    {
        if ( (n == 0) || (this->size() == 0) ) return (*this);
        if (n >= (this->size() * VALUE_DIGITS)) return ((*this) = 0);
        if (n == VALUE_DIGITS)
        {
            this->values_.pop_front();
        }
        else if (n < VALUE_DIGITS)
        {
            for (size_t i = 0; i < this->size() - 1; ++i)
            {
                (this->values_[i] >>= n) |= (this->values_[i + 1] << (VALUE_DIGITS - n));
            }
            this->values_.back() >>= n;
            if (this->values_.back() == 0) this->values_.pop_back();
        }
        else // n > VALUE_DIGITS
        {
            do this->values_.pop_front();
            while ((n -= VALUE_DIGITS) > VALUE_DIGITS);
            this->operator>>=(n);
        }
        return (*this);
    }

    integer& operator<<=(size_t n)
    {
        if ( (n == 0) || (this->size() == 0) ) return (*this);
        if (n == VALUE_DIGITS)
        {
            this->values_.push_front(0);
        }
        else if (n < VALUE_DIGITS)
        {
            arr_t::value_type carry = this->values_.back() >> (VALUE_DIGITS - n);
            for (size_t i = this->size() - 1; i > 0; --i)
            {
                (this->values_[i] <<= n) |= (this->values_[i - 1] >> (VALUE_DIGITS - n));
            }
            this->values_.front() <<= n;
            if (carry != 0) this->values_.push_back(carry);
        }
        else // n > VALUE_DIGITS
        {
            do this->values_.push_front(0);
            while ((n -= VALUE_DIGITS) > VALUE_DIGITS);
            this->operator<<=(n);
        }
        return (*this);
    }

    integer& operator+=(const integer& y)
    {
        if (this->minus_ == y.minus_)
            return unsign_add(y);
        else
        {
            int c = unsign_compare(*this, y);
            if (c == 0) return (*this) = 0;
            if (c  > 0) return unsign_sub(y);
            // (c  < 0)
            integer t = y;
            t.unsign_sub(*this);
            return (*this) = std::move(t);
        }
    }

    integer& operator-=(const integer& y)
    {
        this->minus_ = !this->minus_;
        (*this) += (y);
        this->minus_ = !this->minus_;
        return (*this);
    }

    integer& operator*=(const integer& y)
    {
        unsign_mul(y);
        this->minus_ = (this->minus_ != y.minus_);
        return (*this);
    }

    integer& operator/=(const integer& y)
    {
        unsign_div(y, nullptr);
        this->minus_ = (this->minus_ != y.minus_);
        return (*this);
    }

    integer& operator%=(const integer& y)
    {
        unsign_div(y, this);
        return (*this);
    }

    friend std::ostream& operator<<(std::ostream& s, const integer& me)
    {
        if (me.minus_) s << "-";
        if (me.size() == 0)
        {
            s << 0;
            return s;
        }
        auto ff = s.flags();
        if (ff & std::ios_base::hex)
        {
            auto old_w = s.width();
            auto old_f = s.fill();
            for (auto it = me.values_.rbegin(); it != me.values_.rend(); ++it)
            {
                s << (*it);
                s.width(8);
                s.fill('0');
            }
            s.width(old_w);
            s.fill(old_f);
        }
        else
        {
            static const auto gain_unit = [](uint64_t x)
            {
                uint64_t a = x / 10;
                return static_cast<uint32_t>(x - a * 10);
            };
            integer tmp = me;
            std::vector<uint32_t> r_digits;
            while (integer::compare(tmp, 0))
            {
                std::vector<uint32_t> digits(tmp.values_.size());
                int i = tmp.values_.size() - 1;
                for (auto it = tmp.values_.rbegin(); it != tmp.values_.rend(); ++it, --i)
                {
                    digits[i] = gain_unit(*it);
                }
                uint32_t r = 0;
                for (unsigned k = 0; k < digits.size(); ++k)
                {
                    r += digits[k] * static_cast<uint32_t>(std::pow(6, k));
                }
                r_digits.push_back(gain_unit(r));
                tmp.unsign_div(10, nullptr);
            }
            for (auto it = r_digits.rbegin(); it != r_digits.rend(); ++it)
                s << (*it);
        }
        return s;
    }

    friend bool operator!=(const integer& x, const integer& y) { return !(x == y); }
    friend bool operator> (const integer& x, const integer& y) { return  (y < x); }
    friend bool operator<=(const integer& x, const integer& y) { return !(x > y); }
    friend bool operator>=(const integer& x, const integer& y) { return !(x < y); }

    friend integer& operator++(integer& x)      { return x += 1; }
    friend integer  operator++(integer& x, int) { integer nrv(x); ++x; return std::move(nrv); }
    friend integer& operator--(integer& x)      { return x -= 1; }
    friend integer  operator--(integer& x, int) { integer nrv(x); --x; return std::move(nrv); }

    friend integer operator+(const integer & x, const integer & y) { return std::move(integer(x) += y); }
    friend integer operator+(      integer&& x,       integer&& y) { return std::move(        x  += y); }
    friend integer operator+(      integer&& x, const integer & y) { return std::move(        x  += y); }
    friend integer operator+(const integer & x,       integer&& y) { return std::move(        y  += x); }
    friend integer operator-(const integer & x, const integer & y) { return std::move(integer(x) -= y); }
    friend integer operator-(      integer&& x,       integer&& y) { return std::move(        x  -= y); }
    friend integer operator-(      integer&& x, const integer & y) { return std::move(        x  -= y); }
    friend integer operator-(const integer & x,       integer&& y) { return std::move(        y  -= x); }
    friend integer operator*(const integer & x, const integer & y) { return std::move(integer(x) *= y); }
    friend integer operator*(      integer&& x,       integer&& y) { return std::move(        x  *= y); }
    friend integer operator*(      integer&& x, const integer & y) { return std::move(        x  *= y); }
    friend integer operator*(const integer & x,       integer&& y) { return std::move(        y  *= x); }
    friend integer operator/(const integer & x, const integer & y) { return std::move(integer(x) /= y); }
    friend integer operator/(      integer&& x,       integer&& y) { return std::move(        x  /= y); }
    friend integer operator/(      integer&& x, const integer & y) { return std::move(        x  /= y); }
    friend integer operator/(const integer & x,       integer&& y) { return std::move(        y  /= x); }
    friend integer operator%(const integer & x, const integer & y) { return std::move(integer(x) %= y); }
    friend integer operator%(      integer&& x,       integer&& y) { return std::move(        x  %= y); }
    friend integer operator%(      integer&& x, const integer & y) { return std::move(        x  %= y); }
    friend integer operator%(const integer & x,       integer&& y) { return std::move(        y  %= x); }
};
