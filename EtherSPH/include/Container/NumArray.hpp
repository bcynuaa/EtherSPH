/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-10-24 01:16:51
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: Container: header file for NumArray
 */

#ifndef NUM_ARRAY_HPP_
#define NUM_ARRAY_HPP_

#include <cmath>
#include "Container/Array.hpp"

namespace Container {

template <typename NumType>
struct NumArray : public Array<NumType> {
    using Array<NumType>::Array;
    using Array<NumType>::operator=;
    using Array<NumType>::operator[];
    using Array<NumType>::size;
    using Array<NumType>::begin;
    using Array<NumType>::end;
    using Array<NumType>::push_back;
    using Array<NumType>::pop_back;
    using Array<NumType>::insert;
    using Array<NumType>::erase;
    using Array<NumType>::clear;
    using Array<NumType>::empty;
    using Array<NumType>::resize;
    using Array<NumType>::reserve;
    using Array<NumType>::shrink_to_fit;
    using Array<NumType>::front;
    using Array<NumType>::back;
    using Array<NumType>::data;
    using Array<NumType>::swap;
    using Array<NumType>::emplace;
    using Array<NumType>::emplace_back;
    using Array<NumType>::at;
    using Array<NumType>::capacity;
    using Array<NumType>::max_size;
    using Array<NumType>::assign;
    using Array<NumType>::get_allocator;
    
    using Array<NumType>::operator==;
    using Array<NumType>::operator!=;

    // overload operator=
    NumArray<NumType>& operator=(const NumType& elem) {
        for (int i = 0; i < this->size(); i++) {
            (*this)[i] = elem;
        }
        return *this;
    };

    const bool operator<(const NumArray<NumType>& array) const {
        if (this->size() != array.size()) {
            std::cout << "Error: NumArray::operator<(): size not match" << std::endl;
            exit(1);
        }
        for (int i = 0; i < this->size(); i++) {
            if ((*this)[i] >= array[i]) {
                return false;
            }
        }
        return true;
    };

    const bool operator<(const NumType& elem) const {
        for (int i = 0; i < this->size(); i++) {
            if ((*this)[i] >= elem) {
                return false;
            }
        }
        return true;
    };

    friend const bool operator<(const NumType& elem, const NumArray<NumType>& array) {
        for (int i = 0; i < array.size(); i++) {
            if (elem >= array[i]) {
                return false;
            }
        }
        return true;
    };

    const bool operator<=(const NumArray<NumType>& array) const {
        if (this->size() != array.size()) {
            std::cout << "Error: NumArray::operator<=(): size not match" << std::endl;
            exit(1);
        }
        for (int i = 0; i < this->size(); i++) {
            if ((*this)[i] > array[i]) {
                return false;
            }
        }
        return true;
    };

    const bool operator<=(const NumType& elem) const {
        for (int i = 0; i < this->size(); i++) {
            if ((*this)[i] > elem) {
                return false;
            }
        }
        return true;
    };

    friend const bool operator<=(const NumType& elem, const NumArray<NumType>& array) {
        for (int i = 0; i < array.size(); i++) {
            if (elem > array[i]) {
                return false;
            }
        }
        return true;
    };

    const bool operator>(const NumArray<NumType>& array) const {
        if (this->size() != array.size()) {
            std::cout << "Error: NumArray::operator>(): size not match" << std::endl;
            exit(1);
        }
        for (int i = 0; i < this->size(); i++) {
            if ((*this)[i] <= array[i]) {
                return false;
            }
        }
        return true;
    };

    const bool operator>(const NumType& elem) const {
        for (int i = 0; i < this->size(); i++) {
            if ((*this)[i] <= elem) {
                return false;
            }
        }
        return true;
    };

    friend const bool operator>(const NumType& elem, const NumArray<NumType>& array) {
        for (int i = 0; i < array.size(); i++) {
            if (elem <= array[i]) {
                return false;
            }
        }
        return true;
    };

    const bool operator>=(const NumArray<NumType>& array) const {
        if (this->size() != array.size()) {
            std::cout << "Error: NumArray::operator>=(): size not match" << std::endl;
            exit(1);
        }
        for (int i = 0; i < this->size(); i++) {
            if ((*this)[i] < array[i]) {
                return false;
            }
        }
        return true;
    };

    const bool operator>=(const NumType& elem) const {
        for (int i = 0; i < this->size(); i++) {
            if ((*this)[i] < elem) {
                return false;
            }
        }
        return true;
    };

    friend const bool operator>=(const NumType& elem, const NumArray<NumType>& array) {
        for (int i = 0; i < array.size(); i++) {
            if (elem < array[i]) {
                return false;
            }
        }
        return true;
    };

    const NumArray<NumType> operator+(const NumArray<NumType>& array) const {
        if (this->size() != array.size()) {
            std::cout << "Error: NumArray::operator+(): size not match" << std::endl;
            exit(1);
        }
        NumArray<NumType> result(this->size());
        for (int i = 0; i < this->size(); i++) {
            result[i] = (*this)[i] + array[i];
        }
        return result;
    };

    const NumArray<NumType> operator-(const NumArray<NumType>& array) const {
        if (this->size() != array.size()) {
            std::cout << "Error: NumArray::operator-(): size not match" << std::endl;
            exit(1);
        }
        NumArray<NumType> result(this->size());
        for (int i = 0; i < this->size(); i++) {
            result[i] = (*this)[i] - array[i];
        }
        return result;
    };

    const NumArray<NumType> operator-() const {
        NumArray<NumType> result(this->size());
        for (int i = 0; i < this->size(); i++) {
            result[i] = -(*this)[i];
        }
        return result;
    };

    const NumArray<NumType> operator*(const NumArray<NumType>& array) const {
        if (this->size() != array.size()) {
            std::cout << "Error: NumArray::operator*(): size not match" << std::endl;
            exit(1);
        }
        NumArray<NumType> result(this->size());
        for (int i = 0; i < this->size(); i++) {
            result[i] = (*this)[i] * array[i];
        }
        return result;
    };

    const NumArray<NumType> operator/(const NumArray<NumType>& array) const {
        if (this->size() != array.size()) {
            std::cout << "Error: NumArray::operator/(): size not match" << std::endl;
            exit(1);
        }
        NumArray<NumType> result(this->size());
        for (int i = 0; i < this->size(); i++) {
            result[i] = (*this)[i] / array[i];
        }
        return result;
    };

    const NumArray<NumType> operator+(const NumType& elem) const {
        NumArray<NumType> result(this->size());
        for (int i = 0; i < this->size(); i++) {
            result[i] = (*this)[i] + elem;
        }
        return result;
    };

    const NumArray<NumType> operator-(const NumType& elem) const {
        NumArray<NumType> result(this->size());
        for (int i = 0; i < this->size(); i++) {
            result[i] = (*this)[i] - elem;
        }
        return result;
    };

    const NumArray<NumType> operator*(const NumType& elem) const {
        NumArray<NumType> result(this->size());
        for (int i = 0; i < this->size(); i++) {
            result[i] = (*this)[i] * elem;
        }
        return result;
    };

    const NumArray<NumType> operator/(const NumType& elem) const {
        NumArray<NumType> result(this->size());
        for (int i = 0; i < this->size(); i++) {
            result[i] = (*this)[i] / elem;
        }
        return result;
    };

    friend const NumArray<NumType> operator+(const NumType& elem, const NumArray<NumType>& array) {
        NumArray<NumType> result(array.size());
        for (int i = 0; i < array.size(); i++) {
            result[i] = elem + array[i];
        }
        return result;
    };

    friend const NumArray<NumType> operator-(const NumType& elem, const NumArray<NumType>& array) {
        NumArray<NumType> result(array.size());
        for (int i = 0; i < array.size(); i++) {
            result[i] = elem - array[i];
        }
        return result;
    };

    friend const NumArray<NumType> operator*(const NumType& elem, const NumArray<NumType>& array) {
        NumArray<NumType> result(array.size());
        for (int i = 0; i < array.size(); i++) {
            result[i] = elem * array[i];
        }
        return result;
    };

    friend const NumArray<NumType> operator/(const NumType& elem, const NumArray<NumType>& array) {
        NumArray<NumType> result(array.size());
        for (int i = 0; i < array.size(); i++) {
            result[i] = elem / array[i];
        }
        return result;
    };
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
* define addWith(), addBy(), subtractWith(), subtractBy(), multiplyWith(), multiplyBy(), divideWith(), divideBy()
*/

template <typename NumType>
static void addWith(NumArray<NumType>& array, const NumType& elem) {
    for (int i = 0; i < array.size(); i++) {
        array[i] += elem;
    }
};

template <typename NumType>
static void addBy(NumArray<NumType>& array, const NumType& elem) {
    addWith(array, elem);
};

// overload += operator
template <typename NumType>
static void operator+=(NumArray<NumType>& array, const NumType& elem) {
    addWith(array, elem);
    return;
};

template <typename NumType>
static void addWith(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    if (array1.size() != array2.size()) {
        std::cout << "Error: addWith(): size not match" << std::endl;
        exit(1);
    }
    for (int i = 0; i < array1.size(); i++) {
        array1[i] += array2[i];
    }
};

template <typename NumType>
static void addBy(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    addWith(array1, array2);
};

// overload += operator
template <typename NumType>
static void operator+=(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    addWith(array1, array2);
    return;
};

template <typename NumType>
static void subtractWith(NumArray<NumType>& array, const NumType& elem) {
    for (int i = 0; i < array.size(); i++) {
        array[i] -= elem;
    }
};

template <typename NumType>
static void subtractBy(NumArray<NumType>& array, const NumType& elem) {
    for (int i = 0; i < array.size(); i++) {
        array[i] = elem - array[i];
    }
};

// overload -= operator
template <typename NumType>
static void operator-=(NumArray<NumType>& array, const NumType& elem) {
    subtractWith(array, elem);
    return;
};

template <typename NumType>
static void subtractWith(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    if (array1.size() != array2.size()) {
        std::cout << "Error: subtractWith(): size not match" << std::endl;
        exit(1);
    }
    for (int i = 0; i < array1.size(); i++) {
        array1[i] -= array2[i];
    }
};

template <typename NumType>
static void subtractBy(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    if (array1.size() != array2.size()) {
        std::cout << "Error: subtractBy(): size not match" << std::endl;
        exit(1);
    }
    for (int i = 0; i < array1.size(); i++) {
        array1[i] = array2[i] - array1[i];
    }
};

// overload -= operator
template <typename NumType>
static void operator-=(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    subtractWith(array1, array2);
    return;
};

template <typename NumType>
static void multiplyWith(NumArray<NumType>& array, const NumType& elem) {
    for (int i = 0; i < array.size(); i++) {
        array[i] *= elem;
    }
};

template <typename NumType>
static void multiplyBy(NumArray<NumType>& array, const NumType& elem) {
    multiplyWith(array, elem);
};

// overload *= operator
template <typename NumType>
static void operator*=(NumArray<NumType>& array, const NumType& elem) {
    multiplyWith(array, elem);
    return;
};

template <typename NumType>
static void multiplyWith(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    if (array1.size() != array2.size()) {
        std::cout << "Error: multiplyWith(): size not match" << std::endl;
        exit(1);
    }
    for (int i = 0; i < array1.size(); i++) {
        array1[i] *= array2[i];
    }
};

template <typename NumType>
static void multiplyBy(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    multiplyWith(array1, array2);
};

// overload *= operator
template <typename NumType>
static void operator*=(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    multiplyWith(array1, array2);
    return;
};

template <typename NumType>
static void divideWith(NumArray<NumType>& array, const NumType& elem) {
    for (int i = 0; i < array.size(); i++) {
        array[i] /= elem;
    }
};

template <typename NumType>
static void divideBy(NumArray<NumType>& array, const NumType& elem) {
    for (int i = 0; i < array.size(); i++) {
        array[i] = elem / array[i];
    }
};

// overload /= operator
template <typename NumType>
static void operator/=(NumArray<NumType>& array, const NumType& elem) {
    divideWith(array, elem);
    return;
};

template <typename NumType>
static void divideWith(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    if (array1.size() != array2.size()) {
        std::cout << "Error: divideWith(): size not match" << std::endl;
        exit(1);
    }
    for (int i = 0; i < array1.size(); i++) {
        array1[i] /= array2[i];
    }
};

template <typename NumType>
static void divideBy(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    if (array1.size() != array2.size()) {
        std::cout << "Error: divideBy(): size not match" << std::endl;
        exit(1);
    }
    for (int i = 0; i < array1.size(); i++) {
        array1[i] = array2[i] / array1[i];
    }
};

// overload /= operator
template <typename NumType>
static void operator/=(NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    divideWith(array1, array2);
    return;
};

/*
* end define addWith(), addBy(), subtractWith(), subtractBy(), multiplyWith(), multiplyBy(), divideWith(), divideBy()
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// function library

template <typename NumType>
NumType sum(const NumArray<NumType>& array) {
    NumType result = 0;
    for (int i = 0; i < array.size(); ++i) {
        result += array[i];
    }
    return result;
};

template <typename NumType>
NumType prod(const NumArray<NumType>& array) {
    NumType result = 1;
    for (int i = 0; i < array.size(); ++i) {
        result *= array[i];
    }
    return result;
};

template <typename NumType>
double mean(const NumArray<NumType>& array) {
    return static_cast<double>(sum(array)) / array.size();
};

template <typename NumType>
NumType minimum(const NumArray<NumType>& array) {
    NumType result = array[0];
    for (int i = 1; i < array.size(); ++i) {
        if (array[i] < result) {
            result = array[i];
        }
    }
    return result;
};

template <typename NumType>
NumType maximum(const NumArray<NumType>& array) {
    NumType result = array[0];
    for (int i = 1; i < array.size(); ++i) {
        if (array[i] > result) {
            result = array[i];
        }
    }
    return result;
};

template <typename NumType>
NumArray<NumType> zeros(const int num) {
    NumArray<NumType> array(num);
    return array;
};

template <typename NumType>
NumArray<NumType> ones(const int num) {
    NumArray<NumType> array(num);
    for (int i = 0; i < num; ++i) {
        array[i] = 1;
    }
    return array;
};

template <typename NumType>
NumArray<NumType> linspace(const NumType& start, const NumType& end, const int num) {
    NumArray<NumType> array(num);
    for (int i = 0; i < num; ++i) {
        array[i] = start + i * (end - start) / (num - 1);
    }
    return array;
};

template <typename NumType>
NumArray<NumType> arange(const NumType& start, const NumType& end, const NumType& step) {
    const int num = static_cast<int>((end - start) / step);
    NumArray<NumType> array(num);
    for (int i = 0; i < num; ++i) {
        array[i] = start + i * step;
    };
    return array;
};

template <typename NumType>
NumType dot(const NumArray<NumType>& array1, const NumArray<NumType>& array2) {
    if (array1.size() != array2.size()) {
        std::cout << "Error: dot(): size not match" << std::endl;
        exit(1);
    }
    NumType result = 0;
    for (int i = 0; i < array1.size(); i++) {
        result += array1[i] * array2[i];
    }
    return result;
};

template <typename NumType>
double norm(const NumArray<NumType>& array) {
    return std::sqrt(static_cast<double>(dot(array, array)));
};

// function library
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// * declare IntArray and DoubleArray

using IntArray = NumArray<int>;
using DoubleArray = NumArray<double>;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}; // namespace Container

#endif // NUM_ARRAY_HPP_