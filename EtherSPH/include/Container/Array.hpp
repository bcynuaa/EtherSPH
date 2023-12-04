/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-10-20 18:49:56
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: Container: header file for Array
 */

#ifndef ARRAY_HPP_
#define ARRAY_HPP_

#include <iostream>
#include <vector>

namespace Container {

template <typename ArrayElemType>
struct Array : public std::vector<ArrayElemType> {
    using std::vector<ArrayElemType>::vector;
    using std::vector<ArrayElemType>::operator=;
    using std::vector<ArrayElemType>::operator[];
    using std::vector<ArrayElemType>::size;
    using std::vector<ArrayElemType>::begin;
    using std::vector<ArrayElemType>::end;
    using std::vector<ArrayElemType>::push_back;
    using std::vector<ArrayElemType>::pop_back;
    using std::vector<ArrayElemType>::insert;
    using std::vector<ArrayElemType>::erase;
    using std::vector<ArrayElemType>::clear;
    using std::vector<ArrayElemType>::empty;
    using std::vector<ArrayElemType>::resize;
    using std::vector<ArrayElemType>::reserve;
    using std::vector<ArrayElemType>::shrink_to_fit;
    using std::vector<ArrayElemType>::front;
    using std::vector<ArrayElemType>::back;
    using std::vector<ArrayElemType>::data;
    using std::vector<ArrayElemType>::swap;
    using std::vector<ArrayElemType>::emplace;
    using std::vector<ArrayElemType>::emplace_back;
    using std::vector<ArrayElemType>::at;
    using std::vector<ArrayElemType>::capacity;
    using std::vector<ArrayElemType>::max_size;
    using std::vector<ArrayElemType>::assign;
    using std::vector<ArrayElemType>::get_allocator;

    // overload operator =
    Array<ArrayElemType>& operator=(const Array<ArrayElemType>& array) {
        if (this != &array) {
            this->resize(array.size());
            for (int i = 0; i < array.size(); i++) {
                (*this)[i] = array[i];
            }
        }
        else;
        return *this;
    };

    // operator == array
    bool operator==(const Array<ArrayElemType>& array) const {
        if (this->size() != array.size()) {
            return false;
        }
        else {
            for (int i = 0; i < this->size(); i++) {
                if ((*this)[i] != array[i]) {
                    return false;
                }
            }
            return true;
        }
    };

    // operator == elem
    bool operator==(const ArrayElemType& elem) const {
        for (int i = 0; i < this->size(); i++) {
            if ((*this)[i] != elem) {
                return false;
            }
        }
        return true;
    };

    // operator != array
    bool operator!=(const Array<ArrayElemType>& array) const {
        if (this->size() != array.size()) {
            return true;
        }
        else {
            for (int i = 0; i < this->size(); i++) {
                if ((*this)[i] != array[i]) {
                    return true;
                }
            }
            return false;
        }
    };

    // operator != elem
    bool operator!=(const ArrayElemType& elem) const {
        for (int i = 0; i < this->size(); i++) {
            if ((*this)[i] != elem) {
                return true;
            }
        }
        return false;
    };

    // operator <<
    friend std::ostream& operator<<(std::ostream& os, const Array<ArrayElemType>& array) {
        os << "[";
        for (int i = 0; i < array.size(); i++) {
            os << array[i];
            if (i != array.size() - 1) {
                os << ", ";
            }
            else;
        }
        os << "]";
        return os;
    };

    // operator >>
    friend std::istream& operator>>(std::istream& is, Array<ArrayElemType>& array) {
        for (int i = 0; i < array.size(); i++) {
            is >> array[i];
        }
        return is;
    };

};

template <typename ArrayElemType>
static void allocate(Array<ArrayElemType>& array, int size) {
    if (size < 1) {
        std::cout << "Error: Array::allocate(): size < 1" << std::endl;
        exit(1);
    }
    else if (array.size() != size) {
        array.resize(size);
    }
    else;
    return;
};

template <typename ArrayElemType>
static void allocate(Array<ArrayElemType>& array, int size, const ArrayElemType& init_elem) {
    if (size < 1) {
        std::cout << "Error: Array::allocate(): size < 1" << std::endl;
        exit(1);
    }
    else if (array.size() != size) {
        array.resize(size, init_elem);
    }
    else;
    return;
};

template <typename ArrayElemType>
static void allocate(Array<ArrayElemType>& array, ArrayElemType* ptr, int size) {
    if (size < 1) {
        std::cout << "Error: Array::allocate(): size < 1" << std::endl;
        exit(1);
    }
    else if (array.size() != size) {
        array.resize(size);
        for (int i = 0; i < size; i++) {
            array[i] = ptr[i];
        }
    }
    else;
    return;
};

template <typename ArrayElemType>
static void allocate(Array<ArrayElemType>& array, const Array<ArrayElemType>& array2) {
    if (array.size() != array2.size()) {
        array.resize(array2.size());
    }
    for (int i = 0; i < array2.size(); i++) {
        array[i] = array2[i];
    }
    return;
};

template <typename ArrayElemType>
static void fill(Array<ArrayElemType>& array, const ArrayElemType& elem) {
    for (int i = 0; i < array.size(); i++) {
        array[i] = elem;
    }
    return;
};

// apply1: apply a function to each element of an array
template <typename ArrayElemType>
Array<ArrayElemType> apply(const Array<ArrayElemType>& array, ArrayElemType (*func)(ArrayElemType)) {
    Array<ArrayElemType> result(array.size());
    for (int i = 0; i < array.size(); i++) {
        result[i] = func(array[i]);
    }
    return result;
};

// apply2: apply a function to each element of an array
template <typename ArrayElemType>
Array<ArrayElemType> apply(const Array<ArrayElemType>& array, ArrayElemType (*func)(const ArrayElemType&)) {
    Array<ArrayElemType> result(array.size());
    for (int i = 0; i < array.size(); i++) {
        result[i] = func(array[i]);
    }
    return result;
};

// apply3: apply a function to each element of an array
template <typename ArrayElemType>
Array<ArrayElemType> apply(const Array<ArrayElemType>& array, void (*func)(ArrayElemType&)) {
    Array<ArrayElemType> result(array);
    for (int i = 0; i < array.size(); i++) {
        func(result[i]);
    }
    return result;
};

// applyOn1: apply a function to each element of an array
template <typename ArrayElemType>
static void applyOn(Array<ArrayElemType>& array, ArrayElemType (*func)(ArrayElemType)) {
    for (int i = 0; i < array.size(); i++) {
        array[i] = func(array[i]);
    }
    return;
};

template <typename ArrayElemType>
static void applyOn(Array<ArrayElemType>& array, ArrayElemType (*func)(const ArrayElemType&)) {
    for (int i = 0; i < array.size(); i++) {
        array[i] = func(array[i]);
    }
    return;
};

template <typename ArrayElemType>
static void applyOn(Array<ArrayElemType>& array, void (*func)(ArrayElemType&)) {
    for (int i = 0; i < array.size(); i++) {
        func(array[i]);
    }
    return;
};

template <typename ArrayElemType>
static void append(Array<ArrayElemType>& array, const ArrayElemType& elem) {
    array.push_back(elem);
    return;
};

template <typename ArrayElemType>
static void append(Array<ArrayElemType>& array, const Array<ArrayElemType>& array2) {
    for (int i = 0; i < array2.size(); i++) {
        array.push_back(array2[i]);
    }
    return;
};

template <typename ArrayElemType>
static void pop(Array<ArrayElemType>& array) {
    if (array.size() == 0) {
        std::cout << "Warning: Array::pop(): array.size() == 0" << std::endl;
        return;
    }
    else {
        array.pop_back();
    }
    return;
};

template <typename ArrayElemType>

Array<ArrayElemType> connect(const Array<ArrayElemType>& array1, const Array<ArrayElemType>& array2) {
    Array<ArrayElemType> result(array1.size() + array2.size());
    for (int i = 0; i < array1.size(); i++) {
        result[i] = array1[i];
    }
    for (int i = 0; i < array2.size(); i++) {
        result[array1.size() + i] = array2[i];
    }
    return result;
};

template <typename ArrayElemType>
Array<ArrayElemType> slice(const Array<ArrayElemType>& array, const Array<int>& index_list) {
    Array<ArrayElemType> result(index_list.size());
    for (int i = 0; i < index_list.size(); i++) {
        result[i] = array[index_list[i]];
    }
    return result;
};

}; // namespace Container

#endif // ARRAY_HPP_