/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-10-24 19:42:25
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: Container: test for Array
 */

#include <iostream>

#include "Container/Container.hpp"

// suppoerts functions declaration as follows

int add1(const int& a) {
    return a + 1;
}

int add2(int b) {
    return b + 2;
}

void add3(int& c) {
    c += 3;
    return;
}

int main(int argc, char* argv[]) {
    Container::Array<int> array0;
    Container::allocate(array0, 10);
    std::cout << "- constructor0, array0" << std::endl;
    std::cout << "  use allocate(array0, 10)" << std::endl;
    std::cout << array0 << std::endl;
    Container::Array<int> array1(10);
    Container::fill(array1, 1);
    std::cout << "- constructor1, array1" << std::endl;
    std::cout << "  use fill(array1, 1)" << std::endl;
    std::cout << array1 << std::endl;
    Container::Array<int> array2(10, 2);
    std::cout << "- constructor2, array2" << std::endl;
    std::cout << array2 << std::endl;
    Container::Array<int> array3(array2);
    std::cout << "- constructor3, array3" << std::endl;
    std::cout << "  use array2" << std::endl;
    std::cout << array3 << std::endl;
    Container::Array<int> array4 = {1, 2, 3};
    std::cout << "- constructor4, array4" << std::endl;
    std::cout << "  use int p[3] = {1, 2, 3}" << std::endl;
    std::cout << array4 << std::endl;
    Container::Array<int> array5 = array4;
    std::cout << "- constructor5, array5" << std::endl;
    std::cout << "  use = array4" << std::endl;
    std::cout << array5 << std::endl;

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "array1 == array2: " << (array1 == array2) << std::endl;
    std::cout << "array5 == array4: " << (array1 == array1) << std::endl;

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Container::allocate(array0, 10, 100)" << std::endl;
    Container::allocate(array0, 10, 100);
    std::cout << array0 << std::endl;
    std::cout << "Container::fill(array0, 8)" << std::endl;
    Container::fill(array0, 8);
    std::cout << array0 << std::endl;

    std::cout << "----------------------------------------" << std::endl;
    Container::Array<int> array6 = Container::apply(array0, add1);
    std::cout << "Container::apply(array0, add1)" << std::endl;
    std::cout << array6 << std::endl;
    std::cout << "Container::apply(array0, add2)" << std::endl;
    std::cout << Container::apply(array0, add2) << std::endl;
    std::cout << "Container::apply(array0, add3)" << std::endl;
    Container::apply(array0, add3);
    std::cout << array0 << std::endl;

    std::cout << "----------------------------------------" << std::endl;
    Container::applyOn(array0, add1);
    std::cout << "Container::applyOn(array0, add1)" << std::endl;
    std::cout << array0 << std::endl;
    Container::applyOn(array0, add2);
    std::cout << "Container::applyOn(array0, add2)" << std::endl;
    std::cout << array0 << std::endl;
    Container::applyOn(array0, add3);
    std::cout << "Container::applyOn(array0, add3)" << std::endl;
    std::cout << array0 << std::endl;

    std::cout << "----------------------------------------" << std::endl;
    Container::append(array0, 1);
    std::cout << "Container::append(array0, 1)" << std::endl;
    std::cout << array0 << std::endl;
    Container::append(array0, array1);
    std::cout << "Container::append(array0, array1)" << std::endl;
    std::cout << array0 << std::endl;
    Container::pop(array0);
    std::cout << "Container::pop(array0)" << std::endl;
    std::cout << array0 << std::endl;
    Container::Array<int> array7 = Container::connect(array0, array1);
    std::cout << "Container::connect(array0, array1)" << std::endl;
    std::cout << array7 << std::endl;
    Container::Array<int> index_list(7);
    for (int i = 0; i < index_list.size(); i++) {
        index_list[i] = 2*i + 1;
    }
    Container::Array<int> array8 = Container::slice(array7, index_list);
    std::cout << "Container::slice(array7, index_list)" << std::endl;
    std::cout << array8 << std::endl;
    return 0;
};