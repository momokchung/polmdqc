// Author: Moses KJ Chung
// Year:   2024

#pragma once

namespace polmdqc
{
//////////////////////////////////////////////////
//                                              //
//  darray  --  allocate and deallocate arrays  //
//                                              //
//////////////////////////////////////////////////

// "darray" defines functions that allocates and
// deallocates device dynamic array

// MDQCArray     array of dimension m
// MDQCArray2D   array of dimention m x n, const n
// MDQCArray3D   array of dimention m x n x o, const n and o

template <typename T>
struct MDQCArray
{
private:
    int len;
    T *array;

public:
    // overload [] operator for array access
    T& operator[](int index)
    {
        return array[index];
    }
    const T& operator[](int index) const
    {
        return array[index];
    }

    // return len
    int size() {
        return len;
    }

    // Return pointer to the underlying array
    T* ptr() {
        return array;
    }

    // allocate array
    void allocate(int m)
    {
        if (array) {
            if (len != m) {
                delete[] array;
                len = m;
                array = new T[m];
            }
        }
        else {
            len = m;
            array = new T[m];
        }
    }

    // deallocate array
    void deallocate()
    {
        if (array) {
            delete[] array;
            array = nullptr;
            len = 0;
        }
    }
};

template <typename T, const int n>
struct MDQCArray2D
{
private:
    int len;
    int len2;
    T (*array)[n];

public:
    // overload [] operator for array access
    T* operator[](int index)
    {
        return array[index];
    }
    const T* operator[](int index) const
    {
        return array[index];
    }

    // return len
    int size() {
        return len;
    }
    int size2() {
        return len2;
    }

    // Return pointer to the underlying array
    T* ptr() {
        return array;
    }

    // allocate array
    void allocate(int m)
    {
        if (array) {
            if (len != m) {
                delete[] array;
                len = m;
                len2 = n;
                array = new T[m][n];
            }
        }
        else {
            len = m;
            len2 = n;
            array = new T[m][n];
        }
    }

    // deallocate array
    void deallocate()
    {
        if (array) {
            delete[] array;
            array = nullptr;
            len = 0;
            len2 = 0;
        }
    }
};

template <typename T, const int n, const int o>
struct MDQCArray3D
{
private:
    int len;
    int len2;
    int len3;
    T (*array)[n][o];

public:
    // overload [] operator for array access
    T* operator[](int index)
    {
        return array[index];
    }
    const T* operator[](int index) const
    {
        return array[index];
    }

    // return size
    int size() {
        return len;
    }
    int size2() {
        return len2;
    }
    int size3() {
        return len3;
    }

    // Return pointer to the underlying array
    T* ptr() {
        return array;
    }

    // allocate array
    void allocate(int m)
    {
        if (array) {
            if (len != m) {
                delete[] array;
                len = m;
                len2 = n;
                len3 = o;
                array = new T[m][n][o];
            }
        }
        else {
            len = m;
            len2 = n;
            len3 = o;
            array = new T[m][n][o];
        }
    }

    // deallocate array
    void deallocate()
    {
        if (array) {
            delete[] array;
            array = nullptr;
            len = 0;
            len2 = 0;
            len3 = 0;
        }
    }
};
}
