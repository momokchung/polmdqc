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

// MDQCArray      array of dimension m
// MDQCArray2D    array of dimention m x n, const n
// MDQCArray3D    array of dimention m x n x o, const n and o
// MDQCVector2D   array of dimention m x n, variable n

template <typename T>
struct MDQCArray
{
private:
    int len = 0;
    T *array;
    void alloc (int m)
    {
        len = m;
        array = new T[m];
    }

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
                deallocate();
                alloc(m);
            }
        }
        else {
            alloc(m);
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
    int len = 0;
    int len2 = 0;
    T (*array)[n];
    void alloc(int m)
    {
        len = m;
        len2 = n;
        array = new T[m][n];
    }

public:
    // overload [] operator for array access
    T (&operator[](int index))[n]
    {
        return array[index];
    }
    const T (&operator[](int index) const)[n]
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
    T (*ptr())[n] {
        return array;
    }

    // allocate array
    void allocate(int m)
    {
        if (array) {
            if (len != m) {
                deallocate();
                alloc(m);
            }
        }
        else {
            alloc(m);
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
    int len = 0;
    int len2 = 0;
    int len3 = 0;
    T (*array)[n][o];
    void alloc(int m)
    {
        len = m;
        len2 = n;
        len3 = o;
        array = new T[m][n][o];
    }

public:
    // overload [] operator for array access
    T (&operator[](int index))[n][o]
    {
        return array[index];
    }
    const T (&operator[](int index) const)[n][o]
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
    T (*ptr())[n][o] {
        return array;
    }

    // allocate array
    void allocate(int m)
    {
        if (array) {
            if (len != m) {
                deallocate();
                alloc(m);
            }
        }
        else {
            alloc(m);
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

template <typename T>
struct MDQCVector2D
{
private:
    int len = 0;
    int len2 = 0;
    T** array;
    void alloc(int m, int n)
    {
        len = m;
        len2 = n;
        array = new T*[len];
        for (int i = 0; i < len; i++) {
            array[i] = new T[len2];
        }
    }

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
    T** ptr() {
        return array;
    }

    // allocate array
    void allocate(int m, int n)
    {
        if (array) {
            if (len != m or len2 != n) {
                deallocate();
                alloc(m,n);
            }
        }
        else {
            alloc(m,n);
        }
    }

    // deallocate array
    void deallocate()
    {
        if (array) {
            for (int i = 0; i < len; i++) {
                delete[] array[i];
            }
            delete[] array;
            array = nullptr;
            len = 0;
            len2 = 0;
        }
    }
};

template <typename T>
struct MDQCVector3D
{
private:
    int len = 0;
    int len2 = 0;
    int len3 = 0;
    T*** array;
    void alloc(int m, int n, int o)
    {
        len = m;
        len2 = n;
        len3 = o;
        array = new T**[len];
        for (int i = 0; i < len; i++) {
            array[i] = new T*[len2];
            for (int j = 0; j < len2; j++) {
                array[i][j] = new T[len3];
            }
        }
    }

public:
    // overload [] operator for array access
    T** operator[](int index)
    {
        return array[index];
    }
    const T** operator[](int index) const
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
    int size3() {
        return len3;
    }

    // Return pointer to the underlying array
    T*** ptr() {
        return array;
    }

    // allocate array
    void allocate(int m, int n, int o)
    {
        if (array) {
            if (len != m or len2 != n or len3 != o) {
                deallocate();
                alloc(m,n,o);
            }
        }
        else {
            alloc(m,n,o);
        }
    }

    // deallocate array
    void deallocate()
    {
        if (array) {
            for (int i = 0; i < len; i++) {
                for (int j = 0; j < len2; j++) {
                    delete[] array[i][j];
                }
                delete[] array[i];
            }
            delete[] array;
            array = nullptr;
            len = 0;
            len2 = 0;
            len3 = 0;
        }
    }
};

template <typename T>
struct MDQCVector4D
{
private:
    int len = 0;
    int len2 = 0;
    int len3 = 0;
    int len4 = 0;
    T**** array;
    void alloc(int m, int n, int o, int p)
    {
        len = m;
        len2 = n;
        len3 = o;
        len4 = p;
        array = new T***[len];
        for (int i = 0; i < len; i++) {
            array[i] = new T**[len2];
            for (int j = 0; j < len2; j++) {
                array[i][j] = new T*[len3];
                for (int k = 0; k < len3; k++) {
                    array[i][j][k] = new T[len4];
                }
            }
        }
    }

public:
    // overload [] operator for array access
    T*** operator[](int index)
    {
        return array[index];
    }
    const T*** operator[](int index) const
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
    int size3() {
        return len3;
    }
    int size4() {
        return len4;
    }

    // Return pointer to the underlying array
    T**** ptr() {
        return array;
    }

    // allocate array
    void allocate(int m, int n, int o, int p)
    {
        if (array) {
            if (len != m or len2 != n or len3 != o or len4 != p) {
                deallocate();
                alloc(m,n,o,p);
            }
        }
        else {
            alloc(m,n,o,p);
        }
    }

    // deallocate array
    void deallocate()
    {
        if (array) {
            for (int i = 0; i < len; i++) {
                for (int j = 0; j < len2; j++) {
                    for (int k = 0; k < len3; k++) {
                        delete[] array[i][j][k];
                    }
                    delete[] array[i][j];
                }
                delete[] array[i];
            }
            delete[] array;
            array = nullptr;
            len = 0;
            len2 = 0;
            len3 = 0;
            len4 = 0;
        }
    }
};
}
