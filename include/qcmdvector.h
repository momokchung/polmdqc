//////////////////////////////////////////////////////////////
//                                                          //
//  qcmdvector.h  --  define vectors of various dimensions  //
//                                                          //
//////////////////////////////////////////////////////////////


#pragma once

template <typename T>
class Vector1D 
{
private:
    T* vector1;
    size_t length;

public:
    Vector1D() : vector1(nullptr), length(0) {}

    Vector1D(size_t initialLength) : vector1(nullptr), length(0)
    {
        allocate(initialLength);
    }

    ~Vector1D() 
    {
        delete[] vector1;
    }

    void allocate(size_t newLength)
    {
        if (newLength == length) return;

        if (vector1) delete[] vector1;
        length = newLength;
        vector1 = new T[length];
    }

    size_t size() const
    {
        return length;
    }

    T& operator[](size_t index)
    {
        return vector1[index];
    }

    const T& operator[](size_t index) const {
        return vector1[index];
    }
};
