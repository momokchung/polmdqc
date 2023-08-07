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
        deallocate();
    }

    void allocate(size_t newLength)
    {
        if (newLength == length) return;

        deallocate();
        length = newLength;
        vector1 = new T[length];
    }

    void deallocate() {
        if (vector1) {
            delete[] vector1;
            vector1 = nullptr;
            length = 0;
        }
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

template <typename T>
class Vector2D {
private:
    T** vector2;
    size_t rows;
    size_t cols;

public:
    Vector2D() : vector2(nullptr), rows(0), cols(0) {}

    Vector2D(size_t numRows, size_t numCols) : vector2(nullptr), rows(0), cols(0) {
        allocate(numRows, numCols);
    }

    ~Vector2D() {
        deallocate();
    }

    void allocate(size_t numRows, size_t numCols) {
        if (numRows == rows && numCols == cols) return;

        deallocate();

        rows = numRows;
        cols = numCols;

        vector2 = new T*[rows];
        for (size_t i = 0; i < rows; ++i) {
            vector2[i] = new T[cols];
        }
    }

    void deallocate() {
        if (vector2) {
            for (size_t i = 0; i < rows; ++i) {
                delete[] vector2[i];
            }
            delete[] vector2;
            vector2 = nullptr;
            rows = 0;
            cols = 0;
        }
    }

    size_t size() const {
        return rows;
    }

    size_t numRows() const {
        return rows;
    }

    size_t numCols() const {
        return cols;
    }

    T*& operator[](size_t rowIndex) {
        return vector2[rowIndex];
    }

    const T*& operator[](size_t rowIndex) const {
        return vector2[rowIndex];
    }
};

template <typename T>
class Vector3D {
private:
    T*** vector3;
    size_t xSize;
    size_t ySize;
    size_t zSize;

public:
    Vector3D() : vector3(nullptr), xSize(0), ySize(0), zSize(0) {}

    Vector3D(size_t sizeX, size_t sizeY, size_t sizeZ)
        : vector3(nullptr), xSize(0), ySize(0), zSize(0) {
        allocate(sizeX, sizeY, sizeZ);
    }

    ~Vector3D() {
        deallocate();
    }

    void allocate(size_t sizeX, size_t sizeY, size_t sizeZ) {
        if (sizeX == xSize && sizeY == ySize && sizeZ == zSize) return;

        deallocate();

        xSize = sizeX;
        ySize = sizeY;
        zSize = sizeZ;

        vector3 = new T**[xSize];
        for (size_t i = 0; i < xSize; ++i) {
            vector3[i] = new T*[ySize];
            for (size_t j = 0; j < ySize; ++j) {
                vector3[i][j] = new T[zSize];
            }
        }
    }

    void deallocate() {
        if (vector3) {
            for (size_t i = 0; i < xSize; ++i) {
                for (size_t j = 0; j < ySize; ++j) {
                    delete[] vector3[i][j];
                }
                delete[] vector3[i];
            }
            delete[] vector3;
            vector3 = nullptr;
            xSize = 0;
            ySize = 0;
            zSize = 0;
        }
    }

    size_t size() const {
        return xSize;
    }

    size_t sizeX() const {
        return xSize;
    }

    size_t sizeY() const {
        return ySize;
    }

    size_t sizeZ() const {
        return zSize;
    }

    T** operator[](size_t xIndex) {
        return vector3[xIndex];
    }

    const T** operator[](size_t xIndex) const {
        return vector3[xIndex];
    }
};

template <typename T>
class Vector4D {
private:
    T**** vector4;
    size_t wSize;
    size_t xSize;
    size_t ySize;
    size_t zSize;

public:
    Vector4D() : vector4(nullptr), wSize(0), xSize(0), ySize(0), zSize(0) {}

    Vector4D(size_t sizeW, size_t sizeX, size_t sizeY, size_t sizeZ)
        : vector4(nullptr), wSize(0), xSize(0), ySize(0), zSize(0) {
        allocate(sizeW, sizeX, sizeY, sizeZ);
    }

    ~Vector4D() {
        deallocate();
    }

    void allocate(size_t sizeW, size_t sizeX, size_t sizeY, size_t sizeZ) {
        if (sizeW == wSize && sizeX == xSize && sizeY == ySize && sizeZ == zSize) return;

        deallocate();

        wSize = sizeW;
        xSize = sizeX;
        ySize = sizeY;
        zSize = sizeZ;

        vector4 = new T***[wSize];
        for (size_t i = 0; i < wSize; ++i) {
            vector4[i] = new T**[xSize];
            for (size_t j = 0; j < xSize; ++j) {
                vector4[i][j] = new T*[ySize];
                for (size_t k = 0; k < ySize; ++k) {
                    vector4[i][j][k] = new T[zSize];
                }
            }
        }
    }

    void deallocate() {
        if (vector4) {
            for (size_t i = 0; i < wSize; ++i) {
                for (size_t j = 0; j < xSize; ++j) {
                    for (size_t k = 0; k < ySize; ++k) {
                        delete[] vector4[i][j][k];
                    }
                    delete[] vector4[i][j];
                }
                delete[] vector4[i];
            }
            delete[] vector4;
            vector4 = nullptr;
            wSize = 0;
            xSize = 0;
            ySize = 0;
            zSize = 0;
        }
    }

    size_t size() const {
        return wSize;
    }

    size_t sizeW() const {
        return wSize;
    }

    size_t sizeX() const {
        return xSize;
    }

    size_t sizeY() const {
        return ySize;
    }

    size_t sizeZ() const {
        return zSize;
    }

    T*** operator[](size_t wIndex) {
        return vector4[wIndex];
    }

    const T*** operator[](size_t wIndex) const {
        return vector4[wIndex];
    }
};
