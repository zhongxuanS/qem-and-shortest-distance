#ifndef METRICS_H
#define METRICS_H
class Matrix {
private:
    double data[16];
public:
    Matrix() {
        for(int i = 0; i < 16; i++)
            data[i] = 0.0;
    }
    Matrix(const Matrix& m) {
        for(int i = 0; i < 16; i++)
            data[i] = m.data[i];
    }
    Matrix(double value) {
        for(int i = 0; i < 16; i++)
            data[i] = value;
    }
    Matrix(	double m11, double m12, double m13, double m14,
            double m21, double m22, double m23, double m24,
            double m31, double m32, double m33, double m34,
            double m41, double m42, double m43, double m44) {
        data[0] = m11;
        data[1] = m12;
        data[2] = m13;
        data[3] = m14;
        data[4] = m21;
        data[5] = m22;
        data[6] = m23;
        data[7] = m24;
        data[8] = m31;
        data[9] = m32;
        data[10] = m33;
        data[11] = m34;
        data[12] = m41;
        data[13] = m42;
        data[14] = m43;
        data[15] = m44;
    }

    /*
    *	用一个列向量plane去初始化一个矩阵M = plane × planeT
    */
    Matrix(double plane[]) {
        double a = plane[0];
        double b = plane[1];
        double c = plane[2];
        double d = plane[3];
        data[0] = a*a;
        data[1] = a*b;
        data[2] = a*c;
        data[3] = a*d;
        data[4] = a*b;
        data[5] = b*b;
        data[6] = b*c;
        data[7] = b*d;
        data[8] = a*c;
        data[9] = b*c;
        data[10] = c*c;
        data[11] = c*d;
        data[12] = a*d;
        data[13] = b*d;
        data[14] = c*d;
        data[15] = d*d;
    }
    /*
     *  计算矩阵的3*3子矩阵行列式
     *	参数为矩阵的索引号
     */
    double det(	int a11, int a12, int a13,
                int a21, int a22, int a23,
                int a31, int a32, int a33) {
        double det = data[a11]*data[a22]*data[a33] + data[a13]*data[a21]*data[a32] +
                     data[a12]*data[a23]*data[a31]
                     - data[a13]*data[a22]*data[a31] - data[a11]*data[a23]*data[a32] - data[a12]*data[a21]*data[a33];
        return det;
    }


    double operator[](int i) const {
        return data[i];
    }
    Matrix operator=(const Matrix& m) {
        for(int i = 0; i < 16; i++)
            data[i] = m.data[i];
        return *this;
    }
    Matrix operator+(const Matrix& m) {
        return Matrix(data[0] + m.data[0], data[1] + m.data[1],data[2] + m.data[2],data[3] + m.data[3],
                      data[4] + m.data[4], data[5] + m.data[5],data[6] + m.data[6],data[7] + m.data[7],
                      data[8] + m.data[8],data[9] + m.data[9],data[10] + m.data[10],data[11] + m.data[11],
                      data[12] + m.data[12],data[13] + m.data[13],data[14] + m.data[14],data[15] + m.data[15]);
    }
    Matrix operator+=(const Matrix& m) {
        for(int i = 0; i < 16; i++)
            data[i] += m[i];
        return *this;
    }
    void print() {
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                printf("%.8lf  ", data[i * 4 + j]);
            }
            printf("\n");
        }
    }
    const double* getData() const {
        return data;
    }
};
#endif