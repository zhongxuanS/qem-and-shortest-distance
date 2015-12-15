#ifndef QUADRIC_H
#define QUADRIC_H
#include <map>
#include <deque>
#include <string>
#include <algorithm>
#include "Matrix.h"
#include <queue>
#include <set>

class Vertex {
public:
    double x;
    double y;
    double z;
    /* 用在最短边中，表示该顶点的度量值 */
    double cost;
    /* 用在最短边中，表示该顶点预计收缩顶点 */
    int collapseTarget;
    /* 该顶点的邻点 */
    std::set<int> neighbors;
};

class Face {
public:
    int idVertex[3];//顶点索引
    double plane[4];//面方程的4个参数ax+by+cz+d=0
};

typedef std::map<int, Vertex> Vertices;
typedef std::deque<Face> Faces;
typedef std::pair<int, int> Pair;//v1和v2的索引号
typedef std::map<Pair, double> Errors;
//typedef std::deque<Vsplit> Vsplits;
typedef std::map<int, Matrix> Matrices;

class Quadrics {
private:
    Vertices vertices;// 顶点集合
    Faces faces;// 面集合
    Errors errors;// 误差集合
    Matrices QMatrices;// Q矩阵集合
    double vertixError(double x, double y, double z, Matrix q) {
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[5]*y*y
               + 2*q[6]*y*z + 2*q[7]*y + q[10]*z*z + 2*q[11]*z + q[15];
    }
public:
    Quadrics();
    ~Quadrics();
    //读取obj文件，初始化vertices和faces
    void readObj(std::string filePath);
    //解析obj文件
    void parse(FILE* file);
    //初始化每个顶点的QMatrix
    void initQMatrices(void);
    //选择所有合法的pair，同时初始化errors
    void selectPairs();
    //已知Vi和Vj，计算该pair收缩以后的新坐标V，返回deltV
    double calculateErrors(int idVi, int idVj, double* vx = 0, double* vy = 0, double* vz = 0);
    //每次从errors中选择delta最小的pair，进行收缩处理，并且更新相关内容
    void constructContract(int targetFaceNum);
    void write_smf(const char* filename);
    const int getFacesNum() const {
        return faces.size();
    }
    const int getVerticesNum() const {
        return vertices.size();
    }
    // 计算简化以后和原始网格的相似度
    double calculateSimilarity(Quadrics& originQuadrics);
    const Vertices getVertices() const {
        return vertices;
    }
    double minDisSquaredToAllPlance(Vertex vertex);
};
#endif