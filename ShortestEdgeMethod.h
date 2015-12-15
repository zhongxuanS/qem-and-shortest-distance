#pragma once
#include "Quadric.h"
#include <float.h>

typedef std::map<int, double> VertexCost;
typedef VertexCost::value_type VertexCostVal;
class ShortestEdgeMethod {
private:
    /* 顶点集合 */
    Vertices vertices;
    /* 面集合 */
    Faces faces;
    /* 记录每个顶点的cost */
    VertexCost vertexCost;
private:
    //解析obj文件
    void parse(FILE* file);
    // 初始化所有顶点cost
    void initCost();
    // 初始化点的邻点
    void initNeighbors();
    // 计算某一个顶点的cost和收缩对应点
    float calcVertexCollapseCost(int vertexId);
    // 收缩边操作,给定vertexId，按照该点记录的target顶点收缩
    void collapse(int vertexId);
    // 执行收缩以后的面关系更新
    void refreshFaces(int fromVertexId, int toVertexId);
    // 判断三角面中顶点是否合法
    bool isCorrectFace(Face face);
    // 刷新顶点邻点的cost
    void refreshNeighborsCost(std::set<int> neighbors);
    // 刷新邻点的邻点
    void refreshNeighbors(std::set<int> neighbors, int fromVertexId);
public:
    ShortestEdgeMethod(void);
    ~ShortestEdgeMethod(void);

    //读取obj文件，初始化vertices和faces
    void readObj(std::string filePath);

    void write_smf(const char* filename);

    // 简化
    void simplifyMesh(int targetVertices);
    double calculateEdgeLengthSqur(Vertex v0, Vertex v1);

};

