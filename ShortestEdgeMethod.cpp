#include "StdAfx.h"
#include "ShortestEdgeMethod.h"


ShortestEdgeMethod::ShortestEdgeMethod(void) {}


ShortestEdgeMethod::~ShortestEdgeMethod(void) {}

//读取filePath的obj文件，并且利用该obj文件中的内容初始化
//vertices和faces字段，obj文件要求只包含点和面的信息
void ShortestEdgeMethod::readObj(std::string filePath) {
    if(filePath.empty() == true) {
        printf("readObj参数非法！\n");
        return;
    }
    FILE *file;
    errno_t err = fopen_s(&file, &filePath[0], "r");
    if (err != 0) {
        fprintf(stderr, "read_smf() failed: can't open data file \"%s\".\n", &filePath[0]);
        system("PAUSE");
        exit(-1);
    }
    parse(file);
    fclose(file);
}
//解析obj文件中的内容
void ShortestEdgeMethod::parse(FILE* file) {
    if(file == NULL)
        return;

    char ch;
    char buf[1024];
    Vertex v;
    Face f;
    int local_num_vertices = 0;
    double x0, y0, z0;        /* ax + by + cz = 0 */
    double x1, y1, z1;
    double x2, y2, z2;
    double a, b, c, M;

    while ( fscanf(file, "%c", &ch) != EOF ) {
        switch(ch) {
        case ' ' :         /* blanks */
        case '\t':
        case '\n':
            continue;//continue对while生效，switch中只有break生效
        case '#':/* comment */
        case 'g':
            fgets(buf, sizeof(buf), file);//fgets读到换行符为止，或者buf满了
            break;
        case 'v':          /* vertex */
        case 'V':
            local_num_vertices++;							/* vertex index starts from 1 */
            fscanf_s(file, "%lf %lf %lf",	&v.x, &v.y, &v.z);
            vertices.insert(Vertices::value_type(local_num_vertices, v));
            break;
        case 'f':          /* face */
        case 'F':
            fscanf_s(file, "%d %d %d", &f.idVertex[0], &f.idVertex[1], &f.idVertex[2]);
            x0 = vertices[f.idVertex[0]].x;
            y0 = vertices[f.idVertex[0]].y;
            z0 = vertices[f.idVertex[0]].z;
            x1 = vertices[f.idVertex[1]].x;
            y1 = vertices[f.idVertex[1]].y;
            z1 = vertices[f.idVertex[1]].z;
            x2 = vertices[f.idVertex[2]].x;
            y2 = vertices[f.idVertex[2]].y;
            z2 = vertices[f.idVertex[2]].z;
            a = (y1-y0)*(z2-z0) - (z1-z0)*(y2-y0);   /* a1*b2 - a2*b1;        */
            b = (z1-z0)*(x2-x0) - (x1-x0)*(z2-z0);   /* a2*b0 - a0*b2;        */
            c = (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);   /* a0*b1 - a1*b0;        */
            M = sqrt(a*a + b*b + c*c);
            a = a/M;
            b = b/M;
            c = c/M;
            f.plane[0] = a;
            f.plane[1] = b;
            f.plane[2] = c;
            f.plane[3] = -1*(a*x0 + b*y0 + c*z0);				/* -1*(a*x + b*y + c*z); */
            faces.push_back(f);
            break;
        default:          /* invalid commands */
            fgets(buf, sizeof(buf), file);
            fprintf(stderr, "Parse() failed: invalid attributes: \"%c\".\n", ch);
            system("PAUSE");
            exit(-2);
        }
    }
}

void ShortestEdgeMethod::write_smf(const char* filename) {
    FILE *file;
    fopen_s(&file, filename, "w");
    if (file == NULL) {
        fprintf(stderr, "write_smf() failed: can't write data file \"%s\".\n", filename);
        system("PAUSE");
        exit(-1);
    }

    /* print header info */
    fprintf(file, "#$PM 0.1\n");
    fprintf(file, "#$vertices %d\n", vertices.size());
    fprintf(file, "#$faces %d\n", faces.size());
    fprintf(file, "#\n");

    /* print vertices */
    int ii = 1;
    std::map<int, int> imap;
    for (Vertices::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
        fprintf(file, "v %lf %lf %lf\n", /*iter->first, */iter->second.x, iter->second.y, iter->second.z);
        imap.insert(std::map<int, int>::value_type(iter->first, ii));
        ii++;
    }
    int v1, v2, v3;
    for (int i = 0; i < static_cast<int>(faces.size()); i++) {
        v1 = imap[faces[i].idVertex[0]];
        v2 = imap[faces[i].idVertex[1]];
        v3 = imap[faces[i].idVertex[2]];
        fprintf(file, "f %d %d %d\n", v1, v2, v3);
    }
    fclose(file);
}
void ShortestEdgeMethod::initNeighbors() {
    for(Faces::iterator iter = faces.begin(); iter != faces.end(); ++iter) {
        Face face = *iter;
        int v0 = face.idVertex[0];
        int v1 = face.idVertex[1];
        int v2 = face.idVertex[2];

        vertices[v0].neighbors.insert(v1);
        vertices[v0].neighbors.insert(v2);

        vertices[v1].neighbors.insert(v0);
        vertices[v1].neighbors.insert(v2);

        vertices[v2].neighbors.insert(v0);
        vertices[v2].neighbors.insert(v1);
    }
}
float ShortestEdgeMethod::calcVertexCollapseCost(int vertexId) {
    if(vertices.find(vertexId) == vertices.end()) {
        fprintf(stderr, "calcVertexCollapseCost vertexId:%d error!\n", vertexId);
        return FLT_MAX;
    }

    std::set<int>& neighbors = vertices[vertexId].neighbors;
    float minLength = FLT_MAX;
    int collapseTarget = -1;
    for(std::set<int>::iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        float length = (float)calculateEdgeLengthSqur(vertices[vertexId], vertices[*iter]);
        if(length < minLength) {
            collapseTarget = *iter;
            minLength = length;
        }
    }
    vertices[vertexId].cost = minLength;
    vertices[vertexId].collapseTarget = collapseTarget;
    return minLength;
}


// 初始化所有顶点cost
void ShortestEdgeMethod::initCost() {
    for(Vertices::iterator iter = vertices.begin(); iter != vertices.end(); ++iter) {
        int vertexId = iter->first;
        float cost = calcVertexCollapseCost(vertexId);
        vertexCost.insert(VertexCostVal(vertexId, cost));
    }
}

bool ShortestEdgeMethod::isCorrectFace(Face face) {
    if(face.idVertex[0] == face.idVertex[1]) return false;
    else if(face.idVertex[0] == face.idVertex[2]) return false;
    else if(face.idVertex[1] == face.idVertex[2]) return false;
    else return true;
}
// 刷新面关系
void ShortestEdgeMethod::refreshFaces(int fromVertexId, int toVertexId) {
    bool changeFromVertex = false;
    for(Faces::iterator iter = faces.begin(); iter != faces.end(); ) {
        for(int i = 0; i < 3; ++i) {
            if(iter->idVertex[i] == fromVertexId) {
                iter->idVertex[i] = toVertexId;
                changeFromVertex = true;
                break;
            }
        }
        // 将from变成to以后，会出现一个面有两个相同的顶点，删除
        if(changeFromVertex && !isCorrectFace(*iter)) {
            iter = faces.erase(iter);
            changeFromVertex = false;
        } else {
            ++iter;
        }
    }
}
void ShortestEdgeMethod::refreshNeighborsCost(std::set<int> neighbors) {
    for(std::set<int>::iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        int neighborId = *iter;
        float cost = calcVertexCollapseCost(neighborId);
        if(vertices.find(neighborId) == vertices.end()) {
            fprintf(stderr, "refreshNeighborsCost failed! neighborId:%d\n", neighborId);
            return;
        }
        vertexCost[neighborId] = cost;
    }
}
void ShortestEdgeMethod::refreshNeighbors(std::set<int> neighbors, int fromVertexId) {
    // 将neighbors中顶点的邻点fromVertexId删掉
    for(std::set<int>::iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        Vertex& v = vertices[*iter];
        std::set<int>& n = v.neighbors;
        for(std::set<int>::iterator iterN = n.begin(); iterN != n.end();) {
            if(*iterN == fromVertexId) {
                iterN = n.erase(iterN);
            } else {
                ++iterN;
            }
        }
    }
    // 重新初始化所有顶点的邻点
    initNeighbors();
}
void ShortestEdgeMethod::collapse(int vertexId) {
    if(vertices.find(vertexId) == vertices.end()) {
        fprintf(stderr, "collapse function error,vertexId:%d\n", vertexId);
        return;
    }
    //1、拿到顶点的邻点
    std::set<int> neighbors = vertices[vertexId].neighbors;
    //2、更新面关系
    int collapseTarget = vertices[vertexId].collapseTarget;
    // 若该点没有要收缩的对象
    if(collapseTarget == -1) {
        vertices.erase(vertexId);
        return;
    }
    refreshFaces(vertexId, collapseTarget);
    //3、根据原来邻点更新所有邻点的cost和target，还有vertexCost
    refreshNeighbors(neighbors, vertexId);
    refreshNeighborsCost(neighbors);
    //4、删除vertexId
    vertices.erase(vertexId);
}
void ShortestEdgeMethod::simplifyMesh(int targetVertices) {
    if(targetVertices <= 0) {
        fprintf(stderr, "targetVertices must bigger than zero\n");
        return;
    }
    initNeighbors();
    initCost();
    while(vertices.size() > targetVertices) {
        VertexCost::iterator deleteIter = vertexCost.begin();
        VertexCost::iterator minCostIter = vertexCost.begin();
        float minCost = minCostIter->second;
        int minId = minCostIter->first;
        for(; minCostIter != vertexCost.end(); ++minCostIter) {
            if(minCostIter->second < minCost) {
                minCost = minCostIter->second;
                minId = minCostIter->first;
                deleteIter = minCostIter;
            }
        }
        collapse(minId);
        vertexCost.erase(deleteIter);
    }
}

double ShortestEdgeMethod::calculateEdgeLengthSqur(Vertex v0, Vertex v1) {
    return (v0.x - v1.x) * (v0.x - v1.x)
           + (v0.y - v1.y) * (v0.y - v1.y)
           + (v0.z - v1.z) * (v0.z - v1.z);
}