// qem.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Matrix.h"
#include "Quadric.h"
#include <string>
#include <iostream>
#include "ShortestEdgeMethod.h"
#define _CRTDBG_MAP_ALLOC

int main(int argc, char* argv[]) {
    _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

    if(argc != 2) {
        fprintf(stderr,
                "输入参数不正确，例子：exe名 简化算法名(QEM, Shortest distance)\n");
        return -1;
    }

    if(strcmp(argv[1], "QEM") == 0) {
        Quadrics quadrics;

        Quadrics originQuadrics;
        printf("请输入文件路径:\n");
        std::string path;
        std::cin>>path;
        if(path.empty() == true) return -1;
        quadrics.readObj(path);
        originQuadrics.readObj(path);
        fprintf(stderr, "共计%d顶点，%d面\n", quadrics.getVerticesNum(), quadrics.getFacesNum());
        quadrics.initQMatrices();
        quadrics.selectPairs();
        printf("请输入目标模型的顶点数\n");

        int targetNums;
        scanf("%d", &targetNums);
        printf("开始简化\n");
        quadrics.constructContract(targetNums);
        std::string pathOutput = path.substr(0, path.find_last_of('.')) + "_after.obj";
        quadrics.write_smf(pathOutput.c_str());
        fprintf(stderr, "简化成功\n");

        double e = quadrics.calculateSimilarity(quadrics);
        printf("相似度为%e\n", e);
    } else if(strcmp(argv[1], "Shortest distance") == 0) {
        ShortestEdgeMethod shortestEdgeMethod;
        fprintf(stdout, "请输入原始模型路径\n");
        std::string path;
        std::cin>>path;
        if(path.empty()) return -1;
        shortestEdgeMethod.readObj(path);
        fprintf(stdout, "请输入模型目标顶点数：\n");
        int targetNums;
        scanf("%d", &targetNums);
        fprintf(stdout, "开始简化:\n");
        shortestEdgeMethod.simplifyMesh(targetNums);
        std::string pathOutput = path.substr(0, path.find_last_of('.')) + "_after.obj";
        shortestEdgeMethod.write_smf(pathOutput.c_str());
    }
    return 0;
}

