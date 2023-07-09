#pragma once


#include "ColorTransform/constants.h"
#include <string>
#include <vector>

using namespace std;
#include "NSGA2\NSGA2.h"
class CController;
class StyleExtract;
class MAPSTYLETRANSFORM_EXPORT styleTransferManager
{
	CController* m_pSom;
	double RWidth = 100;
	double RHeight = 100;


	string linePicturepath; //线路径
	string linePicturepath2; //线路径2-粗
	string pointPicturepath;//
	string backPicturepath;

	float Optimizing_wc=1;  //习惯用色权重  //面积大时，提高权重
	float Optimizing_wse = 1;//语义权重
	float Optimizing_wsi = 1; //相似性权重
	float Optimizing_wha = 1;//调和权重
	float Optimizing_whi = 1;//视觉层次权重
	float Optimizing_wseq = 1;//序列权重
	float Optimizing_wdi = 1.0;

	std::vector<ctCategory> m_ctAllCategory2;  //输入地图信息
	std::vector<std::vector<int>> clusters;

	std::vector<bool> IsCustoms;
	std::vector<bool> IsTranfers;

	bool judgeStyleExist(const std::string &mlinePicturepath, const std::string &mlinePicturepath2,
		const std::string & mpointPicturepath, const std::string &mbackPicturepath);
	void calIamgeStyle(const std::string &msrcPicturepath, const std::string &mlinePicturepath, const std::string &mlinePicturepath2,
		const std::string & mpointPicturepath, const std::string &mbackPicturepath);
	vector<vector<double>> resultScores;   //转移结果得分 0是功能得分，1是美得分
public:
	styleTransferManager();
	~styleTransferManager();
	void SetBrightLimit(double is);
	std::vector<ctColor> getBlob();
	std::vector<ctColor> getFill();
	std::vector<ctColor> getStroke();
	std::vector<ctColor> getBackground();

	std::vector<string> Imagefilepath; //输入图片路径

	string mapfilepath;

	double minColor_Dis = 1; //训练颜色最小分离值
	double minColorBack2 = 20; //背景颜色的最小分离值
	int m_edNSGAPopNum = 200; //种群个数
	int m_edIterationTime = 1000; //迭代次数
	double BrightLimit = 10;
	bool isBrightLimit = true;
	int m_ResultNum = 3;// 输出结果个数
	double ValueFeaLimitation = 5;//约束限制值

	float mColorDis2 = 40;//前后背景颜色最小距离
	float mcolorSDis = 15;//颜色间最小距离
	float mColorDis3 = 10;//习惯颜色最小距离

	float backgroundAdjust = 1.0f;  //背景亮颜色占比增强因子

	double Mulit_Picture = 0; //图片个数

	int  StyleLearningState = 0;  //风格学习完成状态1面提取完成，2线,3点，4背景

	std::vector<std::vector<ctColor>> resG;  //转移结果
	vector<vector<double>>  getResultScores(); //转移结果得分 0是功能得分，1是美得分
	string getResultPathByIndex(int index);
	
	StyleExtract* mStyleExtract = NULL;

	bool styleExtract();
	bool StyleRead();
	bool SampleLearning();
	bool getMapfile();
	bool StyleTransfer();

	void Output(const population& pop);
	
	void SetMinBlob(float );
	void SetMaxBlob(float );
	
	void SetlineWidth(float);

	void MapSemanticExtract(const string& mapfile, float eps=15.0f );
	void saveMapRelationshipfile(const string &mapfile, const std::vector<string> &layernames, const vector<vector<int>>& clustersAll
	, const std::vector<bool>& iscus, const std::vector<bool>& isLock);
	void saveMapfile(const string &mapfile, const std::vector<string> &layernames,
		const std::vector<int> &types, const std::vector<vector<double>>& rgbs, const std::vector<double>& pros);

	float calMapH(const vector<vector<double>>&);
};

