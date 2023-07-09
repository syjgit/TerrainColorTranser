#pragma once

#include <algorithm>
#include <vector>

#include "../ColorTransform/constants.h"

#include "../ColorTransform/ColorSapceTransfer.h"

#include <vector> 
#include <iostream>
#include <iostream>                                                                                                                                
#include <fstream>
#include <math.h>
#include <algorithm>
#include<omp.h>

#include <fstream>
#include "../base.h"

using namespace std;




#define URAND (rand()/(RAND_MAX+1.0))//产生随机数


class NSGA2
{
public:
	
	std::vector<int> temp1;//临时数组
	std::vector<int> temp2;//临时数组
	std::vector<int> mark;//标记数组

	int Dimension2;//基因维数，在这里即ZDT1问题xi的i的最大值
	int FeaLimitation;

	vector<int>type2;
	std::vector<ctColor> m_colorRGBDataP;
	std::vector<ctColor> m_colorRGBDataLine;
	std::vector<ctColor> m_colorRGBDataPoint;
	std::vector<ctColor> m_colorRGBDataBack;
	vector<ctCategory> m_dPictureInformation;

	float proportionAll;
	double m_Mapproportion2;

	vector<int> constLayer;
	vector<vector<int>> Clusters;
	std::vector<vector<double>> sourceLabs;;
	int backLayerID;

	double wc;   //习惯用色权重  //面积大时，提高权重
	double wse;   //语义权重
	double wsi;  //相似性权重
	double wha; //调和权重
	double whi; //视觉层次权重
	double wseq; //序列权重
	double wdi; //差异权重

	bool IsHarmony;

	int ubP;
	int ubL;
	int ubPoint;
	int ubBack;


	int NSGA2popsize;//种群大小
	int NSGA2generation;//繁衍代数
	int picnum;

	float delta1 = 15;   //颜色间色差 -语义
	float delta12 = 30;   //颜色间色差 -语义 序列间
	float delta13 = 25;   //颜色间色差 -语义 习惯用色
	float delta3 = 40;   //视觉层次色差
	float delta4 = 10;   //视觉层次亮度差
	float delta2 = 10;   //习惯用色色差

	float LdH = 10;   //H值基本不变
	float LdH2 = 10;   //H值基本不变
	float LdS = 10;
	float LdV = 8;
	float LdV2 = 10;

	std::vector<bool> IsTranfers;
	float AsscoDis = 15;
	float DiverH = 45;  //两端型中间H差
	float DiverV = 10; //两端型中间亮度差

	bool isBrightLimit = true;


	vector<vector<int>> colorRelationShips;  //颜色关系  
	vector<bool> isColorConst;    //颜色是否为习惯颜色图层
	int NSGA2::getColorRelationship(const int& index1, const int &index2);
	void iniParameters(const std::vector<ctColor>& colorP, const std::vector<ctColor>& colorL, const std::vector<ctColor>& colorPoint, const std::vector<ctColor>& colorBack,
		const std::vector<ctCategory>& InformationEntropy, const vector<vector<int>> &clusters, const vector<bool> &IsCustoms, const vector<bool> &IsTranfers);


	void iniColorDis(const float& mcolorSDis, const float& mColorDis2, const float& mColorDis3);
	static bool IsConstLayer(const string &name);
	int getdiffub(const int & index);
	std::vector<ctColor>&  getdiffcolors(const int&index);

	void setOptimizingWeights(const float& w1, const float& w2, const float& w3, const float& w4, const float& w5, const float& w6, const float& w7);

	double rand_real(double low, double high);

	int rand_int(int low, int high);


	int getGroupID(int index);

	

	//static ofstream colorF1;
};

class individual
{
private:
	struct colorXY
	{
		int x;
		int y;
	};

	vector<float> GroupdHs;

	vector<int> GenerateDiffNumber(int min, int max, int num);
	void getManifoldCenter(const int &type, const std::vector<ctColor> &colors, vector<int>&);
	void calColorDisNet();
	bool IsInconstLayer(const int& index1);
	float calSequenceLoss(const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab);
	bool IsInSequence(const int& index1, const int& index2, const vector<vector<int>>& clusters);
	bool IsInDiverging(const int& index1, const int& index2, vector<vector<int>>& clusters);
	
	double getColorSemanticObjective(const int& i, const int& j, const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab);
	//double getColorSemanticObjective2(const int& i, const int& j, const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab);
	double getColorConObjective(const int& index, const vector<vector<double>> &colorsLab);

	void getOrderGroupsDH(const vector<vector<double>> &colorsHSV, vector<float> & colorsdH);

	void getFunction2(const vector<ctColor>& colors, vector<double> &f2);
	void getFunction1_2(const vector<ctColor>& colors, const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab
		, vector<double>& f1);

	double getFunction2(const vector<ctColor>& colors);
	float calHarmony(const vector<ctColor>& colors);
	float calHarmony2(const vector<ctColor>& colors);
	double getFunction1(const vector<ctColor>& colors, const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab);
	double getFunction1_2(const vector<ctColor>& colors, const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab);
	double getFunction3(const vector<vector<double>> &colorsLab);
public:
	NSGA2* mNSGA;
	std::vector<colorXY> value;//xi的值
	vector<int> sp;// [2 * NSGA2popsize];
	//被支配个体集合SP。该量是可行解空间中所有被个体p支配的个体组成的集合。
	int np;
	//支配个数np。该量是在可行解空间中可以支配个体p的所以个体的数量。
	int is_dominated;//集合sp的个数
	void init();//初始化个体
	int rank;//优先级，Pareto级别为当前最高级
	double crowding_distance;//拥挤距离
	double fvalue[2];//ZDT1问题目标函数的值
	void f_count();//计算fvalue的值

	void calObjective(vector<double>& resultscore);

	double IsFeasible = 0.0;	
};




class population
{
public:

	NSGA2* mNSGA;
	individual **F; // [2 * NSGA2popsize][2 * NSGA2popsize];
	population(NSGA2*);//类初始化
	vector<individual> P;// [NSGA2popsize];
	vector<individual> Q;// [NSGA2popsize];
	vector<individual> R;// [2 * NSGA2popsize];
	vector<individual> tempP;
	vector<int> len;// [2 * NSGA2popsize];//各个变异交叉后的群体Fi的长度的集合
	int len_f;//整个群体rank值
	//随机产生一个初始父代P，在此基础上采用二元锦标赛选择、
	//交叉和变异操作产生子代Q。P和Q群体规模均为NSGA2popsize
	//将Pt和Qt并入到Rt中（初始时t=0），对Rt进行快速非支配解排序，
	//构造其所有不同等级的非支配解集F1、F2........
	void set_p_q();
	int Rnum;
	int Pnum;
	int Qnum;
	//P,Q,R中元素的个数
	void make_new_pop();//产生新的子代
	void make_new_pop2();//产生新的子代
	void fast_nondominated_sort();//快速非支配排序
	void calu_crowding_distance(int i);//拥挤距离计算
	void f_sort(int i);//对拥挤距离降序排列
	void maincal();//主要操作
	int choice(int a, int b);
	std::vector<ctColor> getResult(const float&index, vector<double>& resultscore);
	//两个个体属于不同等级的非支配解集，优先考虑等级序号较小的
	//若两个个体属于同一等级的非支配解集，优先考虑拥挤距离较大的


	static string doubleToString(const double &dbNum)
	{
		char *chCode;
		chCode = new(std::nothrow)char[20];
		sprintf(chCode, "%.2lf", dbNum); // .2 是控制输出精度的bai，两位小数
		string strCode(chCode);
		delete[]chCode;
		return strCode;
	}
	bool e_is_dominated(const individual &a, const individual &b);
private:
	int getMinIndex(const int & findex);


	static int cmp1(const void *a, const void *b);

	static int cmp2(const void *a, const void *b);

	static int cmp_c_d(const void *a, const void *b);


	

};









