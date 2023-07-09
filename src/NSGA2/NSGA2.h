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




#define URAND (rand()/(RAND_MAX+1.0))//���������


class NSGA2
{
public:
	
	std::vector<int> temp1;//��ʱ����
	std::vector<int> temp2;//��ʱ����
	std::vector<int> mark;//�������

	int Dimension2;//����ά���������ＴZDT1����xi��i�����ֵ
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

	double wc;   //ϰ����ɫȨ��  //�����ʱ�����Ȩ��
	double wse;   //����Ȩ��
	double wsi;  //������Ȩ��
	double wha; //����Ȩ��
	double whi; //�Ӿ����Ȩ��
	double wseq; //����Ȩ��
	double wdi; //����Ȩ��

	bool IsHarmony;

	int ubP;
	int ubL;
	int ubPoint;
	int ubBack;


	int NSGA2popsize;//��Ⱥ��С
	int NSGA2generation;//���ܴ���
	int picnum;

	float delta1 = 15;   //��ɫ��ɫ�� -����
	float delta12 = 30;   //��ɫ��ɫ�� -���� ���м�
	float delta13 = 25;   //��ɫ��ɫ�� -���� ϰ����ɫ
	float delta3 = 40;   //�Ӿ����ɫ��
	float delta4 = 10;   //�Ӿ�������Ȳ�
	float delta2 = 10;   //ϰ����ɫɫ��

	float LdH = 10;   //Hֵ��������
	float LdH2 = 10;   //Hֵ��������
	float LdS = 10;
	float LdV = 8;
	float LdV2 = 10;

	std::vector<bool> IsTranfers;
	float AsscoDis = 15;
	float DiverH = 45;  //�������м�H��
	float DiverV = 10; //�������м����Ȳ�

	bool isBrightLimit = true;


	vector<vector<int>> colorRelationShips;  //��ɫ��ϵ  
	vector<bool> isColorConst;    //��ɫ�Ƿ�Ϊϰ����ɫͼ��
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
	std::vector<colorXY> value;//xi��ֵ
	vector<int> sp;// [2 * NSGA2popsize];
	//��֧����弯��SP�������ǿ��н�ռ������б�����p֧��ĸ�����ɵļ��ϡ�
	int np;
	//֧�����np���������ڿ��н�ռ��п���֧�����p�����Ը����������
	int is_dominated;//����sp�ĸ���
	void init();//��ʼ������
	int rank;//���ȼ���Pareto����Ϊ��ǰ��߼�
	double crowding_distance;//ӵ������
	double fvalue[2];//ZDT1����Ŀ�꺯����ֵ
	void f_count();//����fvalue��ֵ

	void calObjective(vector<double>& resultscore);

	double IsFeasible = 0.0;	
};




class population
{
public:

	NSGA2* mNSGA;
	individual **F; // [2 * NSGA2popsize][2 * NSGA2popsize];
	population(NSGA2*);//���ʼ��
	vector<individual> P;// [NSGA2popsize];
	vector<individual> Q;// [NSGA2popsize];
	vector<individual> R;// [2 * NSGA2popsize];
	vector<individual> tempP;
	vector<int> len;// [2 * NSGA2popsize];//�������콻����Ⱥ��Fi�ĳ��ȵļ���
	int len_f;//����Ⱥ��rankֵ
	//�������һ����ʼ����P���ڴ˻����ϲ��ö�Ԫ������ѡ��
	//����ͱ�����������Ӵ�Q��P��QȺ���ģ��ΪNSGA2popsize
	//��Pt��Qt���뵽Rt�У���ʼʱt=0������Rt���п��ٷ�֧�������
	//���������в�ͬ�ȼ��ķ�֧��⼯F1��F2........
	void set_p_q();
	int Rnum;
	int Pnum;
	int Qnum;
	//P,Q,R��Ԫ�صĸ���
	void make_new_pop();//�����µ��Ӵ�
	void make_new_pop2();//�����µ��Ӵ�
	void fast_nondominated_sort();//���ٷ�֧������
	void calu_crowding_distance(int i);//ӵ���������
	void f_sort(int i);//��ӵ�����뽵������
	void maincal();//��Ҫ����
	int choice(int a, int b);
	std::vector<ctColor> getResult(const float&index, vector<double>& resultscore);
	//�����������ڲ�ͬ�ȼ��ķ�֧��⼯�����ȿ��ǵȼ���Ž�С��
	//��������������ͬһ�ȼ��ķ�֧��⼯�����ȿ���ӵ������ϴ��


	static string doubleToString(const double &dbNum)
	{
		char *chCode;
		chCode = new(std::nothrow)char[20];
		sprintf(chCode, "%.2lf", dbNum); // .2 �ǿ���������ȵ�bai����λС��
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









