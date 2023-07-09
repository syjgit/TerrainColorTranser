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


	string linePicturepath; //��·��
	string linePicturepath2; //��·��2-��
	string pointPicturepath;//
	string backPicturepath;

	float Optimizing_wc=1;  //ϰ����ɫȨ��  //�����ʱ�����Ȩ��
	float Optimizing_wse = 1;//����Ȩ��
	float Optimizing_wsi = 1; //������Ȩ��
	float Optimizing_wha = 1;//����Ȩ��
	float Optimizing_whi = 1;//�Ӿ����Ȩ��
	float Optimizing_wseq = 1;//����Ȩ��
	float Optimizing_wdi = 1.0;

	std::vector<ctCategory> m_ctAllCategory2;  //�����ͼ��Ϣ
	std::vector<std::vector<int>> clusters;

	std::vector<bool> IsCustoms;
	std::vector<bool> IsTranfers;

	bool judgeStyleExist(const std::string &mlinePicturepath, const std::string &mlinePicturepath2,
		const std::string & mpointPicturepath, const std::string &mbackPicturepath);
	void calIamgeStyle(const std::string &msrcPicturepath, const std::string &mlinePicturepath, const std::string &mlinePicturepath2,
		const std::string & mpointPicturepath, const std::string &mbackPicturepath);
	vector<vector<double>> resultScores;   //ת�ƽ���÷� 0�ǹ��ܵ÷֣�1�����÷�
public:
	styleTransferManager();
	~styleTransferManager();
	void SetBrightLimit(double is);
	std::vector<ctColor> getBlob();
	std::vector<ctColor> getFill();
	std::vector<ctColor> getStroke();
	std::vector<ctColor> getBackground();

	std::vector<string> Imagefilepath; //����ͼƬ·��

	string mapfilepath;

	double minColor_Dis = 1; //ѵ����ɫ��С����ֵ
	double minColorBack2 = 20; //������ɫ����С����ֵ
	int m_edNSGAPopNum = 200; //��Ⱥ����
	int m_edIterationTime = 1000; //��������
	double BrightLimit = 10;
	bool isBrightLimit = true;
	int m_ResultNum = 3;// ����������
	double ValueFeaLimitation = 5;//Լ������ֵ

	float mColorDis2 = 40;//ǰ�󱳾���ɫ��С����
	float mcolorSDis = 15;//��ɫ����С����
	float mColorDis3 = 10;//ϰ����ɫ��С����

	float backgroundAdjust = 1.0f;  //��������ɫռ����ǿ����

	double Mulit_Picture = 0; //ͼƬ����

	int  StyleLearningState = 0;  //���ѧϰ���״̬1����ȡ��ɣ�2��,3�㣬4����

	std::vector<std::vector<ctColor>> resG;  //ת�ƽ��
	vector<vector<double>>  getResultScores(); //ת�ƽ���÷� 0�ǹ��ܵ÷֣�1�����÷�
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

