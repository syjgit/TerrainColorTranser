#include <stb\stb_image.h>

#include "CController.h"
#include <iostream>
#include <Gdiplus.h>

#include "ColorSapceTransfer.h"
#include "../ColorTransform/colorHistogram.h"

using namespace std;
#include <cmath>
#include <vector>
#include <fstream>
#pragma comment(lib, "gdiplus.lib")
using namespace Gdiplus;
//---------------------------- Render ------------------------------------
//
//------------------------------------------------------------------------
void CController::Render(HDC surface)
{
	m_pSOM->Render(surface);
}

void CController::Render(HDC surface, Csom* SOM)
{
	SOM->Render(surface);
}


std::vector<int> CController::Poi2Color(int x, int y)
{
	return m_pSOM->Poi2Color(x, y);
}

std::vector<int> CController::Color2Poi(std::vector<int>& color)
{
	return m_pSOM->Color2Poi(color);

}




std::vector<ctColor> CController::getColor()
{
	////概率值计算
	std::vector<ctColor> color = m_pSOM->getColors();
	
	m_pSOM->setColors(color);

	ofstream out("../output/somColorP.csv");

	//std::sort(color.begin(), color.end(), [ ](const ctColor& x,const ctColor& y) {
	//	return x.proportion > y.proportion;
	//});

	for (const auto& v : color)
	{
		out << v.r << "," << v.g << "," << v.b << "," << v.proportion << std::endl;
	}
	out.close();

	return color;
}


std::vector<ctColor> CController::getColorLine()
{
	////概率值计算
	std::vector<ctColor> color = m_pSOMLine->getColors();
	for (int i = 0; i < color.size(); i++)
	{
		color[i].LineLength = color[i].proportion;
	}

	m_pSOMLine->setColors(color);

	ofstream out("../output/somColorLine.csv");

	//std::sort(color.begin(), color.end(), [ ](const ctColor& x,const ctColor& y) {
	//	return x.proportion > y.proportion;
	//});

	for (const auto& v : color)
	{
		out << v.r << "," << v.g << "," << v.b << "," << v.LineLength << std::endl;
	}
	out.close();

	return color;
}


std::vector<ctColor> CController::getColorBack()
{
	ofstream out("../output/somColorBack.csv");
	for (const auto& v : m_ctInutPictureBack)
	{
		out << v.r << "," << v.g << "," << v.b << "," << v.proportion << std::endl;
	}
	out.close();

	return m_ctInutPictureBack;
}


std::vector<ctColor> CController::getColorPoint()
{
	////概率值计算
	std::vector<ctColor> color = m_pSOMPoint->getColors();
	

	m_pSOMPoint->setColors(color);

	ofstream out("../output/somColorPoint.csv");

	//std::sort(color.begin(), color.end(), [ ](const ctColor& x,const ctColor& y) {
	//	return x.proportion > y.proportion;
	//});

	for (const auto& v : color)
	{
		out << v.r << "," << v.g << "," << v.b << "," << v.proportion << std::endl;
	}
	out.close();

	return color;
}


void CController::matchingColor(Csom* m_pSOM,const vector<int>& matchingindex, const vector<vector<double> >&TrainingSet)
{
	for (int i = 0; i < matchingindex.size(); i++)
	{
		if (matchingindex[i] > -1)
		{
			m_pSOM->m_SOM[i].selected = true;

			m_pSOM->m_SOM[i].m_dWeights[0] = TrainingSet[matchingindex[i]][0];
			m_pSOM->m_SOM[i].m_dWeights[1] = TrainingSet[matchingindex[i]][1];
			m_pSOM->m_SOM[i].m_dWeights[2] = TrainingSet[matchingindex[i]][2];
			m_pSOM->m_SOM[i].m_dWeights[3] = TrainingSet[matchingindex[i]][3];
		}
		else
		{
			m_pSOM->m_SOM[i].m_dWeights[0] = 0.0;
			m_pSOM->m_SOM[i].m_dWeights[1] = 0.0;
			m_pSOM->m_SOM[i].m_dWeights[2] = 0.0;
		}
	}
}


//--------------------------- Train --------------------------------------
//
// trains the network given a std::vector of input vectors
//------------------------------------------------------------------------
bool CController::Train()
{
	
	if (!m_pSOM->FinishedTraining())
	{
		if (!m_pSOM->Epoch(m_TrainingSet, ijs))
		{
			vector<int> indexs;
			vector<int> indexs2;

			matchingindex.resize(m_pSOM->numCellsAcross*m_pSOM->numCellsUp);
			memset(matchingindex.data(), -1, m_pSOM->numCellsAcross*m_pSOM->numCellsUp *4);

			for (int i = ijs.size() - 1; i >= 0; i--)
			{
				if (find(indexs.begin(), indexs.end(), ijs[i][0]) == indexs.end())
				{
					if (find(indexs2.begin(), indexs2.end(), ijs[i][1]) == indexs2.end())
					{		
						matchingindex[ijs[i][0]] = ijs[i][1];
						indexs.push_back(ijs[i][0]);
						indexs2.push_back(ijs[i][1]);
					}
				}

			}
			vector<vector<int>>().swap(ijs);
			vector<int>().swap(indexs);
			vector<int>().swap(indexs2);


			/*ofstream infile("color-matching.txt");
			for (int i = 0; i<ijs.size(); i++)
			{
				infile << ijs[i][0] << "," << ijs[i][1] << endl;
			}
			infile.close();*/

			ProximityMatch(matchingindex,m_TrainingSet);
			matchingColor(m_pSOM, matchingindex,m_TrainingSet);

			return false;
		}
	}

	return true;
}


bool CController::Train2()
{
	
	if (!m_pSOMLine->FinishedTraining())
	{
		if (!m_pSOMLine->Epoch(m_TrainingSetLine, ijs))
		{
			vector<int> indexs;
			vector<int> indexs2;

			matchingindex.resize(m_pSOMLine->numCellsUp*m_pSOMLine->numCellsAcross);
			memset(matchingindex.data(), -1, m_pSOMLine->numCellsUp*m_pSOMLine->numCellsAcross * 4);

			
			for (int i = ijs.size() - 1; i >= 0; i--)
			{
				if (find(indexs.begin(), indexs.end(), ijs[i][0]) == indexs.end())
				{
					if (find(indexs2.begin(), indexs2.end(), ijs[i][1]) == indexs2.end())
					{
						matchingindex[ijs[i][0]] = ijs[i][1];
						indexs.push_back(ijs[i][0]);
						indexs2.push_back(ijs[i][1]);
					}
				}
			}
			
			vector<vector<int>>().swap(ijs);
			vector<int>().swap(indexs);
			vector<int>().swap(indexs2);

			ProximityMatch(matchingindex, m_TrainingSetLine);

			/*ofstream infile("color-matchingLine.txt");
			for(int i=0;i<ijs.size();i++)
			{
				infile << ijs[i][0] << "," << ijs[i][1] << endl;
			}
			infile.close();*/
			matchingColor(m_pSOMLine,matchingindex,m_TrainingSetLine);

			return false;
		}
	}

	return true;
}

bool CController::Train3()
{

	if (!m_pSOMPoint->FinishedTraining())
	{
		if (!m_pSOMPoint->Epoch(m_TrainingSetPoint, ijs))
		{
			vector<int> indexs;
			vector<int> indexs2;

			matchingindex.resize(m_pSOMPoint->numCellsUp*m_pSOMPoint->numCellsAcross);
			memset(matchingindex.data(), -1, m_pSOMPoint->numCellsUp*m_pSOMPoint->numCellsAcross * 4);


			for (int i = ijs.size() - 1; i >= 0; i--)
			{
				if (find(indexs.begin(), indexs.end(), ijs[i][0]) == indexs.end())
				{
					if (find(indexs2.begin(), indexs2.end(), ijs[i][1]) == indexs2.end())
					{
						matchingindex[ijs[i][0]] = ijs[i][1];
						indexs.push_back(ijs[i][0]);
						indexs2.push_back(ijs[i][1]);
					}
				}
			}

			vector<vector<int>>().swap(ijs);
			vector<int>().swap(indexs);
			vector<int>().swap(indexs2);

			ProximityMatch(matchingindex, m_TrainingSetPoint);

			/*ofstream infile("color-matchingLine.txt");
			for (int i = 0; i<ijs.size(); i++)
			{
				infile << ijs[i][0] << "," << ijs[i][1] << endl;
			}
			infile.close();*/
			matchingColor(m_pSOMPoint, matchingindex, m_TrainingSetPoint);

			return false;
		}
	}

	return true;
}



bool CController::Train4()
{

	if (!m_pSOMBack->FinishedTraining())
	{
		if (!m_pSOMBack->Epoch(m_TrainingSetBack, ijs))
		{
			vector<int> indexs;
			vector<int> indexs2;

			matchingindex.resize(m_pSOMBack->numCellsUp*m_pSOMBack->numCellsAcross);
			memset(matchingindex.data(), -1, m_pSOMBack->numCellsUp*m_pSOMBack->numCellsAcross * 4);


			for (int i = ijs.size() - 1; i >= 0; i--)
			{
				if (find(indexs.begin(), indexs.end(), ijs[i][0]) == indexs.end())
				{
					if (find(indexs2.begin(), indexs2.end(), ijs[i][1]) == indexs2.end())
					{
						matchingindex[ijs[i][0]] = ijs[i][1];
						indexs.push_back(ijs[i][0]);
						indexs2.push_back(ijs[i][1]);
					}
				}
			}

			vector<vector<int>>().swap(ijs);
			vector<int>().swap(indexs);
			vector<int>().swap(indexs2);

			ProximityMatch(matchingindex, m_TrainingSetBack);

			matchingColor(m_pSOMBack, matchingindex, m_TrainingSetBack);

			return false;
		}
	}

	return true;
}


float colordistance(tuple<double, double, double> c1, tuple<double, double, double>c2)
{
	return sqrt(pow(get<0>(c1) - get<0>(c2), 2) + pow(get<1>(c1) - get<1>(c2), 2) + pow(get<2>(c1) - get<2>(c2), 2));

}

void CController::ProximityMatch(vector<int>& matchingindex,const vector<vector<double> >&TrainingSet)
{
	int plength = TrainingSet.size();
	vector<tuple<double, double, double> > imgLab(plength);
	vector<int> lastindex;
	for (int j = 0; j <plength; j++)
	{
		imgLab[j] = ColorSpaceTransfer::RGB2Lab1(TrainingSet[j][0], TrainingSet[j][1], TrainingSet[j][2]);

		if (find(matchingindex.begin(), matchingindex.end(), j) == matchingindex.end())
			lastindex.emplace_back(j);
	}


	vector<tuple<double, double, double>> SOMvalue;
	tuple<double, double, double> lab;
	for (auto& weight : m_pSOM->m_SOM)
	{
		lab = ColorSpaceTransfer::RGB2Lab1(weight.m_dWeights[0], weight.m_dWeights[1], weight.m_dWeights[2]);
		SOMvalue.emplace_back(lab);
	}

	int length = matchingindex.size();
	for (int i = 0; i < length; i++)
	{
		double mindis = 999999999999999;

		if (matchingindex[i] == -1)
		{
			vector<int>::iterator lindex;
			tuple<double, double, double> &Slab = SOMvalue[i];
			int index;
			for (vector<int>::iterator iindex = lastindex.begin(); iindex < lastindex.end(); iindex++)
			{
				index = *iindex;
				tuple<double, double, double> &Plab = imgLab[index];
				float dis = colordistance(Slab, Plab);
				if (dis < mindis)
				{
					mindis = dis;
					lindex = iindex;
				}
			}
			index = *lindex;
			matchingindex[i] = *lindex;
			lastindex.erase(lindex);
		}

	}
	vector<tuple<double, double, double> >().swap(SOMvalue);
	vector<tuple<double, double, double> >().swap(imgLab);
}

void CController::getPictureLine(const std::string& linepath, int **&picture)
{
	cv::Mat linemat = cv::imread(linepath);
	int heigth = linemat.cols;
	int width =  linemat.rows;
	picture = new int*[width];
	for (size_t i = 0; i < width; i++)
	{
		picture[i]= new int[heigth];
		for (size_t j = 0; j < heigth; j++)
		{
			picture[i][j] = linemat.at<uchar>(i, j);
		}
	}
}

std::vector<ctColor> CController::calPersonalPictureProportion(const std::string& img_path)
{
	std::vector<ctColor> iMapColor;
	colorHist::compute_colorHist(img_path, iMapColor, 50);


	std::sort(iMapColor.begin(), iMapColor.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});

	std::vector<ctColor> ctMapColor;
	int index = 0, pindex = 0;
	for (auto& v : iMapColor)
	{
		//计算色差
		index = 0; 
		double minD = 999999.0;
		std::vector<double> v_lab = ColorSpaceTransfer::RGB2Lab(v.r, v.g, v.b);
		for (size_t i=0;i<ctMapColor.size();i++)
		{
			std::vector<double> ct_lab = ColorSpaceTransfer::RGB2Lab(ctMapColor[i].r, ctMapColor[i].g, ctMapColor[i].b);
			
			double color_dis = sqrt(pow((ct_lab[0] - v_lab[0]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0));
			if (color_dis < m_dColorDis)
			{
				if (color_dis < minD)
				{
					minD = color_dis;
					pindex = i;
				}	
				index++;
			}
		}
		if (index > 0)
		{
			ctMapColor[pindex].proportion += v.proportion;
			continue;
		}
			
		ctMapColor.emplace_back(v);
		if (ctMapColor.size() > 150)
		{
			message = "色差过小";
			return{};
		}

	}


	if (ctMapColor.size() < 100)
	{
		message = "色差过大";
		return{};
	}
	ofstream mapp("inputP.csv");
	for (auto& v : ctMapColor)
	{
		mapp << v.r << "," << v.g << "," << v.b << "," << v.proportion << endl;
	}
	mapp.close();

	std::vector<ctColor>().swap(iMapColor);
	std::sort(ctMapColor.begin(), ctMapColor.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});


	return ctMapColor;

}

std::vector<ctColor> CController::calPersonalPictureProportion(const std::string& img_path, const std::string& linepath, std::vector<ctColor>& ctMapColorLine)
{
	int limitedV = 39;
	std::vector<ctColor> iMapColor;
	std::vector<ctColor> iMapColorL;
	colorHist::compute_colorHist(img_path, linepath, iMapColor, iMapColorL, 50, limitedV);


	std::vector<ctColor> iMapColor_all;
	colorHist::compute_colorHist(img_path, iMapColor_all, 50);

	std::sort(iMapColor.begin(), iMapColor.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});

	std::sort(iMapColorL.begin(), iMapColorL.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});

	std::sort(iMapColor_all.begin(), iMapColor_all.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});

	std::vector<ctColor> ctMapColor;

	//ofstream colorCount("colorCount.txt", ios_base::app);;
	int index = 0, pindex = 0;
	for (auto& v : iMapColor)
	{
		//计算色差
		index = 0;
		double minD = 999999.0;
		std::vector<double> v_lab = ColorSpaceTransfer::RGB2Lab(v.r, v.g, v.b);
		for (size_t i = 0; i<ctMapColor.size(); i++)
		{
			std::vector<double> ct_lab = ColorSpaceTransfer::RGB2Lab(ctMapColor[i].r, ctMapColor[i].g, ctMapColor[i].b);

			double color_dis = sqrt(pow((ct_lab[0] - v_lab[0]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0));
			if (color_dis < m_dColorDis)
			{
				if (color_dis < minD)
				{
					minD = color_dis;
					pindex = i;
				}
				index++;
			}
		}
		if (index > 0)
		{
			ctMapColor[pindex].proportion += v.proportion;
			continue;
		}
		ctMapColor.emplace_back(v);
		if (ctMapColor.size() > constNumCellsAcross*constNumCellsDown*2.5)
		{
			message = "色差过小";
			//colorCount << m_dColorDis<<","<< ctMapColor.size() << "," << ctMapColorLine.size() << endl;
			return{};
		}
	}

	if (ctMapColor.size() < constNumCellsAcross*constNumCellsDown)
	{
		message = "色差过大";
		//colorCount << m_dColorDis << "," << ctMapColor.size() << "," << ctMapColorLine.size() << endl;
		return{};
	}

	index = 0; pindex = 0;
	for (auto& v : iMapColorL)
	{
		//计算色差
		index = 0;
		double minD = 999999.0;
		std::vector<double> v_lab = ColorSpaceTransfer::RGB2Lab(v.r, v.g, v.b);
		for (size_t i = 0; i<ctMapColorLine.size(); i++)
		{
			std::vector<double> ct_lab = ColorSpaceTransfer::RGB2Lab(ctMapColorLine[i].r, ctMapColorLine[i].g, ctMapColorLine[i].b);

			double color_dis = sqrt(pow((ct_lab[0] - v_lab[0]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0));
			if (color_dis < m_dColorDis)
			{
				if (color_dis < minD)
				{
					minD = color_dis;
					pindex = i;
				}
				index++;
			}
		}
		if (index > 0)
		{
			ctMapColorLine[pindex].proportion += v.proportion;
			continue;
		}

		ctMapColorLine.emplace_back(v);
		if (ctMapColorLine.size() > 2.0*constNumCellsAcrossLine*constNumCellsDownLine)
		{
			message = "色差过小";
			//colorCount << m_dColorDis << "," << ctMapColor.size() << "," << ctMapColorLine.size() << endl;
			return{};
		}
	}
	
	if (ctMapColor.size() < constNumCellsAcross*constNumCellsDown || ctMapColorLine.size()<constNumCellsAcrossLine*constNumCellsDownLine)
	{
		message = "色差过大";
		//colorCount << m_dColorDis << "," << ctMapColor.size() << "," << ctMapColorLine.size() << endl;
		return{};
	}


	std::vector<ctColor>().swap(iMapColor);
	std::vector<ctColor>().swap(iMapColorL);

	std::sort(ctMapColor.begin(), ctMapColor.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});
	std::sort(ctMapColorLine.begin(), ctMapColorLine.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});

	ofstream mapp("inputP.csv");
	for (auto& v : ctMapColor)
	{
		mapp << v.r << "," << v.g << "," << v.b << "," << v.proportion << endl;
	}
	mapp.close();
	ofstream mapL("inputL.csv");
	for (int i=0;i<ctMapColorLine.size();i++)
	{
		mapL << ctMapColorLine[i].r << "," << ctMapColorLine[i].g << "," << ctMapColorLine[i].b << "," << ctMapColorLine[i].proportion << endl;
	}
	mapL.close();

	return ctMapColor;
}

std::vector<ctColor> CController::calPersonalPictureProportion(const std::string& img_path, const std::string& linepath,
	const std::string& pointpath, const std::string& backpath, std::vector<ctColor>& ctMapColorLine, std::vector<ctColor>& ctMapColorPoint,std::vector<ctColor>& ctMapColorBack)
{
	int limitedV = 39;
	std::vector<ctColor> iMapColor;
	std::vector<ctColor> iMapColorL;
	std::vector<ctColor> iMapColorPoint;
	std::vector<ctColor> iMapColorBack;
	colorHist::compute_colorHist(img_path, linepath, pointpath, backpath, iMapColor, iMapColorL, iMapColorPoint, iMapColorBack,30, limitedV,1);


	std::sort(iMapColorPoint.begin(), iMapColorPoint.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});

	std::sort(iMapColor.begin(), iMapColor.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});

	std::sort(iMapColorL.begin(), iMapColorL.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});

	

	std::vector<ctColor> ctMapColor;

	ofstream colorCount("colorCount.txt", ios_base::app);;
	int index = 0, pindex = 0;

	//计算背景比例
	ctMapColorBack = iMapColorBack;


	//计算面
	for (auto& v : iMapColor)
	{
		//计算色差
		index = 0;
		double minD = 999999.0;
		std::vector<double> v_lab = ColorSpaceTransfer::RGB2Lab(v.r, v.g, v.b);
		for (size_t i = 0; i<ctMapColor.size(); i++)
		{
			std::vector<double> ct_lab = ColorSpaceTransfer::RGB2Lab(ctMapColor[i].r, ctMapColor[i].g, ctMapColor[i].b);

			double color_dis = sqrt(pow((ct_lab[0] - v_lab[0]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0));
			if (color_dis < m_dColorDis)
			{
				if (color_dis < minD)
				{
					minD = color_dis;
					pindex = i;
				}
				index++;
			}
		}
		if (index > 0)
		{
			ctMapColor[pindex].proportion += v.proportion;
			continue;
		}
		ctMapColor.emplace_back(v);
		if (ctMapColor.size() > constNumCellsAcross*constNumCellsDown*2.5)
		{
			message = "色差过小";
			colorCount << m_dColorDis << "," << ctMapColor.size() << "," << ctMapColorLine.size() << "," << ctMapColorPoint.size() << endl;
			return{};
		}
	}

	if (ctMapColor.size() < constNumCellsAcross*constNumCellsDown)
	{
		message = "色差过大";
		colorCount << m_dColorDis << "," << ctMapColor.size() << "," << ctMapColorLine.size() << "," << ctMapColorPoint.size() << endl;
		return{};
	}

	index = 0; pindex = 0;
	//计算线
	for (auto& v : iMapColorL)
	{
		//计算色差
		index = 0;
		double minD = 999999.0;
		std::vector<double> v_lab = ColorSpaceTransfer::RGB2Lab(v.r, v.g, v.b);
		for (size_t i = 0; i<ctMapColorLine.size(); i++)
		{
			std::vector<double> ct_lab = ColorSpaceTransfer::RGB2Lab(ctMapColorLine[i].r, ctMapColorLine[i].g, ctMapColorLine[i].b);

			double color_dis = sqrt(pow((ct_lab[0] - v_lab[0]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0));
			if (color_dis < m_dColorDis)
			{
				if (color_dis < minD)
				{
					minD = color_dis;
					pindex = i;
				}
				index++;
			}
		}
		if (index > 0)
		{
			ctMapColorLine[pindex].proportion += v.proportion;
			continue;
		}

		ctMapColorLine.emplace_back(v);
		if (ctMapColorLine.size() > 2.0*constNumCellsAcrossLine*constNumCellsDownLine)
		{
			message = "色差过小";
			colorCount << m_dColorDis << "," << ctMapColor.size() << "," << ctMapColorLine.size() << "," << ctMapColorPoint.size() << endl;
			return{};
		}
	}

	index = 0; pindex = 0;
	//计算点
	for (auto& v : iMapColorPoint)
	{
		//计算色差
		index = 0;
		double minD = 999999.0;
		std::vector<double> v_lab = ColorSpaceTransfer::RGB2Lab(v.r, v.g, v.b);
		for (size_t i = 0; i<ctMapColorPoint.size(); i++)
		{
			std::vector<double> ct_lab = ColorSpaceTransfer::RGB2Lab(ctMapColorPoint[i].r, ctMapColorPoint[i].g, ctMapColorPoint[i].b);

			double color_dis = sqrt(pow((ct_lab[0] - v_lab[0]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0));
			if (color_dis < m_dColorDis)
			{
				if (color_dis < minD)
				{
					minD = color_dis;
					pindex = i;
				}
				index++;
			}
		}
		if (index > 0)
		{
			ctMapColorPoint[pindex].proportion += v.proportion;
			continue;
		}

		ctMapColorPoint.emplace_back(v);
		if (ctMapColorPoint.size() > 2.0*constNumCellsAcrossPoint*constNumCellsDownPoint)
		{
			message = "色差过小";
			colorCount << m_dColorDis << "," << ctMapColor.size() << "," << ctMapColorLine.size() << "," << ctMapColorPoint.size() << endl;
			return{};
		}
	}



	if (ctMapColor.size() < constNumCellsAcross*constNumCellsDown || ctMapColorLine.size()<constNumCellsAcrossLine*constNumCellsDownLine || ctMapColorPoint.size()<constNumCellsAcrossPoint*constNumCellsDownPoint)
	{
		message = "色差过大";
		colorCount << m_dColorDis << "," << ctMapColor.size() << "," << ctMapColorLine.size() << "," << ctMapColorPoint.size() << endl;
		return{};
	}


	std::vector<ctColor>().swap(iMapColor);
	std::vector<ctColor>().swap(iMapColorL);
	std::vector<ctColor>().swap(iMapColorPoint);

	std::sort(ctMapColor.begin(), ctMapColor.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});
	std::sort(ctMapColorLine.begin(), ctMapColorLine.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});
	std::sort(ctMapColorPoint.begin(), ctMapColorPoint.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});

	/*ofstream mapp("inputP.csv");
	for (auto& v : ctMapColor)
	{
		mapp << v.r << "," << v.g << "," << v.b << "," << v.proportion << endl;
	}
	mapp.close();
	ofstream mapL("inputL.csv");
	for (int i = 0; i<ctMapColorLine.size(); i++)
	{
		mapL << ctMapColorLine[i].r << "," << ctMapColorLine[i].g << "," << ctMapColorLine[i].b << "," << ctMapColorLine[i].proportion << endl;
	}
	mapL.close();
	ofstream mapPoint("inputPoint.csv");
	for (int i = 0; i<ctMapColorPoint.size(); i++)
	{
		mapL << ctMapColorPoint[i].r << "," << ctMapColorPoint[i].g << "," << ctMapColorPoint[i].b << "," << ctMapColorPoint[i].proportion << endl;
	}
	mapPoint.close();*/

	return ctMapColor;
}


void CController::colorAggregation(std::vector<ctColor>& iMapColor,const int& mlimitation)
{
	vector<ctColor>::iterator c1;
	vector<ctColor>::iterator c2;
	for (c1 = iMapColor.begin(); c1 != iMapColor.end(); c1++)
	{
		std::vector<double> v_lab = ColorSpaceTransfer::RGB2Lab((*c1).r, (*c1).g, (*c1).b);
		for (c2 = c1 + 1; c2 != iMapColor.end(); c2++)
		{
			std::vector<double> ct_lab = ColorSpaceTransfer::RGB2Lab((*c2).r, (*c2).g, (*c2).b);

			double color_dis = sqrt(pow((ct_lab[0] - v_lab[0]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0) + pow((ct_lab[1] - v_lab[1]), 2.0));
			if (color_dis < mlimitation)
			{
				double p0 = (*c1).proportion;
				(*c1).proportion += (*c2).proportion;
				(*c2).proportion += p0;
			}
		}
	}
}

std::vector<ctColor> CController::calPersonalPictureProportion(const std::vector<ctColor>& iMapColor, const std::vector<ctColor>&  iMapColorL,
	const std::vector<ctColor>&  iMapColorPoint, const std::vector<ctColor>&  iMapColorBack, std::vector<ctColor>& ctMapColorLine, std::vector<ctColor>& ctMapColorPoint, std::vector<ctColor>& ctMapColorBack)
{
	std::vector<ctColor> ctMapColor;

	//ofstream colorCount("colorCount.txt", ios_base::app);;
	int index = 0, pindex = 0;

	//计算背景比例
	ctMapColorBack = iMapColorBack;

	ctMapColorPoint = iMapColorPoint;
	ctMapColorLine = iMapColorL;
	ctMapColor = iMapColor;

	/*double limition = 5;
	colorAggregation(ctMapColor, limition);
	colorAggregation(ctMapColorPoint, limition);
	colorAggregation(ctMapColorLine, limition);
	colorAggregation(ctMapColorBack, limition);
	*/

	std::sort(ctMapColor.begin(), ctMapColor.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});
	std::sort(ctMapColorLine.begin(), ctMapColorLine.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});
	std::sort(ctMapColorPoint.begin(), ctMapColorPoint.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});
	std::sort(ctMapColorBack.begin(), ctMapColorBack.end(), [](const ctColor& x, const ctColor& y) {
		return x.proportion > y.proportion;
	});

	/*ofstream mapp("inputP.csv");
	for (auto& v : ctMapColor)
	{
		mapp << v.r << "," << v.g << "," << v.b << "," << v.proportion << endl;
	}
	mapp.close();
	ofstream mapL("inputL.csv");
	for (int i = 0; i<ctMapColorLine.size(); i++)
	{
		mapL << ctMapColorLine[i].r << "," << ctMapColorLine[i].g << "," << ctMapColorLine[i].b << "," << ctMapColorLine[i].proportion << endl;
	}
	mapL.close();
	ofstream mapPoint("inputPoint.csv");
	for (int i = 0; i<ctMapColorPoint.size(); i++)
	{
		mapPoint << ctMapColorPoint[i].r << "," << ctMapColorPoint[i].g << "," << ctMapColorPoint[i].b << "," << ctMapColorPoint[i].proportion << endl;
	}
	mapPoint.close();*/

	return ctMapColor;
}


std::vector<ctColor> CController::calPersonalPictureProportion(std::vector<int>& ctMapPixels, const int& w, const int& h, int comp)
{

	comp = 3;
	int ctMapColorSize = w*h;

	//统计颜色数量，只适用于三通道颜色
	std::vector<ctColor> ctMapColor;
	std::vector<vector<float>> MapcolorP;
	vector<float> color(4);
	for (int i = 0; i < ctMapColorSize * comp; i += comp)
	{
		color[0] = ctMapPixels[i + 0];
		color[1] = ctMapPixels[i + 1];
		color[2] = ctMapPixels[i + 2];
		color[3] = 1;
		for (auto& v : MapcolorP)
		{
			if (v[0] == color[0] && v[1] == color[1] && v[2] == color[2])
			{
				v[3]++;
				continue;
			}
		}
		MapcolorP.emplace_back(color);
	}


	for (int i = 0; i < ctMapColorSize * comp; i += comp)
	{
		int r = ctMapPixels[i + 0];
		int g = ctMapPixels[i + 1];
		int b = ctMapPixels[i + 2];

		bool bFind = false;
		for (auto& v : ctMapColor)
		{
			//计算色差
			double color_dis = sqrt(pow((r - v.r), 2.0) + pow((g - v.g), 2.0) + pow((b - v.b), 2.0));
			if (color_dis < m_dColorDis)
			{
				bFind = true;
				v.count++;

				break;
			}

		}

		if (!bFind)
		{
			ctColor t;
			t.r = r;
			t.g = g;
			t.b = b;
			ctMapColor.emplace_back(t);

			if (ctMapColor.size() > 200)
			{
				message = "色差过小";
				return{};
			}

		}

	}

	//计算占比
	for (auto& v : ctMapColor)
	{
		v.proportion = (double)v.count / (double)ctMapColorSize;
	}

	if (ctMapColor.size() < 100)
	{
		message = "色差过大";
		return{};
	}

	return ctMapColor;
}




//-------------------------- CreateDataSet -------------------------------
//
//------------------------------------------------------------------------
void CController::CreateDataSet(const std::string& img_path)
{

#ifndef RANDOM_TRAINING_SETS

	m_ctInutPicture = calPersonalPictureProportion(img_path);

	if (m_ctInutPicture.size() < 1 || message != "")
	{
		return;
	}

	std::sort(m_ctInutPicture.begin(), m_ctInutPicture.end(), [](const ctColor& x, const ctColor y)
	{
		return x.proportion > y.proportion;
	});


	
	for (const auto& v : m_ctInutPicture)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSet.emplace_back(tmp);

	}





#else

	//choose a random number of training sets
	int NumSets = RandInt(constMinNumTrainingSets, constMaxNumTrainingSets);

	for (int s = 0; s<NumSets; ++s)
	{

		vector<double> set;

		set.push_back(RandFloat());
		set.push_back(RandFloat());
		set.push_back(RandFloat());

		m_TrainingSet.push_back(set);
	}

#endif

}

void CController::CreateDataSet(const std::string& img_path, const std::string& linepath)
{

#ifndef RANDOM_TRAINING_SETS

	m_ctInutPicture = calPersonalPictureProportion(img_path, linepath, m_ctInutPictureLine);

	if (m_ctInutPicture.size() < 1 || message != "")
	{
		return;
	}


	for (const auto& v : m_ctInutPicture)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSet.emplace_back(tmp);

	}

	for (const auto& v : m_ctInutPictureLine)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSetLine.emplace_back(tmp);

	}

#else

	//choose a random number of training sets
	int NumSets = RandInt(constMinNumTrainingSets, constMaxNumTrainingSets);

	for (int s = 0; s<NumSets; ++s)
	{

		vector<double> set;

		set.push_back(RandFloat());
		set.push_back(RandFloat());
		set.push_back(RandFloat());

		m_TrainingSet.push_back(set);
	}

#endif

}

void CController::CreateDataSet(const std::string& img_path, const std::string& linepath, 
	const std::string& pointpath, const std::string& backpath)
{

#ifndef RANDOM_TRAINING_SETS

	m_ctInutPicture = calPersonalPictureProportion(img_path, linepath,pointpath, backpath, m_ctInutPictureLine, m_ctInutPicturePoint, m_ctInutPictureBack);

	if (m_ctInutPicture.size() < 1 || message != "")
	{
		return;
	}


	for (const auto& v : m_ctInutPicture)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSet.emplace_back(tmp);

	}

	for (const auto& v : m_ctInutPictureLine)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSetLine.emplace_back(tmp);

	}

	for (const auto& v : m_ctInutPicturePoint)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSetPoint.emplace_back(tmp);

	}

	for (const auto& v : m_ctInutPictureBack)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSetBack.emplace_back(tmp);

	}

#else

	//choose a random number of training sets
	int NumSets = RandInt(constMinNumTrainingSets, constMaxNumTrainingSets);

	for (int s = 0; s<NumSets; ++s)
	{

		vector<double> set;

		set.push_back(RandFloat());
		set.push_back(RandFloat());
		set.push_back(RandFloat());

		m_TrainingSet.push_back(set);
	}

#endif

}

void CController::CreateDataSet(const std::vector<ctColor>& img_path, const std::vector<ctColor>& linepath,
	const std::vector<ctColor>& pointpath, const std::vector<ctColor>& backpath)
{

#ifndef RANDOM_TRAINING_SETS

	m_ctInutPicture = calPersonalPictureProportion(img_path, linepath, pointpath, backpath, m_ctInutPictureLine, m_ctInutPicturePoint, m_ctInutPictureBack);

	if (m_ctInutPicture.size() < 1 || message != "")
	{
		return;
	}


	for (const auto& v : m_ctInutPicture)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSet.emplace_back(tmp);

	}

	for (const auto& v : m_ctInutPictureLine)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSetLine.emplace_back(tmp);

	}

	for (const auto& v : m_ctInutPicturePoint)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSetPoint.emplace_back(tmp);

	}

	for (const auto& v : m_ctInutPictureBack)
	{
		std::vector<double> tmp = { v.r / 255.0 ,v.g / 255.0 ,v.b / 255.0 ,v.proportion };
		m_TrainingSetBack.emplace_back(tmp);

	}

#else

	//choose a random number of training sets
	int NumSets = RandInt(constMinNumTrainingSets, constMaxNumTrainingSets);

	for (int s = 0; s<NumSets; ++s)
	{

		vector<double> set;

		set.push_back(RandFloat());
		set.push_back(RandFloat());
		set.push_back(RandFloat());

		m_TrainingSet.push_back(set);
	}

#endif

}