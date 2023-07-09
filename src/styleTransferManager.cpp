
#include "styleTransferManager.h"
#include "NSGA2/NSGA2.h"
#include "ColorHarmony/MapHarmony.h"
#include "ColorTransform/pictureOp.h"
#include "ColorTransform/constants.h"
#include "ColorTransform/colorHistogram.h"
#include "ColorTransform/CController.h"
#include "StyleExtract/StyleExtract.h"
#include "SemanticExtract\MapRelationship.h"


styleTransferManager::styleTransferManager()
{
	mStyleExtract = new StyleExtract();

}


styleTransferManager::~styleTransferManager()
{
	if (mStyleExtract)
	{
		delete mStyleExtract;
		mStyleExtract = NULL;
	}
}


void styleTransferManager::SetMinBlob(float minblob)
{
	if (mStyleExtract == NULL)
		mStyleExtract = new StyleExtract();
	mStyleExtract->minBlob = minblob;
}
void styleTransferManager::SetMaxBlob(float maxblob)
{
	if (mStyleExtract == NULL)
		mStyleExtract = new StyleExtract();
	mStyleExtract->maxBlob = maxblob;
}

void styleTransferManager::SetlineWidth(float linewidth)
{
	if (mStyleExtract == NULL)
		mStyleExtract = new StyleExtract();
	mStyleExtract->lineWidth = linewidth;
}


void styleTransferManager::MapSemanticExtract(const string& mapfile,float eps)
{
	MapRelationship mp;
	mp.setDBSCANEps(eps);
	std::vector<std::vector<double>> rgbs;
	std::vector<double> color(3);
	std::vector<string> layernames;
	std::vector<bool> isCus;
	std::vector<bool> isTra;
	std::vector<bool> isLock;
	std::ifstream inFile(mapfile, std::ios::in);//inFile来自fstream,ifstream为输入文件流(从文件读入)
	std::string lineStr;

	while (getline(inFile, lineStr)) //getline来自sstream
	{
		std::stringstream ss(lineStr);//来自sstream
		std::string str;
		std::vector<std::string> lineArray;
		//按照逗号分隔
		while (getline(ss, str, ','))
			lineArray.push_back(str);//一行数据以vector保存
									 //cout<<lineArray[0]<<endl;
		
		color[0] = atoi(lineArray[2].c_str());
		color[1] = atoi(lineArray[3].c_str());
		color[2] = atoi(lineArray[4].c_str());
		
		layernames.emplace_back(lineArray[0]);
		isCus.emplace_back(NSGA2::IsConstLayer(lineArray[0]));
		isTra.emplace_back(true);
		isLock.emplace_back(false);
		rgbs.emplace_back(color);
	}
	inFile.close();
	mp.calCluster(rgbs);

	saveMapRelationshipfile(mapfile, layernames, mp.clustersAll, isCus,isLock);
}

void styleTransferManager::saveMapRelationshipfile(const string &mapfile, const std::vector<string> &layernames,const vector<vector<int>>& clustersAll,
	const std::vector<bool>& iscus, const std::vector<bool>& isLock)
{
	string rpath;
	PictureOP::getRelationshipPath(mapfile, rpath);
	ofstream result(rpath);
	for (int i = 0; i < layernames.size(); i++)
	{
		result << layernames[i] << "," << clustersAll[i][0] << "," << iscus[i] << "," << isLock[i] << endl;
	}
	result.close();
}

void styleTransferManager::saveMapfile(const string &mapfile, const std::vector<string> &layernames, 
	const std::vector<int> &types, const std::vector<vector<double>>& rgbs, const std::vector<double>& pros)
{
	ofstream result(mapfile);
	for (int i = 0; i < layernames.size(); i++)
	{
		result << layernames[i] << "," << types[i]<< "," << rgbs[i][0] << "," << rgbs[i][1] << "," << rgbs[i][2] << "," << pros[i]  << endl;
	}
	result.close();
}


bool styleTransferManager::styleExtract()
{
	int picnum = Imagefilepath.size();
	for (int i = 0; i < picnum; i++)
	{
		PictureOP::getLinePicture(Imagefilepath[i], linePicturepath, linePicturepath2, pointPicturepath, backPicturepath);
		calIamgeStyle(Imagefilepath[i], linePicturepath, linePicturepath2, pointPicturepath, backPicturepath);
		
	}
	return true;
}

void styleTransferManager::calIamgeStyle(const std::string &msrcPicturepath, const std::string &mlinePicturepath, const std::string &mlinePicturepath2,
	const std::string & mpointPicturepath, const std::string &mbackPicturepath)
{
	if (!PictureOP::exists(msrcPicturepath))
		return;
	mStyleExtract->iniParameter(msrcPicturepath);
	mStyleExtract->extractBackground(msrcPicturepath);
	mStyleExtract->extractStroke(mlinePicturepath);
	mStyleExtract->extractEdge2(mlinePicturepath2);
	mStyleExtract->extractBlob(mpointPicturepath);
}

bool styleTransferManager::judgeStyleExist(const std::string &mlinePicturepath, const std::string &mlinePicturepath2,
	const std::string & mpointPicturepath, const std::string &mbackPicturepath)
{
	if (!PictureOP::exists(mlinePicturepath))
		return false;
	if (!PictureOP::exists(mlinePicturepath2))
		return false;
	if (!PictureOP::exists(mpointPicturepath))
		return false;
	if (!PictureOP::exists(mbackPicturepath))
		return false;
	return true;
}


bool styleTransferManager::StyleRead()
{
	string m_strAnalysePictruePath;
	
	int picnum = Imagefilepath.size();
	if (picnum > 1)
	{
		vector<ctColor>iMapColor;
		vector<ctColor>iMapColorL;
		vector<ctColor>iMapColorPoint;
		vector<ctColor>iMapColorBack;

		vector<ctColor>iMapColor0;
		vector<ctColor>iMapColorL0;
		vector<ctColor>iMapColorPoint0;
		vector<ctColor>iMapColorBack0;
		for (int i = 0; i < picnum; i++)
		{
			m_strAnalysePictruePath = Imagefilepath[i];
			PictureOP::getLinePicture(m_strAnalysePictruePath, linePicturepath, linePicturepath2, pointPicturepath, backPicturepath);
			if (!judgeStyleExist(linePicturepath, linePicturepath2, pointPicturepath, backPicturepath))
				calIamgeStyle(m_strAnalysePictruePath, linePicturepath, linePicturepath2, pointPicturepath, backPicturepath);
			int limitedV = 39;
			colorHist::colorDis = minColor_Dis;
			colorHist::compute_colorHist2(m_strAnalysePictruePath, linePicturepath, linePicturepath2, pointPicturepath, backPicturepath, iMapColor0, iMapColorL0, iMapColorPoint0, iMapColorBack0, 30, limitedV, minColorBack2);
			iMapColor.insert(iMapColor.end(), iMapColor0.begin(), iMapColor0.end());
			iMapColorL.insert(iMapColorL.end(), iMapColorL0.begin(), iMapColorL0.end());
			iMapColorPoint.insert(iMapColorPoint.end(), iMapColorPoint0.begin(), iMapColorPoint0.end());
			iMapColorBack.insert(iMapColorBack.end(), iMapColorBack0.begin(), iMapColorBack0.end());
		}

		if (backgroundAdjust != 0.0)
		{
			colorHist::colorTransfer(iMapColorBack);
			colorHist::BrightnessAdjust(iMapColorBack, backgroundAdjust);
		}

		

		m_pSom = new CController(iMapColor, iMapColorL, iMapColorPoint, iMapColorBack,
			RWidth,
			RHeight,
			sqrt(picnum),
			constNumIterations);

		Mulit_Picture = picnum;
	}
	else if (picnum == 1)
	{
		Mulit_Picture = 1;
		m_strAnalysePictruePath = Imagefilepath[0];
		PictureOP::getLinePicture(m_strAnalysePictruePath, linePicturepath, linePicturepath2, pointPicturepath, backPicturepath);
		if (!judgeStyleExist(linePicturepath, linePicturepath2, pointPicturepath, backPicturepath))
			calIamgeStyle(m_strAnalysePictruePath, linePicturepath, linePicturepath2, pointPicturepath, backPicturepath);
		int limitedV = 39;
		std::vector<ctColor> iMapColor;
		std::vector<ctColor> iMapColorL;
		std::vector<ctColor> iMapColorPoint;
		std::vector<ctColor> iMapColorBack;
		colorHist::colorDis = minColor_Dis;
		colorHist::compute_colorHist2(m_strAnalysePictruePath, linePicturepath, linePicturepath2, pointPicturepath, backPicturepath, iMapColor, iMapColorL, iMapColorPoint, iMapColorBack, 30, limitedV, minColorBack2);
		
		if (backgroundAdjust != 0.0)
		{
			colorHist::colorTransfer(iMapColorBack);
			colorHist::BrightnessAdjust(iMapColorBack, backgroundAdjust);
		}
		
		std::sort(iMapColorPoint.begin(), iMapColorPoint.end(), [](const ctColor& x, const ctColor& y) {
			return x.proportion > y.proportion;
		});

		std::sort(iMapColor.begin(), iMapColor.end(), [](const ctColor& x, const ctColor& y) {
			return x.proportion > y.proportion;
		});

		std::sort(iMapColorL.begin(), iMapColorL.end(), [](const ctColor& x, const ctColor& y) {
			return x.proportion > y.proportion;
		});


		double tetsDis = 1.0;

		if (linePicturepath != "" && pointPicturepath != "" && backPicturepath != "")
		{
			m_pSom = new CController(iMapColor, iMapColorL, iMapColorPoint, iMapColorBack,
				tetsDis, RWidth,
				RHeight,
				constNumCellsAcross,
				constNumCellsDown,
				constNumIterations);	
		}
	}

	return true;
}

bool styleTransferManager::SampleLearning()
{
	int time = 0;
	while (1)
	{
		//Sleep(1000);
		bool bDone = false;
		if (!m_pSom->Finished())
		{
			if (!m_pSom->Train())
			{
				break;
			}		
		}
		else
		{
			break;
		}
		StyleLearningState = 1;
	}
	
	if (linePicturepath != "")
	{
		m_pSom->ijs.clear();
		m_pSom->matchingindex.clear();
		while (1)
		{
			if (!m_pSom->Finished2())
			{

				if (!m_pSom->Train2())
				{
					//bDone = true; //quit if there is a problem
					break;
				}				
			}
			else
			{
				break;
			}

		}
		StyleLearningState = 2;
	}

	if (pointPicturepath != "")
	{
		m_pSom->ijs.clear();
		m_pSom->matchingindex.clear();
		while (1)
		{
			if (!m_pSom->Finished3())
			{

				if (!m_pSom->Train3())
				{
					//bDone = true; //quit if there is a problem
					break;
				}
				
			}
			else
			{
				break;
			}

		}
		StyleLearningState = 3;
	}

	if (backPicturepath != "")
	{
		m_pSom->ijs.clear();
		m_pSom->matchingindex.clear();
		while (1)
		{
			if (!m_pSom->Finished4())
			{

				if (!m_pSom->Train4())
				{
					//bDone = true; //quit if there is a problem
					break;
				}
				
			}
			else
			{
				break;
			}

		}
		StyleLearningState = 4;
	}

	return true;
}


bool styleTransferManager::getMapfile()
{
	m_ctAllCategory2.clear();
	clusters.clear();
	IsCustoms.clear();
	IsTranfers.clear();

	std::ifstream inFile(mapfilepath, std::ios::in);//inFile来自fstream,ifstream为输入文件流(从文件读入)
	std::string lineStr;
	ctCategory layer0;
	layer0.list.resize(1);
	ctColor c0;
	while (getline(inFile, lineStr)) //getline来自sstream
	{
		std::stringstream ss(lineStr);//来自sstream
		std::string str;
		std::vector<std::string> lineArray;
		//按照逗号分隔
		while (getline(ss, str, ','))
			lineArray.push_back(str);//一行数据以vector保存
									 //cout<<lineArray[0]<<endl;
		layer0.id = lineArray[0];
		layer0.type2 = atoi(lineArray[1].c_str());
		c0.r = atoi(lineArray[2].c_str());
		c0.g = atoi(lineArray[3].c_str());
		c0.b = atoi(lineArray[4].c_str());
		c0.proportion = atof(lineArray[5].c_str());
		c0.type2 = layer0.type2;
		colorHist::colorTransfer(c0);
		layer0.list[0] = c0;
		layer0.proportion = c0.proportion;
		m_ctAllCategory2.emplace_back(layer0);
	}
	inFile.close();
	if (m_ctAllCategory2.size() == 0)
		return false;
	string relationshipFile;
	map<int, vector<int>> relationships;
	PictureOP::getRelationshipPath(mapfilepath, relationshipFile);
	std::ifstream inFile2(relationshipFile, std::ios::in);
	int index = 0;

	while (getline(inFile2, lineStr)) //getline来自sstream
	{
		std::stringstream ss(lineStr);//来自sstream
		std::string str;
		std::vector<std::string> lineArray;
		//按照逗号分隔
		while (getline(ss, str, ','))
			lineArray.push_back(str);//一行数据以vector保存
									 //cout<<lineArray[0]<<endl;
		if (atoi(lineArray[1].c_str()) == -2)
		{

		}
		else
		{
			int id = atoi(lineArray[1].c_str());
			map<int, vector<int>>::iterator iter = relationships.find(id);
			if (iter != relationships.end()) //查找到关系
			{
				iter->second.emplace_back(index);
			}
			else  //新关系
			{
				vector<int> list;				
				list.emplace_back(index);
				relationships.insert(pair<int, vector<int>>(atoi(lineArray[1].c_str()), list));
			}
		}

		if (lineArray.size()>3)
		{
			IsCustoms.emplace_back(atoi(lineArray[2].c_str()));
			IsTranfers.emplace_back(!atoi(lineArray[3].c_str()));
		}
		else
		{
			IsCustoms.emplace_back(true);
			IsTranfers.emplace_back(true);
		}
		


		index++;
	}
	inFile2.close();
	for (auto r : relationships)
		clusters.emplace_back(r.second);

	return true;
}


bool styleTransferManager::StyleTransfer()
{
	if (m_ctAllCategory2.size() < 1)
	{	
		return false;
	}
	std::vector<ctColor> colors = m_pSom->getColor();
	{
		std::sort(colors.begin(), colors.end(), [](const ctColor& x, const  ctColor& y)
		{
			return x.proportion > y.proportion;
		});
	}

	std::vector<ctColor> colorsPoint = m_pSom->getColorPoint();
	{
		std::sort(colorsPoint.begin(), colorsPoint.end(), [](const ctColor& x, const  ctColor& y)
		{
			return x.proportion > y.proportion;
		});
	}

	std::vector<ctColor> colorsLine = m_pSom->getColorLine();
	{
		std::sort(colorsLine.begin(), colorsLine.end(), [](const ctColor& x, const  ctColor& y)
		{
			return x.proportion > y.proportion;
		});
	}

	std::vector<ctColor> colorsBack = m_pSom->getColorBack();
	{
		std::sort(colorsBack.begin(), colorsBack.end(), [](const ctColor& x, const  ctColor& y)
		{
			return x.proportion > y.proportion;
		});
	}

	colorHist::colorTransfer(colors);
	colorHist::colorTransfer(colorsPoint);
	colorHist::colorTransfer(colorsLine);
	colorHist::colorTransfer(colorsBack);

	if (m_edNSGAPopNum % 2 == 1)
		m_edNSGAPopNum++;

	NSGA2* mnsga = new NSGA2();
	 
	mnsga->IsHarmony = false;
	mnsga->isBrightLimit = isBrightLimit;
	mnsga->delta4=BrightLimit;
	mnsga->NSGA2popsize = m_edNSGAPopNum;
	mnsga->FeaLimitation = ValueFeaLimitation;
	mnsga->NSGA2generation = m_edIterationTime;
	mnsga->picnum = Mulit_Picture;
	mnsga->setOptimizingWeights(Optimizing_wc, Optimizing_wse, Optimizing_wsi, Optimizing_wha, Optimizing_whi, Optimizing_wseq, Optimizing_wdi);
	mnsga->iniColorDis(mcolorSDis, mColorDis2, mColorDis3);
	mnsga->iniParameters(colors, colorsLine, colorsPoint, colorsBack, m_ctAllCategory2, clusters, IsCustoms, IsTranfers);
	population pop(mnsga);
	pop.maincal();

	resG.resize(m_ResultNum);
	resultScores.resize(m_ResultNum);
	for (int i = 0; i < m_ResultNum; i++)
	{
		resG[i] = pop.getResult(i*1.0 / (m_ResultNum - 1), resultScores[i]);
	}
	Output(pop);
	return true;
}

string styleTransferManager::getResultPathByIndex(int index)
{
	if (index < 0)
		index = 0;
	if (index >= m_ResultNum)
		index = m_ResultNum - 1;
	
	float rindex = index*1.0 / (m_ResultNum - 1);
	return  "../output/result_" + population::doubleToString(rindex) + ".csv";
}

void styleTransferManager::SetBrightLimit(double is)
{
	if (is < 0)
	{
		isBrightLimit = false;
		BrightLimit = 0;
	}
	else
	{
		isBrightLimit = true;
		BrightLimit = is;
	}
	
}

std::vector<ctColor> styleTransferManager::getBlob()
{
	return m_pSom->getColorPoint();
}

std::vector<ctColor> styleTransferManager::getFill()
{
	return m_pSom->getColor();
}

std::vector<ctColor> styleTransferManager::getStroke()
{
	return m_pSom->getColorLine();
}

std::vector<ctColor> styleTransferManager::getBackground()
{
	return m_pSom->getColorBack();
}




vector<vector<double>>styleTransferManager::getResultScores()
{
	return resultScores;
}

void styleTransferManager::Output(const population& pop)
{
	if (!PictureOP::exists("../output"))
	{
		int ret = _mkdir("../output");
		if (ret != 0)
			return;
	}
	for (int j = 0; j < m_ResultNum; j++)
	{
		string  outp = "../output/MapColor" + population::doubleToString((j*1.0 / (m_ResultNum - 1))) + ".csv";
		auto &res = resG[j];
		ofstream xxx(outp);
		for (int i = 0; i < m_ctAllCategory2.size(); i++)
		{
			xxx << m_ctAllCategory2[i].id << "," << m_ctAllCategory2[i].type2 << "," << m_ctAllCategory2[i].list[0].r << "," << m_ctAllCategory2[i].list[0].g << "," << m_ctAllCategory2[i].list[0].b << ","
				<< res[i].r << "," << res[i].g << "," << res[i].b << endl;
		}
		xxx.close();
	}

	FILE *p;
	p = fopen("../output/My_NSGA2.csv", "w+");
	for (int i = 0; i<m_edNSGAPopNum; i++)
	{
		fprintf(p, "%f,%f,%f\n", 1 - pop.P[i].fvalue[0], 1 - pop.P[i].fvalue[1], pop.P[i].IsFeasible);
	}
	fclose(p);
}


float styleTransferManager::calMapH(const vector<vector<double>>& rgbs0)
{
	int index = 0;
	int D = rgbs0.size();
	cv::Mat rgbs(3, D, CV_64F);
	hueProbs hueProb;
	FeatureCompute::HInit(hueProb);
	MapHarmony::hueProb = hueProb;
	for (int i = 0; i < D; i++)
	{
		{
			rgbs.at<double>(0, index) = rgbs0[i][0];
			rgbs.at<double>(1, index) = rgbs0[i][1];
			rgbs.at<double>(2, index) = rgbs0[i][2];

			index++;
		};
	}
	rgbs /= 255;
	return MapHarmony::glmnetPredict2(rgbs) / 5;
}