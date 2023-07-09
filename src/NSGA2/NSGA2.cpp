#include "NSGA2.h"
#include "opencv2/core/core.hpp"
#include "../ColorHarmony/MapHarmony.h"


void individual::calColorDisNet()
{
	int allNum = mNSGA->ubP*mNSGA->ubP + mNSGA->ubL*mNSGA->ubL + mNSGA->ubBack*mNSGA->ubBack + mNSGA->ubPoint*mNSGA->ubPoint;
	vector<ctColor> allColor;
	allColor.insert(allColor.end(), mNSGA->m_colorRGBDataPoint.begin(), mNSGA->m_colorRGBDataPoint.end());
	allColor.insert(allColor.end(), mNSGA->m_colorRGBDataLine.begin(), mNSGA->m_colorRGBDataLine.end());
	allColor.insert(allColor.end(), mNSGA->m_colorRGBDataP.begin(), mNSGA->m_colorRGBDataP.end());
	allColor.insert(allColor.end(), mNSGA->m_colorRGBDataBack.begin(), mNSGA->m_colorRGBDataBack.end());
	vector<vector<float>> disNet(allNum);
	vector<vector<double>> colorsLabs(allNum);

	for (int i = 0; i < allNum; i++)
		colorsLabs[i] = ColorSpaceTransfer::RGB2Lab(allColor[i].r, allColor[i].g, allColor[i].b);


	for (int i = 0; i < allNum - 1; i++)
	{
		disNet[i].resize(allNum, 0.0f);

		for (int j = i + 1; j < allNum; j++)
		{
			disNet[i][j] = sqrt(pow(colorsLabs[i][0] - colorsLabs[j][0], 2) + pow(colorsLabs[i][1] - colorsLabs[j][1], 2) + pow(colorsLabs[i][2] - colorsLabs[j][2], 2));
		}
	}

	for (int i = 0; i < allNum - 1; i++)
	{

	}

}


void individual::getManifoldCenter(const int &type,const std::vector<ctColor> &Mainfoldcolors,vector<int>& colorCenters)
{
	cv::Mat centers;
	cv::Mat labels;
	int size = Mainfoldcolors.size();
	cv::Mat colorpoint(size, 3, CV_32F, cv::Scalar(10));

	int numCluster = count(mNSGA->type2.begin(), mNSGA->type2.end(), type);
	if (numCluster == 0)
		return;
	//vector<vector<double>> mcolorlabs(Mainfoldcolors.size());
	for (size_t i = 0; i < size; i++)
	{
		colorpoint.at<float>(i, 0) = (float)Mainfoldcolors[i].L;
		colorpoint.at<float>(i, 1) = (float)Mainfoldcolors[i].A;
		colorpoint.at<float>(i, 2) = (float)Mainfoldcolors[i].B;
	}
	cv::kmeans(colorpoint, numCluster, labels,
	cv::TermCriteria(cv::TermCriteria::EPS + cv::TermCriteria::COUNT, 10, 0.1), 3, cv::KMEANS_PP_CENTERS, centers);

	colorCenters.resize(numCluster);
	vector<double> mindists(numCluster,INT32_MAX);
	for (int j = 0; j < size; j++)
	{
		int index = labels.at<int>(j);
		double dist = sqrt(pow(colorpoint.at<float>(j, 0) - centers.at<float>(index, 0), 2)
						+ pow(colorpoint.at<float>(j, 1) - centers.at<float>(index, 1), 2)
						+pow(colorpoint.at<float>(j, 2) - centers.at<float>(index, 2), 2));
		if (dist < mindists[index])
		{
			mindists[index] = dist;
			colorCenters[index] = j;
		}
	}

}


vector<int> individual::GenerateDiffNumber(int min, int max, int num)
{
	int rnd;
	vector<int> diff;
	vector<int> tmp;//存储剩余的数
					//初始化
	for (int i = min; i < max + 1; i++)
	{
		tmp.push_back(i);
	}
	for (int i = 0; i < num; i++)
	{
		do {
			rnd = min + rand() % (max - min + 1);

		} while (tmp.at(rnd - min) == -1);
		diff.push_back(rnd);
		tmp.at(rnd - min) = -1;
	}
	return diff;
}

 


void individual::init()
{
	value.resize(mNSGA->Dimension2);

	vector<int> indexs(4, 0);
	vector<vector<int>> colorCenters(indexs.size());
	vector<vector<int>> indexsAll(indexs.size());
	for (int i = 0; i < indexs.size(); i++)
	{		
		getManifoldCenter(i,mNSGA->getdiffcolors(i), colorCenters[i]);
		indexsAll[i] = GenerateDiffNumber(0, colorCenters[i].size()-1, colorCenters[i].size());
	}

	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		int ub = mNSGA->getdiffub(i);
		value[i].x = colorCenters[mNSGA->type2[i]][indexsAll[mNSGA->type2[i]][indexs[mNSGA->type2[i]]]] / (ub + 1);
		value[i].y = colorCenters[mNSGA->type2[i]][indexsAll[mNSGA->type2[i]][indexs[mNSGA->type2[i]]]] % (ub + 1);
		indexs[mNSGA->type2[i]]++;
	}
}

float individual::calHarmony(const vector<ctColor>& colors)
{
	/*ColorData rgbs;
	vector<float> hs;
	float h = 0.0f;
	int num = 0;
	float colorp = 0.0;
	int D = mNSGA->Dimension2;
	while (num + 5 < D)
	{
		colorp = 0.0;
		for (int i = num; i < num + 5; i++)
		{
			rgbs.colors.emplace_back(HColor{ 1.0*colors[i].r / 255,1.0*colors[i].g / 255,1.0*colors[i].b / 255 });
			colorp += mNSGA->m_dPictureInformation[i].proportion;
		}
		num += 5;
		hs.emplace_back(colorp*colorH::glmnetPredict(rgbs));

		rgbs.colors.clear();
	}
	if (num < D - 1)
	{
		colorp = 0.0;
		for (int i = 0; i < 5; i++)
		{
			rgbs.colors.emplace_back(HColor{ 1.0*colors[D - i - 1].r / 255,1.0*colors[D - i - 1].g / 255,1.0*colors[D - i - 1].b / 255 });
			colorp += mNSGA->m_dPictureInformation[D - i - 1].proportion;
		}
		hs.emplace_back(colorp*colorH::glmnetPredict(rgbs));
		rgbs.colors.clear();
	}

	for (auto &hh : hs)
	{
		h += hh / 5;
	}
	h /= hs.size();
	h /= 3;
	return -log(h);*/

	return 1;
}

float individual::calHarmony2(const vector<ctColor>& colors)
{
	float h = 0.0f;
	int D = mNSGA->Dimension2;

	
	int index = 0;
	cv::Mat rgbs(3, D, CV_64F);

	for (int i = 0; i < D; i++)
	{
		{
			rgbs.at<double>(0, index) = colors[i].r;
			rgbs.at<double>(1, index) = colors[i].g;
			rgbs.at<double>(2, index) = colors[i].b;

			index++;
		};
	}	
	rgbs /= 255;
	h = MapHarmony::glmnetPredict2(rgbs)/5;
	
	return h;
}


double individual::getFunction2(const vector<ctColor>& colors)
{
	double harmony = 1.0;
	double similarity = 0.0;

	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		similarity += 1 - abs(colors[i].proportion - mNSGA->m_dPictureInformation[i].proportion) / max(colors[i].proportion, mNSGA->m_dPictureInformation[i].proportion);
	}
	
	if (mNSGA->IsHarmony)
		harmony = calHarmony2(colors);
		
	//similarity /= mNSGA->m_Mapproportion2;
	similarity /= mNSGA->Dimension2;
	//mNSGA->colorF1 << mNSGA->wsi*similarity << "," << mNSGA->wha*harmony;
	return 1 - similarity*harmony;
}

void individual::getFunction2(const vector<ctColor>& colors, vector<double>& f2)
{
	f2.resize(2);
	double harmony = 1.0;
	double similarity = 0.0;

	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		similarity += 1 - abs(colors[i].proportion - mNSGA->m_dPictureInformation[i].proportion) / max(colors[i].proportion, mNSGA->m_dPictureInformation[i].proportion);
	}

	if (mNSGA->IsHarmony)
		harmony = calHarmony2(colors);

	similarity /= mNSGA->Dimension2;
	f2[0] = similarity;
	f2[1] = harmony;
}


float individual::calSequenceLoss(const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab)
{
	float SequenceLossAll = 0.0f;
	for (auto& cluster : mNSGA->Clusters)
	{
		double SequenceLoss = 0.0f;
		double p = 0.0f;
		if (cluster.size() > 2)   //order
		{
			//级别型  S递增，V递减
			if (cluster[0] != -1)
			{
				float dH = 0.0f;
				for (int i = 0; i < cluster.size() - 1; i++)
				{
					p += mNSGA->m_dPictureInformation[cluster[i]].proportion;

					dH = ColorSpaceTransfer::getColorDH(colorsHSV[cluster[i]][0], colorsHSV[cluster[i + 1]][0]);

					if (dH >= mNSGA->LdH2)
						SequenceLoss += ((dH - mNSGA->LdH2));
						//SequenceLoss += 1;

					{
						if (colorsHSV[cluster[i + 1]][1] < colorsHSV[cluster[i]][1])
							SequenceLoss += (((colorsHSV[cluster[i]][1]) - colorsHSV[cluster[i + 1]][1]));
							//SequenceLoss += 1;

						if (colorsHSV[cluster[i]][2] < colorsHSV[cluster[i + 1]][2])
							SequenceLoss += (((colorsHSV[cluster[i + 1]][2]) - colorsHSV[cluster[i]][2]));
							//SequenceLoss += 1;
					}
				}
			}
			//两端型  S先递减后递增  V先递增后递减
			else
			{
				float dH = 0.0f;

				for (int i = 1; i < cluster.size() - 1; i++)
				{
					p = mNSGA->m_dPictureInformation[cluster[i]].proportion + mNSGA->m_dPictureInformation[cluster[i + 1]].proportion;
					dH = ColorSpaceTransfer::getColorDH(colorsHSV[cluster[i]][0], colorsHSV[cluster[i + 1]][0]);
				

					if (i < cluster.size() / 2)
					{
						if (dH > mNSGA->LdH2)
							SequenceLoss += p*(dH - mNSGA->LdH2);
						
						float dS = colorsHSV[cluster[i]][1] - colorsHSV[cluster[i + 1]][1];
						float dV = colorsHSV[cluster[i + 1]][2] - colorsHSV[cluster[i]][2];
						if (dS <= mNSGA->LdS)
							SequenceLoss += p*(mNSGA->LdS - dS);
						if (dV < mNSGA->LdV)
							SequenceLoss += p*(mNSGA->LdV - dV);

					}
					if (i == cluster.size() / 2)
					{
						if (dH <= mNSGA->DiverH)
							SequenceLoss += p*((mNSGA->DiverH - dH));
						
						float dV = abs(colorsHSV[cluster[i]][2] - colorsHSV[cluster[i + 1]][2]);
						if (dV >= mNSGA->DiverV)
							SequenceLoss += p*((dV -mNSGA->DiverV));
					}
					if (i > cluster.size() / 2)
					{
						if (dH > mNSGA->LdH2)
							SequenceLoss += p*(dH - mNSGA->LdH2);
						
						float dS = colorsHSV[cluster[i + 1]][1] - colorsHSV[cluster[i]][1];
						float dV = colorsHSV[cluster[i]][2] - colorsHSV[cluster[i + 1]][2];
						if (dS <= mNSGA->LdS)
							SequenceLoss += p*(mNSGA->LdS - dS);
						if (dV < mNSGA->LdV)
							SequenceLoss += p*(mNSGA->LdV - dV);

					}
				}
				p = mNSGA->m_dPictureInformation[cluster[1]].proportion + mNSGA->m_dPictureInformation[cluster[cluster.size()-1]].proportion;
				dH = ColorSpaceTransfer::getColorDH(colorsHSV[cluster[1]][0], colorsHSV[cluster[cluster.size()-1]][0]);
				if (dH <= mNSGA->DiverH)
					SequenceLoss += p*((mNSGA->DiverH - dH));
			}
		}
		//association
		else
		{
			for (int i = 0; i < cluster.size() - 1; i++)
			{
				float dis = sqrt(pow(colorsLab[cluster[i]][0] - colorsLab[cluster[i + 1]][0], 2) + pow(colorsLab[cluster[i]][1] - colorsLab[cluster[i + 1]][1], 2) + pow(colorsLab[cluster[i]][2] - colorsLab[cluster[i + 1]][2], 2));
				p = mNSGA->m_dPictureInformation[cluster[i]].proportion + mNSGA->m_dPictureInformation[cluster[i + 1]].proportion;


				if (dis < mNSGA->AsscoDis)
					SequenceLoss += 0.0;
				else
					SequenceLoss += p*((dis - mNSGA->AsscoDis));
				//SequenceLoss += 1;
			}
		}
		SequenceLossAll += SequenceLoss;
		//SequenceLossAll += SequenceLoss;
	}
	return SequenceLossAll;
}

bool individual::IsInDiverging(const int& index1, const int& index2, vector<vector<int>>& clusters)
{
	for (auto &cluster : clusters)
	{
		if (cluster[0] == -1)
		{
			vector<int>::iterator  ite1 = find(cluster.begin(), cluster.end(), index1);
			vector<int>::iterator ite2 = find(cluster.begin(), cluster.end(), index2);
			if (ite1 != cluster.end() && ite2 != cluster.end())
			{
				if(std::distance(cluster.begin(), ite1)<cluster.size()/2 && std::distance(cluster.begin(), ite2)<cluster.size() / 2
					|| std::distance(cluster.begin(), ite1)>cluster.size() / 2 && std::distance(cluster.begin(), ite2)>cluster.size() / 2)
					return true;
			}
		}
			
	}
	return false;
}


bool individual::IsInSequence(const int& index1, const int& index2, const vector<vector<int>>& clusters)
{
	for (auto &cluster : clusters)
	{
		if (cluster[0] != -1)
			if (find(cluster.begin(), cluster.end(), index1) != cluster.end()
				&& find(cluster.begin(), cluster.end(), index2) != cluster.end())
				return true;
	}
	return false;
}


bool individual::IsInconstLayer(const int& index1)
{
	if (find(mNSGA->constLayer.begin(), mNSGA->constLayer.end(), index1) != mNSGA->constLayer.end())
		return true;
	return false;
}

int NSGA2::getGroupID(int index)
{
	int num = 0;
	for (auto &cluster : Clusters)
	{
		if (find(cluster.begin(), cluster.end(), index) != cluster.end())
		{
			return num;
		}
		num++;
	}
	return -1;
}

int NSGA2::getColorRelationship(const int& index1,const int &index2)
{
	for (auto &cluster : Clusters)
	{
		if (find(cluster.begin(), cluster.end(), index1) != cluster.end()
			&& find(cluster.begin(), cluster.end(), index2) != cluster.end())
		{
			if (cluster.size() == 2)
				return 1;
			else
				return 2;
		}
	}
	return -1;
}

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
double individual::getColorConObjective(const int& index, const vector<vector<double>> &colorsLab)
{
	if (!mNSGA->isColorConst[index])
		return 1.0;
	else
	{
		float dis = sqrt(pow(colorsLab[index][0] - mNSGA->sourceLabs[index][0], 2) + pow(colorsLab[index][1] - mNSGA->sourceLabs[index][1], 2) + pow(colorsLab[index][2] - mNSGA->sourceLabs[index][2], 2));
		
		return min(1.0f, mNSGA->delta2 / dis);
		//return max(1.0f, dis - mNSGA->delta2);
	}
}



double individual::getColorSemanticObjective(const int& i, const int& j, const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab)
{

	double dist = ColorSpaceTransfer::calColorDist(colorsLab[i], colorsLab[j]);
	float dH = ColorSpaceTransfer::getColorDH(colorsHSV[i][0], colorsHSV[j][0]);
	float dS = colorsHSV[j][1] - colorsHSV[i][1];
	float dV = colorsHSV[j][2] - colorsHSV[i][2];
	int classNum = j - i;

	float H1,H2;
	float S1,S2;
	float V1,V2;

	
	if (mNSGA->colorRelationShips[i][j] == 1)
	{
		return min(1.0, mNSGA->AsscoDis / dist);
	}
	else if(mNSGA->colorRelationShips[i][j]==2)
	{
		
		int groupID = mNSGA->getGroupID(i);
		float GroupdH = GroupdHs[groupID];
		if (GroupdH > 0)
		{
			if (dH > 0)
				return min(min(1.0f, dH / (classNum*mNSGA->LdH2 - mNSGA->LdH)), (classNum*mNSGA->LdH2 + mNSGA->LdH) / dH)
				*min(min(1.0, colorsHSV[j][1] / (classNum*mNSGA->LdV2 + colorsHSV[i][1])), ((classNum + 1)*mNSGA->LdV2 + colorsHSV[i][1]) / colorsHSV[j][1])
				*min(min(1.0, colorsHSV[i][2] / (classNum*mNSGA->LdV2 + colorsHSV[j][2])), ((classNum + 1)*mNSGA->LdV2 + colorsHSV[j][2]) / colorsHSV[i][2]);
			else
				return 0.0;
		}
		else
		{
			if (dH<0)
				return min(min(1.0f, -dH / (classNum*mNSGA->LdH2 - mNSGA->LdH)), (classNum*mNSGA->LdH2 + mNSGA->LdH) / -dH)
				*min(min(1.0, colorsHSV[j][1] / (classNum*mNSGA->LdV2 + colorsHSV[i][1])), ((classNum + 1)*mNSGA->LdV2 + colorsHSV[i][1]) / colorsHSV[j][1])
				*min(min(1.0, colorsHSV[i][2] / (classNum*mNSGA->LdV2 + colorsHSV[j][2])), ((classNum + 1)*mNSGA->LdV2 + colorsHSV[j][2]) / colorsHSV[i][2]);
			else
				return 0.0;
		}

		return min(1.0f, mNSGA->LdH / abs(dH))
			*min(min(1.0, colorsHSV[j][1] / (classNum*mNSGA->LdV2 + colorsHSV[i][1])), ((classNum + 1)*mNSGA->LdV2 + colorsHSV[i][1]) / colorsHSV[j][1])
			*min(min(1.0, colorsHSV[i][2] / (classNum*mNSGA->LdV2 + colorsHSV[j][2])), ((classNum + 1)*mNSGA->LdV2 + colorsHSV[j][2]) / colorsHSV[i][2]);
	}	
	else
	{
		return min(1.0, dist / mNSGA->AsscoDis);
	}
}

double individual::getFunction1_2(const vector<ctColor>& colors, const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab)
{
	vector<vector<float>> Semantics(mNSGA->Dimension2, vector<float>(mNSGA->Dimension2, 0.0));
	double Semantic = 0.0;
	double sem = 0.0;
	double convention = 0.0;
	double objective1 = 0.0;

	float p = 1.0f;
	int i;
	for (i = 0; i < mNSGA->Dimension2-1; i++)
	{
		for (int j = i + 1; j < mNSGA->Dimension2; j++)
		{
			if (!mNSGA->IsTranfers[i] && !mNSGA->IsTranfers[j])
			{
				Semantics[i][j] = 1;
				Semantics[j][i] = 1;
				continue;
			}
			sem = getColorSemanticObjective(i, j, colorsHSV, colorsLab);
			Semantics[i][j] = sem;
			Semantics[j][i] = sem;
			
			//Semantic += p* getColorSemanticObjective(i, j, colorsHSV, colorsLab);			
		}
		Semantic += accumulate(Semantics[i].begin(), Semantics[i].end(), 0.0f) / (mNSGA->Dimension2 - 1);
		convention += getColorConObjective(i, colorsLab);
	}
	Semantic += accumulate(Semantics[i].begin(), Semantics[i].end(), 0.0f) / (mNSGA->Dimension2 - 1);
	convention += getColorConObjective(i, colorsLab);
	objective1 = Semantic / mNSGA->Dimension2*convention / mNSGA->Dimension2;
	// return (exp(objective1)- exp(-objective1))/ (exp(objective1) + exp(-objective1));
	return 1 - objective1;
}


void individual::getFunction1_2(const vector<ctColor>& colors, const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab
								, vector<double>& f1)
{
	f1.resize(2);
	vector<vector<float>> Semantics(mNSGA->Dimension2, vector<float>(mNSGA->Dimension2, 0.0));
	double Semantic = 0.0;
	double sem = 0.0;
	double convention = 0.0;
	double objective1 = 0.0;
	float p = 1.0f;
	int i;
	for (i = 0; i < mNSGA->Dimension2 - 1; i++)
	{
		for (int j = i + 1; j < mNSGA->Dimension2; j++)
		{
			sem = getColorSemanticObjective(i, j, colorsHSV, colorsLab);
			Semantics[i][j] = sem;
			Semantics[j][i] = sem;

			//Semantic += p* getColorSemanticObjective(i, j, colorsHSV, colorsLab);			
		}
		Semantic += accumulate(Semantics[i].begin(), Semantics[i].end(), 0.0f) / (mNSGA->Dimension2 - 1);
		convention += getColorConObjective(i, colorsLab);


	}
	Semantic += accumulate(Semantics[i].begin(), Semantics[i].end(), 0.0f) / (mNSGA->Dimension2 - 1);
	convention += getColorConObjective(i, colorsLab);
	

	f1[0] = Semantic / mNSGA->Dimension2;
	f1[1]= convention/ mNSGA->Dimension2;
}


double individual::getFunction1(const vector<ctColor>& colors, const vector<vector<double>>& colorsHSV, const vector<vector<double>> &colorsLab)
{
	//double semantic = 0.0;
	double convention1 = 0.0;
	double convention2 = 0.0;
	double convention = 0.0;
	//double hierarchy = 0.0;
	double SequenceLoss = 0.0;
	double differentiation = 0.0;

	for (int i = 0; i < mNSGA->constLayer.size(); i++)
	{
		int index = mNSGA->constLayer[i];
		float dis = sqrt(pow(colorsLab[index][0] - mNSGA->sourceLabs[index][0], 2) + pow(colorsLab[index][1] - mNSGA->sourceLabs[index][1], 2) + pow(colorsLab[index][2] - mNSGA->sourceLabs[index][2], 2));
		if (dis < mNSGA->delta2)
			convention1 += 0.0;
		else
			//convention += mNSGA->m_dPictureInformation[index].proportion;
			convention1 += ((dis - mNSGA->delta2))*mNSGA->m_dPictureInformation[index].proportion;
	}
	for (int i = 0; i < mNSGA->constLayer.size(); i++)
	{
		int index = mNSGA->constLayer[i];
		for (int j = 0; j < mNSGA->Dimension2; j++)
		{
			if (index == j || IsInSequence(index, j, mNSGA->Clusters))
				continue;
			else
			{
				float dis = sqrt(pow(colorsLab[index][0] - colorsLab[j][0], 2) + pow(colorsLab[index][1] - colorsLab[j][1], 2) + pow(colorsLab[index][2] - colorsLab[j][2], 2));
				if (dis < 0.01)
					convention2 += INT_MAX;
				else if (dis > mNSGA->delta13)
					convention2 += 0.0;
				else
					//convention2 += (mNSGA->m_dPictureInformation[index].proportion + mNSGA->m_dPictureInformation[j].proportion)/2;
					convention2 += ((mNSGA->delta13 - dis))*(mNSGA->m_dPictureInformation[j].proportion);
			}
		}
	}
	convention = convention1 + convention2/(mNSGA->Dimension2 - mNSGA->constLayer.size())/3;
	//convention = convention1;

	SequenceLoss = calSequenceLoss(colorsHSV, colorsLab);

	for (int i = 0; i < mNSGA->Dimension2 - 1; i++)
	{
		for (int j = i + 1; j < mNSGA->Dimension2; j++)
		{
			if (IsInSequence(i, j, mNSGA->Clusters) || IsInDiverging(i, j, mNSGA->Clusters))
				continue;

			float dis = sqrt(pow(colorsLab[i][0] - colorsLab[j][0], 2) + pow(colorsLab[i][1] - colorsLab[j][1], 2) + pow(colorsLab[i][2] - colorsLab[j][2], 2));
			if (dis < 0.01)
				differentiation += INT_MAX;
			else if (dis > mNSGA->delta12)
				differentiation += 0.0;
			else
				//differentiation += (mNSGA->m_dPictureInformation[i].proportion + mNSGA->m_dPictureInformation[j].proportion)/2;
				differentiation += ((mNSGA->delta12 - dis))*(mNSGA->m_dPictureInformation[i].proportion + mNSGA->m_dPictureInformation[j].proportion);
		}
	}
	/*for (int i = 0; i < mNSGA->Dimension2 - 1; i++)
	{
		for (auto cluster : mNSGA->Clusters)
			if (find(cluster.begin(), cluster.end(), i) != cluster.end())
				continue;
		for (int j = i + 1; j < mNSGA->Dimension2; j++)
		{
			float dis = sqrt(pow(colorsLab[i][0] - colorsLab[j][0], 2) + pow(colorsLab[i][1] - colorsLab[j][1], 2) + pow(colorsLab[i][2] - colorsLab[j][2], 2));
			if (dis < 0.01)
				differentiation += INT_MAX;
			else if (dis > mNSGA->delta12)
				differentiation += 0.0;
			else
				//differentiation += (mNSGA->m_dPictureInformation[i].proportion + mNSGA->m_dPictureInformation[j].proportion)/2;
				differentiation += ((mNSGA->delta12 - dis))*(mNSGA->m_dPictureInformation[i].proportion + mNSGA->m_dPictureInformation[j].proportion);
		}
	}*/


	//mNSGA->colorF1 << mNSGA->wc*(convention) << "," << mNSGA->wseq*SequenceLoss << ",";
	return mNSGA->wc*(convention)+mNSGA->wseq*SequenceLoss + mNSGA->wdi*differentiation;
}

double individual::getFunction3(const vector<vector<double>> &colorsLab)
{
	double semantic = 0.0;
	double hierarchy = 0.0;

	for (int i = 0; i < mNSGA->Dimension2 - 1; i++)
	{
		for (int j = i + 1; j < mNSGA->Dimension2; j++)
		{
			if (!mNSGA->IsTranfers[i] && !mNSGA->IsTranfers[j])
			{				
				continue;
			}


			float dis = sqrt(pow(colorsLab[i][0] - colorsLab[j][0], 2) + pow(colorsLab[i][1] - colorsLab[j][1], 2) + pow(colorsLab[i][2] - colorsLab[j][2], 2));
			if (dis < 5)
				semantic += INT_MAX;
			else if (dis > mNSGA->delta1)
				semantic += 0.0;
			else
				semantic += ((mNSGA->delta1 - dis));
		}
	}

	for (int i = 0; i < mNSGA->Dimension2 - 1; i++)
	{
		if (mNSGA->type2[i] != 2)
		{
			float dis = sqrt(pow(colorsLab[i][0] - colorsLab[mNSGA->backLayerID][0], 2) +
				pow(colorsLab[i][1] - colorsLab[mNSGA->backLayerID][1], 2) + pow(colorsLab[i][2] - colorsLab[mNSGA->backLayerID][2], 2));
			if (dis < 5)
				hierarchy += INT_MAX;
			else if (dis > mNSGA->delta3)
				hierarchy += 0.0;
			else
				hierarchy += ((mNSGA->delta3 - dis));

			if (mNSGA->isBrightLimit)
			{
				float dL = abs(colorsLab[i][0] - colorsLab[mNSGA->backLayerID][0]);
				if (dL < mNSGA->delta4)
					hierarchy += (mNSGA->delta4 - dL);
			}
			
			/*if (mNSGA->type2[i] == 1)
			{
				if (dis > mNSGA->delta3 + 10)
					hierarchy += 0.0;
				else
					hierarchy += ((mNSGA->delta3 + 10 - dis));
			}*/
				
		}
	}
	return hierarchy + semantic;
}


void individual::getOrderGroupsDH(const vector<vector<double>> &colorsHSV,vector<float> & colorsdH)
{
	for (auto &cluster : mNSGA->Clusters)
	{
		float dh = 0.0f;

		if (cluster.size() >= 3)
		{
			for (int i = 0; i < cluster.size() - 1; i++)
			{
				for (int j = i + 1; j < cluster.size(); j++)
				{
					dh += ColorSpaceTransfer::getColorDH(colorsHSV[i][0],colorsHSV[j][0]) / (j - i);
				}		
			}	
			colorsdH.emplace_back(dh);
		}
		else
			colorsdH.emplace_back(0.0f);
	}
}

void individual::f_count()
{
	ctColor color;
	vector<ctColor> colors(mNSGA->Dimension2);
	vector<vector<double>> colorsLab(mNSGA->Dimension2);
	vector<vector<double>> colorsHSV(mNSGA->Dimension2);
	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		int x = value[i].x;
		int y = value[i].y;

		if (mNSGA->type2[i] == 0)
			color = mNSGA->m_colorRGBDataP[x * (mNSGA->ubP+1) + y];
		else if (mNSGA->type2[i] == 1)
			color = mNSGA->m_colorRGBDataLine[x * (mNSGA->ubL + 1) + y];
		else if (mNSGA->type2[i] == 2)
			color = mNSGA->m_colorRGBDataBack[x * (mNSGA->ubBack + 1) + y];
		else if (mNSGA->type2[i] == 3)
			color = mNSGA->m_colorRGBDataPoint[x * (mNSGA->ubPoint + 1) + y];

		if (!mNSGA->IsTranfers[i])
			color = mNSGA->m_dPictureInformation[i].list[0];

		colorsLab[i] = { color.L,color.A,color.B };
		colorsHSV[i] = { color.H,color.S,color.V };
		colors[i] = color;
	}

	
	getOrderGroupsDH(colorsHSV, GroupdHs);

	fvalue[0] = getFunction1_2(colors, colorsHSV, colorsLab);
	fvalue[1] = getFunction2(colors);
	//mNSGA->colorF1 << endl;
	IsFeasible = getFunction3(colorsLab);
	//fvalue[0] += IsFeasible;
	//fvalue[1] += IsFeasible;
}


void population::f_sort(int i)
{
	int n;
	n = len[i];
	qsort(F[i], n, sizeof(individual), cmp_c_d);
}

population::population(NSGA2* nsga)
{
	mNSGA = nsga;
	F = new individual*[2 * mNSGA->NSGA2popsize];
	for (int i = 0; i < 2 * mNSGA->NSGA2popsize; i++) {
		F[i] = new individual[2 * mNSGA->NSGA2popsize];
	}

	P.resize(mNSGA->NSGA2popsize);
	Q.resize(mNSGA->NSGA2popsize);
	R.resize(2 * mNSGA->NSGA2popsize);
	len.resize(2 * mNSGA->NSGA2popsize);
	int i;

	for (i = 0; i < mNSGA->NSGA2popsize; i++)
	{
		P[i].mNSGA = mNSGA;
		P[i].init();
		P[i].sp.resize(2 * mNSGA->NSGA2popsize);
	}
//#pragma omp parallel
	for (i = 0; i < mNSGA->NSGA2popsize; i++)
	{
		P[i].f_count();
	}
	for (int i = 0; i < mNSGA->NSGA2popsize; i++)
	{
		Q[i].value.resize(mNSGA->Dimension2);
		Q[i].sp.resize(2 * mNSGA->NSGA2popsize);
		R[i].value.resize(mNSGA->Dimension2);
		R[mNSGA->NSGA2popsize + i].value.resize(mNSGA->Dimension2);
		R[i].sp.resize(2 * mNSGA->NSGA2popsize);
		R[mNSGA->NSGA2popsize + i].sp.resize(2 * mNSGA->NSGA2popsize);
	}
	Pnum = mNSGA->NSGA2popsize;
	Qnum = 0;
	Rnum = 0;
}


void population::make_new_pop()
{
	int i, j, x, y, t1, t2, t3;
	double  u, b;
	//mNSGA->mark.assign(mNSGA->NSGA2popsize, 0);	
	//t3 = 0;
	//while (t3<mNSGA->NSGA2popsize)
	//{
	//	while (t1 = t2 = mNSGA->rand_int(0, mNSGA->NSGA2popsize - 1), mNSGA->mark[t1]);
	//	while (t1 == t2 || mNSGA->mark[t2])
	//	{
	//		t2 = mNSGA->rand_int(0, mNSGA->NSGA2popsize - 1);
	//	}
	//	t1 = choice(t1, t2);
	//	mNSGA->temp1[t3] = t1;
	//	//mNSGA->mark[t1] = 1;
	//}
	for (i = 0; i < mNSGA->NSGA2popsize; i += 2)
	{
		while (t1 = t2 = mNSGA->rand_int(0, mNSGA->NSGA2popsize - 1), mNSGA->mark[t1]);
		while (t1 == t2 || mNSGA->mark[t2])
		{
			t2 = mNSGA->rand_int(0, mNSGA->NSGA2popsize - 1);
		}
		t1 = choice(t1, t2);

		while (t1 = t2 = mNSGA->rand_int(0, mNSGA->NSGA2popsize - 1), mNSGA->mark[t1]);
		while (t1 == t2 || mNSGA->mark[t2])
		{
			t2 = mNSGA->rand_int(0, mNSGA->NSGA2popsize - 1);
		}
		t2 = choice(t1, t2);
		//if (s <= 0.9)
		{
			for (j = 0; j < mNSGA->Dimension2; j++)
			{
				u = mNSGA->rand_real((0.0 + 1e-6), (1.0 - 1e-6));
				if (u <= 0.5)
					b = pow(2 * u, 1.0 / 21);
				else
					b = 1.0 / pow(2 * (1 - u), 1.0 / 21);
				/*x = y = mNSGA->rand_int(0, mNSGA->NSGA2popsize / 2 - 1);
				while (x == y)
					y = mNSGA->rand_int(0, mNSGA->NSGA2popsize / 2 - 1);*/
				Q[i].value[j].x = 1.0 / 2 * ((1 - b)*P[t1].value[j].x + (1 + b)*P[t2].value[j].x);
				Q[i].value[j].y = 1.0 / 2 * ((1 - b)*P[t1].value[j].y + (1 + b)*P[t2].value[j].y);

				int ub = mNSGA->getdiffub(j);
				if (Q[i].value[j].x < 0)
					Q[i].value[j].x = 0;
				else if (Q[i].value[j].x > ub)
					Q[i].value[j].x = ub;
				if (Q[i].value[j].y < 0)
					Q[i].value[j].y = 0;
				else if (Q[i].value[j].y > ub)
					Q[i].value[j].y = ub;

				u = mNSGA->rand_real((0.0 + 1e-6), (1.0 - 1e-6));
				if (u <= 0.5)
					b = pow(2 * u, 1.0 / 21);
				else
					b = 1.0 / pow(2 * (1 - u), 1.0 / 21);
				if (i + 1 < mNSGA->NSGA2popsize)
				{
					Q[i + 1].value[j].x = 1.0 / 2 * ((1 + b)*P[t1].value[j].x + (1 - b)*P[t2].value[j].x);
					Q[i + 1].value[j].y = 1.0 / 2 * ((1 + b)*P[t1].value[j].y + (1 - b)*P[t2].value[j].y);
					if (Q[i + 1].value[j].x < 0)
						Q[i + 1].value[j].x = 0;
					else if (Q[i + 1].value[j].x > ub)
						Q[i + 1].value[j].x = ub;
					if (Q[i + 1].value[j].y < 0)
						Q[i + 1].value[j].y = 0;
					else if (Q[i + 1].value[j].y > ub)
						Q[i + 1].value[j].y = ub;
				}
			}

		}
		//else
		{
			for (j = 0; j < mNSGA->Dimension2; j++)
			{
				int um = 10;
				int ub = mNSGA->getdiffub(j);
				//x = mNSGA->rand_int(0, mNSGA->NSGA2popsize / 2 - 1);
				u = mNSGA->rand_real(0.0 + (1e-6), 1.0 - (1e-6));
				if (u < 0.5)
					u = pow(2 * u, 1.0 / (1 + um)) - 1;
				else
					u = 1 - pow(2 * (1 - u), 1.0 / (1 + um));
				Q[i].value[j].x = Q[i].value[j].x + ub*u;
				Q[i].value[j].y = Q[i].value[j].y + ub*u;
				if (Q[i].value[j].x < 0)
					Q[i].value[j].x = 0;
				else if (Q[i].value[j].x > ub)
					Q[i].value[j].x = ub;
				if (Q[i].value[j].y < 0)
					Q[i].value[j].y = 0;
				else if (Q[i].value[j].y > ub)
					Q[i].value[j].y = ub;


				//x = mNSGA->rand_int(0, mNSGA->NSGA2popsize / 2 - 1);
				u = mNSGA->rand_real(0.0 + (1e-6), 1.0 - (1e-6));
				if (u < 0.5)
					u = pow(2 * u, 1.0 / (1 + um)) - 1;
				else
					u = 1 - pow(2 * (1 - u), 1.0 / (1 + um));
				Q[i + 1].value[j].x = Q[i + 1].value[j].x + ub*u;
				Q[i + 1].value[j].y = Q[i + 1].value[j].y + ub*u;
				if (Q[i + 1].value[j].x < 0)
					Q[i + 1].value[j].x = 0;
				else if (Q[i + 1].value[j].x > ub)
					Q[i + 1].value[j].x = ub;
				if (Q[i + 1].value[j].y < 0)
					Q[i + 1].value[j].y = 0;
				else if (Q[i + 1].value[j].y > ub)
					Q[i + 1].value[j].y = ub;
			}


		}
	}
	Qnum = mNSGA->NSGA2popsize;
	for (i = 0; i < mNSGA->NSGA2popsize; i++)
	{
		Q[i].mNSGA = mNSGA;
	}
//    #pragma omp parallel
	for (i = 0; i < mNSGA->NSGA2popsize; i++)
	{
		Q[i].f_count();
	}
		
}

void population::set_p_q()
{
	Rnum = 0;
	Qnum = mNSGA->NSGA2popsize;
	int i;
	for (i = 0; i < Pnum; i++)
		R[Rnum++] = P[i];
	for (i = 0; i < Qnum; i++)
		R[Rnum++] = Q[i];
	/*for (i = 0; i<2 * mNSGA->NSGA2popsize; i++)
		R[i].f_count();*/
}

void population::fast_nondominated_sort()
{
	int i, j, k;
	vector<individual> H(2 * mNSGA->NSGA2popsize);

	int h_len = 0;
	for (i = 0; i < 2 * mNSGA->NSGA2popsize; i++)
	{
		R[i].np = 0;
		R[i].is_dominated = 0;
		len[i] = 0;
	}
	for (i = 0; i < 2 * mNSGA->NSGA2popsize; i++)
	{
		for (j = 0; j < 2 * mNSGA->NSGA2popsize; j++)
		{
			if (i != j)
			{
				if (e_is_dominated(R[i], R[j]))
					R[i].sp[R[i].is_dominated++] = j;
				else if (e_is_dominated(R[j], R[i]))
					R[i].np += 1;
			}
		}
		if (R[i].np == 0)
		{
			len_f = 1;
			F[0][len[0]++] = R[i];
		}

	}
	i = 0;
	while (len[i] != 0)
	{
		h_len = 0;
		for (j = 0; j < len[i]; j++)
		{
			for (k = 0; k < F[i][j].is_dominated; k++)
			{
				R[F[i][j].sp[k]].np--;
				if (R[F[i][j].sp[k]].np == 0)
				{
					H[h_len++] = R[F[i][j].sp[k]];
					R[F[i][j].sp[k]].rank = i + 2;
				}
			}
		}
		++i;
		if (i >= 2 * mNSGA->NSGA2popsize)
			break;
		len[i] = h_len;
		if (h_len != 0)
		{
			len_f++;
			for (j = 0; j < len[i]; j++)
				F[i][j] = H[j];
		}
	}
}

void population::calu_crowding_distance(int i)
{
	int n = len[i];
	double m_max, m_min;
	int j;
	for (j = 0; j < n; j++)
		F[i][j].crowding_distance = 0;

	qsort(F[i], n, sizeof(individual), cmp1);
	F[i][0].crowding_distance = F[i][n - 1].crowding_distance = 0xffffff;
	m_max = -0xfffff;
	m_min = 0xfffff;
	for (j = 0; j < n; j++)
	{
		if (m_max < F[i][j].fvalue[0])
			m_max = F[i][j].fvalue[0];
		if (m_min > F[i][j].fvalue[0])
			m_min = F[i][j].fvalue[0];
	}
	for (j = 1; j < n - 1; j++)
		F[i][j].crowding_distance += (F[i][j + 1].fvalue[0] - F[i][j - 1].fvalue[0]) / (m_max - m_min);

	qsort(F[i], n, sizeof(individual), cmp2);
	F[i][0].crowding_distance = F[i][n - 1].crowding_distance = 0xffffff;
	m_max = -0xfffff;
	m_min = 0xfffff;
	for (j = 0; j < n; j++)
	{
		if (m_max < F[i][j].fvalue[1])
			m_max = F[i][j].fvalue[1];
		if (m_min > F[i][j].fvalue[1])
			m_min = F[i][j].fvalue[1];
	}
	for (j = 1; j < n - 1; j++)
		F[i][j].crowding_distance += (F[i][j + 1].fvalue[1] - F[i][j - 1].fvalue[1]) / (m_max - m_min);
}

int population::choice(int a, int b)
{
	
	{
		if (P[a].rank < P[b].rank)
			return a;
		else if (P[a].rank == P[b].rank)
		{
			if (P[a].crowding_distance > P[b].crowding_distance)
				return a;
			else
				return b;
		}
		else
			return b;
	}
}

int population::getMinIndex(const int & findex)
{
	double minF =INT_MAX;
	int index = -1;
	for (int i = 0; i < len[0]; i++)
	{
		if (F[0][i].fvalue[findex] < minF)
		{
			minF = F[0][i].fvalue[findex];
			index = minF;
		}
	}
	return index;
}

void population::maincal()
{
	int s, i, j;
	s = mNSGA->NSGA2generation;
	make_new_pop();
	ofstream it("../output/迭代.csv");
	while (s--)
	{
		
		//std::cout << ("The %d NSGA2generation\n", s) << std::endl;
		set_p_q();
		fast_nondominated_sort();
		Pnum = 0;
		i = 0;

		while (Pnum + len[i] <= mNSGA->NSGA2popsize)
		{
			calu_crowding_distance(i);
			for (j = 0; j < len[i]; j++)
				P[Pnum++] = F[i][j];
			i++;
			if (i >= len_f)break;
		}

		if (i < len_f)
		{
			calu_crowding_distance(i);
			f_sort(i);
		}
		for (j = 0; j < mNSGA->NSGA2popsize - Pnum; j++)
			P[Pnum++] = F[i][j];

		if (s % 10 == 0)
		{
			int f1index = getMinIndex(0);
			int f2index = getMinIndex(1);
			it << mNSGA->NSGA2generation - s << "," << F[0][f1index].fvalue[0] << "," << F[0][f1index].IsFeasible << "," << F[0][f2index].fvalue[1] << "," << F[0][f2index].IsFeasible << endl;
		}
		make_new_pop();
	}
	it.close();
}

void individual::calObjective( vector<double>& resultscore)
{
	ctColor color;
	vector<ctColor> colors(mNSGA->Dimension2);
	vector<vector<double>> colorsLab(mNSGA->Dimension2);
	vector<vector<double>> colorsHSV(mNSGA->Dimension2);
	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		int x = value[i].x;
		int y = value[i].y;

		if (mNSGA->type2[i] == 0)
			color = mNSGA->m_colorRGBDataP[x * (mNSGA->ubP + 1) + y];
		else if (mNSGA->type2[i] == 1)
			color = mNSGA->m_colorRGBDataLine[x * (mNSGA->ubL + 1) + y];
		else if (mNSGA->type2[i] == 2)
			color = mNSGA->m_colorRGBDataBack[x * (mNSGA->ubBack + 1) + y];
		else if (mNSGA->type2[i] == 3)
			color = mNSGA->m_colorRGBDataPoint[x * (mNSGA->ubPoint + 1) + y];

		if (!mNSGA->IsTranfers[i])
			color = mNSGA->m_dPictureInformation[i].list[0];

		colorsLab[i] = { color.L,color.A,color.B };
		colorsHSV[i] = { color.H,color.S,color.V };
		colors[i] = color;
	}

	getOrderGroupsDH(colorsHSV, GroupdHs);

	vector<double> f1;
	vector<double> f2;

	getFunction1_2(colors, colorsHSV, colorsLab, f1);
	getFunction2(colors, f2);

	resultscore[2] = f1[0];
	resultscore[3] = f1[1];
	resultscore[4] = f2[0];
	resultscore[5] = f2[1];

}

std::vector<ctColor> population::getResult(const float &Rindex,vector<double>& resultscore)
{
	std::vector<ctColor> res;
	int Pindex = 0;

	if (Rindex == 0)
		qsort(&P[0], mNSGA->NSGA2popsize, sizeof(individual),cmp1);
	else if (Rindex == 1)
		qsort(&P[0], mNSGA->NSGA2popsize, sizeof(individual), cmp2);
	else
	{
		qsort(&F[0][0], len[0], sizeof(individual), cmp1);
		for (int i = 0; i < mNSGA->NSGA2popsize; i++)
		{
			if (P[i].fvalue[0] == F[0][int(len[0] * Rindex)].fvalue[0] && P[i].fvalue[1] == F[0][int(len[0] * Rindex)].fvalue[1])
			{
				Pindex = i;
				break;
			}
		}
	}

	ofstream sceq("../output/序列_" + doubleToString(Rindex) + ".csv");
	for (auto & cluster : mNSGA->Clusters)
	{
		for (int i = 0; i < cluster.size(); i++)
		{
			int index = cluster[i];
			if (index < 0)
				continue;
			int x = P[Pindex].value[index].x;
			int y = P[Pindex].value[index].y;

			ctColor color;
			if (mNSGA->type2[index] == 0)
				color = mNSGA->m_colorRGBDataP[x * (mNSGA->ubP+1) + y];
			if (mNSGA->type2[index] == 1)
				color = mNSGA->m_colorRGBDataLine[x * (mNSGA->ubL + 1) + y];
			if (mNSGA->type2[index] == 2)
				color = mNSGA->m_colorRGBDataBack[x * (mNSGA->ubBack + 1) + y];
			if (mNSGA->type2[index] == 3)
				color = mNSGA->m_colorRGBDataPoint[x * (mNSGA->ubPoint + 1) + y];

			if (!mNSGA->IsTranfers[i])
				color = mNSGA->m_dPictureInformation[i].list[0];
			vector<double> hsv = ColorSpaceTransfer::RGB2HSV(color.r, color.g, color.b);

			sceq << mNSGA->m_dPictureInformation[index].id << "," << hsv[0] << "," << hsv[1] << "," << hsv[2] << endl;
		}
		sceq << ".........................." << endl;
	}
	sceq.close();
	resultscore.resize(6);
	resultscore[0] = 1 - P[Pindex].fvalue[0];
	resultscore[1] = 1 - P[Pindex].fvalue[1];

	P[Pindex].calObjective(resultscore);

	ofstream iii("../output/result_" + doubleToString(Rindex) + ".csv");
	iii << 1 - P[Pindex].fvalue[0] << "," << 1 - P[Pindex].fvalue[1] << "," << P[Pindex].IsFeasible << endl;
	for (int i = 0; i < mNSGA->Dimension2; i++)
	{
		int x = P[Pindex].value[i].x;
		int y = P[Pindex].value[i].y;
		//int dub = mNSGA->getdiffub(i);
		if (mNSGA->type2[i] == 0)
			res.emplace_back(mNSGA->m_colorRGBDataP[x * (mNSGA->ubP + 1) + y]);
		if (mNSGA->type2[i] == 1)
			res.emplace_back(mNSGA->m_colorRGBDataLine[x * (mNSGA->ubL + 1) + y]);
		if (mNSGA->type2[i] == 2)
			res.emplace_back(mNSGA->m_colorRGBDataBack[x * (mNSGA->ubBack + 1) + y]);
		if (mNSGA->type2[i] == 3)
			res.emplace_back(mNSGA->m_colorRGBDataPoint[x * (mNSGA->ubPoint + 1) + y]);

		if (!mNSGA->IsTranfers[i])
			res[i] = mNSGA->m_dPictureInformation[i].list[0];
		auto vvv = res[i];
		iii << mNSGA->m_dPictureInformation[i].id << "," << mNSGA->type2[i] << "," << vvv.r << "," << vvv.g << "," << vvv.b << "," << mNSGA->m_dPictureInformation[i].proportion << ","<<  vvv.proportion << std::endl;

	}
	iii.close();


	return res;
}




void NSGA2::iniColorDis(const float& mcolorSDis,const float& mColorDis2, const float& mColorDis3)
{
	delta1 = mcolorSDis;   //颜色间色差 -语义
	delta3 = mColorDis2;   //视觉层次色差
	delta2 = mColorDis3;  //习惯用色色差
}

void NSGA2::iniParameters(const std::vector<ctColor>& colorP, const std::vector<ctColor>& colorL, const std::vector<ctColor>& colorPoint, const std::vector<ctColor>& colorBack,
	const std::vector<ctCategory>& InformationEntropy, const vector<vector<int>> &clusters, const vector<bool> &IsCustoms, const vector<bool> &IsTranfers0)
{
	Dimension2 = InformationEntropy.size();

	//初始化数据
	m_colorRGBDataP = colorP;
	m_colorRGBDataLine = colorL;
	m_colorRGBDataPoint = colorPoint;
	m_colorRGBDataBack = colorBack;
	Clusters = clusters;
	m_dPictureInformation = InformationEntropy;

	std::vector<vector<double>>().swap(sourceLabs);
	ofstream inputMapColor("../output/inputMapColor.csv");
	proportionAll = 0.0f;
	m_Mapproportion2 = 0.0;
	for (auto &MapColor : m_dPictureInformation)
	{
		m_Mapproportion2 += pow(MapColor.proportion, 2);
		proportionAll += MapColor.proportion;
		if (MapColor.list.size() > 0)
		{
			inputMapColor << MapColor.id << "," << MapColor.type2 << "," << MapColor.list[0].r << "," << MapColor.list[0].g << "," <<
				MapColor.list[0].b << "," << MapColor.list[0].proportion << endl;
		}
		sourceLabs.emplace_back(ColorSpaceTransfer::RGB2Lab(MapColor.list[0].r,
			MapColor.list[0].g, MapColor.list[0].b));
	}
	inputMapColor.close();


	vector<int>().swap(type2);
	float AllMaplineLength;
	{
		ofstream lineLength;
		lineLength.open("../output/linelength.csv");
		AllMaplineLength = 0;
		for (auto &layer : m_dPictureInformation)
		{
			type2.emplace_back(layer.list[0].type2);
			if (layer.list[0].type2 == 2)
			{
				backLayerID = type2.size() - 1;
			}
			if (layer.list[0].type2 == 1)
			{
				AllMaplineLength += layer.list[0].LineLength;
				lineLength << layer.id << "," << layer.list[0].LineLength << endl;
			}
		}
		lineLength << AllMaplineLength << endl;
		lineLength.close();
	}

	if (AllMaplineLength > 0.0)
		for (int i = 0; i < Dimension2; i++)
		{
			if (type2[i] == 1)
			{
				m_dPictureInformation[i].proportion = m_dPictureInformation[i].list[0].LineLength / AllMaplineLength;
			}
		}

	ubP = sqrt(picnum)* constNumCellsAcross - 1;
	ubL = sqrt(picnum)*constNumCellsAcrossLine - 1;
	ubPoint = sqrt(picnum)*constNumCellsAcrossPoint - 1;
	ubBack = sqrt(picnum)*constNumCellsAcrossBack - 1;


	constLayer.clear();
	ofstream constLayerName("../output/constLayer.csv");
	isColorConst.resize(m_dPictureInformation.size(),false);
	for (int i = 0; i < m_dPictureInformation.size(); i++)
	{
		if (IsCustoms[i])
		{
			isColorConst[i] = true;
			constLayer.emplace_back(i);
			constLayerName << i << "," << m_dPictureInformation[i].id << endl;
		}
	}
	constLayerName.close();

	IsTranfers = IsTranfers0;

	

	Dimension2 = InformationEntropy.size();
	

	mark.resize(NSGA2popsize);//A
	temp1.resize(NSGA2popsize);//临时数组
	temp2.resize(NSGA2popsize);//临时数组

	hueProbs hueProb;
	FeatureCompute::HInit(hueProb);
	MapHarmony::hueProb= hueProb;

	ofstream relationships("../output/relationship.csv");
	colorRelationShips.resize(Dimension2,vector<int>(Dimension2,-1));
	for (int i = 0; i < Dimension2-1; i++)
	{
		for (int j = i+1; j < Dimension2; j++)
		{
			colorRelationShips[i][j] = getColorRelationship(i, j);
			colorRelationShips[j][i] = colorRelationShips[i][j];

			relationships << i << "," << j << "," << m_dPictureInformation[i].id << "," 
				<< m_dPictureInformation[j].id << "," << colorRelationShips[i][j] << endl;
		}
	}
}

bool NSGA2::IsConstLayer(const string &name)
{
	if (name.find("Water") != std::string::npos)
		return true;
	if (name.find("Green") != std::string::npos)
		return true;
	if (name.find("water") != std::string::npos)
		return true;
	if (name.find("grass") != std::string::npos)
		return true;
	if (name.find("wood") != std::string::npos)
		return true;
	if (name.find("green") != std::string::npos)
		return true;
	if (name.find("forest") != std::string::npos)
		return true;
	if (name.find("river") != std::string::npos)
		return true;
	if (name.find("ocean") != std::string::npos)
		return true;
	if (name.find("lake") != std::string::npos)
		return true;

	if (name.find("0 - 200") != std::string::npos)
		return true;

	if (name.find("湖") != std::string::npos)
		return true;
	if (name.find("水体") != std::string::npos)
		return true;
	if (name.find("绿地") != std::string::npos)
		return true;
	return false;
}

std::vector<ctColor>&  NSGA2::getdiffcolors(const int&index)
{
	if (index == 0)
		return m_colorRGBDataP;
	else if (index == 1)
		return m_colorRGBDataLine;
	else if (index == 3)
		return m_colorRGBDataPoint;
	else if (index == 2)
		return m_colorRGBDataBack;
}

int NSGA2::getdiffub(const int & index)
{
	if (type2.size() == 0)
		return ubP;
	else if (type2[index] == 0)
		return ubP;
	else if (type2[index] == 1)
		return ubL;
	else if (type2[index] == 3)
		return ubPoint;
	else if (type2[index] == 2)
		return ubBack;

}



void NSGA2::setOptimizingWeights(const float& w1, const float& w2, const float& w3, const float& w4, const float& w5, const float& w6, const float& w7)
{
	wc = w1;
	wse = w2;   //语义权重
	wsi = w3;  //相似性权重
	wha = w4; //调和权重
	whi = w5; //视觉层次权重
	wseq = w6; //序列权重

	wdi = w7;

	ofstream ww("../output/权重.csv");
	ww << "习惯用色" << "," << "语义" << "," << "相似性" << "," << "调和" << "视觉层次,序列" << endl;
	ww << wc << "," << wse << "," << wsi << "," << wha << "," << whi << "," << wseq << endl;
	ww.close();
}

double NSGA2::rand_real(double low, double high)
//产生随机实数
{
	double h;
	h = (high - low)*URAND + low + 0.001;
	if (h >= high)
		h = high - 0.001;
	return h;
}

int NSGA2::rand_int(int low, int high)
//产生随机整数
{
	return int((high - low + 1)*URAND) + low;
}

int population::cmp1(const void *a, const void *b)
//目标函数f1的升序排序
{
	const individual *e = (const individual *)a;
	const individual *f = (const individual *)b;
	if (e->fvalue[0] == f->fvalue[0])
		return 0;
	else if (e->fvalue[0] < f->fvalue[0])
		return -1;
	else return 1;
}

int population::cmp2(const void *a, const void *b)
//目标函数f2的升序排序
{
	const individual *e = (const individual *)a;
	const individual *f = (const individual *)b;
	if (e->fvalue[1] == f->fvalue[1])
		return 0;
	else if (e->fvalue[1] < f->fvalue[1])
		return -1;
	else return 1;
}

int population::cmp_c_d(const void *a, const void *b)
//对拥挤距离降序排序
{
	const individual *e = (const individual *)a;
	const individual *f = (const individual *)b;
	if (e->crowding_distance == f->crowding_distance)
		return 0;
	else if (e->crowding_distance < f->crowding_distance)
		return 1;
	else
		return -1;
}

bool population::e_is_dominated(const individual &a, const individual &b)
{
	if (a.IsFeasible >mNSGA->FeaLimitation || b.IsFeasible > mNSGA->FeaLimitation) {   //照片55  绘画25
		if (a.IsFeasible < b.IsFeasible)
			return true;
		if (b.IsFeasible < a.IsFeasible)
			return false;
	}

	if ((a.fvalue[0] <= b.fvalue[0]) && (a.fvalue[1] <= b.fvalue[1]))
	{
		if ((a.fvalue[0] == b.fvalue[0]) && a.fvalue[1] == b.fvalue[1])
			return false;
		else
			return true;
	}
	else
		return false;
}

