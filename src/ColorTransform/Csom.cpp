#include "Csom.h"
#include<iostream>
#include<fstream>

using namespace std;

#include "ColorSapceTransfer.h"

void Csom::Create(int cxClient,
                  int cyClient,
                  int CellsUp,
                  int CellsAcross,
                  int NumIterations,
				  vector<vector<double> > m_TrainingSet)
{
	numCellsUp = CellsUp;
	numCellsAcross = CellsAcross;
  m_dCellWidth  = (double)cxClient / (double)CellsAcross;

  m_dCellHeight = (double)cyClient / (double)CellsUp;

  m_iNumIterations = NumIterations;
  getIniColor(iniColor, m_TrainingSet);
  //create all the nodes
  for (int row=0; row<CellsUp; ++row)
  {
    for (int col=0; col<CellsAcross; ++col)
    {
      m_SOM.push_back(CNode(col*m_dCellWidth,           //left
                            (col+1)*m_dCellWidth,       //right
                            row*m_dCellHeight,          //top
                            (row+1)*m_dCellHeight,      //bottom
							constSizeOfInputVector, iniColor[row*CellsAcross + col]));   //num weights
    }
  }

  //this is the topological 'radius' of the feature map
  m_dMapRadius = max(constWindowWidth, constWindowHeight)/2;

   //used in the calculation of the neighbourhood width of m_dInfluence
  m_dTimeConstant = m_iNumIterations/log(m_dMapRadius);
}  


void  Csom::getIniColor(vector<vector<double>>& Colors, const vector<vector<double> >&m_TrainingSet)
{
	srand(10); /*随便一个数字，只要是不变的*/
	int length = numCellsAcross *numCellsUp;
	vector<double> c;
	for (int i = 0; i <length; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			c.push_back(Random(0, 1));
		}
		c.push_back(Random(0, 1));

		Colors.push_back(c);
		vector<double>().swap(c);
	}

}


//---------------------------- Render ------------------------------------
//
//------------------------------------------------------------------------
void Csom::Render(HDC surface)
{
  //render all the cells
  for (int nd=0; nd<m_SOM.size(); ++nd)
  {
    m_SOM[nd].Render(surface);
    
  }

  SetBkMode(surface, TRANSPARENT);
  SetTextColor(surface, RGB(255,255,255));
  
 /* string s = "迭代: " + itos(m_iIterationCount);
 
  TextOut(surface, 5, 5, s.c_str() , s.size());*/

 // TextOut(surface, 100, 330, TEXT("Press 'R' to retrain"), strlen("Press 'R' to retrain"));

/* 	SOMDemo.exe!CController::Render(HDC__ * surface) 行 14	C++

  s = "Learning: " + ftos(m_dLearningRate);
  TextOut(surface, 5, 20, s.c_str(), s.size());

  s = "Radius: " + ftos(m_dNeighbourhoodRadius);
  TextOut(surface, 5, 40, s.c_str(), s.size());
*/
}
          

std::vector<int> Csom::Poi2Color(int x, int y)
{
	//for (int nd = 0; nd<m_SOM.size(); ++nd)
	{
	


		int cx = x / m_dCellWidth;
		int cy = y / m_dCellHeight;
		
		cx = cx < 0 ? 0 : cx; cy = cy < 0 ? 0 : cy;
		cx = cx > numCellsAcross - 1 ? (numCellsAcross - 1) : cx;
		cy = cy > numCellsAcross - 1 ? (numCellsAcross - 1) : cy;

		int index = cy * constNumCellsAcross + cx;
		CNode node = m_SOM[index];

		//if (node.m_iLeft <= x && node.m_iRight >= x && node.m_iTop <= x && node.m_iBottom >= y)
		{
			int red = (int)(node.m_dWeights[0] * 255);
			int green = (int)(node.m_dWeights[1] * 255);
			int blue = (int)(node.m_dWeights[2] * 255);
			return{ red ,green ,blue };
		}

	}
}

std::vector<int> Csom::Color2Poi(std::vector<int>& color)
{
	std::vector<double> res = ColorSpaceTransfer::RGB2Lab(color[0],color[1],color[2]);
	double MIN_COLOR_DIS = 100000;
	std::vector<int> retruns;

	for (int nd = 0; nd < m_SOM.size(); ++nd)
	{
		CNode node = m_SOM[nd];

		int red = (int)(node.m_dWeights[0] * 255);
		int green = (int)(node.m_dWeights[1] * 255);
		int blue = (int)(node.m_dWeights[2] * 255);

		std::vector<double> tmp = ColorSpaceTransfer::RGB2Lab(red, green, blue);

		//计算色差
		double color_dis = sqrt(pow((res[0] - tmp[0]), 2.0) + pow((res[1] - tmp[1]), 2.0) + pow((res[2] - tmp[2]), 2.0));
		if (color_dis < MIN_COLOR_DIS)
		{
			MIN_COLOR_DIS = color_dis;

			retruns = { (int)((node.m_iLeft + node.m_iRight)/2.0),(int)((node.m_iTop+ node.m_iBottom)/2.0) };
		}

	}

	return retruns;
}
      

void Csom::setColors(std::vector<ctColor> colors)
{


	for (int nd = 0; nd<m_SOM.size(); ++nd)
	{

		ctColor color;

		CNode& node = m_SOM[nd];
		{
			node.m_dWeights[0] = colors[nd].r / 255.0;
			node.m_dWeights[1] = colors[nd].g / 255.0;
			node.m_dWeights[2] = colors[nd].b / 255.0;
		}

		
	}
}

std::vector<ctColor>  Csom::getColors()
{
	std::vector<ctColor> results;
	for (int nd = 0; nd<m_SOM.size(); ++nd)
	{
		ctColor color;

		CNode node = m_SOM[nd];
		{
			int red = (int)(node.m_dWeights[0] * 255);
			int green = (int)(node.m_dWeights[1] * 255);
			int blue = (int)(node.m_dWeights[2] * 255);

			color.r = red;
			color.g = green;
			color.b = blue;
			
			color.selected = node.selected;
			color.proportion = node.m_dWeights[3];
			if (node.selected == false)
			{
				color.proportion = 0.0;
			}
		}

		results.emplace_back(color);
	}

	return results;
}











//--------------------------- Epoch --------------------------------------
//
//  Given a std::vector of input vectors this method choses one at random
//  and runs the network through one training epoch
//------------------------------------------------------------------------
bool Csom::Epoch(const vector<vector<double> > &data)
{
  //make sure the size of the input vector matches the size of each node's 
  //weight vector
  if (data[0].size() != constSizeOfInputVector) return false;

  //return if the training is complete
  if (m_bDone) return true;
 
  
  //enter the training loop
  if (--m_iNumIterations > 0)
  {

	  //the input vectors are presented to the network at random
	 int ThisVector = RandInt(0, data.size()-1);



    //present the vector to each node and determine the BMU
    m_pWinningNode = FindBestMatchingNode(data[ThisVector]);

    //calculate the width of the neighbourhood for this timestep
    m_dNeighbourhoodRadius = m_dMapRadius * exp(-(double)m_iIterationCount/m_dTimeConstant);
	//m_dNeighbourhoodRadius = 50.0;

    //Now to adjust the weight vector of the BMU and its
    //neighbours

    //For each node calculate the m_dInfluence (Theta from equation 6 in
    //the tutorial. If it is greater than zero adjust the node's weights
    //accordingly
    for (int n=0; n<m_SOM.size(); ++n)
    {
      //calculate the Euclidean distance (squared) to this node from the
      //BMU
      double DistToNodeSq = (m_pWinningNode->X()-m_SOM[n].X()) *
                            (m_pWinningNode->X()-m_SOM[n].X()) +
                            (m_pWinningNode->Y()-m_SOM[n].Y()) *
                            (m_pWinningNode->Y()-m_SOM[n].Y()) ;

      double WidthSq = m_dNeighbourhoodRadius * m_dNeighbourhoodRadius;
	  //double WidthSq = m_dNeighbourhoodRadius;

      //if within the neighbourhood adjust its weights
      if (DistToNodeSq < (m_dNeighbourhoodRadius * m_dNeighbourhoodRadius))
      {

        //calculate by how much its weights are adjusted
        m_dInfluence = exp(-(DistToNodeSq) / (2*WidthSq));
		 // m_dInfluence = exp(1.0/(-(DistToNodeSq) / (2 * WidthSq)));

        m_SOM[n].AdjustWeights(data[ThisVector],
                               m_dLearningRate,
                               m_dInfluence);
      }

    }//next node


    //reduce the learning rate
    m_dLearningRate = constStartLearningRate * exp(-(double)m_iIterationCount/m_iNumIterations);
    
    ++m_iIterationCount;

  }

  else
  {
    m_bDone = true;
	//return false;
	ofstream output;
	output.open("SOMResult.txt");

	vector<vector < double >> value;
	for (auto& weight : m_SOM)
	{
		value.emplace_back(weight.m_dWeights);
		//output << weight.m_dWeights[0] << " " << weight.m_dWeights[1] << " " << weight.m_dWeights[2] << " " << weight.m_dWeights[3] << std::endl;
	}
	output.close();
  }

  return true;
}


bool Csom::Epoch(const vector<vector<double> > &data, vector<vector<int>>&ijs)
{
	//make sure the size of the input vector matches the size of each node's 
	//weight vector
	if (data[0].size() != constSizeOfInputVector) return false;

	//return if the training is complete
	if (m_bDone) return true;


	//enter the training loop
	vector<int> ij(2);
	if (--m_iNumIterations > 0)
	{

		//the input vectors are presented to the network at random
		int ThisVector = RandInt(0, data.size() - 1);

		//int ThisVector = m_icur;
		++m_icur;

		ij[1] = ThisVector;
		//present the vector to each node and determine the BMU
		m_pWinningNode = FindBestMatchingNode(data[ThisVector], ij[0]);

		ijs.push_back(ij);

		//calculate the width of the neighbourhood for this timestep
		m_dNeighbourhoodRadius = m_dMapRadius * exp(-(double)m_iIterationCount / m_dTimeConstant);

		//Now to adjust the weight vector of the BMU and its
		//neighbours

		//For each node calculate the m_dInfluence (Theta from equation 6 in
		//the tutorial. If it is greater than zero adjust the node's weights
		//accordingly
		for (int n = 0; n<m_SOM.size(); ++n)
		{
			//calculate the Euclidean distance (squared) to this node from the
			//BMU
			double DistToNodeSq = (m_pWinningNode->X() - m_SOM[n].X()) *
				(m_pWinningNode->X() - m_SOM[n].X()) +
				(m_pWinningNode->Y() - m_SOM[n].Y()) *
				(m_pWinningNode->Y() - m_SOM[n].Y());

			double WidthSq = m_dNeighbourhoodRadius * m_dNeighbourhoodRadius;

			//if within the neighbourhood adjust its weights
			if (DistToNodeSq < (m_dNeighbourhoodRadius * m_dNeighbourhoodRadius))
			{

				//calculate by how much its weights are adjusted
				m_dInfluence = exp(-(DistToNodeSq) / (2 * WidthSq));

				m_SOM[n].AdjustWeights(data[ThisVector],
					m_dLearningRate,
					m_dInfluence);
			}

		}//next node


		 //reduce the learning rate
		m_dLearningRate = constStartLearningRate * exp(-(double)m_iIterationCount / m_iNumIterations);

		++m_iIterationCount;

	}

	else
	{
		m_bDone = true;
		return false;
		ofstream output;
		output.open("SOMResult.txt");

		vector<vector < double >> value;
		for (auto& weight : m_SOM)
		{
			value.emplace_back(weight.m_dWeights);
			output << weight.m_dWeights[0] << " " << weight.m_dWeights[1] << " " << weight.m_dWeights[2] << " " << weight.m_dWeights[3] << std::endl;
		}
		output.close();
	}

	return true;
}

CNode* Csom::FindBestMatchingNode(const vector<double> &vec, int & index)
{
	CNode* winner = NULL;

	double LowestDistance = 999999;
	for (int n = 0; n<m_SOM.size(); ++n)
	{


		double dist = m_SOM[n].CalculateDistance(vec);

		if (dist < LowestDistance)
		{
			LowestDistance = dist;

			winner = &m_SOM[n];
			index = n;
		}
	}

	return winner;
}

//--------------------- CalculateBestMatchingNode ------------------------
//
//  this function presents an input vector to each node in the network
//  and calculates the Euclidean distance between the vectors for each
//  node. It returns a pointer to the best performer
//------------------------------------------------------------------------
CNode* Csom::FindBestMatchingNode(const vector<double> &vec)
{
  CNode* winner = NULL;

  double LowestDistance = 999999;
 
  for (int n=0; n<m_SOM.size(); ++n)
  {
    double dist = m_SOM[n].CalculateDistance(vec);

    if (dist < LowestDistance)
    {
      LowestDistance = dist;

      winner = &m_SOM[n];
    }
  }

  return winner;
}

