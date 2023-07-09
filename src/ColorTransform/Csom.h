#ifndef CSOM_H
#define CSOM_H

//------------------------------------------------------------------------
//
//  Name:   Csom.h
//
//  Desc:   class to define a Kohonen Self Organising Feature Map
//
//  Author: Mat Buckland 2002 (fup@ai-junkie.com)
//
//------------------------------------------------------------------------

#include <windows.h>
#include <vector>



using namespace std;
#include "constants.h"
#include "CNode.h"



class MAPSTYLETRANSFORM_EXPORT Csom
{

private:

	

  //this holds the address of the winning node from the current iteration
  CNode*              m_pWinningNode; 

   //this is the topological 'radius' of the feature map
  double              m_dMapRadius;

   //used in the calculation of the neighbourhood width of influence
  double              m_dTimeConstant;

  //the number of training iterations
  int                 m_iNumIterations;
  int                m_icur=0;

  //keeps track of what iteration the epoch method has reached
  int                 m_iIterationCount;

  //the current width of the winning node's area of influence
  double              m_dNeighbourhoodRadius;

  //how much the learning rate is adjusted for nodes within
  //the area of influence
  double              m_dInfluence;

  double              m_dLearningRate;

  //set true when training is finished
  bool                m_bDone;

  //the height and width of the cells that the nodes occupy when 
  //rendered into 2D space.
  double              m_dCellWidth;
  double              m_dCellHeight;

  vector<vector<double>> iniColor;  //³õÊ¼ÑÕÉ«

  CNode*    FindBestMatchingNode(const vector<double> &vec);
  CNode*  FindBestMatchingNode(const vector<double> &vec, int&);

  inline    double Gaussian(const double dist, const double sigma);


public:

	//the neurons representing the Self Organizing Map
	vector<CNode>       m_SOM;
	int numCellsUp;
	int numCellsAcross;
  Csom():m_dCellWidth(0),
         m_dCellHeight(0),
         m_pWinningNode(NULL),
         m_iIterationCount(1),
         m_iNumIterations(0),
         m_dTimeConstant(0),
         m_dMapRadius(0),
         m_dNeighbourhoodRadius(0),
         m_dInfluence(0),
         m_dLearningRate(constStartLearningRate),
         m_bDone(false),
		 m_icur(0)
               
  {}

  void Create(int cxClient, int cyClient, int CellsUp, int CellsAcross, int NumIterations, vector<vector<double> > m_TrainingSet);

  void Render(HDC surface);

  bool Epoch(const vector<vector<double> > &data);

  bool Epoch(const vector<vector<double> > &data, vector<vector<int>>&);

  bool FinishedTraining()const{return m_bDone;}

  void getIniColor(vector<vector<double>>&, const vector<vector<double> >&);
  std::vector<int> Poi2Color(int x, int y);
  std::vector<int> Color2Poi(std::vector<int>&);
  std::vector<ctColor>  getColors();
  void setColors(std::vector<ctColor>);
};


#endif