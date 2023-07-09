#ifndef CCONTROLLER_H
#define CCONTROLLER_H


//------------------------------------------------------------------------
//
//  Name:   CController.h
//
//  Desc:   controller class for the SOM demo
//
//  Author: Mat Buckland 2002 (fup@ai-junkie.com)
//
//------------------------------------------------------------------------

#include <windows.h>
#include <vector>

using namespace std;

#include "CNode.h"
#include "Csom.h"
#include "constants.h"

class MAPSTYLETRANSFORM_EXPORT CController
{

private:


  //the data for the training
  vector<vector<double> > m_TrainingSet;
  vector<vector<double> > m_TrainingSetLine;
  vector<vector<double> > m_TrainingSetPoint;
  vector<vector<double> > m_TrainingSetBack;

  //this method creates a small data set of color values
  //that are used to train the network with
  void CreateDataSet() {};

  void CreateDataSet(const std::string&);
  void CreateDataSet(const std::string&,const std::string&);
  void CreateDataSet(const std::string&, const std::string&, const std::string&, const std::string&);
  void CreateDataSet(const std::vector<ctColor>&, const std::vector<ctColor>&, const std::vector<ctColor>&, const std::vector<ctColor>&);
  void getPictureLine(const std::string& linepath, int ** &picture);

  void colorAggregation(std::vector<ctColor>& iMapColor,const int& mlimitation);

  std::vector<ctColor>  calPersonalPictureProportion(const std::string& img_path);
  std::vector<ctColor>  calPersonalPictureProportion(const std::string& img_path, const std::string&, std::vector<ctColor>&);
  std::vector<ctColor>  calPersonalPictureProportion(const std::string& img_path, const std::string&, const std::string&, const std::string&, std::vector<ctColor>&, std::vector<ctColor>&, std::vector<ctColor>&);
  std::vector<ctColor> CController::calPersonalPictureProportion(const std::vector<ctColor>& iMapColor, const std::vector<ctColor>&  iMapColorL,
	  const std::vector<ctColor>&  iMapColorPoint, const std::vector<ctColor>&  iMapColorBack, std::vector<ctColor>& ctMapColorLine, std::vector<ctColor>& ctMapColorPoint, std::vector<ctColor>& ctMapColorBack);

  // CColorTransformAnalysis ¶Ô»°¿ò
  std::vector<ctColor> calPersonalPictureProportion(std::vector<int>& ctMapPixels, const int& w, const int& h, int comp);
  std::vector<ctColor> m_ctInutPicture;
  std::vector<ctColor> m_ctInutPictureLine;
  std::vector<ctColor> m_ctInutPicturePoint;
  std::vector<ctColor> m_ctInutPictureBack;

  std::string message;
  double m_dColorDis;
public:
	//pointer to a Self Organising Map
	Csom*   m_pSOM;
	Csom*   m_pSOMLine;
	Csom*   m_pSOMPoint;
	Csom*   m_pSOMBack;
  CController(int cxClient, int cyClient, int CellsUp, int CellsAcross, int NumIterations)
  {
    //create the SOM
    m_pSOM = new Csom();

	//create the training set
	CreateDataSet();

	m_pSOM->Create(cxClient, cyClient, CellsUp, CellsAcross, NumIterations, m_TrainingSet);

  }


  CController(const std::string& img_path, const double& dis, int cxClient, int cyClient, int CellsUp, int CellsAcross, int NumIterations)
  {
	  message = ""; m_dColorDis = dis;
	  //create the SOM
	  m_pSOM = new Csom();
	  //create the training set
	  CreateDataSet(img_path);

	  if (message == "")
	  {
		  m_pSOM->Create(cxClient, cyClient, CellsUp, CellsAcross, NumIterations, m_TrainingSet);
	  }
  }

  CController(const std::string& img_path, 
	  const std::string& linepath, const double& dis,
	  int cxClient, int cyClient, int CellsUp, int CellsAcross, int NumIterations)
  {
	  message = ""; m_dColorDis = dis;
	  //create the SOM
	  m_pSOM = new Csom();
	  m_pSOMLine = new Csom();
	  //create the training set
	  CreateDataSet(img_path, linepath);

	  if (message == "")
	  {
		  m_pSOM->Create(cxClient, cyClient, CellsUp, CellsAcross, NumIterations, m_TrainingSet);
		  m_pSOMLine->Create(cxClient, cyClient, constNumCellsDownLine, constNumCellsAcrossLine, NumIterations, m_TrainingSetLine);
	  }
  }
      

  CController(const std::string& img_path, const std::string& linepath,
	  const std::string& pointpath, const std::string& backpath,
	  const double& dis, int cxClient, int cyClient, int CellsUp, int CellsAcross, int NumIterations)
  {
	  message = ""; m_dColorDis = dis;
	  //create the SOM
	  m_pSOM = new Csom();
	  m_pSOMLine = new Csom();
	  m_pSOMPoint = new Csom();
	  //create the training set
	  CreateDataSet(img_path, linepath, pointpath, backpath);

	  if (message == "")
	  {
		  m_pSOM->Create(cxClient, cyClient, CellsUp, CellsAcross, NumIterations, m_TrainingSet);
		  m_pSOMLine->Create(cxClient, cyClient, constNumCellsDownLine, constNumCellsAcrossLine, NumIterations, m_TrainingSetLine);
		  m_pSOMPoint->Create(cxClient, cyClient, constNumCellsDownPoint, constNumCellsAcrossPoint, NumIterations, m_TrainingSetPoint);
	  }
  }

  CController(const std::vector<ctColor>& img_path, const std::vector<ctColor>& linepath,
	  const std::vector<ctColor>& pointpath, const std::vector<ctColor>& backpath,
	  const double& dis, int cxClient, int cyClient, int CellsUp, int CellsAcross, int NumIterations)
  {
	  message = ""; m_dColorDis = dis;
	  //create the SOM
	  m_pSOM = new Csom();
	  m_pSOMLine = new Csom();
	  m_pSOMPoint = new Csom();
	  m_pSOMBack = new Csom();
	  //create the training set
	  CreateDataSet(img_path, linepath,pointpath, backpath);

	  if (message == "")
	  {
		  m_pSOM->Create(cxClient, cyClient, CellsUp, CellsAcross, NumIterations, m_TrainingSet);
		  m_pSOMLine->Create(cxClient, cyClient, constNumCellsDownLine, constNumCellsAcrossLine, NumIterations, m_TrainingSetLine);
		  m_pSOMPoint->Create(cxClient, cyClient, constNumCellsDownPoint, constNumCellsAcrossPoint, NumIterations, m_TrainingSetPoint);
		  m_pSOMBack->Create(cxClient, cyClient, constNumCellsAcrossBack, constNumCellsAcrossBack, NumIterations, m_TrainingSetBack);
	  }
  }


  CController(const std::vector<ctColor>& img_path, const std::vector<ctColor>& linepath,
	  const vector<ctColor>& pointpath, const std::vector<ctColor>& backpath,
	 int cxClient, int cyClient, float picnum,int NumIterations)
  {
	  message = ""; 
	  //create the SOM
	  m_pSOM = new Csom();
	  m_pSOMLine = new Csom();
	  m_pSOMPoint = new Csom();
	  m_pSOMBack = new Csom();
	  //create the training set


	  CreateDataSet(img_path, linepath, pointpath, backpath);

	  if (message == "")
	  {
		  m_pSOM->Create(cxClient, cyClient, constNumCellsDown*picnum, constNumCellsAcross*picnum, NumIterations, m_TrainingSet);
		  m_pSOMLine->Create(cxClient, cyClient, constNumCellsDownLine*picnum, constNumCellsAcrossLine*picnum, NumIterations, m_TrainingSetLine);
		  m_pSOMPoint->Create(cxClient, cyClient, constNumCellsDownPoint*picnum, constNumCellsAcrossPoint*picnum, NumIterations, m_TrainingSetPoint);
		  m_pSOMBack->Create(cxClient, cyClient, constNumCellsAcrossBack*picnum, constNumCellsAcrossBack*picnum, NumIterations, m_TrainingSetBack);
	  }
  }


  ~CController()
  {
    delete m_pSOM;
  }
  
  std::vector<ctColor> getColor();
  std::vector<ctColor> getColorLine();
  std::vector<ctColor> getColorPoint();
  std::vector<ctColor> getColorBack();
  std::vector<int> Poi2Color(int x, int y);
  std::vector<int> Color2Poi(std::vector<int>& color);
  void Render(HDC surface);
  void Render(HDC surface, Csom* SOM);

  bool Train();
  bool Train2();
  bool Train3();
  bool Train4();
  bool Train(vector<vector<int>>&ijs);


  void CalPro();
  void CalPro(const vector<int>&);
  void ProximityMatch(vector<int>&, const vector<vector<double> >&);
  

  bool Finished()const{return m_pSOM->FinishedTraining();}

  bool Finished2()const { return m_pSOMLine->FinishedTraining(); }

  bool Finished3()const { return m_pSOMPoint->FinishedTraining(); }

  bool Finished4()const { return m_pSOMBack->FinishedTraining(); }

  void runABC();
  vector<vector<int>>ijs;
  vector<int> matchingindex;
  void matchingColor(Csom* ,const vector<int>& matchingindex, const vector<vector<double> >&);
  std::string GetMessageeX() { return message; }

};


#endif