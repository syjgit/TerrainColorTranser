#pragma once

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "constants.h"

class MAPSTYLETRANSFORM_EXPORT colorHist
{
public:
	static double colorDis; //ÑÕÉ«ÈÝ²î

	static void compute_colorHist(const std::string& imgPath, std::vector<ctColor>& ctMapColor, const int h_bins);

	static void compute_colorHist(const std::string& imgPath, int**picture, std::vector<ctColor>&  , std::vector<ctColor> &, const int& h_bins, const int &);
	static void compute_colorHist(const std::string& imgPath, const std::string&, std::vector<ctColor>&, std::vector<ctColor> &, const int& h_bins, const int &);
	static void compute_colorHist(const std::string& imgPath, const std::string& linePath, const std::string& pointpath, const std::string& ,std::vector<ctColor>& ctMapColorP, std::vector<ctColor>& ctMapColorL, std::vector<ctColor>& ctMapColorPoint, std::vector<ctColor>&, const int &h_bins, const int &limitedV, const int &back);
	static void compute_colorHist2(const std::string& imgPath, const std::string& linePath, const std::string& linePath2, const std::string& pointpath, const std::string&, std::vector<ctColor>& ctMapColorP, std::vector<ctColor>& ctMapColorL, std::vector<ctColor>& ctMapColorPoint, std::vector<ctColor>&, const int &h_bins, const int &limitedV, const int &back);
	
	static void colorTransfer(ctColor& color0);
	static void colorTransfer(std::vector<ctColor>& colors);
	static void BrightnessAdjust(std::vector<ctColor>& Colors0, float backgroundadjust);
	
private:
	static void thinner(cv::Mat &image, cv::Mat &result, int );
	static void colorHist::cleanSmallArea(cv::Mat &sourceMat);
	static float getMeanL(const cv::Mat &Lab);
	static float getMeanV(const cv::Mat &Lab);
	void DOG2(cv::Mat &src, cv::Mat &dst, cv::Size wsize, double sigma, double k );
	static void getBlob(const std::string& pointpath, const cv::Mat& mat, cv::Mat& matPoint);
	static int Distance(const int&, const int&, const int&, const int&);
	static void getKmeansColorPoint(const int& num, const int& size, const cv::Mat & color, std::vector<ctColor>& ctMapColor);

	
	
	struct blob
	{
		int x;
		int y;
		int r;
	};

};



