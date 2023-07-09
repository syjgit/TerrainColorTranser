#pragma once
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <string>

#include "../base.h"

class MAPSTYLETRANSFORM_EXPORT StyleExtract
{
	cv::Mat* src;
	cv::Mat* edge;

	int pixnumStroke = 100;
	int pixnumEdge2 = 30;

	bool extractPStroke();

	void DOG2(const cv::Mat &, cv::Mat &dst, cv::Size wsize, double sigma, double k );
public:
	
	~StyleExtract() {
	};

	void iniParameter(std::string);
	
	std::string srcfile; //输入文件

	int lineWidth = 9; //线宽
	float minBlob = 10;
	float maxBlob = 300;

	bool extractBlob(std::string blobfile);
	bool extractStroke(std::string Strokefile);
	bool extractEdge2(std::string Fillfile);
	bool extractBackground(std::string Backgroundfile);
	bool useSP = true;

private:

};
