#include "StyleExtract.h"
#include "../ColorTransform/pictureOp.h"
#include "DoH.h"
#include "BackgroundDetector.h"
#include <opencv2/calib3d.hpp>
#include <fstream>





void StyleExtract::iniParameter(std::string srcfilepath)
{
	src = new cv::Mat();
	edge = new cv::Mat();

	srcfile = srcfilepath;
	*src = cv::imread(srcfile);
	if (src->empty()) {
		return;
	}
	extractPStroke();
}


bool StyleExtract::extractBlob(std::string filepath)
{

	vector<cv::KeyPoint> keyPoints;
	cv::SimpleBlobDetector::Params params;
	//阈值控制
	//params.minThreshold = minBlob;
	//params.maxThreshold = maxBlob;
	params.filterByArea = true;
	params.minArea = 5*minBlob;
	params.maxArea = 5 * maxBlob;

	//Filter by Circularity
	params.filterByCircularity = false;
	params.minCircularity = 0.1;

		// Filter by Convexity
	params.filterByConvexity = false;
	params.minConvexity = 0.87;

	// Filter by Inertia
	params.filterByInertia = false;
	params.minInertiaRatio = 0.01;
	cv::Ptr<cv::SimpleBlobDetector> blobDetect= cv::SimpleBlobDetector::create(params);
	//blobDetect->create();
	blobDetect->detect(*src, keyPoints);

	string ImagefoldName = std::string(filepath).substr(0, std::string(filepath).find_last_of("\\"));
	if (!PictureOP::exists(ImagefoldName))
	{
		int ret = _mkdir(ImagefoldName.c_str());
		if (ret != 0)
			return false;
	}

	
	int count = 0;
	ofstream blob(filepath);
	for (auto b : keyPoints)
	{
		if (count > 70)
			break;
		blob << b.pt.y << "," << b.pt.x << "," << b.size << endl;
		count++;
	}

	blob.close();
	return true;
}

void StyleExtract::DOG2(const cv::Mat &src1, cv::Mat &dst, cv::Size wsize, double sigma, double k = 1.6) {
	cv::Mat gaussian_dst1, gaussian_dst2;
	//高斯滤波
	cv::GaussianBlur(src1, gaussian_dst1, wsize, k*sigma);

	cv::GaussianBlur(src1, gaussian_dst2, wsize, sigma);

	dst = (gaussian_dst1 - gaussian_dst2);
}

bool StyleExtract::extractPStroke()
{	
	DOG2(*src, *edge, cv::Size(lineWidth, lineWidth), 2);
	cvtColor(*edge,*edge, CV_RGB2GRAY);
	return true;
}

bool StyleExtract::extractStroke(std::string filepath)
{
	cv::Mat img_edge, labels, centroids, img_color, stats;

	cv::threshold(*edge, img_edge, 0, 255, cv::THRESH_OTSU);
	int nccomps = connectedComponentsWithStats(img_edge, labels, stats, centroids);
	std::vector<cv::Vec3b>colors(nccomps + 1);;
	colors[0] = cv::Vec3b(0, 0, 0);

	for (int i = 1; i < nccomps; i++)
	{
		colors[i] = cv::Vec3b(0, 0, 0);
		int area = stats.at<int>(i, cv::CC_STAT_AREA);//连通域的面积
		int height = stats.at<int>(i, cv::CC_STAT_HEIGHT);
		int width = stats.at<int>(i, cv::CC_STAT_WIDTH);
		float ratio = min(height, width) / max(height, width);//连通域的长宽比
		
		
		if (ratio<0.80 && area > pixnumStroke * 5)
			colors[i] = cv::Vec3b(255, 255, 255);
		else if (ratio<0.30 &&area>pixnumStroke * 3)
			colors[i] = cv::Vec3b(255, 255, 255);
		else if (ratio<0.10 && area>pixnumStroke)
			colors[i] = cv::Vec3b(255, 255, 255);
	
	}
	img_color = cv::Mat::zeros(edge->size(), CV_8UC3);
	for (int y = 0; y < img_color.rows; y++)
		for (int x = 0; x < img_color.cols; x++)
		{
			int label = labels.at<int>(y, x);
			CV_Assert(0 <= label && label <= nccomps);
			img_color.at<cv::Vec3b>(y, x) = colors[label];
		}

	string ImagefoldName = std::string(filepath).substr(0, std::string(filepath).find_last_of("\\"));
	if (!PictureOP::exists(ImagefoldName))
	{
		int ret = _mkdir(ImagefoldName.c_str());
		if (ret != 0)
			return false;
	}
	imwrite(filepath, img_color);
	return true;
}


bool StyleExtract::extractEdge2(std::string filepath)
{
	cv::Mat img_edge, labels, centroids, img_color, stats;

	cv::threshold(*edge, img_edge, 0, 255, cv::THRESH_OTSU);
	int nccomps = connectedComponentsWithStats(img_edge, labels, stats, centroids);
	std::vector<cv::Vec3b>colors(nccomps + 1);;
	colors[0] = cv::Vec3b(0, 0, 0);

	for (int i = 1; i < nccomps; i++)
	{
		colors[i] = cv::Vec3b(255, 255, 255);
		
		int area = stats.at<int>(i, cv::CC_STAT_AREA);//连通域的面积
		if (area<pixnumEdge2)
			colors[i] = cv::Vec3b(0, 0, 0);
	}

	
	img_color = cv::Mat::zeros(edge->size(), CV_8UC3);
	for (int y = 0; y < img_color.rows; y++)
		for (int x = 0; x < img_color.cols; x++)
		{
			int label = labels.at<int>(y, x);
			CV_Assert(0 <= label && label <= nccomps);
			img_color.at<cv::Vec3b>(y, x) = colors[label];
		}
	cv::Mat kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5, 5));
	dilate(img_color, img_color, kernel);

	string ImagefoldName = std::string(filepath).substr(0, std::string(filepath).find_last_of("\\"));
	if (!PictureOP::exists(ImagefoldName))
	{	
		int ret = _mkdir(ImagefoldName.c_str());
		if (ret != 0)
			return false;
	}

	imwrite(filepath, img_color);
	return true;
}

bool StyleExtract::extractBackground(std::string filepath)
{
	string ImagefoldName,backfile;
	PictureOP::getBackgroundPicture(filepath, backfile);
	ImagefoldName = std::string(backfile).substr(0, std::string(backfile).find_last_of("\\"));
	if (!PictureOP::exists(ImagefoldName))
	{
		int ret = _mkdir(ImagefoldName.c_str());
		if (ret != 0)
			return false;
	}
	BackgroundDetector::Backgrounddetect(filepath.c_str(), backfile);

	return true;
}

