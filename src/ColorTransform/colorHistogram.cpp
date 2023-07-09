#include "colorHistogram.h"
#include "ColorSapceTransfer.h"
#include <iostream>
#include <fstream>
#include <cassert>

double colorHist::colorDis; //颜色容差

float colorHist::getMeanL(const cv::Mat &mat)
{
	cv::Mat Lab;
	cvtColor(mat, Lab, cv::COLOR_BGR2Lab);
	
	CvScalar L_mean = mean(Lab);
	return L_mean.val[0];
}

float colorHist::getMeanV(const cv::Mat &mat)
{
	cv::Mat HSV;
	cvtColor(mat, HSV, cv::COLOR_BGR2HSV);

	CvScalar V_mean = mean(HSV);
	return V_mean.val[2];
}

void colorHist::cleanSmallArea(cv::Mat &sourceMat)
{
	cv::Mat grayImage;
	cv::Mat thresholdImage;
	vector< vector< cv::Point> > contours;
	vector< vector< cv::Point> > contours2; //用于保存面积不足100的轮廓
	vector<cv::Point> tempV;
	cv::cvtColor(sourceMat, grayImage, CV_BGR2GRAY);
	threshold(grayImage, thresholdImage, 0, 255, CV_THRESH_OTSU + CV_THRESH_BINARY);
	findContours(thresholdImage, contours, CV_RETR_LIST, CV_CHAIN_APPROX_NONE);
	//sort(contours.begin(), contours.end(), ascendSort);//升序排序
	std::sort(contours.begin(), contours.end(), [](const vector< cv::Point>& x, const vector< cv::Point>& y) {
		return x.size() < y.size();
	});
	vector<vector<cv::Point> >::iterator itc = contours.begin();
	int i = 0;
	while (itc != contours.end())
	{		
		if (itc->size() < 100)
		{			
			contours2.push_back(*itc);	
		}
		++itc;
	}
	//删除轮廓面积不足100的区域，即用黑色填充轮廓面积不足100的区域：
	cv::drawContours(sourceMat, contours2, -1, cv::Scalar(255, 255, 255), CV_FILLED);
}

void colorHist::compute_colorHist(const std::string& imgPath, std::vector<ctColor>& ctMapColor, const int h_bins) {
	cv::Mat mat = cv::imread(imgPath);

	int histSize[] = { h_bins, h_bins,h_bins };
	// hue varies from 0 to 256, saturation from 0 to 180
	float s_ranges[] = { 0, 256 };
	const float* ranges[] = { s_ranges, s_ranges,s_ranges };
	// Use the o-th and 1-st channels
	int channels[] = { 0,1,2 };
	/// Histograms
	cv::MatND color_features;
	/// Calculate the histograms for the HSV images
	cv::calcHist(&mat, 1, channels, cv::Mat(), color_features, 3, histSize, ranges,
		true, false);
	//cv::normalize(color_features, color_features, 0, 1, cv::NORM_MINMAX, -1, cv::Mat());
	ctColor c1;
	float colorD = 255.0 / h_bins;
	for (int i = 0; i < h_bins; i++) {
		for (int j = 0; j < h_bins; j++) {
			for (int z = 0; z < h_bins; z++)
			{
				if (color_features.at<float>(i, j, z) > 0)
				{
					c1.b = colorD*i;
					c1.g = colorD*j;
					c1.r = colorD*z;
					c1.proportion = color_features.at<float>(i, j, z) / (mat.rows*mat.cols);
					ctMapColor.emplace_back(c1);
				}

			}
		}
	}
}




void colorHist::compute_colorHist(const std::string& imgPath, int** picture, std::vector<ctColor>& ctMapColorP, std::vector<ctColor>& ctMapColorL, const int &h_bins,const int &limitedV) {
	cv::Mat mat = cv::imread(imgPath);
	cv::Mat matP;
	cv::Mat matL;
	int sizeP=0, sizeL=0;
	for (size_t i = 0; i <mat.rows ; i++)
	{
		for (size_t j = 0; j <mat.cols; j++)
		{

			if (picture[i][j] < limitedV)
			{
				sizeP++;
				matP.push_back(mat.at<cv::Vec3b>(i, j));

				//matP.at<cv::Vec3b>(i, j) = cv::Vec3b(255, 255, 255);
			}
				
			else
			{
				//matL.at<cv::Vec3b>(i, j) = cv::Vec3b(255, 255, 255);
				sizeL++;
				matL.push_back(mat.at<cv::Vec3b>(i, j));
			}
				
		}
	}


	{
		int histSize[] = { h_bins, h_bins,h_bins };
		// hue varies from 0 to 256, saturation from 0 to 180
		float s_ranges[] = { 0, 256 };
		const float* ranges[] = { s_ranges, s_ranges,s_ranges };
		// Use the o-th and 1-st channels
		int channels[] = { 0,1,2 };
		/// Histograms
		cv::MatND color_features;
		/// Calculate the histograms for the HSV images
		cv::calcHist(&matL, 1, channels, cv::Mat(), color_features, 3, histSize, ranges,
			true, false);
		//cv::normalize(color_features, color_features, 0, 1, cv::NORM_MINMAX, -1, cv::Mat());
		ctColor c1;
		float colorD = 255.0 / h_bins;
		for (int i = 0; i < h_bins; i++) {
			for (int j = 0; j < h_bins; j++) {
				for (int z = 0; z < h_bins; z++)
				{
					if (color_features.at<float>(i, j, z) > 0)
					{
						c1.b = colorD*i;
						c1.g = colorD*j;
						c1.r = colorD*z;
						c1.proportion = color_features.at<float>(i, j, z) / sizeL;
						ctMapColorL.emplace_back(c1);
					}

				}
			}
		}
	}
	{
		int histSize[] = { h_bins, h_bins,h_bins };
		// hue varies from 0 to 256, saturation from 0 to 180
		float s_ranges[] = { 0, 256 };
		const float* ranges[] = { s_ranges, s_ranges,s_ranges };
		// Use the o-th and 1-st channels
		int channels[] = { 0,1,2 };
		/// Histograms
		cv::MatND color_features;
		/// Calculate the histograms for the HSV images
		cv::calcHist(&matP, 1, channels, cv::Mat(), color_features, 3, histSize, ranges,
			true, false);
		//cv::normalize(color_features, color_features, 0, 1, cv::NORM_MINMAX, -1, cv::Mat());
		ctColor c1;
		float colorD = 255.0 / h_bins;
		for (int i = 0; i < h_bins; i++) {
			for (int j = 0; j < h_bins; j++) {
				for (int z = 0; z < h_bins; z++)
				{
					if (color_features.at<float>(i, j, z) > 0)
					{
						c1.b = colorD*i;
						c1.g = colorD*j;
						c1.r = colorD*z;
						c1.proportion = color_features.at<float>(i, j, z) / sizeP;
						ctMapColorP.emplace_back(c1);
					}

				}
			}
		}
	}
}

void colorHist::compute_colorHist(const std::string& imgPath, const std::string& linePath, std::vector<ctColor>& ctMapColorP, std::vector<ctColor>& ctMapColorL, const int &h_bins, const int &limitedV) {
	cv::Mat mat = cv::imread(imgPath);

	cv::Mat LineP = cv::imread(linePath,0);
	cv::Mat picture;
	cv::Mat binary;

	
	//图像转化HSV颜色空间图像
	

	float meanL = getMeanL(mat);
	float meanV = getMeanV(mat);

	/*medianBlur(LineP, LineP, 3);
	threshold(LineP, binary, 128, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
	threshold(binary, binary, 128, 1, cv::THRESH_BINARY);//这个不是一般的二值，而是根据阈值来设定灰度值
	thinner(binary, picture,-1);
	for (int i = 0; i<picture.rows; i++)
	{
		for (int j = 0; j<picture.cols; j++)
		{
			if (picture.at<uchar>(i, j) == 1)
				picture.at<uchar>(i, j) = 255;
		}
	}*/

	
	picture = LineP;
	cv::Mat matP;
	cv::Mat matL;
	int sizeP = 0, sizeL = 0;
	for (size_t i = 0; i <mat.rows; i++)
	{
		for (size_t j = 0; j <mat.cols; j++)
		{
			double bgr[3];
			
			bgr[0] = mat.at<cv::Vec3b>(i, j)[0];
			bgr[1] = mat.at<cv::Vec3b>(i, j)[1];
			bgr[2] = mat.at<cv::Vec3b>(i, j)[2];
			vector<double>lab = ColorSpaceTransfer::RGB2Lab(bgr[2], bgr[1], bgr[0]);

			//if (picture.at<uchar>(i, j) > limitedV &&lab[0] < 1.0 / 3 * meanL)
			if (picture.at<uchar>(i, j) > limitedV )
			{
				sizeL++;
				matL.push_back(mat.at<cv::Vec3b>(i, j));
			}
			else
			{
				sizeP++;
				matP.push_back(mat.at<cv::Vec3b>(i, j));
			}
		}
	}


	{
		int histSize[] = { h_bins, h_bins,h_bins };
		// hue varies from 0 to 256, saturation from 0 to 180
		float s_ranges[] = { 0, 256 };
		const float* ranges[] = { s_ranges, s_ranges,s_ranges };
		// Use the o-th and 1-st channels
		int channels[] = { 0,1,2 };
		/// Histograms
		cv::MatND color_features;
		/// Calculate the histograms for the HSV images
		cv::calcHist(&matL, 1, channels, cv::Mat(), color_features, 3, histSize, ranges,
			true, false);
		//cv::normalize(color_features, color_features, 0, 1, cv::NORM_MINMAX, -1, cv::Mat());
		ctColor c1;
		float colorD = 255.0 / h_bins;
		for (int i = 0; i < h_bins; i++) {
			for (int j = 0; j < h_bins; j++) {
				for (int z = 0; z < h_bins; z++)
				{
					if (color_features.at<float>(i, j, z) > 0)
					{
						c1.b = colorD*i;
						c1.g = colorD*j;
						c1.r = colorD*z;
						c1.proportion = color_features.at<float>(i, j, z) / sizeL;
						ctMapColorL.emplace_back(c1);
					}

				}
			}
		}
	}
	{
		int histSize[] = { h_bins, h_bins,h_bins };
		// hue varies from 0 to 256, saturation from 0 to 180
		float s_ranges[] = { 0, 256 };
		const float* ranges[] = { s_ranges, s_ranges,s_ranges };
		// Use the o-th and 1-st channels
		int channels[] = { 0,1,2 };
		/// Histograms
		cv::MatND color_features;
		/// Calculate the histograms for the HSV images
		cv::calcHist(&matP, 1, channels, cv::Mat(), color_features, 3, histSize, ranges,
			true, false);
		//cv::normalize(color_features, color_features, 0, 1, cv::NORM_MINMAX, -1, cv::Mat());
		ctColor c1;
		float colorD = 255.0 / h_bins;
		for (int i = 0; i < h_bins; i++) {
			for (int j = 0; j < h_bins; j++) {
				for (int z = 0; z < h_bins; z++)
				{
					if (color_features.at<float>(i, j, z) > 0)
					{
						c1.b = colorD*i;
						c1.g = colorD*j;
						c1.r = colorD*z;
						c1.proportion = color_features.at<float>(i, j, z) / sizeP;
						ctMapColorP.emplace_back(c1);
					}

				}
			}
		}
	}
}

void colorHist::DOG2(cv::Mat &src, cv::Mat &dst, cv::Size wsize, double sigma, double k = 1.6) {
	cv::Mat gaussian_dst1, gaussian_dst2;
	//高斯滤波
	cv::GaussianBlur(src, gaussian_dst1, wsize, k*sigma);

	cv::GaussianBlur(src, gaussian_dst2, wsize, sigma);

	dst = gaussian_dst1 - gaussian_dst2;
	cv::threshold(dst, dst, 0, 255, cv::THRESH_BINARY);
}

int colorHist::Distance(const int& x1, const int& y1, const int& x2, const int& y2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

void colorHist::getBlob(const std::string& pointpath, const cv::Mat& mat, cv::Mat& matPoint)
{
	std::vector<blob> blobs;
	
	std::ifstream inFile(pointpath, std::ios::in);//inFile来自fstream,ifstream为输入文件流(从文件读入)
	std::string lineStr;
	std::vector<std::vector<std::string>> strArray;
	while (getline(inFile, lineStr)) //getline来自sstream
	{
		//打印整行字符串
		//cout<<lineStr<<endl;
		//存成二维表结构
		std::stringstream ss(lineStr);//来自sstream
		std::string str;
		std::vector<std::string> lineArray;
		//按照逗号分隔
		while (getline(ss, str, ','))
			lineArray.push_back(str);//一行数据以vector保存
									 //cout<<lineArray[0]<<endl;
		blob b;
		b.x=atoi(lineArray[0].c_str());
		b.y = atoi(lineArray[1].c_str());
		b.r= atoi(lineArray[2].c_str());
		blobs.emplace_back(b);
		
	}

	for (size_t i = 0; i < mat.rows; i++)
	{
		for (size_t j = 0; j < mat.cols; j++)
		{
			for (blob b : blobs)
			{
				if (Distance(i, j, b.x,b.y) < b.r)
				{
					matPoint.at<cv::Vec3b>(i, j)[0]=mat.at<cv::Vec3b>(i, j)[0];
					matPoint.at<cv::Vec3b>(i, j)[1] = mat.at<cv::Vec3b>(i, j)[1];
					matPoint.at<cv::Vec3b>(i, j)[2] = mat.at<cv::Vec3b>(i, j)[2];
				}
			}
			
		}
	}
}


void colorHist::compute_colorHist(const std::string& imgPath, const std::string& linePath, 
	const std::string& pointpath, const std::string& backpath, std::vector<ctColor>& ctMapColorP, std::vector<ctColor>& ctMapColorL, 
	std::vector<ctColor>& ctMapColorPoint, std::vector<ctColor>& ctMapColorBack, const int &h_bins, const int &limitedV, const int &backlimit) {
	cv::Mat mat = cv::imread(imgPath);

	cv::Mat matback = cv::imread(backpath, CV_LOAD_IMAGE_ANYDEPTH);
	cv::Mat matback2;
	cv::Mat matback3(mat.rows, mat.cols, CV_8UC3, cv::Scalar(255,255,255));
	
	//imreadLinePoint();
	
	cv::Mat matPoint(mat.rows, mat.cols,CV_8UC3,cv::Scalar(0, 0, 0));
	getBlob(pointpath, mat, matPoint);
	cv::Mat matPoint2;
	cv::Mat matPoint3(mat.rows, mat.cols, CV_8UC3, cv::Scalar(255, 255,255));

	cv::Mat LineP = cv::imread(linePath, 0);


	cv::Mat picture;
	cv::Mat binary;


	//图像转化HSV颜色空间图像


	float meanL = getMeanL(mat);
	float meanV = getMeanV(mat);

	/*medianBlur(LineP, LineP, 3);
	threshold(LineP, binary, 128, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
	threshold(binary, binary, 128, 1, cv::THRESH_BINARY);//这个不是一般的二值，而是根据阈值来设定灰度值
	thinner(binary, picture,-1);
	for (int i = 0; i<picture.rows; i++)
	{
	for (int j = 0; j<picture.cols; j++)
	{
	if (picture.at<uchar>(i, j) == 1)
	picture.at<uchar>(i, j) = 255;
	}
	}*/


	picture = LineP;
	cv::Mat matP;
	cv::Mat matL;
	cv::Mat matP2(mat.rows, mat.cols, CV_8UC3, cv::Scalar(255, 255, 255));
	cv::Mat matL2(mat.rows, mat.cols, CV_8UC3, cv::Scalar(255, 255, 255));
	int sizeP = 0, sizeL = 0,sizePoint=0,sizeBack=0;
	for (size_t i = 0; i <mat.rows; i++)
	{
		for (size_t j = 0; j <mat.cols; j++)
		{		
			if (matback.at<uchar>(i, j)<= backlimit)
			{
				matback2.push_back(mat.at<cv::Vec3b>(i, j));
				matback3.at<cv::Vec3b>(i, j)[0] = mat.at<cv::Vec3b>(i, j)[0];
				matback3.at<cv::Vec3b>(i, j)[1] = mat.at<cv::Vec3b>(i, j)[1];
				matback3.at<cv::Vec3b>(i, j)[2] = mat.at<cv::Vec3b>(i, j)[2];
				sizeBack++;
				continue;
			}

			if (matPoint.at<cv::Vec3b>(i, j) != cv::Vec3b(0, 0, 0))
			{
				matPoint2.push_back(mat.at<cv::Vec3b>(i, j));
				matPoint3.at<cv::Vec3b>(i, j)[0] = mat.at<cv::Vec3b>(i, j)[0];
				matPoint3.at<cv::Vec3b>(i, j)[1] = mat.at<cv::Vec3b>(i, j)[1];
				matPoint3.at<cv::Vec3b>(i, j)[2] = mat.at<cv::Vec3b>(i, j)[2];
				sizePoint++;
				continue;
			}
			else
			{
				/*double bgr[3];

				bgr[0] = mat.at<cv::Vec3b>(i, j)[0];
				bgr[1] = mat.at<cv::Vec3b>(i, j)[1];
				bgr[2] = mat.at<cv::Vec3b>(i, j)[2];
				vector<double>lab = ColorSpaceTransfer::RGB2Lab(bgr[2], bgr[1], bgr[0]);
				//limitedV 灰度值
				//if (picture.at<uchar>(i, j) > limitedV &&lab[0] < 1.0 / 3 * meanL)*/
				if (picture.at<uchar>(i, j) > limitedV)
				{
					sizeL++;
					matL.push_back(mat.at<cv::Vec3b>(i, j));
					matL2.at<cv::Vec3b>(i, j)[0] = mat.at<cv::Vec3b>(i, j)[0];
					matL2.at<cv::Vec3b>(i, j)[1] = mat.at<cv::Vec3b>(i, j)[1];
					matL2.at<cv::Vec3b>(i, j)[2] = mat.at<cv::Vec3b>(i, j)[2];
				}
				else
				{
					sizeP++;
					matP.push_back(mat.at<cv::Vec3b>(i, j));
					matP2.at<cv::Vec3b>(i, j)[0] = mat.at<cv::Vec3b>(i, j)[0];
					matP2.at<cv::Vec3b>(i, j)[1] = mat.at<cv::Vec3b>(i, j)[1];
					matP2.at<cv::Vec3b>(i, j)[2] = mat.at<cv::Vec3b>(i, j)[2];
				}
			}
		}
	}

	cv::imwrite("back.jpg", matback3);
	cv::imwrite("point0.jpg", matPoint);
	cv::imwrite("point.jpg", matPoint3);
	cv::imwrite("line.jpg", matL2);
	cv::imwrite("polygon.jpg", matP2);

	//计算背景
	{
		int numCluster = constNumCellsAcrossBack*constNumCellsAcrossBack;
		cv::Mat labels;
		cv::Mat centers;
		cv::Mat PointsBack(sizeBack, 3, CV_32F, cv::Scalar(10));
		for (size_t i = 0; i < sizeBack; i++)
		{
			cv::Vec3b bgr = matback2.at<cv::Vec3b>(i, 0);
			PointsBack.at<float>(i, 0) = static_cast<int>(bgr[0]);
			PointsBack.at<float>(i, 1) = static_cast<int>(bgr[1]);
			PointsBack.at<float>(i, 2) = static_cast<int>(bgr[2]);
		}
		cv::kmeans(PointsBack, numCluster, labels, cv::TermCriteria(cv::TermCriteria::EPS + cv::TermCriteria::COUNT, 10, 0.1), 3, cv::KMEANS_PP_CENTERS, centers);
		ctColor c1;
		std::vector<float> backP(numCluster);
		for (int i = 0; i < sizeBack; i++)
		{
			int index = labels.at<int>(i);
			backP[index]++;
		}
		for (int i = 0; i < numCluster; i++)
		{
			c1.b = centers.at<float>(i, 0);
			c1.g = centers.at<float>(i, 1);
			c1.r = centers.at<float>(i, 2);
			c1.proportion = backP[i] / sizeBack;
			ctMapColorBack.emplace_back(c1);
		}
	}
	
	


	//计算点
	{
		int histSize[] = { h_bins, h_bins,h_bins };
		// hue varies from 0 to 256, saturation from 0 to 180
		float s_ranges[] = { 0, 256 };
		const float* ranges[] = { s_ranges, s_ranges,s_ranges };
		// Use the o-th and 1-st channels
		int channels[] = { 0,1,2 };
		/// Histograms
		cv::MatND color_features;
		/// Calculate the histograms for the HSV images
		cv::calcHist(&matPoint, 1, channels, cv::Mat(), color_features, 3, histSize, ranges,
			true, false);
		//cv::normalize(color_features, color_features, 0, 1, cv::NORM_MINMAX, -1, cv::Mat());
		ctColor c1;
		float colorD = 255.0 / h_bins;
		for (int i = 0; i < h_bins; i++) {
			for (int j = 0; j < h_bins; j++) {
				for (int z = 0; z < h_bins; z++)
				{
					if (color_features.at<float>(i, j, z) > 0)
					{
						c1.b = colorD*i;
						c1.g = colorD*j;
						c1.r = colorD*z;
						c1.proportion = color_features.at<float>(i, j, z) / sizePoint;
						ctMapColorPoint.emplace_back(c1);
					}

				}
			}
		}
	}

	//计算线
	{
		int histSize[] = { h_bins, h_bins,h_bins };
		// hue varies from 0 to 256, saturation from 0 to 180
		float s_ranges[] = { 0, 256 };
		const float* ranges[] = { s_ranges, s_ranges,s_ranges };
		// Use the o-th and 1-st channels
		int channels[] = { 0,1,2 };
		/// Histograms
		cv::MatND color_features;
		/// Calculate the histograms for the HSV images
		cv::calcHist(&matL, 1, channels, cv::Mat(), color_features, 3, histSize, ranges,
			true, false);
		//cv::normalize(color_features, color_features, 0, 1, cv::NORM_MINMAX, -1, cv::Mat());
		ctColor c1;
		float colorD = 255.0 / h_bins;
		for (int i = 0; i < h_bins; i++) {
			for (int j = 0; j < h_bins; j++) {
				for (int z = 0; z < h_bins; z++)
				{
					if (color_features.at<float>(i, j, z) > 0)
					{
						c1.b = colorD*i;
						c1.g = colorD*j;
						c1.r = colorD*z;
						c1.proportion = color_features.at<float>(i, j, z) / sizeL;
						ctMapColorL.emplace_back(c1);
					}

				}
			}
		}
	}

	//计算面
	{
		int histSize[] = { h_bins, h_bins,h_bins };
		// hue varies from 0 to 256, saturation from 0 to 180
		float s_ranges[] = { 0, 256 };
		const float* ranges[] = { s_ranges, s_ranges,s_ranges };
		// Use the o-th and 1-st channels
		int channels[] = { 0,1,2 };
		/// Histograms
		cv::MatND color_features;
		/// Calculate the histograms for the HSV images
		cv::calcHist(&matP, 1, channels, cv::Mat(), color_features, 3, histSize, ranges,
			true, false);
		//cv::normalize(color_features, color_features, 0, 1, cv::NORM_MINMAX, -1, cv::Mat());
		ctColor c1;
		float colorD = 255.0 / h_bins;
		for (int i = 0; i < h_bins; i++) {
			for (int j = 0; j < h_bins; j++) {
				for (int z = 0; z < h_bins; z++)
				{
					if (color_features.at<float>(i, j, z) > 0)
					{
						c1.b = colorD*i;
						c1.g = colorD*j;
						c1.r = colorD*z;
						c1.proportion = color_features.at<float>(i, j, z) / sizeP;
						ctMapColorP.emplace_back(c1);
					}

				}
			}
		}
	}
}

void colorHist::BrightnessAdjust(vector<ctColor>& Colors0,float backgroundadjust)
{
	std::sort(Colors0.begin(), Colors0.end(), [](const ctColor& x, const ctColor& y) {
		return x.L > y.L;
	});

	int num = Colors0.size()*0.2;
	int index = 0;
	double maxL = Colors0[0].L;

	for (auto &c : Colors0)
	{
		//c.proportion = c.proportion*c.L / maxL* backgroundadjust;
		c.L += (maxL - c.L)*backgroundadjust;
		auto rgb = ColorSpaceTransfer::LabToRgb(c.L, c.A, c.B);
		c.r = rgb[0];
		c.g = rgb[1];
		c.b = rgb[2];
		auto hsv = ColorSpaceTransfer::RGB2HSV(c.r, c.g, c.b);
		c.H = hsv[0];
		c.S = hsv[1];
		c.V = hsv[2];
		index++;
	}
}

void colorHist::colorTransfer(std::vector<ctColor>& colors)
{
	for (auto &color0 : colors)
	{
		auto lab = ColorSpaceTransfer::RGB2Lab(color0.r, color0.g, color0.b);
		auto hsv = ColorSpaceTransfer::RGB2HSV(color0.r, color0.g, color0.b);
		color0.L = lab[0];
		color0.A = lab[1];
		color0.B = lab[2];
		color0.H = hsv[0];
		color0.S = hsv[1];
		color0.V = hsv[2];
	}
}

void colorHist::colorTransfer(ctColor& color0)
{
	auto lab = ColorSpaceTransfer::RGB2Lab(color0.r, color0.g, color0.b);
	auto hsv = ColorSpaceTransfer::RGB2HSV(color0.r, color0.g, color0.b);
	color0.L = lab[0];
	color0.A = lab[1];
	color0.B = lab[2];
	color0.H = hsv[0];
	color0.S = hsv[1];
	color0.V = hsv[2];
}



void colorHist::getKmeansColorPoint(const int& numCluster, const int& size, const cv::Mat & matcolor, std::vector<ctColor>& ctMapColor)
{
	cv::Mat labels;
	cv::Mat centers;
	cv::Mat PointsBack(size, 3, CV_32F, cv::Scalar(10));
	for (size_t i = 0; i < size; i++)
	{
		cv::Vec3b bgr = matcolor.at<cv::Vec3b>(i, 0);
		PointsBack.at<float>(i, 0) = static_cast<int>(bgr[0]);
		PointsBack.at<float>(i, 1) = static_cast<int>(bgr[1]);
		PointsBack.at<float>(i, 2) = static_cast<int>(bgr[2]);
	}
	cv::kmeans(PointsBack, numCluster, labels, cv::TermCriteria(cv::TermCriteria::EPS + cv::TermCriteria::COUNT, 10, 0.1), 3, cv::KMEANS_PP_CENTERS, centers);
	ctColor c1;
	std::vector<float> backPoint(numCluster);
	for (int i = 0; i < size; i++)
	{
		int index = labels.at<int>(i);
		backPoint[index]++;
	}

	vector<double> tempP(numCluster);
	ctMapColor.resize(numCluster);
	for (int i = 0; i < numCluster; i++)
	{
		c1.b = centers.at<float>(i, 0);
		c1.g = centers.at<float>(i, 1);
		c1.r = centers.at<float>(i, 2);
		c1.proportion = backPoint[i] / size;
		ctMapColor[i] = (c1);
		tempP[i] = c1.proportion;

	}

	for (int i = 0; i < numCluster - 1; i++)
	{
		vector<double> lab1 = ColorSpaceTransfer::RGB2Lab(ctMapColor[i].r, ctMapColor[i].g, ctMapColor[i].b);
		for (int j = i + 1; j < numCluster; j++)
		{
			vector<double> lab2 = ColorSpaceTransfer::RGB2Lab(ctMapColor[j].r, ctMapColor[j].g, ctMapColor[j].b);
			double dis = sqrt(pow(lab1[0] - lab2[0], 2) + pow(lab1[1] - lab2[1], 2) + pow(lab1[2] - lab2[2], 2));
			if (dis< colorDis)
			{
				ctMapColor[i].proportion += tempP[j];
				ctMapColor[j].proportion += tempP[i];
			}
		}
	}
}

void colorHist::compute_colorHist2(const std::string& imgPath, const std::string& linePath, const std::string& linePath2,
	const std::string& pointpath, const std::string& backpath, std::vector<ctColor>& ctMapColorP, std::vector<ctColor>& ctMapColorL,
	std::vector<ctColor>& ctMapColorPoint, std::vector<ctColor>& ctMapColorBack, const int &h_bins, const int &limitedV, const int &backlimit) {
	cv::Mat mat = cv::imread(imgPath);

	cv::Mat matback = cv::imread(backpath, CV_LOAD_IMAGE_ANYDEPTH);
	cv::Mat matback2;
	cv::Mat matback3(mat.rows, mat.cols, CV_8UC3, cv::Scalar(255, 255, 255));

	//imreadLinePoint();

	cv::Mat matPoint(mat.rows, mat.cols, CV_8UC3, cv::Scalar(255,255,255));
	getBlob(pointpath, mat, matPoint);
	cv::Mat matPoint2;
	cv::Mat matPoint3(mat.rows, mat.cols, CV_8UC3, cv::Scalar(255, 255, 255));

	cv::Mat LineP = cv::imread(linePath, 0);

	cv::Mat LineP2 = cv::imread(linePath2, 0);

	cv::Mat picture;
	cv::Mat binary;


	//图像转化HSV颜色空间图像


	float meanL = getMeanL(mat);
	float meanV = getMeanV(mat);

	
	picture = LineP;
	cv::Mat matP;
	cv::Mat matL;
	cv::Mat matP2(mat.rows, mat.cols, CV_8UC3, cv::Scalar(255, 255, 255));  //平面图像
	cv::Mat matL2(mat.rows, mat.cols, CV_8UC3, cv::Scalar(255, 255, 255));
	int sizeP = 0, sizeL = 0, sizePoint = 0, sizeBack = 0;
	for (size_t i = 0; i <mat.rows; i++)
	{
		for (size_t j = 0; j <mat.cols; j++)
		{
			if (matPoint.at<cv::Vec3b>(i, j) != cv::Vec3b(255, 255, 255))
			{
				matPoint2.push_back(mat.at<cv::Vec3b>(i, j));
				matPoint3.at<cv::Vec3b>(i, j)[0] = mat.at<cv::Vec3b>(i, j)[0];
				matPoint3.at<cv::Vec3b>(i, j)[1] = mat.at<cv::Vec3b>(i, j)[1];
				matPoint3.at<cv::Vec3b>(i, j)[2] = mat.at<cv::Vec3b>(i, j)[2];
				sizePoint++;
				continue;
			}
			else if (picture.at<uchar>(i, j) > limitedV)
			{
				sizeL++;
				matL.push_back(mat.at<cv::Vec3b>(i, j));
				matL2.at<cv::Vec3b>(i, j)[0] = mat.at<cv::Vec3b>(i, j)[0];
				matL2.at<cv::Vec3b>(i, j)[1] = mat.at<cv::Vec3b>(i, j)[1];
				matL2.at<cv::Vec3b>(i, j)[2] = mat.at<cv::Vec3b>(i, j)[2];
			}
			else if (matback.at<uchar>(i, j) <= backlimit)
			{
				matback2.push_back(mat.at<cv::Vec3b>(i, j));
				matback3.at<cv::Vec3b>(i, j)[0] = mat.at<cv::Vec3b>(i, j)[0];
				matback3.at<cv::Vec3b>(i, j)[1] = mat.at<cv::Vec3b>(i, j)[1];
				matback3.at<cv::Vec3b>(i, j)[2] = mat.at<cv::Vec3b>(i, j)[2];
				sizeBack++;
				continue;
			}
			else
			{
				if (!LineP2.empty())
					if (LineP2.at<uchar>(i, j) > 1)
						continue;
				
				matP2.at<cv::Vec3b>(i, j)[0] = mat.at<cv::Vec3b>(i, j)[0];
				matP2.at<cv::Vec3b>(i, j)[1] = mat.at<cv::Vec3b>(i, j)[1];
				matP2.at<cv::Vec3b>(i, j)[2] = mat.at<cv::Vec3b>(i, j)[2];
			}	
		}
	}
	
	//cleanSmallArea(matP2);
	for (size_t i = 0; i <mat.rows; i++)
	{
		for (size_t j = 0; j < mat.cols; j++)
		{
			if (matP2.at<cv::Vec3b>(i, j) != cv::Vec3b(255, 255, 255))
			{
				sizeP++;
				matP.push_back(matP2.at<cv::Vec3b>(i, j));
			}
		}
	}

	cv::imwrite("../output/back.jpg", matback3);
	//cv::imwrite("point0.jpg", matPoint);
	cv::imwrite("../output/point.jpg", matPoint3);
	cv::imwrite("../output/line.jpg", matL2);
	cv::imwrite("../output/polygon.jpg", matP2);

	getKmeansColorPoint(constNumCellsAcrossBack*constNumCellsAcrossBack, sizeBack, matback2, ctMapColorBack);

	getKmeansColorPoint(constNumCellsAcrossPoint*constNumCellsAcrossPoint, sizePoint, matPoint2, ctMapColorPoint);

	getKmeansColorPoint(constNumCellsAcrossLine*constNumCellsAcrossLine, sizeL, matL, ctMapColorL);

	getKmeansColorPoint(constNumCellsAcross*constNumCellsAcross, sizeP, matP, ctMapColorP);

}


void colorHist::thinner(cv::Mat &image, cv::Mat &result, int maxitrator = -1)
{
	image.copyTo(result);
	int count = 0;//记录迭代次数
	while (1)
	{
		count++;
		if (maxitrator != -1 && count>maxitrator)
			break;
		std::vector<cv::Point>flagPoint;//用于标记要删除的点//vector<pair<int,int>>flagPoint
							   //进行点标记
		for (int i = 0; i<image.rows; i++)
		{
			for (int j = 0; j<image.cols; j++)
			{
				//满足以下四个条件，进行标记
				//p9 p2 p3
				//p8 p1 p4
				//p7 p6 p5
				int p1 = result.at<uchar>(i, j);
				int p2 = (i == 0) ? 0 : result.at<uchar>(i - 1, j);
				int p3 = (i == 0 || j == result.cols - 1) ? 0 : result.at<uchar>(i - 1, j + 1);
				int p4 = (j == result.cols - 1) ? 0 : result.at<uchar>(i, j + 1);
				int p5 = (i == result.rows - 1 || j == image.cols - 1) ? 0 : result.at<uchar>(i + 1, j + 1);
				int p6 = (i == result.rows - 1) ? 0 : result.at<uchar>(i + 1, j);
				int p7 = (i == result.rows - 1 || j == 0) ? 0 : result.at<uchar>(i + 1, j - 1);
				int p8 = (j == 0) ? 0 : result.at<uchar>(i, j - 1);
				int p9 = (i == 0 || j == 0) ? 0 : result.at<uchar>(i - 1, j - 1);
				if ((p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) >= 2 && (p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) <= 6)
				{
					int ap = 0;
					if (p2 == 0 && p3 == 1) ++ap;
					if (p3 == 0 && p4 == 1) ++ap;
					if (p4 == 0 && p5 == 1) ++ap;
					if (p5 == 0 && p6 == 1) ++ap;
					if (p6 == 0 && p7 == 1) ++ap;
					if (p7 == 0 && p8 == 1) ++ap;
					if (p8 == 0 && p9 == 1) ++ap;
					if (p9 == 0 && p2 == 1) ++ap;
					if (ap == 1)
					{
						if (p2*p4*p6 == 0)
						{
							if (p4*p6*p8 == 0)
								//记入标记点
								flagPoint.push_back(cv::Point(i, j));
						}
					}
				}
			}
		}
		//将标记点删除
		//vector<Vec2b>::const_iterator it=flagPoint.begin();
		for (std::vector<cv::Point>::const_iterator it = flagPoint.begin(); it != flagPoint.end(); it++)
			result.at<uchar>((*it).x, (*it).y) = 0;
		//直到没有点满足，算法结束
		if (flagPoint.size() == 0)
			break;
		else
			flagPoint.clear();
		//再对点进行标记
		for (int i = 0; i<image.rows; i++)
		{
			for (int j = 0; j<image.cols; j++)
			{
				//满足以下四个条件，进行标记
				//p9 p2 p3
				//p8 p1 p4
				//p7 p6 p5
				int p1 = result.at<uchar>(i, j);
				if (p1 != 1) continue;
				int p2 = (i == 0) ? 0 : result.at<uchar>(i - 1, j);
				int p3 = (i == 0 || j == result.cols - 1) ? 0 : result.at<uchar>(i - 1, j + 1);
				int p4 = (j == result.cols - 1) ? 0 : result.at<uchar>(i, j + 1);
				int p5 = (i == result.rows - 1 || j == image.cols - 1) ? 0 : result.at<uchar>(i + 1, j + 1);
				int p6 = (i == result.rows - 1) ? 0 : result.at<uchar>(i + 1, j);
				int p7 = (i == result.rows - 1 || j == 0) ? 0 : result.at<uchar>(i + 1, j - 1);
				int p8 = (j == 0) ? 0 : result.at<uchar>(i, j - 1);
				int p9 = (i == 0 || j == 0) ? 0 : result.at<uchar>(i - 1, j - 1);
				if ((p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) >= 2 && (p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) <= 6)
				{
					int ap = 0;
					if (p2 == 0 && p3 == 1) ++ap;
					if (p3 == 0 && p4 == 1) ++ap;
					if (p4 == 0 && p5 == 1) ++ap;
					if (p5 == 0 && p6 == 1) ++ap;
					if (p6 == 0 && p7 == 1) ++ap;
					if (p7 == 0 && p8 == 1) ++ap;
					if (p8 == 0 && p9 == 1) ++ap;
					if (p9 == 0 && p2 == 1) ++ap;
					if (ap == 1)
					{
						if (p2*p4*p8 == 0)
						{
							if (p2*p6*p8 == 0)
								//记入标记点
								flagPoint.push_back(cv::Point(i, j));
						}
					}
				}
			}
		}
		//将标记点删除
		for (std::vector<cv::Point>::const_iterator it = flagPoint.begin(); it != flagPoint.end(); it++)
			result.at<uchar>((*it).x, (*it).y) = 0;
		//直到没有点满足，算法结束
		if (flagPoint.size() == 0)
			break;
		else
			flagPoint.clear();
	}
}

