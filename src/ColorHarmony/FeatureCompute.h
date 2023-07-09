#pragma once
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
                                                                                                                               
#include <fstream>
#define M_PI    3.14159265358979323846
class FeatureCompute
{

	static double stdToSStd(const double &input, const int& num)
	{
		return sqrt(pow(input, 2)*num / (num - 1));
	}

	static double median(const cv::Mat& input)
	{
		cv::Mat v = input.clone();
		if (v.cols % 2 != 0)
		{
			std::nth_element(v.begin<double>(), v.begin<double>() + v.cols / 2, v.end<double>());
			return v.at<double>(v.cols / 2);
		}
		else
		{
			std::nth_element(v.begin<double>(), v.begin<double>() + v.cols / 2, v.end<double>());
			return (v.at<double>(v.cols / 2) + v.at<double>(v.cols / 2) - 1) / 2.0;
		}

	}

	static cv::Mat TRound(const cv::Mat& input)
	{
		cv::Mat matix_reduced_precision;
		input.convertTo(matix_reduced_precision, CV_32S);
		return matix_reduced_precision;
	}

	static cv::Mat getBasicStats(const std::vector<double>& x, const int& addLog = 1)
	{

		cv::Mat result = cv::Mat::zeros(1, 8, CV_64F);
		if (x.size() > 0)
		{
			cv::Mat log_x;
			cv::log(cv::Mat(x).t() + 0.000001, log_x);

			/////x
			cv::Mat _mean ;
			cv::Mat _stdev;
			cv::meanStdDev(x, _mean, _stdev);
			double _min, _max;
			cv::minMaxLoc(x, &_min, &_max);
			

			/////log_x
			cv::Mat _mean_log_x;
			cv::Mat _stdev_log_x;
			cv::meanStdDev(log_x, _mean_log_x, _stdev_log_x);


			double _min_log_x, _max_log_x;
			cv::minMaxLoc(log_x, &_min_log_x, &_max_log_x);

			result=(cv::Mat_<double>(1, 8) << _mean.at<double>(0,0), stdToSStd(_stdev.at<double>(0, 0), x.size()), _min, _max, _mean_log_x.at<double>(0, 0), stdToSStd(_stdev_log_x.at<double>(0, 0),log_x.cols), _min_log_x, _max_log_x );
		}
		//MatrixOpera::nan(result);

		return result;

	}


	static cv::Mat cvCosMat(const cv::Mat& a)
	{
		int rows = a.rows;
		int cols = a.cols;
		cv::Mat out = cv::Mat::ones(rows, cols, CV_64F);
		for (int i = 0; i<rows; i++)
		{
			double* ptra = (double*)(a.data + i*a.step);
			double* ptrout = (double*)(out.data + i*out.step);
			for (int j = 0; j<cols; j++)
			{
				*ptrout = cos(*ptra);
				ptra++;
				ptrout++;
			}
		}
		return out;
	}

	static  cv::Mat circ_vmpdf(cv::Mat alpha, const double& thetahat, const double& kappa)
	{

		double C = 0.001827088328087;

		cv::Mat temp = cvCosMat(alpha - thetahat)*kappa;

		cv::exp(temp, temp);
		return temp*C;
	}


	//mode==1 计算每列的均值， mode==2 计算每行的均值
	static cv::Mat mean(const cv::Mat& mat, const int& mode)
	{

		cv::Mat result;
		if (mode == 1)
		{
			result= cv::Mat::zeros(1, mat.cols, CV_64F);
			for (int i = 0; i < mat.cols; i++)
			{
				result.at<double>(0, i) = cv::mean(mat.col(i))[0];
			}
			
		}
		else if (mode == 2)
		{

			result = cv::Mat::zeros(mat.rows, 1, CV_64F);
			for (int i = 0; i < mat.rows; i++)
			{
				result.at<double>(i, 0) = cv::mean(mat.row(i))[0];
			}
		}

		return result;
	}

	static std::vector<cv::Mat> pca2(const cv::Mat& data)
	{
		int M = data.cols;
		int N = data.rows;

		cv::Mat mn= mean(data, 2);
		cv::Mat data1 = data - cv::repeat(mn, 1, M);

		cv::Mat Y = data1.t() / std::sqrt(M - 1);

		int k = 3;
		
		std::vector<cv::Mat> result;

		cv::Mat S, U, PC; //u is left,v is left

		cv::SVDecomp(Y, S, U, PC);



		//_matrix<double> mat_S(S.size(), 1);
		//mat_S.data = S;

		cv::Mat mat_V = S.mul(S);

		PC *= -1;

		cv::Mat mat_signals;
		mat_signals = (PC* data1);

		return{ mat_signals, PC.t(), mat_V };
	}

	static cv::Mat extractByMat(const cv::Mat& data, const cv::Mat& indexs)
	{
		cv::Mat result(indexs.size(), CV_64F);
		for (int i = 0; i < indexs.rows; i++)
		{
			int index = indexs.at<uchar>(i, 0);
			result.at<double>(i,0) = data.at<int>(0, index);
		}
		return result;
	}


	static cv::Mat thresholdIndexs(const cv::Mat& data, double value)
	{
		std::vector<uchar> indexs;
		for (int i = 0; i < data.cols; i++)
		{
			if (data.at<double>(0, i) >= value)
				indexs.push_back(i);

		}

		return cv::Mat(indexs).clone();
	}

	static cv::Mat Tmin(const cv::Mat& data)
	{
		cv::Mat temp1 = data.row(1);
		cv::Mat temp2 = data.row(2);
		return cv::min(temp1, temp2);
	}

	


public:

	struct PlaneFeatures
	{
		cv::Mat _normal;
		cv::Mat pctExplained;
		cv::Mat meanX;
		double sse;
	};


	struct HarmonyFeature
	{
		cv::Mat math_compute;
		cv::Mat Plane;
		cv::Mat HueProb;

		double CH;
		double Grads;

	};


	static void HInit(hueProbs& hue)
	{
		std::ifstream hueJoint("..\\data\\hueProbs.hueJoint.txt");
		if (hueJoint.fail())
		{
			return;
		}
		double value = 0;
		while (!hueJoint.eof())
		{
			hueJoint >> value;
			hue.hueJoint.push_back(value);
		}
		hue.hueJoint.pop_back();
		hueJoint.close();

		std::ifstream hueAdjacency("..\\data\\hueProbs.hueAdjacency.txt");
		if (hueAdjacency.fail())
		{
			return;
		}
		while (!hueAdjacency.eof())
		{
			double value = 0;
			hueAdjacency >> value;
			hue.hueAdjacency.push_back(value);
		}
		hue.hueAdjacency.pop_back();
		hueAdjacency.close();
		
	}



	//n*6
	static cv::Mat mathcompute(const cv::Mat& colors)
	{
		cv::Mat result(colors.rows, 6, CV_64F, cv::Scalar(0));


		for (int i = 0; i < colors.rows; i++)
		{
			cv::Mat row = colors.rowRange(i, i + 1).clone();

			//min max
			double maxvalue = 0, minvalue = 0;
			cv::minMaxLoc(row, &minvalue, &maxvalue);
			result.at<double>(i, 4) = minvalue;
			result.at<double>(i, 3) = maxvalue;

			//mean stdDev
			cv::Mat mean;
			cv::Mat stdDev;
			cv::meanStdDev(row, mean, stdDev);
			result.at<double>(i, 0) = mean.at<double>(0, 0);
			result.at<double>(i, 1) = stdToSStd(stdDev.at<double>(0, 0), row.cols);

			result.at<double>(i, 2) = median(row);

			//max-min
			result.at<double>(i, 5) = result.at<double>(i, 3) - result.at<double>(i, 4);
		
		}

		return result;
		
	}

	static cv::Mat getHueProbFeatures(const cv::Mat& hsv,const double& satValThresh, const hueProbs& hueProbs,int theme=1)
	{
		
		cv::Mat hueFeatures = cv::Mat::ones(theme, 25, CV_64F)*-99;;
		
		cv::Mat temp =Tmin(hsv);

		cv::Mat selectColors = thresholdIndexs(temp, satValThresh);
		
		cv::Mat mo = (cv::Mat_<double>(3, 1) << 359, 100, 100);
		
		temp = cv::repeat(mo, 1, hsv.cols);
		cv::Mat hsv2 = TRound(hsv.mul(temp));
		hsv2 = hsv2 + 1;
		

		cv::Mat visHues = extractByMat(hsv2.row(0), selectColors);

		std::vector<double> hueJointList;
		for (int h1 = 0; h1 < visHues.rows; h1++)
		{
			for (int h2 = h1; h2 < visHues.rows; h2++)
			{
				int index = (visHues.at<double>(h2,0) - 1) * 360 + (visHues.at<double>(h1,0) - 1);
				hueJointList.emplace_back(hueProbs.hueJoint.at<double>(index, 0));
			}
		}

		std::vector<double> hueAdjList;
		for (int h1 = 0; h1 < visHues.rows - 1; h1++)
		{
			int index = (visHues.at<double>( h1,0) - 1) * 360 + (visHues.at<double>( h1 + 1,0) - 1);
			hueAdjList.emplace_back(hueProbs.hueAdjacency.at<double>(index, 0));
		}

		std::vector<double> hueProbList;

		for (int h1 = 0; h1 < visHues.rows; h1++)
		{
			hueProbList.emplace_back(hueProbs.hueProb.at<double>(visHues.at<double>(h1, 0) - 1, 0));
		}

		cv::Mat hueProbFeatures = getBasicStats(hueProbList, 1);
		cv::Mat hueJointProbFeatures = getBasicStats(hueJointList, 1);
		cv::Mat hueAdjProbFeatures = getBasicStats(hueAdjList, 1);

		cv::Mat alpha = m_alpha.clone();
		cv::Mat pMix = 0.001*cv::Mat::ones(alpha.size(), CV_64F);

		for (int j = 0; j < visHues.rows; j++)
		{
			pMix = pMix + circ_vmpdf(alpha, (visHues.at<double>(j, 0)) * 2 * M_PI, 2 * M_PI);
		}

		double sum = cv::sum(pMix)[0];
		pMix = pMix / sum;

		double entropy;
		if (visHues.rows != 0)
		{
			
			log(pMix,temp);
			entropy = -cv::sum(pMix.mul(temp))[0];
		}
		else
		{
			entropy = 5.9;
		}


		cv::hconcat(hueProbFeatures, hueJointProbFeatures, hueFeatures);
		cv::hconcat(hueFeatures, hueAdjProbFeatures, hueFeatures);
		cv::hconcat(hueFeatures, entropy, hueFeatures);
		
		//hueFeatures = (cv::Mat_<double>(theme, 25) << hueProbFeatures, hueJointProbFeatures, hueAdjProbFeatures, entropy);

		/*HueFeatures hueFeatures;
		hueFeatures.hueProbFeatures = hueProbFeatures;
		hueFeatures.hueJointProbFeatures = hueJointProbFeatures;
		hueFeatures.hueAdjProbFeatures = hueAdjProbFeatures;
		hueFeatures.entropy = entropy;
	*/
		//测试完成
		return hueFeatures;
	}

	static PlaneFeatures getPlaneFeatures(cv::Mat& X)
	{
		cv::Mat eee = X.t();
		std::vector<cv::Mat> res = pca2(eee);

		cv::Mat signals = res[0];
		cv::Mat coeff = res[1];
		cv::Mat roots = res[2];

		cv::Mat _normal(coeff.rows, 1, CV_64F);

		int index = 0;
		for (int i = 0; i < coeff.rows; i++)
		{
			_normal.at<double>(i, 0) = coeff.at<double>(i, coeff.cols - 1);
			index++;
		}		


		if (_normal.at<double>(0, 0) < 0)
		{
			_normal = _normal* -1;
		}

		cv::Mat pctExplained;
		if (cv::sum(roots)[0] == 0)
		{
			pctExplained = cv::Mat::zeros(3, 1, CV_64F);
		}
		else
		{
			pctExplained = roots.t() / cv::sum(roots)[0];
		}

		int n = X.rows;
		int p = X.cols;

		cv::Mat meanX = mean(X, 1);

		cv::Mat tmp = X - cv::repeat(meanX, n, 1);
		cv::Mat tmp2;
		tmp2 = (tmp*_normal);

		cv::Mat error = cv::abs(tmp2);
		
		double sse = cv::sum(error.mul(error))[0];

		return  PlaneFeatures{ _normal, pctExplained, meanX ,sse };
	}


private:

};


