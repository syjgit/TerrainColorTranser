#pragma once
#include "Hdata.h"
#include "FeatureCompute.h"
#include "colorManager.h"


class MapHarmony
{
	static cv::Mat FeatureToMatrix(std::vector <FeatureCompute::HarmonyFeature> Features)
	{
		cv::Mat result(1,118, CV_64F), tmpMat;

		int index = 0;
		int length=0;
		for (auto& value : Features)
		{
			length = value.math_compute.cols* value.math_compute.rows;
			tmpMat = result.colRange(index, index + length);
			cv::Mat temp = value.math_compute.t();
			temp.reshape(0, 1).copyTo(tmpMat);
			temp = temp.reshape(0, 1);
			index += length;

			if (value.HueProb.rows == 1)
			{
				length = value.HueProb.cols* value.HueProb.rows;
				tmpMat = result.colRange(index, index + length);
				value.HueProb.copyTo(tmpMat);
				index += length;
			}
			if (value.Plane.rows == 1)
			{
				length = value.Plane.cols* value.Plane.rows;
				tmpMat = result.colRange(index, index + length);
				value.Plane.copyTo(tmpMat);
				index += length;
			}
		}
	
		return result;
	}

	static std::vector <FeatureCompute::HarmonyFeature> createFeaturesFromData(cv::Mat rgbs)
	{
		std::vector <FeatureCompute::HarmonyFeature> allFeatures(4);

		cv::Mat plane = cv::Mat::zeros(1, 7, CV_64F);
		cv::Mat  hueProbFeatures;
		double satValThresh = 0.2;

		for (int c = 1; c <= 4; c++)
		{
			cv::Mat col;
			cv::Mat hsv;
			FeatureCompute::HarmonyFeature features;

			//chsv
			if (c == 1)
			{
				col = color_manager::RGB2CHSV(rgbs);
			}
			//"lab"
			else if (c == 2)
			{
				//features.type = eColorType::e_lab;
				col = color_manager::RGB2Lab(rgbs);

				//labs = col;
			}
			//"hsv"
			else if (c == 3)
			{
				//features.type = eColorType::e_hsv;
				col = color_manager::RGB2HSV(rgbs);
				hsv = col;
			}
			// "rgb"
			else if (c == 4)
			{
				//features.type = eColorType::e_rgb;
				col = rgbs;
			}

			features.math_compute = FeatureCompute::mathcompute(col);

			if (c == 3)
			{
				hueProbFeatures = FeatureCompute::getHueProbFeatures(hsv, satValThresh, hueProb);
				for (int i = 1; i < hueProbFeatures.cols; i++)
				{
					if (hueProbFeatures.at<double>(0, i) == -99)
					{
						hueProbFeatures.at<double>(0, i) += 0.0001;
					}
				}

				features.HueProb = hueProbFeatures;
			}
			else
			{
				/////c
				cv::Mat col00 = col.t();
				FeatureCompute::PlaneFeatures _plan = FeatureCompute::getPlaneFeatures(col00);

				////
				for (int plan_i = 0; plan_i < 3; plan_i++)
				{
					plane.at<double>(0, plan_i) = _plan._normal.at<double>(plan_i, 0);
					plane.at<double>(0, plan_i + 3) = _plan.pctExplained.at<double>(0, plan_i);
				}

				plane.at<double>(0, 6) = _plan.sse;

				features.Plane = plane;
			}

			allFeatures[c - 1] = (features);
		}
		return allFeatures;
	}


public:
	
	static hueProbs hueProb;

	static double glmnetPredict2(const cv::Mat& rgbs)
	{
		std::vector <FeatureCompute::HarmonyFeature> datapoints = createFeaturesFromData(rgbs);

		cv::Mat mat = FeatureToMatrix(datapoints);
		
		int index = 0;
		for (int i = 0; i < mat.rows; ++i)
		{
			for (int j = 0; j < mat.cols; ++j)
			{
				if (isnan(mat.at<double>(i, j)))
					mat.at<double>(i, j) = 0.0;   // 通过下标遍历像素
			}
		}
		
		FitData fitdata;

		cv::Mat temp = mat*(fitdata.LassoA);
		return temp.at<double>(0, 0) + fitdata.LassoB;

	}
};

