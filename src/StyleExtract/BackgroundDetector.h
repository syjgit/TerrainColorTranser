#pragma once

#include "Saliency\GMRsaliency.h"

class BackgroundDetector
{
public:
	static void Backgrounddetect(const char* str,string backfile)
	{
		Mat sal, img;
		img = imread(str);
		GMRsaliency GMRsal;
		
		sal = GMRsal.GetSal(img);
		imwrite(backfile, sal * 255);
	}
};