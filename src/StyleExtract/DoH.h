#pragma once
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#define DoHMAXPOINTS 50000

class DoH
{
public:
	
	typedef struct myPoint
	{
		int x;
		int y;
		int radiu;
	}_point;

	
	static void cornerDOH(unsigned char* srcImg, _point* corner, int& count, int height, int width, float minScale, float maxScale, int stepScale)
	{
		int numScale = int((maxScale - minScale) / stepScale);
		int** scaledImg = new int*[numScale];
		for (int i = 0; i<numScale; i++)
			scaledImg[i] = new int[height*width];

		for (int k = 0; k<numScale; k++)
		{
			float scale = minScale + stepScale*k;
			int kernelWin = 3 * scale;
			float *kernelxx = new float[(2 * kernelWin + 1)*(2 * kernelWin + 1)];
			float *kernelyy = new float[(2 * kernelWin + 1)*(2 * kernelWin + 1)];
			float *kernelxy = new float[(2 * kernelWin + 1)*(2 * kernelWin + 1)];
			getDOHkernel(kernelWin, scale, &kernelxx, &kernelyy, &kernelxy);

			int** tmpImg = new int*[3];
			for (int i = 0; i<3; i++)
				tmpImg[i] = new int[height*width];
			converlution(srcImg, &(tmpImg[0]), height, width, kernelxx, kernelWin);
			converlution(srcImg, &(tmpImg[1]), height, width, kernelyy, kernelWin);
			converlution(srcImg, &(tmpImg[2]), height, width, kernelxy, kernelWin);

			detValue(&(scaledImg[k]), tmpImg[0], tmpImg[1], tmpImg[2], scale, height, width);

			for (int i = 0; i<3; i++)
				delete[] tmpImg[i];
			delete[] tmpImg;
			delete[] kernelxx;
			delete[] kernelyy;
			delete[] kernelxy;
		}

		count = 0;
		for (int i = 1; i<height - 1; i++)
		{
			for (int j = 1; j<width - 1; j++)
			{
				for (int k = 1; k<numScale - 1; k++)
				{
					if ((count) >= DoHMAXPOINTS)
					{
						for (int m = 0; m<numScale; m++)
							delete[] scaledImg[m];
						delete[] scaledImg;

						return;
					}

					_point cp;
					cp.x = j;
					cp.y = i;
					float scale = minScale + k*stepScale;
					cp.radiu = int(1.414*scale + 0.5);
					if (isMax(scaledImg[k - 1], scaledImg[k][i*width + j], cp, 1, width, false) &&
						isMax(scaledImg[k], scaledImg[k][i*width + j], cp, 1, width, true) &&
						isMax(scaledImg[k + 1], scaledImg[k][i*width + j], cp, 1, width, false))
					{
						corner[(count)++] = cp;
						break;
					}
				}
			}
		}

		for (int i = 0; i<numScale; i++)
			delete[] scaledImg[i];
		delete[] scaledImg;
	}

private:
	static void converlution(unsigned char* srcImg, int** dstImg, int height, int width, float* kernel, int kernelWin)
	{
		for (int i = 0; i<height; i++)
		{
			for (int j = 0; j<width; j++)
			{
				float sumValue = 0;
				int count = 0;
				for (int m = i - kernelWin; m <= i + kernelWin; m++)
				{
					for (int n = j - kernelWin; n <= j + kernelWin; n++)
					{
						if (m >= 0 && m<height && n >= 0 && n<width)
							sumValue += int(srcImg[m*width + n])*kernel[count];
						count++;
					}
				}
				sumValue *= 100;
				(*dstImg)[i*width + j] = int(sumValue);
			}
		}
	}

	static void getDOHkernel(int halfWin, float sita, float** kernelxx, float** kernelyy, float** kernelxy)
	{
		int winSize = 2 * halfWin + 1;
		float tmp1, tmpxx, tmpyy, tmpxy;
		float sumValuexx = 0, sumValueyy = 0, sumValuexy = 0;
		float powsita = sita*sita;
		for (int i = -halfWin; i <= halfWin; i++)
		{
			for (int j = -halfWin; j <= halfWin; j++)
			{
				tmp1 = -1 * (i*i + j*j) / (2 * powsita)/(2*3.141592653*powsita);
				tmpxx = exp(tmp1)*(j*j / powsita - 1);
				tmpyy = exp(tmp1)*(i*i / powsita - 1);
				tmpxy = exp(tmp1)*i*j / powsita;

				sumValuexx += tmpxx;
				sumValueyy += tmpyy;
				sumValuexy += tmpxy;

				(*kernelxx)[(i + halfWin)*winSize + (j + halfWin)] = tmpxx;
				(*kernelyy)[(i + halfWin)*winSize + (j + halfWin)] = tmpyy;
				(*kernelxy)[(i + halfWin)*winSize + (j + halfWin)] = tmpxy;
			}
		}
		for (int i = 0; i<winSize*winSize; i++)
		{
			(*kernelxx)[i] -= sumValuexx / (winSize*winSize);
			(*kernelyy)[i] -= sumValueyy / (winSize*winSize);
			(*kernelxy)[i] -= sumValuexy / (winSize*winSize);
		}
	}


	static bool isMax(int* compImg, int compValue, _point compPoint, int kerWin, int imgWidth, bool isSameK)
	{
		for (int i = compPoint.y - kerWin; i <= compPoint.y + kerWin; i++)
		{
			for (int j = compPoint.x - kerWin; j <= compPoint.x + kerWin; j++)
			{
				if (isSameK && i == compPoint.y && j == compPoint.x)
					continue;
				if (abs(compValue) <= abs(compImg[i*imgWidth + j]))
					return false;
			}
		}
		return true;
	}


	static void detValue(int** imgscale, int* imgxx, int* imgyy, int* imgxy, float scale, int height, int width)
	{
		for (int i = 0; i<height; i++)
		{
			for (int j = 0; j<width; j++)
			{
				float xx = imgxx[i*width + j] / 100.0;
				float yy = imgyy[i*width + j] / 100.0;
				float xy = imgxy[i*width + j] / 100.0;
				int vv = xx*yy - xy*xy;
				(*imgscale)[i*width + j] = int(vv / (scale*scale));
			}
		}
	}



};
