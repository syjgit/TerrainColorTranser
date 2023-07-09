#pragma once
#include <vector>
#include <string>
#include<io.h>
#include "pictureOp.h"
#include <atltypes.h>

#include "CController.h"

class BatchProcessing
{
public:
	BatchProcessing() {
	};
	CController* m_pSom;
private:
	
	std::string m_strMessage;


	double testMinDis = 1.0;
	double testMaxDis = 20.0;

public:
	void GetFiles(std::string folder_path, std::vector<std::string>& files)
	{

		//文件句柄
		intptr_t hFile = 0;//Win10
						   //long hFile = 0;
						   //文件信息  
		struct _finddata_t fileinfo;
		std::string p;
		try
		{
			if ((hFile = _findfirst(p.assign(folder_path).append("\\*.*").c_str(), &fileinfo)) != -1)
			{
				do
				{
					//如果是目录,迭代之  
					//如果不是,加入列表  
					if ((fileinfo.attrib &  _A_SUBDIR))
					{
						if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
							GetFiles(p.assign(folder_path).append("\\").append(fileinfo.name), files);
					}
					else
					{

						//files.push_back(p.assign(folder_path).append("\\").append(fileinfo.name));
						files.emplace_back(fileinfo.name);
					}
				} while (_findnext(hFile, &fileinfo) == 0);

				_findclose(hFile);
			}
		}
		catch (std::exception e)
		{
		}
	}

	void train(const std::string& linePicturepath)
	{
		//som训练样本
		int time = 0;
		while (true)
		{
			bool bDone = false;
			if (!m_pSom->Finished())
			{
				if (!m_pSom->Train()){break;}
			}
			else{break;}

			time++;
		}

		if (linePicturepath != "")
		{
			m_pSom->ijs.clear();
			m_pSom->matchingindex.clear();
			while (1)
			{
				if (!m_pSom->Finished2())
				{
					if (!m_pSom->Train2()){break;}
				}
				else{break;}
			}
		}
	}


	void Analyse(const std::string path, CRect rect)
	{

		//对每个图片进行流形生成
		//for(const auto& p : files)
		{
			//色差设定
			{
				testMinDis = 1.0;
				testMaxDis = 20.0;
			}
			
			std::string linePicturepath;

			//二分法计算图片色差
			while (true)
			{
				//创建som
				double tetsDis = (testMinDis + testMaxDis) / 2.0;
				std::string msg = createSOM(path, linePicturepath, tetsDis,rect);

				if (msg == "色差过大")
				{
					testMaxDis = tetsDis;
				}
				else if (msg == "色差过小")
				{
					testMinDis = tetsDis;
				}
				else if (msg == "")
				{
					break;
				}

			}

			//som训练
			train(linePicturepath);

			

			
		}

	}

	std::string createSOM(std::string img_path,std::string& line_apth,double dis,const CRect& rect)
	{
		PictureOP::getLinePicture(img_path, line_apth);
		if (line_apth != "")
		{
			m_pSom = new CController(img_path, line_apth, dis, rect.Width(),
				rect.Height(),
				constNumCellsAcross,
				constNumCellsDown,
				constNumIterations);
		}
		else
		{
			m_pSom = new CController(img_path, dis, rect.Width(),
				rect.Height(),
				constNumCellsAcross,
				constNumCellsDown,
				constNumIterations);
		}

		return m_pSom->GetMessageeX();
		
	}



};