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

		//�ļ����
		intptr_t hFile = 0;//Win10
						   //long hFile = 0;
						   //�ļ���Ϣ  
		struct _finddata_t fileinfo;
		std::string p;
		try
		{
			if ((hFile = _findfirst(p.assign(folder_path).append("\\*.*").c_str(), &fileinfo)) != -1)
			{
				do
				{
					//�����Ŀ¼,����֮  
					//�������,�����б�  
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
		//somѵ������
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

		//��ÿ��ͼƬ������������
		//for(const auto& p : files)
		{
			//ɫ���趨
			{
				testMinDis = 1.0;
				testMaxDis = 20.0;
			}
			
			std::string linePicturepath;

			//���ַ�����ͼƬɫ��
			while (true)
			{
				//����som
				double tetsDis = (testMinDis + testMaxDis) / 2.0;
				std::string msg = createSOM(path, linePicturepath, tetsDis,rect);

				if (msg == "ɫ�����")
				{
					testMaxDis = tetsDis;
				}
				else if (msg == "ɫ���С")
				{
					testMinDis = tetsDis;
				}
				else if (msg == "")
				{
					break;
				}

			}

			//somѵ��
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