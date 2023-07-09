#ifndef _PictureOP_H
#define _PictureOP_H
#pragma once

#include <vector>
#include <string>
#include<io.h>
#include<direct.h>

using namespace std;
class PictureOP
{
public:
	static void getLinePicture(const std::string &picturepath, std::string &linepath, std::string &linepath2,
		std::string &pointpath, std::string &backpath)
	{
		char szDrive[_MAX_DRIVE];   //磁盘名
		char szDir[_MAX_DIR];       //路径名
		char szFname[_MAX_FNAME];   //文件名
		char szExt[_MAX_EXT];       //后缀名
		_splitpath(picturepath.c_str(), szDrive, szDir, szFname, szExt);

		string ImagefoldName = std::string(szDir).substr(0, std::string(szDir).find_last_of("\\"));
		ImagefoldName = std::string(szDir).substr(0, ImagefoldName.find_last_of("\\"));
		ImagefoldName = ImagefoldName + "\\edge2\\";

		linepath = std::string(szDrive) + ImagefoldName + std::string(szFname) + "_edges" + std::string(szExt);

		ImagefoldName = std::string(szDir).substr(0, std::string(szDir).find_last_of("\\"));
		ImagefoldName = std::string(szDir).substr(0, ImagefoldName.find_last_of("\\"));
		ImagefoldName = ImagefoldName + "\\edge3\\";

		linepath2 = std::string(szDrive) + ImagefoldName + std::string(szFname) + "_edges" + std::string(szExt);
		
		ImagefoldName = std::string(szDir).substr(0, std::string(szDir).find_last_of("\\"));
		ImagefoldName = std::string(szDir).substr(0, ImagefoldName.find_last_of("\\"));
		ImagefoldName = ImagefoldName + "\\blobDOH\\";

		pointpath = std::string(szDrive) + ImagefoldName + std::string(szFname) + "_blob_10_30.txt" ;

		ImagefoldName = std::string(szDir).substr(0, std::string(szDir).find_last_of("\\"));
		ImagefoldName = std::string(szDir).substr(0, ImagefoldName.find_last_of("\\"));
		ImagefoldName = ImagefoldName + "\\background2\\";

		backpath = std::string(szDrive) + ImagefoldName + std::string(szFname) + "_back" + std::string(szExt);


	}

	static void getBackgroundPicture(const std::string &picturepath,std::string &backpath)
	{
		char szDrive[_MAX_DRIVE];   //磁盘名
		char szDir[_MAX_DIR];       //路径名
		char szFname[_MAX_FNAME];   //文件名
		char szExt[_MAX_EXT];       //后缀名
		_splitpath(picturepath.c_str(), szDrive, szDir, szFname, szExt);

		string ImagefoldName = std::string(szDir).substr(0, std::string(szDir).find_last_of("\\"));
		ImagefoldName = std::string(szDir).substr(0, ImagefoldName.find_last_of("\\"));
		ImagefoldName = ImagefoldName + "\\background2\\";

		backpath = std::string(szDrive) + ImagefoldName + std::string(szFname) + "_back" + std::string(szExt);
	}

	static void getLinePicture(const std::string &picturepath, std::string &linepath)
	{
		char szDrive[_MAX_DRIVE];   //磁盘名
		char szDir[_MAX_DIR];       //路径名
		char szFname[_MAX_FNAME];   //文件名
		char szExt[_MAX_EXT];       //后缀名
		_splitpath(picturepath.c_str(), szDrive, szDir, szFname, szExt);

		string ImagefoldName = std::string(szDir).substr(0, std::string(szDir).find_last_of("\\"));
		ImagefoldName = std::string(szDir).substr(0, ImagefoldName.find_last_of("\\"));
		ImagefoldName = ImagefoldName + "\\edge2\\";

		linepath = std::string(szDrive) + ImagefoldName + std::string(szFname) + "_edges" + std::string(szExt);

	
		if (_access(linepath.c_str(), 0) == 0)
		{

		}
		else
		{
			linepath = "";
		}
		
	}
	
	static bool exists(const std::string& name) {
		return (_access(name.c_str(), 0) != -1);
	}

	static void getRelationshipPath(const std::string &mappath, std::string &Relationshippath)
	{
		char szDrive[_MAX_DRIVE];   //磁盘名
		char szDir[_MAX_DIR];       //路径名
		char szFname[_MAX_FNAME];   //文件名
		char szExt[_MAX_EXT];       //后缀名
		_splitpath(mappath.c_str(), szDrive, szDir, szFname, szExt);

		
		Relationshippath = std::string(szDrive) + std::string(szDir) + std::string(szFname) + "_relationship" + std::string(szExt);
		if (_access(Relationshippath.c_str(), 0) == 0)
		{

		}
		else
		{
			//Relationshippath = "";
		}
	}
};
#endif
