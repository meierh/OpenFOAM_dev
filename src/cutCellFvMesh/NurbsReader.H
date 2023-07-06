#ifndef NURBSREADER_H
#define NURBSREADER_H

#if defined( _MSC_VER )
	#if !defined( _CRT_SECURE_NO_WARNINGS )
		#define _CRT_SECURE_NO_WARNINGS		// This test file is not intended to be secure.
	#endif
#endif

#include "tinyxml2.h"
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <memory>
#include "Nurbs1D.H"
#include <vector>

#if defined( _MSC_VER ) || defined (WIN32)
	#include <crtdbg.h>
	#define WIN32_LEAN_AND_MEAN
	#include <windows.h>
	_CrtMemState startMemState;
	_CrtMemState endMemState;
#else
	#include <sys/stat.h>
	#include <sys/types.h>
#endif
    
#include <dirent.h>

namespace Foam
{
class NurbsData
{
    public:
        DynamicList<scalar> knots;
        DynamicList<vector> coefs;
        DynamicList<scalar> weights;
        label nurbsDegree;
        scalar radius;
};
    
class NurbsReader
{
    public:
        explicit NurbsReader(fileName runDirectory, fileName caseName);
        NurbsReader(const NurbsReader&) = delete;
        std::shared_ptr<std::vector<Nurbs1D>> getNurbsCurves();
        std::shared_ptr<std::vector<NurbsData>> getNurbsData();

        
    private:
        std::shared_ptr<std::vector<Nurbs1D>> readOutNurbsFromXML();
        word getXMLPath();
        word fullXMLPath;
        std::shared_ptr<std::vector<Nurbs1D>> nurbsCurves;
        std::shared_ptr<std::vector<NurbsData>> nurbsData;
        fileName runDirectory;
        fileName caseName;
};
}
#endif
