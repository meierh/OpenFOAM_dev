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
#include "Nurbs.H"
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

using namespace tinyxml2;
using namespace std;

namespace Foam
{
class NurbsReader
{
    public:
        explicit NurbsReader(fileName runDirectory, fileName caseName);
        NurbsReader(const NurbsReader&) = delete;
        std::shared_ptr<std::vector<Nurbs>> getNurbsCurves();  
        
    private:
        word getXMLPath();
        std::shared_ptr<std::vector<Nurbs>> readOutNurbsFromXML();
        
        word fullXMLPath;
        std::shared_ptr<std::vector<Nurbs>> nurbsCurves;
        fileName runDirectory;
        fileName caseName;
};
}
#endif
