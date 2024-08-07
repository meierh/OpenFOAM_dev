#ifndef APPLICATION_H
#define APPLICATION_H
#include "fvCFD.H"

namespace Foam
{
    
class Application
{
    public:
        explicit Application(int argc, char *argv[]);
        virtual ~Application()=default;
        
        virtual void run();

    protected:
        int argc;
        char **argv;
        
        // --- set Root Case Lists --- //
        void setRootCaseLists();
        
        void listOptions();
        void setRootCase();
        void listOutput();

        Foam::argList args;
        bool _listOptions;
        // --- ------------------- --- //
        
        // --- create Time --- //
        void createTime();

        Foam::Time runTime;
        // --- ----------- --- //
        
        // --- create Time --- //
        void createMesh();

        Foam::fvMesh mesh;
        // --- ----------- --- //
        
        // --- create Fields --- //
        virtual void createFields();
        // --- ----------- --- //
};

}
#endif
