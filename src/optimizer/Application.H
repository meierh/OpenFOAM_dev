#ifndef APPLICATION_H
#define APPLICATION_H
#include "fvMesh.H"
#include "argList.H"
#include <memory>

namespace Foam
{
    
class Application
{
    public:
        explicit Application(int argc, char *argv[]);
        virtual ~Application() = default;
        
        virtual void run()=0;

    protected:
        int argc;
        char **argv;
        
        // --- set Root Case Lists --- //
        void setRootCaseLists();
        
        //void listOptions();
        void setRootCase();
        //void listOutput();

        Foam::argList args;
        bool _listOptions;
        // --- ------------------- --- //
        
        // --- create Time --- //
        void createTime();

        Foam::Time runTime;
        // --- ----------- --- //
        
        // --- create Time --- //
        void createMesh();

        std::unique_ptr<Foam::fvMesh> mesh;
        // --- ----------- --- //
        
        // --- create Fields --- //
        virtual void createFields();
        // --- ----------- --- //
};

}
#endif
