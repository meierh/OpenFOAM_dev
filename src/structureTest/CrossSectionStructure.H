#ifndef CROSSSECTIONSTRUCTURE_H
#define CROSSSECTIONSTRUCTURE_H

#include "LineStructure.H"

namespace Foam
{
class CrossSection
{
    public:
        explicit CrossSection
        (
            std::vector<gsNurbs<double>> a_k,
            std::vector<gsNurbs<double>> b_k
        );
        double operator()(gsMatrix<double,1,1> parameter, double rad);
        std::function<double(double)> getEvalOnPoint(gsMatrix<double,1,1> parameter);

        
    private:
        const std::vector<gsNurbs<double>> a_k;
        const std::vector<gsNurbs<double>> b_k;
        const uint numberCoeffs;
};

class CrossSectionStructure : public LineStructure
{
    public:
        explicit CrossSectionStructure
        (
            Time& runTime,
            const dimensionedScalar& alpha,
            const volScalarField& T,
            const volScalarField& p,
            const volVectorField& U,
            dynamicRefineFvMesh& mesh,
            const dimensionedScalar nu,
            std::vector<CrossSection> rodCrossSection
        );
        CrossSectionStructure(const CrossSectionStructure&) = delete;
    private:
        std::vector<CrossSection> rodCrossSection;
};
}
#endif
