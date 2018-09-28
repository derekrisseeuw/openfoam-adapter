#ifndef FSI_FORCE_H
#define FSI_FORCE_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads stress
class Stress : public CouplingDataUser
{

private:

    //- OpenFOAM fvMesh object (we need to access the objects' registry multiple times)
    const Foam::fvMesh& mesh_;

    //- Stress field
    Foam::volVectorField * Stress_;

    //- Stress tensor (see the OpenFOAM "Stress" function object)
    Foam::tmp<Foam::volSymmTensorField> devRhoReff(dimensionedScalar rho) const;

public:

    //- Constructor
    Stress
    (
        const Foam::fvMesh& mesh,
        const fileName& timeName // We create an IOobject and we need the time directory
        /* TODO: We should add any required field names here.
        /  They would need to be vector fields.
        /  See CHT/Temperature.H for details.
        /  Apply any changes also to Stress.C.
        */
    );

    //- Write the displacement values into the buffer
    void write(double * buffer);

    //- Read the displacement values from the buffer
    void read(double * buffer);

    //- Destructor
    ~Stress();

};

}
}

#endif
