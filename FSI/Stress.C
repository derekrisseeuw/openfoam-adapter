#include "Stress.H"

using namespace Foam;

preciceAdapter::FSI::Stress::Stress
(
    const Foam::fvMesh& mesh,
    const fileName& timeName
    /* TODO: We should add any required field names here.
    /  They would need to be vector fields.
    /  See CHT/Temperature.C for details.
    */
)
:
mesh_(mesh)
{
    dataType_ = vector;

    // TODO: Is this ok?
    Stress_ = new volVectorField
    (
        IOobject
        (
            "Stress",
            timeName,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "fdim",
            dimensionSet(1,1,-2,0,0,0,0),
            Foam::vector::zero
        )
    );
}

Foam::tmp<Foam::volSymmTensorField> preciceAdapter::FSI::Stress::devRhoReff(dimensionedScalar rho) const
{
    // TODO: Only works for laminar flows at the moment.
    // See the OpenFOAM Stresss function object, where an extended
    // version of this method is being used.

    // Get the kinematic viscosity from the transportProperties
    // TODO: Get it from the objects' registry directly?
    const dictionary& transportProperties =
        mesh_.lookupObject<IOdictionary>("transportProperties");

    // TODO: In the pimpleDyMFoam tutorial (v1712), this is just a double,
    // which makes it fail.
    // TODO: Make this more general
    dimensionedScalar nu(transportProperties.lookup("nu"));

    // Get the velocity
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    return -rho * nu * dev(twoSymm(fvc::grad(U)));

}

void preciceAdapter::FSI::Stress::write(double * buffer)
{
    /* TODO: Implement
    * We need two nested for-loops for each patch,
    * the outer for the locations and the inner for the dimensions.
    * See the preCICE writeBlockVectorData() implementation.
    */
    // Compute stress. See the Stresss function object.
    // Normal vectors on the boundary, multiplied with the face areas
    // const surfaceVectorField::Boundary& nf (
    //     mesh_.Sf().boundaryField() / mesh_.magSf().boundaryField()
    //     );

    const surfaceVectorField nf(
            mesh_.Sf() / mesh_.magSf()
            );
    const surfaceVectorField::Boundary& nbf = nf.boundaryField();

    // TODO: Extend to cover also compressible solvers
    dimensionedScalar rho = mesh_.lookupObject<IOdictionary>("transportProperties").lookup("rho");

    // Stress tensor
    tmp<volSymmTensorField> tdevRhoReff = devRhoReff(rho);

    // Stress tensor boundary field
    const volSymmTensorField::Boundary& devRhoReffb =
        tdevRhoReff().boundaryField();
   
    // Pressure
    const volScalarField& p =
        mesh_.lookupObject<volScalarField>("p");


    int bufferIndex = 0;

    // Info<<  Stress_->boundaryFieldRef()[patchIDs_[0]]<< endl;
    
    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Pressure stresses
        // TODO: Extend to cover also compressible solvers
        Stress_->boundaryFieldRef()[patchID] =
            nbf[patchID] * p.boundaryField()[patchID]* rho.value();

        // Viscous stresses
        Stress_->boundaryFieldRef()[patchID] +=
            nbf[patchID] & devRhoReffb[patchID];

        // Write the stress to the preCICE buffer
        // For every cell of the patch
        forAll(Stress_->boundaryFieldRef()[patchID], i)
        {
            // Copy the stress into the buffer
            // x-dimension
            buffer[bufferIndex++]
            =
            Stress_->boundaryFieldRef()[patchID][i].x();

            // y-dimension
            buffer[bufferIndex++]
            =
            Stress_->boundaryFieldRef()[patchID][i].y();

            // z-dimension
            buffer[bufferIndex++]
            =
            Stress_->boundaryFieldRef()[patchID][i].z();
        }
    }
}

void preciceAdapter::FSI::Stress::read(double * buffer)
{
    /* TODO: Implement
    * We need two nested for-loops for each patch,
    * the outer for the locations and the inner for the dimensions.
    * See the preCICE readBlockVectorData() implementation.
    */
    FatalErrorInFunction
        << "Reading stresses is not supported."
        << exit(FatalError);
}

preciceAdapter::FSI::Stress::~Stress()
{
    // TODO: Is this enough?
    delete Stress_;
}
