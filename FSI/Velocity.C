#include "Velocity.H"

using namespace Foam;


preciceAdapter::FSI::Velocity::Velocity
(
    const Foam::fvMesh& mesh,
    const Foam::Time& runTime,
    const std::string nameVelocity
    /* TODO: We should add any required field names here.
    /  They would need to be vector fields.
    /  See CHT/Temperature.C for details.
    */
)
:
mesh_(mesh),
runTime_(runTime),
// time_(0.0),
// timeOld_(0.0),
// velocity_(
//     const_cast<volVectorField*>
//     (
//         &mesh.lookupObject<volVectorField>(nameVelocity)
//     )
velocity_(
    const_cast<volVectorField*>
    (
        &mesh.lookupObject<volVectorField>(nameVelocity)
    )
)
{
    dataType_ = vector;
    // const scalar deltaT
    // // Initialize the current and old face displacement as two IOOBJECTS.
    // // Maybe this is not required, and can be solved in another way.
    // faceDisplacement_ = new volVectorField
    // (
    //     IOobject
    //     (
    //         "faceDisplacement",
    //         runTime_.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedVector
    //     (
    //         "lengthdim",
    //         dimensionSet(0,1,0,0,0,0,0),
    //         Foam::vector::zero
    //     )
    // );

    // faceDisplacementOld_ =  new volVectorField
    // (
    //     IOobject
    //     (
    //         "faceDisplacement_old",
    //         runTime_.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedVector
    //     (
    //         "lengthdim",
    //         dimensionSet(0,1,0,0,0,0,0),
    //         Foam::vector::zero
    //     )
    // );
}

void preciceAdapter::FSI::Velocity::write(double * buffer)
{
    /* TODO: Implement
    * FOR NOW ONLY WORKS IF THE DISPLACEMENT FIELD IS ALREADY UPDATED.
    * Make this function not dependent on the buffer, but rather on the faceDisplacement
    * This can be a function in the displacement.C
    * Create velocity interpolation
    */
    FatalErrorInFunction
        << "Writing velocities is not supported."
        << exit(FatalError);
}

void preciceAdapter::FSI::Velocity::read(double * buffer)
{
    /* TODO: Implement
    * FOR NOW ONLY WORKS IF THE DISPLACEMENT FIELD IS ALREADY UPDATED.
    * check $FOAM_SRC/finiteVolume/fields/fvPatchFields/derived/movingWallVelocity
    * Or myMovingWallVelocity from David Blom for a second order interpolation.
    */
    // if (time_!=runTime_.value())
    // {
    //     // check if the function needs to be called.
    //     timeOld_ = time_;
    //     time_ = runTime_.value();

    //     // save the old displaement at the faceCentres.
    //     *faceDisplacementOld_ = *faceDisplacement_;
    // }
    // if ( time_ == 0 )
    // {
    //     timeOld_ = -1.;
    // }


    // const fvMesh& mesh = internalField().mesh();

    if (mesh_.moving())
    {
        // Normal vectors on the boundary, multiplied with the face areas

        const surfaceVectorField nf(
        mesh_.Sf() / mesh_.magSf()
        );

        // const surfaceScalarField magSf_
        // (
        //     mesh_.magSf()
        // );

        // const volVectorField& U_ =
        //     mesh_.lookupObject<volVectorField>("U");

        const surfaceScalarField& phi_ =
            mesh_.lookupObject<surfaceScalarField>("phi");


        for (uint j = 0; j < patchIDs_.size(); j++)
        {

            int patchID = patchIDs_.at(j);

            // const surfaceVectorField::Boundary& n = nf.boundaryField()[patchID];
            // tmp<vectorField> n = nf.boundaryFieldRef()[patchID];
            // boundaryField()[patchID];

            // mesh.boundaryMesh()[patchIDs_.at(j)].faceCentres()

            // find out how to do this.
            const fvPatch& p = mesh_.boundary()[patchID]; //patch();

            const polyPatch& pp = p.patch();
            // const polyPatch& pp = mesh_.boundaryMesh()[patchID];
            const pointField& oldPoints = mesh_.oldPoints();

            vectorField oldFc(pp.size());

            forAll(oldFc, i)
            {
                oldFc[i] = pp[i].centre(oldPoints);
            }
            
            Info << "new points: " << pp.faceCentres() << endl;
            Info << "old points: " << oldFc << endl;
            Info << "new2 points " << p.patch().faceCentres() << endl;
            // const scalar& deltaT = runTime_.deltaT().value(); //runTime_.deltaTValue();

            const vectorField Up((pp.faceCentres() - oldFc)/runTime_.deltaT().value());

// -------------------------------------------------------
//          TRY THIS FIRST
            const pointVectorField& pointDisplacement_ =
                mesh_.lookupObject<pointVectorField>("pointDisplacement");

            fixedValuePointPatchVectorField& pointDisplacementFluidPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    pointDisplacement_->boundaryFieldRef()[patchID]
                );
            
            // Info << nl << "old value: " << max(pointDisplacement_.oldTime().primitiveField()) << endl;
            
            tmp<scalarField> sweptVols = primitiveMesh::movePoints(
                    pointDisplacementFluidPatch,
                    pointDisplacementFluidPatch.oldTime()
                );

            Info << "sweptVols" << sweptVols << endl;


// -------------------------------------------------------
            // THIS IS THE FLUX ON THE SURFACE. 
            // surfaceScalarField phip 
            // (
            //     // phi_.boundaryFieldRef()[patchID]
            //     p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(velocity_))
            // );

            // surfaceScalarField& phip =
            //     refCast<surfaceScalarField>
            //     (
            //         phi_.boundaryFieldRef()[patchID]
            //     );

            // fvsPatchField<double>& phip = 
            //     phi_.boundaryField()[patchID];

            // const surfaceScalarField& phip =
            //     p.patchField<surfaceScalarField, scalar>(phi_);

            // try if it works with thr fvsPatchScalarField insatead of the surfaceScalarField
            const fvsPatchScalarField& phip = phi_.boundaryField()[patchID];


            // Info << "mesh_ update: " << mesh_.moving() << endl;
            // Info << "old faceCentres: " <<  oldFc[15]  << endl;
            // Info << "new faceCentres: " <<  pp.faceCentres()[15] << endl;

            // scalarField& magSfb = 
            //     refCast<scalarField>
            //     (
            //         magSf_.boundaryFieldRef()[patchID]
            //     );

            // // const scalarField& magSf = pp.magSf();

            // tmp<scalarField> Un = phip/(magSfb + VSMALL);
            const vectorField n(p.nf());
            const scalarField& magSf = p.magSf();
            tmp<scalarField> Un = phip/(magSf + VSMALL);

            // volVectorField::Boundary& Ubf = U.boundaryFieldRef();

            // volVectorField::Boundary& velocityBf =
            //         velocity_->boundaryFieldRef();

            // volVectorField& velocityPatch =
            //     refCast<volVectorField>
            //     (
            //         velocity_->boundaryFieldRef()[patchID]
            //     );

            // vectorField::operator=(Up + n*(Un - (n & Up)));
            // fvPatchVectorFieldField::operator=(Up + n*(Un - (n & Up)));
            
            // Info << nl << "velocity Un: " << Un << endl; 
               // Get the pointMotionU
            vectorField& velocityPatch =
               refCast<vectorField>
                (
                    velocity_->boundaryFieldRef()[patchID]
                );


            // velocityBf[patchID] = 
            //     (Up + n*(Un - (n & Up)));

            velocityPatch = (Up + n*(Un - (n & Up)));

            Info << "velocity on boundary: " << velocityPatch << endl;
        }

        // fixedValueFvPatchVectorField::updateCoeffs();
    }
}



// //- Destructor
// preciceAdapter::FSI::Velocity::~Velocity()
// {
//     // TODO: Is this enough?
//     delete faceDisplacement_;
//     delete faceDisplacementOld_;
// }
