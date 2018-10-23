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
time_(-1.),
// timeOld_(0.0),
velocity_(
    const_cast<volVectorField*>
    (
        &mesh.lookupObject<volVectorField>(nameVelocity)
    )
)
{
    dataType_ = vector;
    // oldpoints = new pointField();
    // CREATE COPY OF THE POINTDISPLACEMENT FIELD?
    oldPointDisplacement_ = new pointVectorField
    (
        IOobject
        (
            "oldPointDisplacement",
            runTime_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "displacement",
            dimensionSet(0,1,0,0,0,0,0),
            Foam::vector::zero
        )
    );
          // dimensionedVector("displacement", dimLength, Zero),


     // DisplOld_ = new vectorField(Velocity_->boundaryFieldRef().size()*3, Foam::vector::zero);
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

Foam::tmp<Foam::scalarField> preciceAdapter::FSI::Velocity::sweptVols
(
    const pointField& newPoints,
    const pointField& oldPoints,
    const faceList& f
)
{

    // Create swept volumes
    // const faceList& f = faces();

    tmp<scalarField> tsweptVols(new scalarField(f.size()));
    scalarField& sweptVols = tsweptVols.ref();
    Info << "HALLO1" << endl;
    forAll(f, facei)
    {
        sweptVols[facei] = f[facei].sweptVol(oldPoints, newPoints);
    }
    Info << "HALLO2" << endl;
    // Force recalculation of all geometric data with new points
    // clearGeom();

    return tsweptVols;
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


    // COLLECT ALL THE RELEVANT FIELDS.     
    const pointVectorField& pointDisplacement_ =
        mesh_.lookupObject<pointVectorField>("pointDisplacement");

    // Normal vectors on the boundary, multiplied with the face areas
    const surfaceVectorField nf(
    mesh_.Sf() / mesh_.magSf()
    );


    const surfaceScalarField& phi_ =
        mesh_.lookupObject<surfaceScalarField>("phi");


    for (uint j = 0; j < patchIDs_.size(); j++)
    {

        int patchID = patchIDs_.at(j);

        // TODO: improve this part.
        
        // TRY. Assign the old pointdisplacemetn only on the first iteration. 
        
        if (time_!=runTime_.value())
        {
            // check if the function needs to be called.
            time_ = runTime_.value();
            *oldPointDisplacement_ = *pointDisplacement_;
            // save the old point displacement
            // const pointField& oldpoints = pointDisplacement_.oldTime().boundaryField()[patchID].patchInternalField();
        }

        

        // const surfaceVectorField::Boundary& n = nf.boundaryField()[patchID];
        // tmp<vectorField> n = nf.boundaryFieldRef()[patchID];
        // boundaryField()[patchID];

        // mesh.boundaryMesh()[patchIDs_.at(j)].faceCentres()

        // find out how to do this.
        const fvPatch& p = mesh_.boundary()[patchID]; //patch();

        const polyPatch& pp = p.patch();
        // const polyPatch& pp = mesh_.boundaryMesh()[patchID];
        const pointField& oldPoints = mesh_.oldPoints();
    

        const pointField& newpoints2 = pp.localPoints();

        const pointField& oldpoints2 = newpoints2;
        
        vectorField oldFc(pp.size());

        forAll(oldFc, i)
        {
            // pp[i] is a face, swith 
            oldFc[i] = pp[i].centre(oldPoints);
        }
        

        // const pointField& oldpoints2(newpoints2.size());
        // oldPoints2 = pp.localPoints(oldPoints);

        // Info << "new points: " << pp.faceCentres() << endl;
        // Info << "old points: " << oldFc << endl;
        // Info << "new2 points " << p.patch().faceCentres() << endl;
        // const scalar& deltaT = runTime_.deltaT().value(); //runTime_.deltaTValue();

        const vectorField Up((pp.faceCentres() - oldFc)/runTime_.deltaT().value());

        // -------------------------------------------------------
        //          TRY THIS FIRST

        // const pointField& newpoints = mesh_.points().oldTime().boundaryMesh()[patchID]
        // fixedValuePointPatchVectorField& pointDisplacementFluidPatch =
        //     refCast<fixedValuePointPatchVectorField>
        //     (
        //         pointDisplacement_->boundaryFieldRef()[patchID]
        //     );
        
        //COMPARE
        // const vector Uw1 = U_.boundaryField()[patchi][patchFacei];
        // const vector& Uw0 =  U_.oldTime().boundaryField()[patchi][patchFacei];

        // Info << nl << "old value: " << max(pointDisplacement_.oldTime().primitiveField()) << endl;
        // const scalarField TpOld(T.oldTime().boundaryField()[patch().index()]);
        // Info << "Pointfield:  " << pointDisplacement_.boundaryField()[patchID].patchInternalField() << endl;
        const pointField& newpoints = pointDisplacement_.boundaryField()[patchID].patchInternalField();
        const pointField& oldpoints = oldPointDisplacement_->boundaryField()[patchID].patchInternalField();
        const faceList& f = pp.localFaces();

        Info << "mesh new points " << newpoints[1] << nl << endl;
        Info << "mesh old points " << oldpoints[1] << nl << endl;
        //GET THIS PART WORKING
        // tmp<scalarField> sweptVols = primitiveMesh::movePoints
        // (
        //     pointDisplacement_.boundaryFieldRef()[patchID].patchInternalField(),
        //     pointDisplacement_.oldTime().boundaryFieldRef()[patchID].patchInternalField()
        // );

/*  
        tmp<scalarField> tsweptVols = sweptVols
        (
            newpoints,
            oldpoints,
            f
        );
*/
        // Info << "sweptVols" << tsweptVols << endl;


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

        // velocityPatch = (Up + n*(Un - (n & Up)));

        // Info << "velocity on boundary: " << velocityPatch << endl;
    }


}



// //- Destructor
// preciceAdapter::FSI::Velocity::~Velocity()
// {
//     // TODO: Is this enough?
//     delete faceDisplacement_;
//     delete faceDisplacementOld_;
// }
