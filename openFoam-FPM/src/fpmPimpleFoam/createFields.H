Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// -----------------------------------------------------------------------------
// TODO: Check if this gives the correct results
// NOTE: This part is taken from the createFields.H in potentialFoam!!
// Infer the velocity potential BCs from the pressure

wordList PhiBCTypes
(
    p.boundaryField().size(),
    zeroGradientFvPatchScalarField::typeName
);

forAll(p.boundaryField(), patchi)
{
    if (p.boundaryField()[patchi].fixesValue())
    {
        PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}



Info<< "Constructing velocity potential field Phi\n" << endl;
volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE // TODO: Don't write this
    ),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, Zero),
    PhiBCTypes
);

volScalarField PhiI
(
    IOobject
    (
        "PhiI",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE // TODO: Don't write this
    ),
    mesh
);



label PhiRefCell = 0;
scalar PhiRefValue = 0;
setRefCell
(
    Phi,
    pimple.dict(),
    PhiRefCell,
    PhiRefValue
);
// mesh.setFluxRequired(Phi.name());

// -----------------------------------------------------------------------------

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
