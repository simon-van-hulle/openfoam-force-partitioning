/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "Pstream.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "dimensionedTypes.H"
#include "fvMesh.H"
#include "fvPatch.H"
#include "surfaceFields.H"
#include "fpmForcePart.H"
#include "volFields.H"

#include "fvCFD.H"
#include "fvcCurl.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(fpmForcePart, 0);
        addToRunTimeSelectionTable(functionObject, fpmForcePart, dictionary);
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fpmForcePart::createFiles()
{
    if (writeToFile() && !fpmForcePartFilePtr_)
    {
        fpmForcePartFilePtr_ = createFile("fpmForcePartCoefficient");
        writeFileHeader("Menon and Mittal Force Partitioning Coefficients",
                        fpmForcePartFilePtr_());
    }
}

void Foam::functionObjects::fpmForcePart::writeFileHeader
(
    const word& header,
    Ostream& os
) const
{
    // TODO: Clean up headers when finished with the implementation.
    writeHeader(os, header);
    writeHeader(os, "");

    writeHeader(os, "IMPLEMENTATION NOT FINISHED - ALL RUBBISH");
    writeHeader(os, "Please do not use these results without verification of "
                    "the implementation.");
    writeHeader(os, "");

    writeHeaderValue(os, "magUInf", magUInf_);
    writeHeaderValue(os, "lRef", lRef_);
    writeHeaderValue(os, "liftDir", coordSys_.e2());
    writeHeaderValue(os, "dragDir", coordSys_.e1());
    writeHeaderValue(os, "forcePartDir", fpmForcePartDir_);
    writeHeaderValue(os, "ReynoldsNr", reynoldsNr_);

    writeHeader(os, "");

    writeCommented(os, "Time");
    writeTabbed(os, "Ckinematic");
    writeTabbed(os, "Cvorticity");
    writeTabbed(os, "Cviscous");
    writeTabbed(os, "Cpotential");
    writeTabbed(os, "Cboundary");
    writeTabbed(os, "Csum");
    os << nl;
}

void Foam::functionObjects::fpmForcePart::initialise()
{
    if (initialised_)
    {
        return;
    }

    // Other fieldnames can be added if they need to be checked as well
    wordList checkVecFields = {UName_};
    wordList checkScalarFields = {PhiName_, PhiIName_};

    for (const word& fieldName : checkVecFields)
    {
        if (!foundObject<volVectorField>(fieldName))
        {
            FatalErrorInFunction << "Could not find field " << fieldName
                                 << " in database" << exit(FatalError);
        }
    }
    for (const word& fieldName : checkScalarFields)
    {
        if (!foundObject<volScalarField>(fieldName))
        {
            FatalErrorInFunction << "Could not find field " << fieldName
                                 << " in database" << exit(FatalError);
        }
    }

    initialised_ = true;
}


void Foam::functionObjects::fpmForcePart::calcVivForcePartCoeffs()
{
    // Set coefficients to 0.0
    Ckinematic_ = 0.0;
    Cvorticity_ = 0.0;
    Cviscous_ = 0.0;
    Cpotential_ = 0.0;
    Cboundary_ = 0.0;

    // Read the velocity field
    // const volVectorField& UFull = lookupObject<volVectorField>(UName_);
    // auto U = (UFull / magUInf_).ref(); // FIXME:

    const volVectorField& U = lookupObject<volVectorField>(UName_);

    // Calculate time derivative of velocity, velocity squared and vorticity
    auto dUdt = fvc::ddt(U).cref();
    auto UU = (U & U).cref();
    auto vorticity = fvc::curl(U).cref();

    // Getting velocity potential field Phi
    const volScalarField& Phi = lookupObject<volScalarField>(PhiName_);

    // const volScalarField &PhiI = solvePhiI();
    const volScalarField& PhiI = lookupObject<volScalarField>(PhiIName_);

    // ------------------------------------------------------------------------
    // Helmholtz Decomposition
    // ------------------------------------------------------------------------
    auto uCurlFree = fvc::grad(Phi).cref();
    auto uDivFree = (U - uCurlFree).cref();

    // ========================================================================
    // KINEMATIC COEFFICIENT
    // ========================================================================
    for (const label patchi : patchSet_)
    {
        auto dUdtBoundary = dUdt.boundaryField()[patchi];
        auto UUBoundary = UU.boundaryField()[patchi];
        auto PhiIBoundary = PhiI.boundaryField()[patchi];
        auto SfBoundary = mesh_.Sf().boundaryField()[patchi];

        Ckinematic_ += sum
            (
                rhoRef_ * (SfBoundary & (dUdtBoundary * PhiIBoundary))
                + 0.5 * rhoRef_ * (UUBoundary * fpmForcePartDir_ & SfBoundary)
            );
    }

    // TODO: make sure fpmForcePartDir_ is a unit vector! otherwise this could
    // cause trouble

    // ========================================================================
    // VORTICITY COEFFICIENT
    // ========================================================================

    auto uExpr = (0.5 * (uDivFree & uDivFree) + (uDivFree & uCurlFree)).cref();
    auto divGradExpr = fvc::div(fvc::grad(uExpr) * PhiI).cref();

    // // FIXME: Reynolds ??
    Cvorticity_ += sum
    (
        - rhoRef_ * (
            (fvc::div(vorticity ^ U) * PhiI) 
            + (1 / reynoldsNr_) * divGradExpr
        ) * mesh_.V()
    );

    // ========================================================================
    // VISCOUS COEFFICIENT
    // ========================================================================

    for (const label patchi : patchSet_)
    {
        auto SfBoundary = mesh_.Sf().boundaryField()[patchi];
        auto omegaBoundary = vorticity.boundaryField()[patchi];

        Cviscous_ += sum
        (
            - mu_ * (
                ((omegaBoundary ^ SfBoundary) & fvc::grad(PhiI).cref()) 
                + ((omegaBoundary ^ SfBoundary) & fpmForcePartDir_)
            )
        );
    }

    // ========================================================================
    // POTENTIAL COEFFICIENT
    // ========================================================================

    auto gradTemp = fvc::grad(0.5 * (uCurlFree & uCurlFree)) * PhiI;
    Cpotential_ += sum(- rhoRef_ * fvc::div(gradTemp) * mesh_.V());

    // ========================================================================
    // BOUNDARY COEFFICIENT
    // ========================================================================

    // Cboundary_ += // TODO: How to do this? How to get to the boundaries?

    // ========================================================================
    // SUM OF ALL COEFFICIENT
    // ========================================================================

    // make the forces dimensionless
    scalar dimensions = 0.5 * rhoRef_ * magUInf_ * magUInf_ * lRef_;

    Ckinematic_ /= dimensions;
    Cvorticity_ /= dimensions;
    Cviscous_ /= dimensions;
    Cpotential_ /= dimensions;
    Cboundary_ /= dimensions;

    Csum_ = Ckinematic_ + Cvorticity_ + Cviscous_ + Cpotential_ + Cboundary_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fpmForcePart::fpmForcePart
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool readFields
)
: 
    forces(name, runTime, dict, false), 
    magUInf_(Zero), 
    lRef_(Zero),
    fpmForcePartFilePtr_(), 
    UName_("U"), 
    PhiName_("Phi"), 
    PhiIName_("PhiI"),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict, "liftDir", "dragDir");
        Info << nl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fpmForcePart::read(const dictionary& dict)
{
    // TODO: What is needed from forces? Don't remove.
    forces::read(dict);

    // Free stream velocity magnitude
    dict.readEntry("magUInf", magUInf_);

    // Reference length // TODO: This should be the diameter.
    dict.readEntry("lRef", lRef_);

    // Reference Density
    dict.readEntry("rhoRef", rhoRef_);

    // Read the direction for force partitioning
    dict.readEntry("fpmForcePartDir", fpmForcePartDir_);
    fpmForcePartDir_.normalise();

    // What to do for boundary Coefficient
    patchSet_ = mesh_.boundaryMesh().patchSet(dict.get<wordRes>("patches"));

    // Get Information from transportProperties
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties", 
            mesh_.time().constant(), 
            mesh_,
            IOobject::MUST_READ, 
            IOobject::NO_WRITE
        )
    );

    // Viscosity
    nu_ = dimensionedScalar("nu", dimViscosity, transportProperties).value();

    // Checking if controlDict specifies this force partitioning.
    fpmForcePartEnabled = dict.getOrDefault("enabled", true);

    if (fpmForcePartEnabled)
    {
        Info << "    VIV Force Partitioning Method switched ON" << nl;
    }
    else
    {
        Info << "    VIV Force Partitioning switched OFF - Not including "
             << "additional force partitioning" << nl;
    }

    return true;
}

bool Foam::functionObjects::fpmForcePart::execute()
{
    // Don't do anything if user as specified to switch off.
    if (!fpmForcePartEnabled)
    {
        return true;
    }

    createFiles();
    initialise();

    // Calculate Reynolds number
    reynoldsNr_ = magUInf_ * lRef_ / nu_;
    mu_ = nu_ * rhoRef_;
    Info << "Reynolds number = " << reynoldsNr_ << nl;

    // Calculate the coefficients
    calcVivForcePartCoeffs();

    return true;
}

bool Foam::functionObjects::fpmForcePart::write()
{
    if (log)
    {
        Log << type() << " " << name() << " execute:" << nl
            << "    Coefficients" << nl;

        Log << "        Ckinematic : " << Ckinematic_ << nl
            << "        Cvorticity : " << Cvorticity_ << nl
            << "        Cviscous   : " << Cviscous_ << nl
            << "        Cpotential : " << Cpotential_ << nl
            << "        Cboundary  : " << Cboundary_ << nl
            << "        Csum       : " << Csum_ << nl << nl;
    }

    if (writeToFile())
    {
        writeCurrentTime(fpmForcePartFilePtr_());
        fpmForcePartFilePtr_()
            << tab << Ckinematic_ << tab << Cvorticity_ << tab << Cviscous_
            << tab << Cpotential_ << tab << Cboundary_ << tab << Csum_ << endl;
    }

    return true;
}

// ************************************************************************* //
