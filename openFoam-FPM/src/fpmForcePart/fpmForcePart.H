/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

Class
    Foam::functionObjects::fpmForcePart

Group
    grpForcesFunctionObjects

Description
    Extends the \c forces functionObject by provideing additional force
    partitioning as proposed by Menon and Mittal (2020).
    - Kinematic forces                          (Ckinematic)
    - Vorticity-induced forces                  (Cvorticity)
    - Forces due to viscous effects             (Cviscous)
    - Potential flow field forces               (Cpotential)
    - Forces induced by the boundary boundary   (Cboundary)

Usage
Example of function object specification
\verbatim
forceCoeffs1
    {
        type        forceCoeffs;
        libs        ("VIV_libfpmForcePart.so");
        ...
        log         yes;
        writeFields yes;
        patches     (CYLINDER);

        // input keywords for directions of force/moment coefficients
        // refer below for options, and relations

        magUInf     100;
        lRef        3.5;
        Aref        2.2;

        CofR        (0 0 0); // Centre of rotation
        dragDir     (1 0 0);
        liftDir     (0 0 1);
    }
\endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::timeControl
    Foam::functionObjects::forces::

SourceFiles
    fpmForcePart.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fpmForcePart_H
#define functionObjects_fpmForcePart_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {

        /*--------------------------------------------------------------------*\
                                 Class fpmForcePart Declaration
        \*--------------------------------------------------------------------*/

        class fpmForcePart : public forces
        {
            // Private data

            scalar currentTime_;

            // Data
            //- Free-stream velocity magnitude
            scalar magUInf_;

            //- Reference length [m]
            scalar lRef_;

            //- Kinematic viscosity [m^2/s]
            scalar nu_;

            //- Dynamic viscosity [kg/m/s]
            scalar mu_;

            //- Reynolds number [-]
            scalar reynoldsNr_;

            // Direction of force;
            vector fpmForcePartDir_; // Direction for force partitioning

            // Coefficients
            scalar Ckinematic_;
            scalar Cvorticity_;
            scalar Cviscous_;
            scalar Cpotential_;
            scalar Cboundary_;
            scalar Csum_;

            // File streams

            //- Integrated coefficients
            autoPtr<OFstream> fpmForcePartFilePtr_;

            // Read from dictionary

            //- Patches to integrate forces over
            labelHashSet patchSet_;

            //- Name of velocity field
            word UName_;

            //- Name of velocity Potential
            word PhiName_;

            //- Name of the auxiliary potential (harmonic function)
            word PhiIName_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Coordinate system used when evaluating forces/moments
            coordSystem::cartesian coordSys_;

            //- Flag to include porosity effects
            bool porosity_;

            //- Flag after initialisation
            bool initialised_;

            // Private Member Functions
            // TODO: Do I need other private member functions?

            //- No copy construct
            fpmForcePart(const fpmForcePart&) = delete;

            //- No copy assignment
            void operator=(const fpmForcePart&) = delete;

        protected:
            // Protected Data

            // Read from Dictionary

            bool fpmForcePartEnabled; // Enable/disable force partitioning

            // Protected Member Functions
            // TODO: Add protected member functions

            //- Create output files
            void createFiles();

            //- Write header for the data
            void writeFileHeader(const word& header, Ostream& os) const;

            // Initialise and check if fields present
            void initialise();

            // Solve the auxiliary potential (harmonic function)
            const volScalarField solvePhiI();

            // Calculate coefficients (Mennon, 2020)
            void calcVivForcePartCoeffs();

        public:
            //- Runtime type information
            TypeName("fpmForcePart");

            // TODO: Adjust constructor
            // Constructors

            //- Construct from Time and dictionary
            fpmForcePart(const word& name, const Time& runTime,
                         const dictionary&, const bool readFields = true);

            //- Destructor
            virtual ~fpmForcePart() = default;

            // Member Functions

            //- Read the forces data
            virtual bool read(const dictionary&);

            //- Execute
            virtual bool execute();

            //- Write the forces
            virtual bool write();
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        // *
        // * * //

    } // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
