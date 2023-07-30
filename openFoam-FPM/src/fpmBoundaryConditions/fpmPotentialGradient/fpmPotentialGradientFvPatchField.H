/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fpmPotentialGradientFvPatchField

Group
    grpGenericBoundaryConditions

Description
    Boundary condition for the auxiliary potential that is used for
    the Force Partitioning method described in Mennon(2020).

        \f[
            n\cdot \nabla\phi^{(i)} = n_i
        \f]

    where
    \vartable
        n           | normal vector to the boundary surface
        i           | direction of force partitioning
        \nabla(\phi)| gradient
    \endvartable

Usage
    \table
        Property        | Description                     | Required    | Default value
        forcePartDir    | Direction of force partitioning | yes         | 1
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fpmPotentialGradient;
        forcePartDir    1;
    }
    \endverbatim

SourceFiles
    fpmPotentialGradientFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fpmPotentialGradientFvPatchField_H
#define fpmPotentialGradientFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fpmPotentialGradientFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fpmPotentialGradientFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

        vector gradientVec_;

        Field<Type> gradient_;

public:

    //- Runtime type information
    TypeName("fpmPotentialGradient");


    // Constructors

        //- Construct from patch and internal field
        fpmPotentialGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fpmPotentialGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given fpmPotentialGradientFvPatchField
        //  onto a new patch
        fpmPotentialGradientFvPatchField
        (
            const fpmPotentialGradientFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fpmPotentialGradientFvPatchField
        (
            const fpmPotentialGradientFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new fpmPotentialGradientFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fpmPotentialGradientFvPatchField
        (
            const fpmPotentialGradientFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new fpmPotentialGradientFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Return defining fields

            //- Return gradient at boundary
            virtual Field<Type>& gradient()
            {
                return gradient_;
            }

            virtual const Field<Type>& gradient() const
            {
                return gradient_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            // virtual void rmap
            // (
            //     const fvPatchField<Type>&,
            //     const labelList&
            // );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const
            {
                return gradient_;
            }

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fpmPotentialGradientFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
