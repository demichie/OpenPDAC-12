/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::dragModels::solidSolidDrag

Description
    MFIX solid-solid drag

    Reference:

SourceFiles
    solidSolidDrag.C

\*---------------------------------------------------------------------------*/

#ifndef solidSolidDrag_H
#define solidSolidDrag_H

#include "mathematicalConstants.H"
#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                      Class solidSolidDrag Declaration
\*---------------------------------------------------------------------------*/

class solidSolidDrag
:
    public dragModel
{
    // Private Data

        //- Interface
        const phaseInterface interface_;

        //- Name of the gaseous phase
        const word gasName_;

        //- Name of the solid phase1
        const word solid1Name_;

        //- Name of the solid phase
        const word solid2Name_;

        //- Coefficient of restitution
        const dimensionedScalar E_;

        //- Coefficient of friction
        const dimensionedScalar Cf_;


    // Private Member Functions

        //- Return the momentum transfer coefficient between solid1 and solid2
        virtual tmp<volScalarField> KSolidSolid
        (
            const phaseModel& gas,
            const phaseModel& solid1,
            const phaseModel& solid2
        ) const;


public:

    //- Runtime type information
    TypeName("solidSolidDrag");


    // Constructors

        //- Construct from a dictionary and an interface
        solidSolidDrag
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~solidSolidDrag();


    // Member Functions

        //- The drag coefficient used in the momentum equation
        virtual tmp<volScalarField> K() const;

        //- The drag coefficient used in the face-momentum equations
        virtual tmp<surfaceScalarField> Kf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
