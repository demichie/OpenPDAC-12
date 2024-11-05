/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::heatTransferModels::CoeffModel

Description
    Zhu et al. (2019) correlation for fixed and fluidised beds with Re < 10^5
    and continuous phase fraction between 0.35 and 1.

    Reference:
    \verbatim
        Zhu, Li-Tao, Yuan-Xing Liu, and Zheng-Hong Luo. 
        "An enhanced correlation for gas-particle heat and mass transfer in 
         packed and fluidized bed reactors." 
         Chemical Engineering Journal 374 (2019): 531-544.
    \endverbatim

SourceFiles
    CoeffModel.C

\*---------------------------------------------------------------------------*/

#ifndef CoeffModel_H
#define CoeffModel_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                            Class CoeffModel Declaration
\*---------------------------------------------------------------------------*/

class CoeffModel
:
    public heatTransferModel
{
    // Private Data

        //- Interface
        const dispersedPhaseInterface interface_;
        const dimensionedScalar a10_;
        const dimensionedScalar a11_;
        const dimensionedScalar a12_;
        const dimensionedScalar a20_;
        const dimensionedScalar a21_;
        const dimensionedScalar a22_;


public:

    //- Runtime type information
    TypeName("CoeffModel");


    // Constructors

        //- Construct from a dictionary and an interface
        CoeffModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~CoeffModel();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
