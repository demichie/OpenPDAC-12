/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "mathematicalConstants.H"
#include "PressureGradientPForce.H"
#include "fvcDdt.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PressureGradientPForce<CloudType>::PressureGradientPForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    ParticleForce<CloudType>(owner, mesh, dict, forceType, true),
    PName_(this->coeffs().template lookupOrDefault<word>("P", "P")),
    gradPInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::PressureGradientPForce<CloudType>::PressureGradientPForce
(
    const PressureGradientPForce& pgf
)
:
    ParticleForce<CloudType>(pgf),
    PName_(pgf.PName_),
    gradPInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PressureGradientPForce<CloudType>::~PressureGradientPForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PressureGradientPForce<CloudType>::cacheFields(const bool store)
{
    static word fName("gradP");

    bool fieldExists = this->mesh().template foundObject<volVectorField>(fName);

    if (store)
    {
        if (!fieldExists)
        {
            const volScalarField& Pc = this->mesh().template
                lookupObject<volScalarField>(PName_);

            volVectorField* gradPPtr = new volVectorField
            (
                fName,
                fvc::grad(Pc)
            );

            gradPPtr->store();
        }

        const volVectorField& gradP = this->mesh().template
            lookupObject<volVectorField>(fName);

        gradPInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                this->owner().solution().interpolationSchemes(),
                gradP
            ).ptr()
        );
    }
    else
    {
        gradPInterpPtr_.clear();

        if (fieldExists)
        {
            const volVectorField& gradP = this->mesh().template
                lookupObject<volVectorField>(fName);

            const_cast<volVectorField&>(gradP).checkOut();
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::PressureGradientPForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    const vector gradP =
        gradPInterp().interpolate
        (
            p.coordinates(),
            p.currentTetIndices(td.mesh)
        );

    value.Su() = - pow3(p.d())*mathematical::pi/6.0*gradP;

    return value;
}


template<class CloudType>
Foam::scalar Foam::PressureGradientPForce<CloudType>::massAdd
(
    const typename CloudType::parcelType&,
    const typename CloudType::parcelType::trackingData& td,
    const scalar
) const
{
    return 0.0;
}


// ************************************************************************* //
