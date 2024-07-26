/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "NonSphereDragVarPhiForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NonSphereDragVarPhiForce<CloudType>::NonSphereDragVarPhiForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true)
{}


template<class CloudType>
Foam::NonSphereDragVarPhiForce<CloudType>::NonSphereDragVarPhiForce
(
    const NonSphereDragVarPhiForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NonSphereDragVarPhiForce<CloudType>::~NonSphereDragVarPhiForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::NonSphereDragVarPhiForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{

    const scalar a = exp(2.3288 - 6.4581*p.phi() 
                                + 2.4486*sqr(p.phi()));
    const scalar b = 0.0964 + 0.5565*p.phi();
    const scalar c = exp(4.9050 - 13.8944*p.phi() 
                         + 18.4222*sqr(p.phi()) - 10.2599*pow3(p.phi()));
    const scalar d = exp(1.4681 + 12.2584*p.phi() 
                         - 20.7322*sqr(p.phi()) + 15.8855*pow3(p.phi()));

    const scalar CdRe =
        24*(1 + a*pow(Re, b)) + Re*c/(1 + d/(Re + rootVSmall));

    return forceSuSp(Zero, mass*0.75*muc*CdRe/(p.rho()*sqr(p.d())));
}


// ************************************************************************* //
