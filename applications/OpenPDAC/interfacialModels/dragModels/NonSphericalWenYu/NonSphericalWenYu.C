/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "NonSphericalWenYu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(NonSphericalWenYu, 0);
    addToRunTimeSelectionTable(dragModel, NonSphericalWenYu, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::NonSphericalWenYu::NonSphericalWenYu
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject),
    sphericity_("sphericity", dimless, dict)    
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::NonSphericalWenYu::~NonSphericalWenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::NonSphericalWenYu::CdRe() const
{
    const volScalarField alpha2
    (
        max(1.0 - interface_.dispersed(), interface_.continuous().residualAlpha())
    );

    const volScalarField Res(alpha2*interface_.Re());

    const volScalarField CdsRes
    (
        neg(Res - 1000)*24*(1.0 + 0.15*pow(Res, 0.687))
      + pos0(Res - 1000)*0.44*Res
    );

    return
        CdsRes
       /sphericity_ 
       *pow(alpha2, -3.65)
       *max(interface_.continuous(), interface_.continuous().residualAlpha());
}


// ************************************************************************* //
