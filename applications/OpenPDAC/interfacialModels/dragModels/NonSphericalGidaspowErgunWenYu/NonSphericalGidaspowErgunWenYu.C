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

#include "NonSphericalGidaspowErgunWenYu.H"
#include "Ergun.H"
#include "WenYu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(NonSphericalGidaspowErgunWenYu, 0);
    addToRunTimeSelectionTable(dragModel, NonSphericalGidaspowErgunWenYu, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::NonSphericalGidaspowErgunWenYu::NonSphericalGidaspowErgunWenYu
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject),
    sphericity_("sphericity", dimless, dict),
    Ergun_(dict, interface, false),
    WenYu_(dict, interface, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::NonSphericalGidaspowErgunWenYu::~NonSphericalGidaspowErgunWenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::NonSphericalGidaspowErgunWenYu::CdRe() const
{

    const phaseModel& dispersed = interface_.dispersed();
    const phaseModel& continuous = interface_.continuous();

    const volScalarField & Ergun_CdRe = (4.0/3.0)
       *(
            150
           *max(1.0 - continuous, dispersed.residualAlpha())
           /max(continuous, continuous.residualAlpha())
           / sqr(sphericity_)
          + 1.75*interface_.Re() / sphericity_
        );

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
    
    const volScalarField & WenYu_CdRe = CdsRes / sphericity_
       *pow(alpha2, -3.65)
       *max(interface_.continuous(), interface_.continuous().residualAlpha());

    return
        pos0(interface_.continuous() - 0.8)*WenYu_CdRe
      + neg(interface_.continuous() - 0.8)*Ergun_CdRe;
}


// ************************************************************************* //
