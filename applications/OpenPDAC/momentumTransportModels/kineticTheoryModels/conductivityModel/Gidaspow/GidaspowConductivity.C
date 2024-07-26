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

#include "GidaspowConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace conductivityModels
{
    defineTypeNameAndDebug(Gidaspow, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        Gidaspow,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::Gidaspow::Gidaspow
(
    const dictionary& dict
)
:
    conductivityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::Gidaspow::~Gidaspow()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::conductivityModels::Gidaspow::kappa
(
    const volScalarField& alpha1,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    return rho1*da*sqrt(Theta)*
    (
        2*sqr(alpha1)*g0*(1 + e)/sqrtPi
      + (9.0/8.0)*sqrtPi*g0*0.5*(1 + e)*sqr(alpha1)
      + (15.0/16.0)*sqrtPi*alpha1
      + (25.0/64.0)*sqrtPi/((1 + e)*g0)
    );
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::conductivityModels::Gidaspow::kappa
(
    const volScalarField& alpha1,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& sumAlphaGs0,
    const volScalarField& beta,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    const scalar Pi = constant::mathematical::pi;

    // Eq. B12 MFIX2012
    const dimensionedScalar eta = 0.5*(1.0 + e);

    // Eq. B10 MFIX2012
    const volScalarField kappa = ( 75*rho1*da*sqrtPi*sqrt(Theta) )
                                 / ( 48*eta*(41-33*eta) );
    // Eq. B9 MFIX2012
    const volScalarField kappaStar = ( rho1*alpha1*g0*Theta*kappa )
                                     / ( rho1*sumAlphaGs0*Theta +
                                         ( 6*beta*kappa ) / ( 5*rho1*alpha1 ) );
    // Eq. B8 MFIX2012
    return kappaStar/g0 *
    (
        ( 1+12/5*eta*sumAlphaGs0 ) *
        ( 1+12/5*sqr(eta)*(4*eta-3)*sumAlphaGs0 ) +
        64/(25*Pi)*(41-33*eta)*sqr(eta)*sqr(sumAlphaGs0)
    );
}

// ************************************************************************* //
