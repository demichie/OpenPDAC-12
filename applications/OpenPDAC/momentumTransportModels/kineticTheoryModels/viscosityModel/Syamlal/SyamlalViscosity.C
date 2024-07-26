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

#include "SyamlalViscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace viscosityModels
{
    defineTypeNameAndDebug(Syamlal, 0);
    addToRunTimeSelectionTable(viscosityModel, Syamlal, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::viscosityModels::Syamlal::Syamlal
(
    const dictionary& dict
)
:
    viscosityModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    alfa_
    (
        "alfa",
        dimless,
        coeffDict_.lookupOrDefault<scalar>("alfa", 1.6)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::viscosityModels::Syamlal::~Syamlal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::viscosityModels::Syamlal::nu
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

    return volScalarField::New
    (
        IOobject::groupName
        (
            Foam::typedName<viscosityModel>("nu"),
            Theta.group()
        ),
        da*sqrt(Theta)
       *(
            (4.0/5.0)*alpha1*g0*(1 + e)/sqrtPi
          + (1.0/15.0)*sqrtPi*g0*(1 + e)*(3*e - 1)*alpha1/(3 - e)
          + (1.0/6.0)*sqrtPi/(3 - e)
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::viscosityModels::Syamlal::nu
(
    const volScalarField& alpha1,
    const volScalarField& Theta,
    const dimensionedScalar& ThetaSmall,
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
    
    // Eq. B6 MFIX2012    
    const volScalarField mu = 5.0/96.0*rho1*da*sqrt(Theta)*sqrtPi; 

    // Eq. B7 MFIX2012
    const volScalarField mu_b = 256.0/(5.0*Pi)*mu*alpha1*sumAlphaGs0;

    // Eq. B5 MFIX2012
    const volScalarField muStar = ( rho1*alpha1*g0*Theta*mu ) /
                                  ( rho1*sumAlphaGs0*Theta + 
                                    (2*beta*mu)/(rho1*alpha1) );
                                     
    // Eq. B4 MFIX2012
    return volScalarField::New
    (
        IOobject::groupName
        (
            Foam::typedName<viscosityModel>("nu"),
            Theta.group()
        ),
        (2+alfa_)/3.0*( 
         muStar / (g0*eta*(2-eta))*
         (1+8/5*eta*sumAlphaGs0)*
         (1+8/5*eta*(3*eta-2)*sumAlphaGs0)+
         3/5*eta*mu_b )
    );
}


// ************************************************************************* //
