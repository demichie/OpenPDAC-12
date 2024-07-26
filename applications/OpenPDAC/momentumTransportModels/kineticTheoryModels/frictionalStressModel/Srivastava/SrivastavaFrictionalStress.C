/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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

#include "SrivastavaFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(Srivastava, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        Srivastava,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
Srivastava::Srivastava
(
    const dictionary& dict
)
:
    frictionalStressModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict_),
    eta_("eta", dimless, coeffDict_),
    p_("p", dimless, coeffDict_),
    phi_("phi", dimless, coeffDict_),
    alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict_)
{
    phi_ *= constant::mathematical::pi/180.0;
    sinPhi_ = sin(phi_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
Srivastava::~Srivastava()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
Srivastava::frictionalPressure
(
    const phaseModel& phase,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax
) const
{
    const volScalarField alphas = 1.0 - continuousPhase;

    const volScalarField JohnsonJacksonP = Fr_*pow(max(alphas - alphaMinFriction, scalar(0)), eta_)
       /pow(max(alphasMax - alphas, alphaDeltaMin_), p_);

    
    const volScalarField LinearP = Fr_*(
				pow(alphasMax-alphaDeltaMin_-alphaMinFriction,eta_)/pow(alphaDeltaMin_,p_) +
				eta_ * pow(alphasMax-alphaDeltaMin_-alphaMinFriction,eta_-1)/pow(alphaDeltaMin_,p_) +
				p_ * pow(alphasMax-alphaDeltaMin_-alphaMinFriction,eta_)/pow(alphaDeltaMin_,p_-1) ) *
				(alphas - alphasMax + alphaDeltaMin_);
				    
    return
        neg0( alphas - (alphasMax - alphaDeltaMin_) )*JohnsonJacksonP
       // + pos( alphas - (alphasMax - alphaDeltaMin_) )*JohnsonJacksonP;
       + pos( alphas - (alphasMax - alphaDeltaMin_) )*LinearP;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
Srivastava::frictionalPressurePrime
(
    const phaseModel& phase,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax
) const
{

    const volScalarField alphas = 1.0 - continuousPhase;

    const volScalarField JohnsonJacksonPprime = Fr_ *
            ( eta_*pow(max(alphas - alphaMinFriction, scalar(0)), eta_ - 1)
	       *max(alphasMax - alphas, alphaDeltaMin_)
	      + p_*pow(max(alphas - alphaMinFriction, scalar(0)), eta_)
	      *pos(alphasMax - alphas - alphaDeltaMin_)
	    )/pow(max(alphasMax - alphas, alphaDeltaMin_), p_ + 1);

    const volScalarField LinearPprime =
            Fr_*(
				pow(alphasMax-alphaDeltaMin_-alphaMinFriction,eta_)/pow(alphaDeltaMin_,p_) +
				eta_ * pow(alphasMax-alphaDeltaMin_-alphaMinFriction,eta_-1)/pow(alphaDeltaMin_,p_) +
				p_ * pow(alphasMax-alphaDeltaMin_-alphaMinFriction,eta_)/pow(alphaDeltaMin_,p_-1) );

    return         
        neg0( alphas - (alphasMax - alphaDeltaMin_) )*JohnsonJacksonPprime
       // + pos( alphas - (alphasMax - alphaDeltaMin_) )*JohnsonJacksonPprime;
       + pos( alphas - (alphasMax - alphaDeltaMin_) )*LinearPprime;

}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
Srivastava::nu
(
    const phaseModel& phase,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax,
    const volScalarField& pf,
    const volScalarField& rho,
    const volSymmTensorField& D
) const
{
    const volScalarField alphas = 1.0 - continuousPhase;
    
    tmp<volScalarField> tnu
    (
        volScalarField::New
        (
            IOobject::groupName
            (
                Foam::typedName<frictionalStressModel>("nu"),
                phase.group()
            ),
            phase.mesh(),
            dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
        )
    );

    volScalarField& nuf = tnu.ref();

    forAll(D, celli)
    {
        if (alphas[celli] > alphaMinFriction.value())
        {
            nuf[celli] =
                0.5*pf[celli]/rho[celli]*sin(phi_.value())
               /(
                    sqrt((1.0/3.0)*sqr(tr(D[celli])) - invariantII(D[celli]))
                  + small
                );
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    volScalarField::Boundary& nufBf = nuf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            nufBf[patchi] =
                (
                    pf.boundaryField()[patchi]/rho.boundaryField()[patchi]
                   *sin(phi_.value())
                   /(
                        mag(phase.U()().boundaryField()[patchi].snGrad())
                      + small
                    )
                );
        }
    }

    // Correct coupled BCs
    nuf.correctBoundaryConditions();

    return tnu;
}


bool Foam::kineticTheoryModels::frictionalStressModels::
Srivastava::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    Fr_.read(coeffDict_);
    eta_.read(coeffDict_);
    p_.read(coeffDict_);

    phi_.read(coeffDict_);
    phi_ *= constant::mathematical::pi/180.0;

    alphaDeltaMin_.read(coeffDict_);

    return true;
}


// ************************************************************************* //
