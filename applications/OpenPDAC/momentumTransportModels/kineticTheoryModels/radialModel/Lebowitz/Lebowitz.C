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

#include "Lebowitz.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(Lebowitz, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        Lebowitz,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::Lebowitz::Lebowitz
(
    const dictionary& dict
)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::Lebowitz::~Lebowitz()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Lebowitz::g0
(
    const volScalarField& alpha,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax
) const
{
    return
        1.0/continuousPhase
      + 3*alpha/(2*sqr(continuousPhase));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Lebowitz::g0prime
(
    const volScalarField& alpha,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax
) const
{
    return
        (alpha+5) / (2.0*pow3(continuousPhase));
}


Foam::PtrList<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Lebowitz::g0
(
    const phaseModel& phasei,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax
) const
{
    const volScalarField& alphai = phasei;
    const label& indexi = phasei.index();
    const phaseSystem& fluid = phasei.fluid();

    PtrList<volScalarField> g0_im(fluid.phases().size());
    
    volScalarField const_sum = alphai/phasei.d();

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];
        
        if ((&phase != &continuousPhase) and !(phaseIdx==indexi))
        {
    	    const volScalarField& alpha = phase;
            const_sum += alpha / phase.d();
        }

    } 
    
    forAll(g0_im, iter)
    {
        const phaseModel& phase = fluid.phases()[iter];

        if (&phase != &continuousPhase)
        {
            g0_im.set
            (
            	iter,
            	volScalarField
            	(
            	    "g0_im" + phasei.name() + "_" + phase.name(),
            	    1.0/continuousPhase + 3 * phasei.d() * phase.d() 
            	    / ( sqr(continuousPhase) * phasei.d() + phase.d() ) * const_sum
            	)
            ); 
        }     
    }

    return g0_im;
}


// Foam::tmp<Foam::volScalarField>
Foam::PtrList<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Lebowitz::g0prime
(
    const phaseModel& phasei,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax
) const
{
    const volScalarField& alphai = phasei;
    const label& indexi = phasei.index();
    const phaseSystem& fluid = phasei.fluid();

    PtrList<volScalarField> g0prime_im(fluid.phases().size());
    
    volScalarField const_sum = alphai/phasei.d();

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];
        
        if ((&phase != &continuousPhase) and !(phaseIdx==indexi))
        {
    	    const volScalarField& alpha = phase;
            const_sum += alpha / phase.d();
        }

    } 
    
    forAll(g0prime_im, iter)
    {
        const phaseModel& phase = fluid.phases()[iter];

        if (&phase != &continuousPhase)
        {
            g0prime_im.set
            (
            	iter,
            	volScalarField
            	(
            	    "g0_im" + phasei.name() + "_" + phase.name(),
            	    1.0/sqr(continuousPhase) + 6 * phasei.d() * phase.d() 
            	    / ( pow(continuousPhase,3) * phasei.d() + phase.d() ) * const_sum
            	    + 3 * phasei.d() * phase.d() / phasei.d()
            	    / ( sqr(continuousPhase) * phasei.d() + phase.d() ) 
            	)
            ); 
        }     
    }
    
    return g0prime_im;
}


// ************************************************************************* //
