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

#include "SinclairJackson.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(SinclairJackson, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        SinclairJackson,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::SinclairJackson::SinclairJackson
(
    const dictionary& dict
)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::SinclairJackson::~SinclairJackson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::SinclairJackson::g0
(
    const volScalarField& alpha,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax
) const
{
    return 1.0/(1 - cbrt(min(alpha, alphaMinFriction)/alphasMax));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::SinclairJackson::g0prime
(
    const volScalarField& alpha,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax
) const
{
    volScalarField aByaMax
    (
        cbrt(min(max(alpha, scalar(1e-3)), alphaMinFriction)/alphasMax)
    );

    // TODO: CHECK IF THIS MAKE THE CONVERGENCE WORST OR BETTER
    volScalarField posCoeff
    (
        pos(alphaMinFriction-alpha) * pos(alpha-scalar(1e-3))
    );
    return posCoeff*(1.0/(3*alphasMax))/sqr(aByaMax - sqr(aByaMax));    
    // return (1.0/(3*alphasMax))/sqr(aByaMax - sqr(aByaMax));
}


//Foam::tmp<Foam::volScalarField>
Foam::PtrList<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::SinclairJackson::g0
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

    PtrList<volScalarField> g0_mm(fluid.phases().size());
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
    
    volScalarField alphas = 1.0 - continuousPhase;

    volScalarField g0 = 1.0/(1 - cbrt(min(alphas, alphaMinFriction)/alphasMax)); 

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];
        
        if (&phase != &continuousPhase)
        {
            g0_mm.set
            (
            	phaseIdx,
            	volScalarField
            	(
            	    "g0_mm" + phasei.name() + "_" + phase.name(),
            	    g0 + 0.5*phase.d()*const_sum
            	)
            ); 
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
            	    ( phase.d()*g0_mm[indexi] + phasei.d()*g0_mm[iter] ) 
            	    / ( phasei.d() + phase.d() )
            	)
            ); 
        }     
    }

    return g0_im;
}


Foam::PtrList<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::SinclairJackson::g0prime
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

    PtrList<volScalarField> g0prime_mm(fluid.phases().size());
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

    volScalarField alphas = 1.0 - continuousPhase;

    volScalarField aByaMax
    (
        cbrt(min(max(alphas, scalar(1e-3)), alphaMinFriction)/alphasMax)
    );
    
    volScalarField g0prime = (1.0/(3*alphasMax))/sqr(aByaMax - sqr(aByaMax));

    // TODO: CHECK IF THIS MAKE THE CONVERGENCE WORST OR BETTER
	volScalarField posCoeff
	(
	    pos(alphaMinFriction-alphas) * pos(alphas-scalar(1e-3))
	);	
    g0prime *= posCoeff;
    // CORRECTION ENDS HERE

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];
        
        if (&phase != &continuousPhase)
        {
            g0prime_mm.set
            (
            	phaseIdx,
            	volScalarField
            	(
            	    "g0prime_mm" + phasei.name() + "_" + phase.name(),
            	    g0prime + 0.5 * phasei.d() / phase.d()
            	)
            );            
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
            	    "g0prime_im" + phasei.name() + "_" + phase.name(),
            	    ( phasei.d()*g0prime_mm[indexi] + phase.d()*g0prime_mm[iter] ) 
            	    / ( phasei.d() + phase.d() ) 
            	)
            );

	// Info << "min g0 " << min(g0prime_im[iter]).value() << endl;
        // Info << "max g0 " << max(g0prime_im[iter]).value() << endl;
	            
        }    
    
    }

    return g0prime_im;
}


// ************************************************************************* //
