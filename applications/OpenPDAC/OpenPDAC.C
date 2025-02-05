/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "OpenPDAC.H"
#include "localEulerDdtScheme.H"
#include "surfaceFields.H"
#include "fvcDiv.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"
#include "myHydrostaticInitialisation.H"

#include "IOobjectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(OpenPDAC, 0);
    addToRunTimeSelectionTable(solver, OpenPDAC, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solvers::OpenPDAC::read()
{
    fluidSolver::read();

    predictMomentum =
        pimple.dict().lookupOrDefault<bool>("momentumPredictor", false);

    faceMomentum =
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false);

    dragCorrection =
        pimple.dict().lookupOrDefault<Switch>("dragCorrection", false);

    nEnergyCorrectors =
        pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1);
        
    lowPressureTimestepCorrection =     
        pimple.dict().lookupOrDefault<Switch>("lowPressureTimestepCorrection", false);        

    correctTdispersed =     
        pimple.dict().lookupOrDefault<Switch>("correctTdispersed", false);        

    nonOrthogonalResidual =     
        pimple.dict().lookupOrDefault<scalar>("nonOrthogonalResidual", 0.0);        

    innerResidual =     
        pimple.dict().lookupOrDefault<scalar>("innerResidual", 0.0);        

    if (pimple.dict().found("energyControl"))
            {
                energyControlDict = pimple.dict().subDict("energyControl");
            }


    return true;
}


void Foam::solvers::OpenPDAC::correctCoNum()
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    forAll(movingPhases, movingPhasei)
    {
        sumPhi = max
        (
            sumPhi,
            fvc::surfaceSum(mag(movingPhases[movingPhasei].phi()))()
           .primitiveField()
        );
    }

    if (lowPressureTimestepCorrection)
    {
        volScalarField alphasMax = fluid_.alfasMax();
        const word&continuousPhaseName = fluid.continuousPhaseName();
        volScalarField alfaCont = fluid.phases()[continuousPhaseName];
    
        scalarField alfa_ratio = pow(max(0*alphasMax,alphasMax-alfaCont)/alphasMax,0.5);

        Info<< "p_ratio = " << p_ratio << endl;
        Info<< "alfa_ratio: min = " << min(alfa_ratio) << endl;
        sumPhi /= sqrt(p_ratio);
    }
    

    CoNum_ = 0.5*gMax(sumPhi/mesh.V().primitiveField())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5
       *(gSum(sumPhi)/gSum(mesh.V().primitiveField()))
       *runTime.deltaTValue();

    if (lowPressureTimestepCorrection)
    {
        Info<< "Courant Number mean: " << meanCoNum*sqrt(p_ratio)
            << " max: " << CoNum*sqrt(p_ratio) << endl;
        Info<< "Modified Courant Number mean: " << meanCoNum
            << " max: " << CoNum << endl;
    }
    else
    {
        Info<< "Courant Number mean: " << meanCoNum*p_ratio
            << " max: " << CoNum << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::OpenPDAC::OpenPDAC(fvMesh& mesh)
:
    fluidSolver(mesh),

    predictMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("momentumPredictor", false)
    ),

    faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    ),

    dragCorrection
    (
        pimple.dict().lookupOrDefault<Switch>("dragCorrection", false)
    ),

    nEnergyCorrectors
    (
        pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1)
    ),

    lowPressureTimestepCorrection
    (
        pimple.dict().lookupOrDefault<Switch>("lowPressureTimestepCorrection", false)
    ),
    
    trDeltaT
    (
        LTS
      ? new volScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTName,
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
      : nullptr
    ),

    trDeltaTf
    (
        LTS && faceMomentum
      ? new surfaceScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTfName,
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, 1)
        )
      : nullptr
    ),

    buoyancy(mesh),

    fluidPtr_(phaseSystem::New(mesh)),

    fluid_(fluidPtr_()),

    phases_(fluid_.phases()),

    movingPhases_(fluid_.movingPhases()),

    phi_(fluid_.phi()),

    p_(movingPhases_[0].fluidThermo().p()),

    p_rgh(buoyancy.p_rgh),

    rho
    (
        IOobject
        (
            "rho",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid_.rho()
    ),
    
    carrierIdx(0),

    muC(phases_[carrierIdx].fluidThermo().mu()),

    muMix
    (
        IOobject
        (
            "muMix",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
        
    U
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    // Initialize cloud
    clouds(rho, U, muMix, buoyancy.g),
    
    pressureReference
    (
        p_,
        p_rgh,
        pimple.dict(),
        fluid_.incompressible()
    ),

    MRF(fluid_.MRF()),

    fluid(fluid_),
    phases(phases_),
    movingPhases(movingPhases_),
    p(p_),
    phi(phi_)
{
    // Read the controls
    read();

    mesh.schemes().setFluxRequired(p_rgh.name());

    // create ph_rgh (p_rgh for hydrostatic pressure)
    volScalarField& ph_rgh = regIOobject::store
    (
        new volScalarField
        (
            IOobject
            (
                "ph_rgh",
                "0",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );
    
    // Initialization of hydrostatic pressure profile
    hydrostaticInitialisation
    (
        p_rgh,
        ph_rgh,
        p_,
        buoyancy.g,
        buoyancy.hRef,
        buoyancy.gh,
        buoyancy.ghf,
        fluid_,
        pimple.dict()
    );
    

    // Correct mixture thermodynamics with new pressure    
    fluid_.correctThermo();
    rho = fluid_.rho();

    Info << "hRef " << buoyancy.hRef.value() << endl;

    Info<< "min p " << min(p_).value() <<
  	               " max p " << max(p_).value() << endl;

    p_ratio = min(p_).value() /p_.weightedAverage(mesh_.V()).value();
    	                 	                 	               
    Info<< "min p_rgh " << min(p_rgh).value() <<
   	               " max p_rgh " << max(p_rgh).value() << endl;
    Info<< "min rho " << min(rho).value() <<
   	               " max rho " << max(rho).value() << endl;

    // Carrier phase viscosity
    const word&continuousPhaseName = fluid.continuousPhaseName();
    muC = phases_[continuousPhaseName].fluidThermo().mu();
    
    Info<< "min muC " << min(muC).value() << " max muC " << max(muC).value() << endl;

    volScalarField alphasMax = fluid_.alfasMax();
    Info<< "min alphasMax " << min(alphasMax).value() << " max alphasMax " << max(alphasMax).value() << endl;
   
    // Mixture viscosity
    muMix = muC * pow( 1.0 - ( 1.0 - max(0.0,phases[continuousPhaseName]) ) / alphasMax , -1.55);
    Info<< "min muMix " << min(muMix).value() << " max muMix " << max(muMix).value() << endl;
   
    // Compute mass-weighted mixture velocity
    U = 0.0* phases_[0].U();
    forAll(phases_, phasei)
    {
        phaseModel& phase = phases_[phasei];
        U += phase * phase.rho() * phase.U() / rho;

    }

    clouds.info();
    
    if (pimple.dict().lookupOrDefault<bool>("hydrostaticInitialisation", false))
    {
 
        const Time& runTime = mesh().time();
        scalar startTime_ = runTime.startTime().value();
        scalar deltaT = runTime.deltaT().value();

        // set small value for deltaT to evolve particles    
        const_cast<Time&>(runTime).setDeltaT(1.e-5*deltaT);

        // increase time iterator            
        const_cast<Time&>(runTime)++;
    
        // evolve particle cloud
        clouds.evolve();

        // restore startTime
        const_cast<Time&>(runTime).setTime(startTime_,startTime_);
                
        // restore deltaT            
        const_cast<Time&>(runTime).setDeltaT(deltaT);

        // write everything (including lagrangian)
        const_cast<Time&>(runTime).writeNow();
    }  
    

    if (transient())
    {
        correctCoNum();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::OpenPDAC::~OpenPDAC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::OpenPDAC::preSolve()
{
    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        setRDeltaT();
    }

    // Store divU from the previous mesh so that it can be
    // mapped and used in correctPhi to ensure the corrected phi
    // has the same divergence
    if (correctPhi || mesh.topoChanging())
    {
        // Construct and register divU for mapping
        divU = new volScalarField
        (
            "divU0",
            fvc::div(fvc::absolute(phi, movingPhases[0].U()))
        );
    }

    fvModels().preUpdateMesh();

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::OpenPDAC::prePredictor()
{
    if (pimple.thermophysics() || pimple.flow())
    {
        fluid_.solve(rAs);
        fluid_.correct();
        fluid_.correctContinuityError();
    }

    if (pimple.flow() && pimple.predictTransport())
    {
        fluid_.predictMomentumTransport();
    }
}


void Foam::solvers::OpenPDAC::postCorrector()
{
    if (pimple.flow() && pimple.correctTransport())
    {
        fluid_.correctMomentumTransport();
        fluid_.correctThermophysicalTransport();
    }
}


void Foam::solvers::OpenPDAC::postSolve()
{
    divU.clear();
    
    volScalarField alphasMax = fluid_.alfasMax();
    
    const word&continuousPhaseName = fluid.continuousPhaseName();
    
    muMix = muC * pow( max(0.0, 1.0 - ( 1.0 - max(0.0,phases[continuousPhaseName]) )) / alphasMax , -1.55);

    rho = fluid_.rho();
    
    U *= 0.0;
    forAll(phases_, phasei)
    {
        phaseModel& phase = phases_[phasei];
        U += phase * phase.rho() * phase.U() / rho;

    }

    Info<< "min mu " << min(muMix).value() << " max mu " << max(muMix).value() << endl;

    clouds.evolve();
            
}


// ************************************************************************* //
