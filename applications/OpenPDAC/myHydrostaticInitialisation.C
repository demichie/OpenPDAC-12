/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "myHydrostaticInitialisation.H"

#include "phaseSystem.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "constrainPressure.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::hydrostaticInitialisation
(
    volScalarField& p_rgh,
    volScalarField& ph_rgh,
    volScalarField& p,
    const uniformDimensionedVectorField& g,
    dimensionedScalar& hRef,
    const volScalarField& gh,
    const surfaceScalarField& ghf,
    phaseSystem& fluid,
    const dictionary& dict
)
{
    dimensionedScalar pRef("pRef",dimensionSet(1, -1, -2, 0, 0),scalar(101325.0));

    if (dict.lookupOrDefault<bool>("hydrostaticInitialisation", false))
    {
        const fvMesh& mesh = p_rgh.mesh();
        
        volScalarField rho("rho", fluid.rho());
        volVectorField U("U", fluid.U());

        dimensionedScalar xMin = min(mesh.Cf().component(0));
        dimensionedScalar xMax = max(mesh.Cf().component(0));
        dimensionedScalar yMin = min(mesh.Cf().component(1));
        dimensionedScalar yMax = max(mesh.Cf().component(1));
        dimensionedScalar zMin = min(mesh.Cf().component(2));
        dimensionedScalar zMax = max(mesh.Cf().component(2));
     
        if (!mesh.time().restart())
        {

            volScalarField& ph = regIOobject::store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "ph",
                        "0",
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
               	    mesh
                )
       	    );

	    label nFixValPatches = 0;
	    label bdryPref = 0;

            forAll(ph.boundaryField(), bdryID)     
   	    {    		
	        Info << "magSfBdry" << mesh.Sf().boundaryField()[bdryID][0] << endl;
   	        //if ( ( mag(mesh.Sf().boundaryField()[bdryID][0]^g).value() == 0.0 ) 
   	        if ( ( mag(mesh.Sf().boundaryField()[bdryID][0]^g).value() <= 
   	               1.0e-8*mag(mesh.Sf().boundaryField()[bdryID][0])*mag(g).value() ) 
   			    	&& ( ph.boundaryField()[bdryID].type() == "fixedValue" ) )
       	        {
   	            nFixValPatches +=1;
   	            const fvPatchScalarField& ph_p = ph.boundaryField()[bdryID];
    			    
		    if ( min(ph_p) == max(ph_p) )
		    {
		        pRef.value() = min(ph_p);
		        bdryPref = bdryID;
		        Info << "pRef " << pRef.value() << endl;
		        if ( g.component(0).value() != 0.0 )
		        {
		            hRef = xMax;
		        }
		        else if ( g.component(1).value() != 0.0 )
		        {
		            hRef = yMax;
		        }
		        else
		        {
		            hRef = zMax;
		        }        
		        
		             
		    }					

   	        }
	    }

            ph = pRef;
            
            for (label i=0; i<5; i++)
            {
                p = ph;
                fluid.correctThermo();
                rho = fluid.rho();

                ph = pRef- ( ( g & mesh.C() ) - (-mag(g)*hRef) ) * 0.5* ( min(rho)+ max(rho));
                Info<< "min ph " << min(ph).value() << 
                       " max ph " << max(ph).value() << endl;
            }
            
            hRef.value() = 0.0;
            Info << "hRef " << hRef.value() << endl;
                    
            // the new hydrostatic pressure profile is used to update the density field 		
            p = ph;
            fluid.correctThermo();
            rho = fluid.rho();         
            
            // we initialize the field ph_rgh with the computed pressure and density
            ph_rgh = ph - rho*gh;
       
            // Find the label of the patch that shall be changed
            label patchID = mesh.boundaryMesh().findIndex("top");
        
            // we change the fixed value b.c. of ph_rgh at the top face, in order to be 
            // consistent with the values of ph, rho and gh
            forAll(ph_rgh.boundaryField()[patchID], faceI)
            {
                ph_rgh.boundaryFieldRef()[patchID][faceI] = ph.boundaryField()[bdryPref][faceI] 
                    - rho.boundaryField()[bdryPref][faceI] * gh.boundaryField()[bdryPref][faceI];
            }           
            
            surfaceScalarField rhof("rhof", fvc::interpolate(rho));
            surfaceScalarField phig
            (
                "phig",-rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
            );

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(ph_rgh, rho, U, phig, rhof);
            p = ph_rgh + rho*gh;
            
            fluid.correctThermo();
            rho = fluid.rho();

            const fvPatchScalarField& ph_rgh_top = ph_rgh.boundaryField()[patchID];

            Info<< "min ph_rgh top " << min(ph_rgh_top) <<
   	           "max ph_rgh top " << max(ph_rgh_top) << endl;
            
            
            p = ph_rgh + rho*gh;
            fluid.correctThermo();
            rho = fluid.rho();

            label nCorr
            (
                dict.lookupOrDefault<label>("nHydrostaticCorrectors", 5)
            );

            for (label i=0; i<nCorr; i++)
            {
                surfaceScalarField rhof("rhof", fvc::interpolate(rho));

                surfaceScalarField phig
                (
                    "phig",
                    -rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
                );

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(ph_rgh, rho, U, phig, rhof);

                fvScalarMatrix ph_rghEqn
                (
                    fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
                );

                ph_rghEqn.solve();

                p = ph_rgh + rho*gh;
                fluid.correctThermo();
                rho = fluid.rho();

                Info<< "Hydrostatic pressure variation "
                    << (max(ph_rgh) - min(ph_rgh)).value() << endl;

                Info<< "min p " << min(p).value() <<
                       "max p " << max(p).value() << endl;
                Info<< "min rho " << min(rho).value() <<
                       "max rho " << max(rho).value() << endl;


            }

            ph_rgh.write();
            p.write();

            p_rgh = ph_rgh;
        }
    }
    else
    {
        volScalarField rho("rho", fluid.rho());
        // Force p_rgh to be consistent with p
        p_rgh = p - rho*gh;
    }
}


// ************************************************************************* //
