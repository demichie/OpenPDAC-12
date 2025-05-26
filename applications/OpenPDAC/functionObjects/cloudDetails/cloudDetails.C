/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

#include "cloudDetails.H"
#include "collidingCloud.H"
#include "momentumParcel.H"

#include "parcelCloud.H"
#include "dictionary.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "List.H"
#include "label.H"
#include "scalar.H"
#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudDetails, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        cloudDetails,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::functionObjects::cloudDetails::writeFileHeader(const label i)
{
    /*
    if (logFiles::file(i).good()) // Controlla se lo stream è valido
    {
        writeTabbed(logFiles::file(i), "origId,x,y,z,Ux,Uy,Uz,rho,d");
        logFiles::file(i) << endl;
    }
    else
    {
        WarningIn("Foam::functionObjects::cloudDetails::writeFileHeader")
            << "Stream for file index " << i << " (cloud: " << logFiles::names()[i]
            << ") is not good. Header not written." << endl;
    }
    */
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudDetails::cloudDetails
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    logFiles(obr_, name)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudDetails::~cloudDetails()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::cloudDetails::gatherScalarData(const Foam::List<Foam::scalar>& localData, Foam::List<Foam::scalar>& gatheredData)
{
    const Foam::label nProcs = Foam::Pstream::nProcs();

    // Lista di liste: ogni processo ha una lista locale, quindi creiamo la struttura per raccoglierle tutte
    Foam::List<Foam::List<Foam::scalar>> allData(nProcs);

    // Ogni processo mette la propria lista nella posizione corrispondente
    allData[Foam::Pstream::myProcNo()] = localData;

    // Gather: raccoglie tutte le liste in 'allData' sul master (proc 0)
    Foam::Pstream::gatherList(allData);

    // Solo sul master concateno tutte le sottoliste in un'unica lista piatta
    if (Foam::Pstream::master())
    {
        // Conta quanti elementi totali
        Foam::label totalSize = 0;
        for (Foam::label procI = 0; procI < nProcs; ++procI)
        {
            totalSize += allData[procI].size();
        }

        gatheredData.setSize(totalSize);
        
        // Copia i dati uno dopo l'altro
        Foam::label pos = 0;
        for (Foam::label procI = 0; procI < nProcs; ++procI)
        {
            forAll(allData[procI], i)
            {
                gatheredData[pos++] = allData[procI][i];                
            }
        }
    }
    else
    {
        // Negli altri processi la lista raccolta rimane vuota
        gatheredData.clear();
    }

    return true;
}

bool Foam::functionObjects::cloudDetails::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);

    Info<< type() << " " << name() << ": ";
    if (names().size())
    {
        Info<< "applying to clouds:" << nl;
        forAll(names(), i)
        {
            Info<< "    " << names()[i] << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no clouds to be processed" << nl << endl;
    }

    resetNames(dict.lookup("clouds"));

    return true;
}


bool Foam::functionObjects::cloudDetails::execute()
{
    return true;
}


bool Foam::functionObjects::cloudDetails::write()
{
    // logFiles::write();

    const fvMesh& pMesh_ = refCast<const fvMesh>(obr_);

 

    
    forAll(names(), i)
    {
        const word& cloudName = names()[i];

        const parcelCloud& cloud = obr_.lookupObject<parcelCloud>(cloudName);

        label nParcels = returnReduce(cloud.nParcels(), sumOp<label>());

        const collidingCloud* momCloud = dynamic_cast<const collidingCloud*>(&cloud);
        
        if (momCloud)
        {
        
            Foam::List<Foam::scalar> localData(cloud.nParcels()*9);

            Info << "Writing cloud with " << nParcels << " nParcels" << endl;
            label j = 0;
            forAllConstIter(collidingCloud, (*momCloud), iter)
            {
                const momentumParcel& p = iter();
                const point pos = p.position(pMesh_);
                const vector& U = p.U();
                scalar rho = p.rho();
                scalar d = p.d();
                label origId = p.origId();
                
                label base = j * 9;

                localData[base + 0] = scalar(origId);
                localData[base + 1] = pos.x();
                localData[base + 2] = pos.y();
                localData[base + 3] = pos.z();
                localData[base + 4] = U.x();
                localData[base + 5] = U.y();
                localData[base + 6] = U.z();
                localData[base + 7] = rho;
                localData[base + 8] = d;
                j++;
            }

            Foam::List<Foam::scalar> allData;
            gatherScalarData(localData, allData);
            
            if (Pstream::master())
            {
            
                const scalar currentTime = time_.value();
                Info << "current time " << currentTime << endl;
    
                fileName outputDir = "postProcessing"
                         / this->name()
                         / Foam::word(Foam::name(currentTime));
                         
                Foam::fileName outputFile = outputDir / "output.csv";                 
                         
                Info << outputFile << endl;                 
                if (!isDir(outputDir))
                {
                    Foam::mkDir(outputDir);
                }
            
                Foam::OFstream ofs(outputFile);

                // writeTime(ofs);

                ofs << "origId,x,y,z,Ux,Uy,Uz,rho,d\n" ;
        

                for (label h = 0; h < allData.size(); h++)
                {
                    if (h% 9 == 0)
                    {
                        ofs << label(allData[h]);
                    }
                    else
                    {
                        ofs << allData[h];
                    }
                    
                    if ((h + 1) % 9 == 0)
                    {
                        ofs << "\n";  // fine riga
                    }
                    else
                    {
                        ofs << ",";
                    }
                }
            }
        }
        else
        {
            FatalErrorInFunction << "Cloud is not a collidingCloud!" << exit(FatalError);
        }        
        
    }

    return true;
}


// ************************************************************************* //
