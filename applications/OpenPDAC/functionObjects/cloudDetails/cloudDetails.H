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

Class
    Foam::functionObjects::cloudDetails

Description
    Outputs Lagrangian cloud information to a file.

    The current outputs include:
    - total current number of parcels
    - total current mass of parcels

    Example of function object specification:
    \verbatim
    cloudDetails1
    {
        type        cloudDetails;
        libs        ("liblagrangianFunctionObjects.so");
        ...
        clouds
        (
            cloud1
            thermoCloud1
        );
    }
    \endverbatim


Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: cloudDetails    | yes         |
        clouds       | list of clouds names to process |yes  |
    \endtable

    The output data of each cloud is written to a file named \<cloudName\>.dat

See also
    Foam::functionObject
    Foam::functionObjects::regionFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    cloudDetails.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudDetails_H
#define functionObjects_cloudDetails_H

#include "regionFunctionObject.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class cloudDetails Declaration
\*---------------------------------------------------------------------------*/

class cloudDetails
:
    public regionFunctionObject,
    public logFiles
{
protected:

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        bool gatherScalarData(const Foam::List<Foam::scalar>& localData, Foam::List<Foam::scalar>& gatheredData);

public:

    //- Runtime type information
    TypeName("cloudDetails");


    // Constructors

        //- Construct from Time and dictionary
        cloudDetails
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        cloudDetails(const cloudDetails&) = delete;


    //- Destructor
    virtual ~cloudDetails();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudDetails&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
