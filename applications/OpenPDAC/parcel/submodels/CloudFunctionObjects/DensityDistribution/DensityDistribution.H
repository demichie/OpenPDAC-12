/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::DensityDistribution

Description
    Creates graphs of a cloud's Density distribution

    Example usage:
    \verbatim
    massFlux1
    {
        type        DensityDistribution;
        nPoints     40;
        setFormat   raw;
    }
    \endverbatim

SourceFiles
    DensityDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef DensityDistribution_H
#define DensityDistribution_H

#include "CloudFunctionObject.H"
#include "setWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class DensityDistribution Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DensityDistribution
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        //- Number of points to plot
        const label nPoints_;

        //- Set formatter
        autoPtr<setWriter> formatter_;


protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("DensityDistribution");


    // Constructors

        //- Construct from dictionary
        DensityDistribution
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        DensityDistribution(const DensityDistribution<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new DensityDistribution<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DensityDistribution();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DensityDistribution.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
