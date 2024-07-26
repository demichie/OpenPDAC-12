/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::PressureGradientPForce

Description
    Calculates particle pressure gradient force

SourceFiles
    PressureGradientPForceI.H
    PressureGradientPForce.C

\*---------------------------------------------------------------------------*/

#ifndef PressureGradientPForce_H
#define PressureGradientPForce_H

#include "ParticleForce.H"
#include "volFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PressureGradientPForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PressureGradientPForce
:
    public ParticleForce<CloudType>
{
protected:

    // Protected data

        //- Name of velocity field
        const word PName_;

        //- Rate of change of carrier phase velocity interpolator
        autoPtr<interpolation<vector>> gradPInterpPtr_;


public:

    //- Runtime type information
    TypeName("PressureGradientP");


    // Constructors

        //- Construct from mesh
        PressureGradientPForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        PressureGradientPForce(const PressureGradientPForce& pgf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new PressureGradientPForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PressureGradientPForce();


    // Member Functions

        // Access

            //- Return the rate of change of carrier phase velocity interpolator
            inline const interpolation<vector>& gradPInterp() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            //- Return the added mass
            virtual scalar massAdd
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar mass
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PressureGradientPForceI.H"

#ifdef NoRepository
    #include "PressureGradientPForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
