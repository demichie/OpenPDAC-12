/*---------------------------------------------------------------------------*\
License
    This file is part of OpenPDAC.

    OpenPDAC is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenPDAC is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenPDAC.  If not, see <http://www.gnu.org/licenses/>.

Application
    topoGrid

Description
    Deforms a polyMesh using an ESRI raster ascii file.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "IStringStream.H"
#include "volPointInterpolation.H"
#include "UniformTable2.H"
#include "RectangularMatrix.H"
#include <fstream>
#include <sstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshNoChangers.H"

    // Read the dictionary file (topoGridDict) from the "system" folder
    IOdictionary topoDict
    (
        IOobject
        (
            "topoGridDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read the raster file name from the dictionary
    const word rasterFile = topoDict.lookup<word>("rasterFile");

    // Read the vent center coordinates from the dictionary
    const scalar xVent = topoDict.lookupOrDefault<scalar>("xVent",0.0);
    const scalar yVent = topoDict.lookupOrDefault<scalar>("yVent",0.0);
    const scalar expFactor = topoDict.lookupOrDefault<scalar>("expFactor",1.0);
    const scalar dzVert = topoDict.lookupOrDefault<scalar>("dzVert",0.0);
    const scalar exp_shape = topoDict.lookupOrDefault<scalar>("exp_shape",1.0);

    // Output the file name to the terminal for verification
    Info << "Raster file specified: " << rasterFile << endl;

    // Read the ESRI ASCII Raster file
    std::ifstream file(rasterFile);
    
    if (!file.is_open()) 
    {
        FatalErrorInFunction
            << "Unable to open the raster file: " << rasterFile << exit(FatalError);
    }


    int ncols = 0, nrows = 0;
    double xllcorner = 0.0, yllcorner = 0.0, cellsize = 0.0;
    double NODATA_value = -9999.0;
    std::string line;

    // Read the header
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "ncols")
            iss >> ncols;
        else if (key == "nrows")
            iss >> nrows;
        else if (key == "xllcorner" || key == "xllcenter")
            iss >> xllcorner;
        else if (key == "yllcorner" || key == "yllcenter")
            iss >> yllcorner;
        else if (key == "cellsize")
            iss >> cellsize;
        else if (key == "NODATA_value")
            iss >> NODATA_value;

        if (key == "NODATA_value")
            break;
    }

    xllcorner -= xVent;
    yllcorner -= yVent;

    // Create a RectangularMatrix to store the elevation data
    RectangularMatrix<double> elevation(nrows, ncols, 0.0);

    // Read the elevation data and store it in the RectangularMatrix
    for (int i = 0; i < nrows; ++i)
    {
        std::getline(file, line);
        std::istringstream iss(line);

        for (int j = 0; j < ncols; ++j)
        {
            double value;
            iss >> value;

            if (value == NODATA_value)
                value = 0.0;  // Handle NODATA_value appropriately

            elevation(i, j) = value;
        }
    }

    double maxTopo(max(elevation));
    // double minTopo(max(elevation));

    scalar zVert(maxTopo + dzVert);

    file.close();

    // Get times list
    instantList Times = runTime.times();

    // skip "constant" time
    for (label timeI = 1; timeI < Times.size(); ++timeI)
    {
        runTime.setTime(Times[timeI], timeI);

        Info<< "Time = " << runTime.userTimeName() << endl;

    scalar zMin = min(mesh.Cf().component(2)).value();
    scalar zMax = max(mesh.Cf().component(2)).value();
    
    Info << "zMin = " << zMin << endl;
    Info << "zMax = " << zMax << endl;
    
    scalar z2Rel(0.0);
    scalar zNew(0.0);

    pointField zeroPoints(mesh.points());

    pointField pDeform
    (
        0.0*zeroPoints
    );


    // Loop over all cells in the mesh to interpolate elevation values
    forAll(pDeform,pointi)
    {
        // Get x, y coordinates of the pointi
        scalar x = mesh.points()[pointi].x();
        scalar y = mesh.points()[pointi].y();
        scalar z = mesh.points()[pointi].z();
        
        scalar zRel = min(1.0, (zMax-z)/(zMax-zMin));

        // Calculate row and column indices in the elevation matrix
        int colIndex = (x - xllcorner) / cellsize;
        int rowIndex = (y - yllcorner) / cellsize;

        // Interpolate elevation value
        if (colIndex >= 0 && colIndex <= ncols  && rowIndex >= 0 && rowIndex <= nrows )
        {
            // Bilinear interpolation
            scalar xLerp = (x - (xllcorner + colIndex * cellsize)) / cellsize;
            scalar yLerp = (y - (yllcorner + rowIndex * cellsize)) / cellsize;

            scalar v00 = elevation(rowIndex, colIndex);
            scalar v01 = elevation(rowIndex, colIndex + 1);
            scalar v10 = elevation(rowIndex + 1, colIndex);
            scalar v11 = elevation(rowIndex + 1, colIndex + 1);

            scalar zInterp = 
                v00 * (1 - xLerp) * (1 - yLerp) +
                v01 * xLerp * (1 - yLerp) +
                v10 * (1 - xLerp) * yLerp +
                v11 * xLerp * yLerp;

            // Assign interpolated value to the volScalarField U
            pDeform[pointi].z() = zRel * zInterp;
                
            zNew = z + zRel * zInterp;
                
            if ( z>= 0.0)
            {
                if (dzVert > 0)
                {
                    // enlarge from a fixed height above the maximum
                    // topography and the top, thus from an horizontal
                    // plane to the top
                    z2Rel = max(0, (zNew - zVert) / (zMax - zVert));
                }
                else
                {
                    // enlarge from the topography to the top
                    z2Rel = (zNew - zInterp) / (zMax - zInterp);
                }
                z2Rel = std::pow(z2Rel,exp_shape);
                
                pDeform[pointi].x() = z2Rel*(expFactor-1.0)*x;
                pDeform[pointi].y() = z2Rel*(expFactor-1.0)*y;
            }
            else
            {
                pDeform[pointi].x() = 0.0;
                pDeform[pointi].y() = 0.0;
            }        
        }
        else
        {
            // If outside the raster bounds, set to a default value (e.g., 0)
            pDeform[pointi].z() = 0.0;
            pDeform[pointi].x() = 0.0;
            pDeform[pointi].y() = 0.0;
        }
    }
           
    
    pointField newPoints
    (
        zeroPoints + pDeform
    );

    mesh.setPoints(newPoints);
    mesh.write();

    Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
