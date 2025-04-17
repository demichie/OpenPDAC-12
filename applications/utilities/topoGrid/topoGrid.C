/*---------------------------------------------------------------------------* \
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
#include "vector.H"
#include "pointFields.H"
#include "IStringStream.H"
#include "volPointInterpolation.H"
#include "UniformTable2.H"
#include "RectangularMatrix.H"
#include <fstream>
#include <sstream>
#include "IOstreams.H"
#include <cstring>
#include "Pstream.H"
#include "tetPointRef.H"
#include "OFstream.H"
#include "globalIndex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void generateCroppedDEM(
    const RectangularMatrix<double>& elevation,
    scalar xllcorner, scalar yllcorner, scalar cellsize,
    scalar xVent, scalar yVent, // Translation factors
    scalar xmin, scalar xmax, scalar ymin, scalar ymax,
    const word& outputFileName)
{
    // Adjust the domain bounds to match the DEM coordinate system
    xmin += xVent;
    xmax += xVent;
    ymin += yVent;
    ymax += yVent;

    // Compute the new cell size ensuring an integer number of rows/cols
    int ncols_new = round((xmax - xmin) / cellsize);
    int nrows_new = round((ymax - ymin) / cellsize);
    scalar cellsize_new = (xmax - xmin) / ncols_new; // Ensure exact fit

    // Compute the new lower-left corner (adjusting for cell center convention)
    scalar xllcorner_new = xmin + 0.5 * cellsize_new;
    scalar yllcorner_new = ymin + 0.5 * cellsize_new;

    Info << "Generating cropped DEM with cellsize: " << cellsize_new << " ("
         << ncols_new << " x " << nrows_new << " grid)" << endl;

    // Open output file
    std::ofstream file(outputFileName);
    if (!file)
    {
        FatalErrorInFunction << "Cannot open output file " << outputFileName << exit(FatalError);
    }

    // Write ASCII Raster Header
    file << "ncols " << ncols_new << "\n";
    file << "nrows " << nrows_new << "\n";
    file << "xllcorner " << xllcorner_new << "\n";
    file << "yllcorner " << yllcorner_new << "\n";
    file << "cellsize " << cellsize_new << "\n";
    file << "NODATA_value -9999\n";

    // Bilinear interpolation over the new grid
    for (int row = 0; row < nrows_new; ++row)
    {
        for (int col = 0; col < ncols_new; ++col)
        {
            // Compute world coordinates of the new cell center
            scalar x = xllcorner_new + ( col + 0.5) * cellsize_new;
            scalar y = yllcorner_new + (nrows_new - row - 0.5) * cellsize_new; // Top to bottom
           
            // Compute corresponding indices in the original DEM
            int i = (y - (yllcorner + 0.5 * cellsize)) / cellsize;
            int j = (x - (xllcorner + 0.5 * cellsize)) / cellsize;

            if (i >= 0 && i < elevation.m() - 1 && j >= 0 && j < elevation.n() - 1)
            {
                // Compute interpolation weights
                scalar xLerp = (x - (xllcorner + j * cellsize + 0.5 * cellsize)) / cellsize;
                scalar yLerp = (y - (yllcorner + i * cellsize + 0.5 * cellsize)) / cellsize;

                // Get the four surrounding elevation values
                scalar v00 = elevation(i, j);
                scalar v01 = elevation(i, j + 1);
                scalar v10 = elevation(i + 1, j);
                scalar v11 = elevation(i + 1, j + 1);

                // Bilinear interpolation
                scalar zInterp =
                    v00 * (1 - xLerp) * (1 - yLerp) +
                    v01 * xLerp * (1 - yLerp) +
                    v10 * (1 - xLerp) * yLerp +
                    v11 * xLerp * yLerp;

                file << zInterp << " ";
            }
            else
            {
                file << "-9999 ";
            }
        }
        file << "\n";
    }

    file.close();
    Info << "Cropped DEM saved as " << outputFileName << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// Function to compute the normal vector of a triangle formed by points 
// p1, p2, p3
vector computeNormal(const point& p1, const point& p2, const point& p3)
{
    vector v1 = p2 - p1;  // Edge vector 1
    vector v2 = p3 - p1;  // Edge vector 2

    // Compute the cross product of v1 and v2 to get the normal
    vector normal = Foam::vector
    (
        v1.y() * v2.z() - v1.z() * v2.y(),  // x-component
        v1.z() * v2.x() - v1.x() * v2.z(),  // y-component
        v1.x() * v2.y() - v1.y() * v2.x()   // z-component
    );

    // Normalize the normal vector
    scalar magnitude = mag(normal);
    if (magnitude > SMALL)
    {
        normal /= magnitude;
    }

    return normal;
}

// Function to write a single triangle in binary format
void writeBinaryTriangle(std::ofstream& stlFile, const vector& normal, 
                         const point& p1, const point& p2, const point& p3)
{
    // Write normal vector (12 bytes: 3 floats)
    float nx(normal.x());
    float ny(normal.y());
    float nz(normal.z());
    stlFile.write(reinterpret_cast<const char*>(&nx), 4);
    stlFile.write(reinterpret_cast<const char*>(&ny), 4);
    stlFile.write(reinterpret_cast<const char*>(&nz), 4);

    // Write vertex 1 (12 bytes: 3 floats)
    float P1x(p1.x());
    float P1y(p1.y());
    float P1z(p1.z());
    stlFile.write(reinterpret_cast<const char*>(&P1x), 4);
    stlFile.write(reinterpret_cast<const char*>(&P1y), 4);
    stlFile.write(reinterpret_cast<const char*>(&P1z), 4);

    // Write vertex 2 (12 bytes: 3 floats)
    float P2x(p2.x());
    float P2y(p2.y());
    float P2z(p2.z());
    stlFile.write(reinterpret_cast<const char*>(&P2x), 4);
    stlFile.write(reinterpret_cast<const char*>(&P2y), 4);
    stlFile.write(reinterpret_cast<const char*>(&P2z), 4);

    // Write vertex 3 (12 bytes: 3 floats)
    float P3x(p3.x());
    float P3y(p3.y());
    float P3z(p3.z());
    stlFile.write(reinterpret_cast<const char*>(&P3x), 4);
    stlFile.write(reinterpret_cast<const char*>(&P3y), 4);
    stlFile.write(reinterpret_cast<const char*>(&P3z), 4);

    // Write attribute byte count (2 bytes, set to 0)
    char attribute[2] = "0";
    stlFile.write(attribute,2);
}


// Function to write STL surface in binary format
void writeBinarySTL(const word& stlFileName, 
                    const RectangularMatrix<scalar>& elevation, 
                    scalar xOffset, 
                    scalar yOffset, 
                    scalar cellSize)
{
    std::ofstream stlFile(stlFileName.c_str(), std::ios::binary);
    if (!stlFile)
    {
        FatalErrorInFunction << "Cannot open STL file " 
                             << stlFileName << " for writing" 
                             << exit(FatalError);
    }

    // Write 80-byte header (just fill with 0 or any message)
    char header[80] = "Generated by OpenFOAM";
    stlFile.write(header, sizeof(header));

    // Get dimensions of the elevation grid
    const label numRows = elevation.m();
    const label numCols = elevation.n();
    
    // Write triangle count (4 bytes)
    auto triangleCount = static_cast<uint32_t>(2 * (numRows - 1) * (numCols - 1));
    stlFile.write(reinterpret_cast<const char*>(&triangleCount), sizeof(triangleCount));  
    
    // Loop over each cell in the grid and create two triangles per cell
    for (label i = 0; i < numRows - 1; ++i)
    {
        for (label j = 0; j < numCols - 1; ++j)
        {
            // Get corner points of the cell

            // Top-left corner
            point p1(xOffset + j * cellSize, yOffset + i * cellSize, elevation(i, j));
            
            // Top-right corner
            point p2(xOffset + (j + 1) * cellSize, yOffset + i * cellSize, elevation(i, j + 1));

            // Bottom-left corner
            point p3(xOffset + j * cellSize, yOffset + (i + 1) * cellSize, elevation(i + 1, j));
            
            // Bottom-right corner
            point p4(xOffset + (j + 1) * cellSize, yOffset + (i + 1) * cellSize, elevation(i + 1, j + 1));

            // First triangle (p1, p2, p3)
            vector normal1 = computeNormal(p1, p2, p3);

            // Info << p1 << p2 << p3 << normal1 << endl;
            writeBinaryTriangle(stlFile, normal1, p1, p2, p3);

            // Second triangle (p2, p4, p3)
            vector normal2 = computeNormal(p2, p4, p3);
            writeBinaryTriangle(stlFile, normal2, p2, p4, p3);
        }
    }

    stlFile.close();
    Info << "Binary STL surface written to " << stlFileName << endl;
}

// Function to write STL surface from the elevation grid
void writeSTL(const word& stlFileName, 
              const RectangularMatrix<scalar>& elevation, 
              scalar xOffset, 
              scalar yOffset, 
              scalar cellSize)
{
    std::ofstream stlFile(stlFileName.c_str());
    if (!stlFile)
    {
        FatalErrorInFunction << "Cannot open STL file " << stlFileName << " for writing" << exit(FatalError);
    }

    // Write STL file header
    stlFile << "solid topoSurface" << "\n";

    // Get dimensions of the elevation grid
    const label numRows = elevation.m();
    const label numCols = elevation.n();

    // Loop over each cell in the grid and create two triangles per cell
    for (label i = 0; i < numRows - 1; ++i)
    {
        for (label j = 0; j < numCols - 1; ++j)
        {
            // Get corner points of the cell (elevation grid)

            // Top-left corner
            vector p1(xOffset + j * cellSize, yOffset + i * cellSize, elevation(i, j));

            // Top-right corner
            vector p2(xOffset + (j + 1) * cellSize, yOffset + i * cellSize, elevation(i, j + 1));

            // Bottom-left corner
            vector p3(xOffset + j * cellSize, yOffset + (i + 1) * cellSize, elevation(i + 1, j));

            // Bottom-right corner
            vector p4(xOffset + (j + 1) * cellSize, yOffset + (i + 1) * cellSize, elevation(i + 1, j + 1));
            // First triangle (p1, p2, p3) - Top-left, top-right, bottom-left
            vector normal1 = computeNormal(p1, p2, p3);
            stlFile << "  facet normal " << normal1.x() << " " << normal1.y() << " " << normal1.z() << "\n";
            stlFile << "    outer loop" << "\n";
            stlFile << "      vertex " << p1.x() << " " << p1.y() << " " << p1.z() << "\n";
            stlFile << "      vertex " << p2.x() << " " << p2.y() << " " << p2.z() << "\n";
            stlFile << "      vertex " << p3.x() << " " << p3.y() << " " << p3.z() << "\n";
            stlFile << "    endloop" << "\n";
            stlFile << "  endfacet" << "\n";

            // Second triangle (p2, p4, p3) - Top-right, bottom-right, bottom-left
            vector normal2 = computeNormal(p2, p4, p3);
            stlFile << "  facet normal " << normal2.x() << " " << normal2.y() << " " << normal2.z() << "\n";
            stlFile << "    outer loop" << "\n";
            stlFile << "      vertex " << p2.x() << " " << p2.y() << " " << p2.z() << "\n";
            stlFile << "      vertex " << p4.x() << " " << p4.y() << " " << p4.z() << "\n";
            stlFile << "      vertex " << p3.x() << " " << p3.y() << " " << p3.z() << "\n";
            stlFile << "    endloop" << "\n";
            stlFile << "  endfacet" << "\n";
        }
    }

    // Write STL file footer
    stlFile << "endsolid topoSurface" << "\n";

    stlFile.close();
    Info << "STL surface written to " << stlFileName << endl;
}



scalar minQuality
(
    const polyMesh& mesh,
    const point& cC,
    const label fI,
    const bool isOwner,
    const label faceBasePtI
)
{
    // Does fan decomposition of face (starting at faceBasePti) and determines
    // min quality over all resulting tets.

    const pointField& pPts = mesh.points();
    const face& f = mesh.faces()[fI];
    const point& tetBasePt = pPts[f[faceBasePtI]];

    scalar thisBaseMinTetQuality = vGreat;

    for (label tetPtI = 1; tetPtI < f.size() - 1; tetPtI++)
    {
        label facePtI = (tetPtI + faceBasePtI) % f.size();
        label otherFacePtI = f.fcIndex(facePtI);

        label ptAI = -1;
        label ptBI = -1;

        if (isOwner)
        {
            ptAI = f[facePtI];
            ptBI = f[otherFacePtI];
        }
        else
        {
            ptAI = f[otherFacePtI];
            ptBI = f[facePtI];
        }

        const point& pA = pPts[ptAI];
        const point& pB = pPts[ptBI];

        tetPointRef tet(cC, tetBasePt, pA, pB);

        scalar tetQuality = tet.quality();

        if (tetQuality < thisBaseMinTetQuality)
        {
            thisBaseMinTetQuality = tetQuality;
        }
    }
    return thisBaseMinTetQuality;
}


// Function to compute the face flatness
scalar calculateFlatness(const face& f, const pointField& points)
{
    // Compute an estimate of the centre as the average of the points
    point pAvg = Zero;
    forAll(f, fp)
    {
        pAvg += points[f[fp]];
    }
    pAvg /= f.size();

    // Compute the face area normal and unit normal
    vector sumA = Zero;
    forAll(f, fp)
    {
        const point& p = points[f[fp]];
        const point& pNext = points[f.nextLabel(fp)];

        const vector a = (pNext - p) ^ (pAvg - p);
        sumA += a;
    }
    const vector sumAHat = normalised(sumA);

    // Compute the area-weighted sum of the triangle centres
    scalar sumAn = 0;
    vector sumAnc = Zero;
    forAll(f, fp)
    {
        const point& p = points[f[fp]];
        const point& pNext = points[f.nextLabel(fp)];

        const vector a = (pNext - p) ^ (pAvg - p);
        const vector c = p + pNext + pAvg;

        const scalar an = a & sumAHat;

        sumAn += an;
        sumAnc += an * c;
    }
    point fc = (1.0 / 3.0) * sumAnc / sumAn;

    // Calculate the sum of the magnitude of areas and compare to magnitude 
    // of sum of areas
    scalar summA = 0.0;
    vector sumN = Zero;

    forAll(f, fp)
    {
        const point& thisPoint = points[f[fp]];
        const point& nextPoint = points[f.nextLabel(fp)];

        // Triangle around fc
        const vector n = 0.5 * ((nextPoint - thisPoint) ^ (fc - thisPoint));

        summA += mag(n);
        sumN += n;
    }

    scalar magArea = mag(sumN);
    scalar faceFlatness = magArea / summA;

    return faceFlatness;
}

// Function to compute the interpolation of dz at any mesh point, based on the
// inverse of the distance
point inverseDistanceInterpolationDz(
    const scalar& Ldef, 
    const scalar& alpha, 
    const scalar& coeffVertDeformation,
    const point& internalPoint, 
    const scalarField& boundaryPointsX, 
    const scalarField& boundaryPointsY, 
    const scalarField& boundaryPointsZ, 
    const scalarField& boundaryDz, 
    const scalarField& boundaryDx, 
    const scalarField& boundaryDy, 
    const scalarField& boundaryAreas)
{
    // Initialize variables
    point DeltaInterp;
    const label n = boundaryDz.size();

    // Precompute alpha^5
    const scalar alpha5 = alpha * alpha * alpha * alpha * alpha;

    // Variables for interpolation
    scalar NumZ(0.0);
    scalar NumX(0.0);
    scalar NumY(0.0);
    scalar Den(0.0);
    scalar Den_z(0.0);
    scalar distance_z, LbyD_z, LbyD3_z, weight_z;
    scalar distance, LbyD, LbyD3, weight;
    scalar dist2_xy;
    scalar dist2_z;

    for (label i = 0; i < n; ++i)
    {
    
        dist2_xy = sqr(internalPoint.x() - boundaryPointsX[i]) +
            sqr(internalPoint.y() - boundaryPointsY[i]);
    
        dist2_z = sqr(internalPoint.z() - boundaryPointsZ[i]);
    
        distance = Foam::sqrt(dist2_xy + dist2_z);

        distance_z = Foam::sqrt(dist2_xy + coeffVertDeformation*dist2_z);

        if (distance < 1.e-5)
        {
            DeltaInterp = vector(boundaryDx[i],boundaryDy[i],boundaryDz[i]);
            return DeltaInterp;
        }

        LbyD = Ldef / distance;
        LbyD3 = LbyD * LbyD * LbyD;
        weight = boundaryAreas[i] * (LbyD3 + alpha5 * LbyD3 * LbyD * LbyD);

        NumX += weight * boundaryDx[i];
        NumY += weight * boundaryDy[i];
        Den += weight;        

        if (coeffVertDeformation<1.0)
        {
            LbyD_z = Ldef / distance_z;
            LbyD3_z = LbyD_z * LbyD_z * LbyD_z;
            weight_z = boundaryAreas[i] * (LbyD3_z + alpha5 * LbyD3_z * LbyD_z * LbyD_z);
            
        }
        else
        {
            weight_z = weight;        
        }
        
        NumZ += weight_z * boundaryDz[i];
        Den_z += weight_z; 
    }

    DeltaInterp = vector(NumX/Den,NumY/Den,NumZ/Den_z);
    
    return DeltaInterp;
}

// Function to compute the interpolation of dz at z=0 points, based on the
// inverse of the distance
Tuple2<scalar, scalar> inverseDistanceInterpolationDzBottom(
    const point& internalPoint,
    const scalarField& boundaryPointsX,
    const scalarField& boundaryPointsY,
    const scalarField& boundaryVal1,
    const scalarField& boundaryVal2,
    const scalar& interpRelRadius)
{
    scalar interpolatedVal1(0.0);
    scalar interpolatedVal2(0.0);

    const label n = boundaryVal1.size();
    scalar minValue = GREAT;
    label minIndex = -1;

    // Calculate distances and find the minimum in a single loop
    scalarField distances(n);
    for (label i = 0; i < n; ++i)
    {
        distances[i] = Foam::sqrt(
            Foam::sqr(internalPoint.x() - boundaryPointsX[i]) +
            Foam::sqr(internalPoint.y() - boundaryPointsY[i])
        );

        if (distances[i] < minValue)
        {
            minValue = distances[i];
            minIndex = i;
        }
    }

    // Special case: very close to a boundary point
    if (minValue < 1.e-5)
    {
        interpolatedVal1 = boundaryVal1[minIndex];
        interpolatedVal2 = boundaryVal2[minIndex];
    }
    else
    {
        // General case: weighted interpolation
        scalar NumVal2(0.0), NumVal1(0.0), Den(0.0);

        const scalar radiusThreshold = interpRelRadius * minValue;

        for (label i = 0; i < n; ++i)
        {
            scalar distance = distances[i];
            scalar weight = minValue / (distance*distance);

            // Neglect points outside the relative radius
            if (distance > radiusThreshold)
                weight = 0.0;

            NumVal1 += weight * boundaryVal1[i];
            NumVal2 += weight * boundaryVal2[i];
            Den += weight;
        }

        interpolatedVal1 = NumVal1 / Den;
        interpolatedVal2 = NumVal2 / Den;
    }

    return Tuple2<scalar, scalar>(interpolatedVal1, interpolatedVal2);
}

// Function to calculate the average of a sub-block 
double calculateBlockAverage(
    const RectangularMatrix<double>& elevation, 
    label startRow, 
    label startCol, 
    label blockSize, 
    label maxRows, 
    label maxCols) 
{
    double sum = 0.0;
    label count = 0;

    for (label i = startRow; i < startRow + blockSize && i < maxRows; ++i) {
        for (label j = startCol; j < startCol + blockSize && j < maxCols; ++j) {
            sum += elevation(i, j);
            ++count;
        }
    }
    return sum / count;
}

// Function to subsample a 2D array (used to sub-sample the
// elevation data and save the STL file)
RectangularMatrix<double> subsampleMatrix(
    const RectangularMatrix<double>& elevation, 
    int ncols, 
    int nrows, 
    label blockSize) 
{
    // New dimensions
    label newRows = nrows / blockSize;
    label newCols = ncols / blockSize;

    // Create a new matrix for the subsampled data
    RectangularMatrix<double> subsampled(newRows, newCols, 0.0);

    // Fill the subsampled matrix
    for (label i = 0; i < newRows; ++i) {
        for (label j = 0; j < newCols; ++j) {
            subsampled(i, j) = calculateBlockAverage(
                elevation, 
                i * blockSize, 
                j * blockSize, 
                blockSize, 
                nrows, 
                ncols);
        }
    }
    return subsampled;
}


Tuple2<scalar, scalar> interpolateNegDeformation(
    scalar z,
    const bool useNegDeformation,
    const scalarList& zNeg,
    const scalarList& dxNeg,
    const scalarList& dyNeg)
{
    if (!useNegDeformation) return Tuple2<scalar, scalar>(0.0, 0.0);

    if (z >= 0.0) return Tuple2<scalar, scalar>(0.0, 0.0); // Above or at z=0 → No deformation

    if (z > zNeg[0]) // Interpolate between (0,0) at z=0 and (dxNeg[0], dyNeg[0]) at zNeg[0]
    {
        scalar w = z / zNeg[0]; // Weight factor (z=0 → w=0, z=zNeg[0] → w=1)
        scalar interpDx = w * dxNeg[0];
        scalar interpDy = w * dyNeg[0];
        return Tuple2<scalar, scalar>(interpDx, interpDy);
    }

    if (z <= zNeg.last()) // Below lowest zNeg → Constant deformation
    {
        return Tuple2<scalar, scalar>(dxNeg.last(), dyNeg.last());
    }

    // Find the two closest points for linear interpolation
    for (label i = 0; i < zNeg.size() - 1; ++i)
    {
        if (zNeg[i] >= z && z > zNeg[i + 1])
        {
            scalar w = (z - zNeg[i + 1]) / (zNeg[i] - zNeg[i + 1]); // Interpolation weight
            scalar interpDx = w * dxNeg[i] + (1 - w) * dxNeg[i + 1];
            scalar interpDy = w * dyNeg[i] + (1 - w) * dyNeg[i + 1];
            return Tuple2<scalar, scalar>(interpDx, interpDy);
        }
    }

    return Tuple2<scalar, scalar>(0.0, 0.0); // Should never reach this point
}

//--------------------------------------------------------------

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
  
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

// Read the expansion factor for the top of the mesh
const scalar expFactor = topoDict.lookupOrDefault<scalar>("expFactor",1.0);

// Read the expansion shape parameter for the expansion of the top of the mesh
const scalar exp_shape = topoDict.lookupOrDefault<scalar>("exp_shape",1.0);

// Read the parameter for the starting elevation of the expansion
const scalar dzVert = topoDict.lookupOrDefault<scalar>("dzVert",0.0);

// Read the relative distance for the smoothing kernel of the topography
const scalar interpRelRadius = topoDict.lookupOrDefault<scalar>("interpRelRadius",4.0);

// Read the switch to save the subsampled topo as STL
const Switch saveSTL = topoDict.lookupOrDefault<Switch>("saveSTL", false);

// Read the switch to save the subsampled topo as binary STL
const Switch saveBinary = topoDict.lookupOrDefault<Switch>("saveBinary", false);

// Read the swtich to perform some mesh quality check at the end 
const Switch checkMesh = topoDict.lookupOrDefault<Switch>("checkMesh", false);

// Read the swtich to raise the top of the mesh by the max elev of the topo
const Switch raiseTop = topoDict.lookupOrDefault<Switch>("raiseTop", true);

const Switch orthogonalCorrection = topoDict.lookupOrDefault<Switch>("orthogonalCorrection", false);
const scalar dist_rel1 = topoDict.lookupOrDefault<scalar>("dist_rel1", 0.1);
const scalar dist_rel2 = topoDict.lookupOrDefault<scalar>("dist_rel2", 0.2);

const scalar distC1 = topoDict.lookupOrDefault<scalar>("distC1", 0.0);
const scalar distC2 = topoDict.lookupOrDefault<scalar>("distC2", 0.0);

    
const scalar coeffVertDeformation = topoDict.lookupOrDefault<scalar>("coeffVertDeformation", 1.0);
 
// Initialize empty lists
scalarList zNeg, dxNeg, dyNeg;
bool useNegDeformation = true;

if (topoDict.found("zNeg") && topoDict.found("dxNeg") && topoDict.found("dyNeg"))
{
    zNeg = topoDict.lookup<scalarList>("zNeg");
    dxNeg = topoDict.lookup<scalarList>("dxNeg");
    dyNeg = topoDict.lookup<scalarList>("dyNeg");

    // Ensure zNeg is sorted in decreasing order
    for (label i = 0; i < zNeg.size() - 1; ++i)
    {
        if (zNeg[i] < zNeg[i + 1]) // Should be decreasing
        {
            FatalErrorInFunction << "zNeg list must be sorted in decreasing order (less negative first)" << exit(FatalError);
        }
    }

    if (zNeg.size() != dxNeg.size() || zNeg.size() != dyNeg.size())
    {
        FatalErrorInFunction << "zNeg, dxNeg, and dyNeg must have the same size" << exit(FatalError);
    }

    Info << "Read " << zNeg.size() << " negative deformation levels." << endl;
}
else
{
    Info << "Missing zNeg, dxNeg, or dyNeg in topoGridDict. Horizontal deformation will be set to zero." << endl;
    useNegDeformation = false;
}
 
 
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

        elevation(nrows-1-i, j) = value;
    }
}

if (saveSTL)
{    
    // Subsample the matrix with a factor of 4
    label factor = 2;
    RectangularMatrix<double> elevationSubsampled = subsampleMatrix(elevation, ncols, nrows, factor);

    scalar xllSubsampled(xllcorner+(0.5*factor)*cellsize);
    scalar yllSubsampled(yllcorner+(0.5*factor)*cellsize);

    scalar cellsizeSubsampled(factor*cellsize);
    
    // Create the output STL file name based on the input raster file
    word stlFileName(rasterFile);
    stlFileName.replace(".asc", ".stl");
    Info << "Saving STL file: " << stlFileName << endl;

    // Write the STL surface to a file
        
    if (saveBinary)
    {        
        writeBinarySTL(stlFileName, elevationSubsampled, xllSubsampled, yllSubsampled, cellsizeSubsampled);
    }
    else
    {
        writeSTL(stlFileName, elevationSubsampled, xllSubsampled, yllSubsampled, cellsizeSubsampled);        
    }
    Info << "Saving completed" << endl;
}
    
file.close();

scalar xMin = min(mesh.Cf().component(0)).value();
scalar xMax = max(mesh.Cf().component(0)).value();

reduce(xMin, minOp<scalar>());
reduce(xMax, maxOp<scalar>());

Info << "xMin = " << xMin << endl;
Info << "xMax = " << xMax << endl;

scalar yMin = min(mesh.Cf().component(1)).value();
scalar yMax = max(mesh.Cf().component(1)).value();

reduce(yMin, minOp<scalar>());
reduce(yMax, maxOp<scalar>());

Info << "yMin = " << yMin << endl;
Info << "yMax = " << yMax << endl;

scalar zMin = min(mesh.Cf().component(2)).value();
scalar zMax = max(mesh.Cf().component(2)).value();

reduce(zMin, minOp<scalar>());
reduce(zMax, maxOp<scalar>());
  
Info << "zMin = " << zMin << endl;
Info << "zMax = " << zMax << endl;

word croppedDEMFile = "DEMcropped.asc";
generateCroppedDEM(elevation, xllcorner+xVent, yllcorner+yVent, cellsize, xVent, yVent, xMin, xMax, yMin, yMax, croppedDEMFile);

 
// Approximation of the maximum distance of any mesh node 
// from the mesh centroid (Sen et al, 2017) 
// scalar Ldef(0.5*std::sqrt( sqr(xMax-xMin) + sqr(yMax-yMin) + sqr(zMax-zMin) ));
scalar Ldef(0.5*std::sqrt( sqr(xMax-xMin) + sqr(yMax-yMin) ));

Info << "Ldef = " << Ldef << endl;
 
scalar noDeformLevel(0.5*Ldef);
Info << "noDeformLevel = " << noDeformLevel << endl << endl;
    
scalar z2Rel(0.0);
scalar zNew(0.0);

const vectorField& faceAreas = mesh.faceAreas();
const vectorField& faceCentres = mesh.faceCentres();
const scalarField magFaceAreas(mag(faceAreas));
const vectorField faceNormals = mesh.faceAreas()/mag(faceAreas);
const faceList& faces = mesh.faces();
    
// List of indexes of faces with z=0
labelList z0FaceIndices;

// Pupolate list of faces at z=0, by iterating over all the mesh faces
forAll(faces, faceI)
{
    // Check z of face
    if (mag(faceCentres[faceI].z()) < 1.e-3) // Usa SMALL per tolleranza numerica
    {
        // Add the face index
        z0FaceIndices.append(faceI);
    }
}

Sout << "Proc" << Pstream::myProcNo() << " z=0 faces " << z0FaceIndices.size() << endl;

// Create fieds for face centres coords, areas and dz
scalarField bottomCentresX(z0FaceIndices.size());
scalarField bottomCentresY(z0FaceIndices.size());
scalarField bottomCentresZ(z0FaceIndices.size());
scalarField bottomAreas(z0FaceIndices.size());
scalarField bottomCentresDz(z0FaceIndices.size());

// Loop through each face in the list and compute dz with bilinear interpolation
forAll(z0FaceIndices, facei)
{
    point pCentre = faceCentres[z0FaceIndices[facei]];

    // Get x, y coordinates of the pointi
    scalar x = pCentre.x(); 
    scalar y = pCentre.y();
        
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

        bottomCentresX[facei] = faceCentres[z0FaceIndices[facei]].x();
        bottomCentresY[facei] = faceCentres[z0FaceIndices[facei]].y();
        bottomCentresZ[facei] = faceCentres[z0FaceIndices[facei]].z();
        bottomAreas[facei] =  magFaceAreas[z0FaceIndices[facei]];
        bottomCentresDz[facei] = zInterp-faceCentres[z0FaceIndices[facei]].z();            
    }
    else
    {
        FatalErrorInFunction << "Check asc size" << exit(FatalError);
    }

}
  
// Create list with nproc fields
List<scalarField> CentresX(Pstream::nProcs());
CentresX[Pstream::myProcNo()].setSize(bottomCentresDz.size(), Pstream::myProcNo());

List<scalarField> CentresY(Pstream::nProcs());
CentresY[Pstream::myProcNo()].setSize(bottomCentresDz.size(), Pstream::myProcNo());

List<scalarField> Dz(Pstream::nProcs());
Dz[Pstream::myProcNo()].setSize(bottomCentresDz.size(), Pstream::myProcNo());

List<scalarField> Areas(Pstream::nProcs());
Areas[Pstream::myProcNo()].setSize(bottomCentresDz.size(), Pstream::myProcNo());

// Each processor populate its field in the list of fields
for (label i = 0; i < bottomCentresX.size(); ++i)
{
    CentresX[Pstream::myProcNo()][i] = bottomCentresX[i];
    CentresY[Pstream::myProcNo()][i] = bottomCentresY[i];
    Areas[Pstream::myProcNo()][i] = bottomAreas[i];
    Dz[Pstream::myProcNo()][i] = bottomCentresDz[i];
}

// Use gather and scatter to have the full lists for all the processors    
Pstream::gatherList<scalarField>(CentresX);
Pstream::scatterList<scalarField>(CentresX);

Pstream::gatherList<scalarField>(CentresY);
Pstream::scatterList<scalarField>(CentresY);

Pstream::gatherList<scalarField>(Areas);
Pstream::scatterList<scalarField>(Areas);

Pstream::gatherList<scalarField>(Dz);
Pstream::scatterList<scalarField>(Dz);

// Create the global fields
scalarField globalBottomCentresX;
scalarField globalBottomCentresY;
scalarField globalBottomCentresDz;
scalarField globalBottomCentresAreas;
    
for (label i = 0; i < Pstream::nProcs(); ++i)
{
    globalBottomCentresX.append(CentresX[i]);
    globalBottomCentresY.append(CentresY[i]);
    globalBottomCentresDz.append(Dz[i]);
    globalBottomCentresAreas.append(Areas[i]);
}
   
double maxTopo;
if (raiseTop)
{
    maxTopo = max(globalBottomCentresDz);
}
else
{
    maxTopo = 0.0;
}

scalar zVert(maxTopo + dzVert);

Info << "z=0 faces " << globalBottomCentresAreas.size() << endl;

pointField zeroPoints(mesh.points());
pointField pDeform(0.0*zeroPoints);

const globalIndex globalPoints(mesh.nPoints());

// Local number of points and cells
label localNumPoints = mesh.points().size();

// Output the global number of points and cells
Info << "Local number of points: " << localNumPoints << endl;
    
reduce(localNumPoints, sumOp<scalar>());  // global sum   
    
label globalNumPoints(localNumPoints);
    
// Output the global number of points and cells
Info << "Global number of points: " << globalNumPoints << endl << endl;
    
// Lists for the mesh points at z=0 and for the area and deformation at these points
// These lists are created for each processor
scalarList bottomPointsX;
scalarList bottomPointsY;
scalarList bottomPointsZ;
scalarList bottomPointsArea;
scalarList bottomPointsDz;
labelList  localIdx;

point pEval(zeroPoints[0]);
    
// Loop over the points with z=0 to compute the deformation from the face centers
forAll(pDeform,pointi)
    {
    pEval = mesh.points()[pointi];

    if ( mag(pEval.z()) < 1e-3 )
    {    
        Tuple2<scalar, scalar> result;

        result = inverseDistanceInterpolationDzBottom(pEval, globalBottomCentresX, 
                 globalBottomCentresY, globalBottomCentresDz, 
                 globalBottomCentresAreas, interpRelRadius);

        scalar interpDz = result.first();
        scalar interpArea = result.second();   
            
        bottomPointsX.append( pEval.x() );
        bottomPointsY.append( pEval.y() );
        bottomPointsZ.append( 0.0 );
        bottomPointsArea.append( interpArea );
        bottomPointsDz.append( interpDz );
        localIdx.append( pointi );
        
        zeroPoints[pointi].z() = interpDz;        
    }
} 


//---------------------------------------------------------------------------//

scalarField dxBottom(z0FaceIndices.size());
scalarField dyBottom(z0FaceIndices.size());

if (orthogonalCorrection)
{
    // Loop over the faces to compute the x,y components of normal unit vector
    forAll(z0FaceIndices, facei)
    {
        const face& f = faces[z0FaceIndices[facei]];
        // Compute an estimate of the centre as the average of the points
        point pAvg = Zero;
        forAll(f, fp)
        {
            pAvg += zeroPoints[f[fp]];
        }
        pAvg /= f.size();

        // Compute the face area normal and unit normal
        vector sumA = Zero;
        forAll(f, fp)
        {
            const point& p = zeroPoints[f[fp]];
            const point& pNext = zeroPoints[f.nextLabel(fp)];
            const vector a = (pNext - p) ^ (pAvg - p);
            sumA += a;
        }
        const vector sumAHat = normalised(sumA);

        if ( sumAHat.z() < 0.0 )
        {
            dxBottom[facei] = -sumAHat.x();
            dyBottom[facei] = -sumAHat.y();
        }
        else
        {
            dxBottom[facei] = sumAHat.x();
            dyBottom[facei] = sumAHat.y();
        }
   }
}
else
{
    dxBottom = 0.0;
    dyBottom = 0.0;
}

List<scalarField> Dx(Pstream::nProcs());
Dx[Pstream::myProcNo()].setSize(dxBottom.size(), Pstream::myProcNo());

List<scalarField> Dy(Pstream::nProcs());
Dy[Pstream::myProcNo()].setSize(dyBottom.size(), Pstream::myProcNo());

// Each processor populate its field in the list of fields
for (label i = 0; i < dyBottom.size(); ++i)
{
    Dx[Pstream::myProcNo()][i] = dxBottom[i];
    Dy[Pstream::myProcNo()][i] = dyBottom[i];
}

Pstream::gatherList<scalarField>(Dx);
Pstream::scatterList<scalarField>(Dx);

Pstream::gatherList<scalarField>(Dy);
Pstream::scatterList<scalarField>(Dy);

// Create the global fields
scalarField globalBottomCentresDx;
scalarField globalBottomCentresDy;
    
for (label i = 0; i < Pstream::nProcs(); ++i)
{
    globalBottomCentresDx.append(Dx[i]);
    globalBottomCentresDy.append(Dy[i]);
}


scalarList bottomPointsDx;
scalarList bottomPointsDy;

// Loop over the points with z=0 to interpolate the x,y components of normal 
// unit vector from the face centers
forAll(pDeform,pointi)
    {
    pEval = mesh.points()[pointi];

    if ( mag(pEval.z()) < 1e-3 )
    {    
        Tuple2<scalar, scalar> result;

        result = inverseDistanceInterpolationDzBottom(pEval, globalBottomCentresX, 
                 globalBottomCentresY, globalBottomCentresDx, 
                 globalBottomCentresDy, interpRelRadius);

        bottomPointsDx.append( result.first() );
        bottomPointsDy.append( result.second() );        
    }
} 

//---------------------------------------------------------------------------//

// Create list of labels in the original global mesh
labelList globalIdx(localIdx.size());
forAll(localIdx, pointi)
{
    globalIdx[pointi] = globalPoints.toGlobal(pointi);
}

syncTools::syncPointList
(
    mesh,
    localIdx,
    globalIdx,
    minEqOp<label>(),
    labelMax
);

Sout << "Proc" << Pstream::myProcNo() << " z=0 points " << bottomPointsArea.size() << endl;

// Local number of points and cells
label globalZ0Points = bottomPointsArea.size();
    
reduce(globalZ0Points, sumOp<scalar>());  // global sum   

Info << "Total z=0 points " << globalZ0Points << endl;
      
// Start the computation of mesh deformation for top face centres  
word patchName = "top";  

// Find the ID# associated with the patchName by iterating through boundaryMesh
label patchID = -1;
forAll(mesh.boundaryMesh(), patchi)
{
    if (mesh.boundaryMesh()[patchi].name() == patchName)
    {
        patchID = patchi;
        break;
    }
}

if (patchID == -1)
{
    FatalErrorInFunction << "Patch " << patchName << " not found in mesh." << exit(FatalError);
}

// Access the patch
const fvPatch& patchTop = mesh.boundary()[patchID];

Sout << "Proc" << Pstream::myProcNo() << " zTop faces/points " << patchTop.size() << endl;

// Local number of faces/points 
label globalZtopPoints = patchTop.size();
    
// Global sum  
reduce(globalZtopPoints, sumOp<scalar>());  

Info << "Total z=zTop faces/points " << globalZtopPoints << endl;

label nGlobalPoints(globalZ0Points+globalZtopPoints); 

Info << "Total interpolation points (including duplicated points) " << nGlobalPoints << endl;

// Local lists for top interpolation points (centres of top faces)
scalarField topCentresX(patchTop.size());
scalarField topCentresY(patchTop.size());
scalarField topCentresZ(patchTop.size());
scalarField topCentresAreas(patchTop.size());
scalarField topCentresDz(patchTop.size());
scalarField topCentresDx(patchTop.size());
scalarField topCentresDy(patchTop.size());
labelField globalIdxTop(patchTop.size());

forAll(patchTop, facei)
{
    topCentresX[facei] = faceCentres[patchTop.start() + facei].x();
    topCentresY[facei] = faceCentres[patchTop.start() + facei].y();
    topCentresZ[facei] = faceCentres[patchTop.start() + facei].z();
    topCentresAreas[facei] =  magFaceAreas[patchTop.start() + facei];
    topCentresDz[facei] = maxTopo; 
    topCentresDx[facei] = 0.0; 
    topCentresDy[facei] = 0.0; 
    globalIdxTop[facei] = -1;       
}

// Create lists for sharing processors values for
// both the z=0 and z=zMax interpolation points
    
List<scalarField> concatenatedPointsX(Pstream::nProcs());
concatenatedPointsX[Pstream::myProcNo()].setSize(bottomPointsDz.size() + topCentresDz.size(), Pstream::myProcNo());

List<scalarField> concatenatedPointsY(Pstream::nProcs());
concatenatedPointsY[Pstream::myProcNo()].setSize(bottomPointsDz.size() + topCentresDz.size(), Pstream::myProcNo());

List<scalarField> concatenatedPointsZ(Pstream::nProcs());
concatenatedPointsZ[Pstream::myProcNo()].setSize(bottomPointsDz.size() + topCentresDz.size(), Pstream::myProcNo());

List<scalarField> concatenatedDz(Pstream::nProcs());
concatenatedDz[Pstream::myProcNo()].setSize(bottomPointsDz.size() + topCentresDz.size(), Pstream::myProcNo());

List<scalarField> concatenatedDx(Pstream::nProcs());
concatenatedDx[Pstream::myProcNo()].setSize(bottomPointsDx.size() + topCentresDx.size(), Pstream::myProcNo());

List<scalarField> concatenatedDy(Pstream::nProcs());
concatenatedDy[Pstream::myProcNo()].setSize(bottomPointsDy.size() + topCentresDy.size(), Pstream::myProcNo());

List<scalarField> concatenatedAreas(Pstream::nProcs());
concatenatedAreas[Pstream::myProcNo()].setSize(bottomPointsDz.size() + topCentresDz.size(), Pstream::myProcNo());

List<labelField> concatenatedGlobalIndex(Pstream::nProcs());
concatenatedGlobalIndex[Pstream::myProcNo()].setSize(bottomPointsDz.size() + topCentresDz.size(), Pstream::myProcNo());

// Copy the z=0 interpolation points into the new field
for (label i = 0; i < bottomPointsDz.size(); ++i)
{
    concatenatedPointsX[Pstream::myProcNo()][i] = bottomPointsX[i];
    concatenatedPointsY[Pstream::myProcNo()][i] = bottomPointsY[i];
    concatenatedPointsZ[Pstream::myProcNo()][i] = bottomPointsZ[i];
    concatenatedAreas[Pstream::myProcNo()][i] = bottomPointsArea[i];
    concatenatedDz[Pstream::myProcNo()][i] = bottomPointsDz[i];
    concatenatedDx[Pstream::myProcNo()][i] = bottomPointsDx[i];
    concatenatedDy[Pstream::myProcNo()][i] = bottomPointsDy[i];
    concatenatedGlobalIndex[Pstream::myProcNo()][i] = globalIdx[i];
}

// Copy the z=zMax interpolation points into the new field, 
// starting after the end of the first
for (label i = 0; i < topCentresX.size(); ++i)
{
    concatenatedPointsX[Pstream::myProcNo()][i + bottomPointsX.size()] = topCentresX[i];
    concatenatedPointsY[Pstream::myProcNo()][i + bottomPointsY.size()] = topCentresY[i];
    // concatenatedPointsZ[Pstream::myProcNo()][i + bottomPointsZ.size()] = topCentresZ[i];
    concatenatedPointsZ[Pstream::myProcNo()][i + bottomPointsZ.size()] = noDeformLevel;
    concatenatedAreas[Pstream::myProcNo()][i + bottomPointsArea.size()] = topCentresAreas[i];
    concatenatedDx[Pstream::myProcNo()][i + bottomPointsDx.size()] = topCentresDx[i];
    concatenatedDy[Pstream::myProcNo()][i + bottomPointsDy.size()] = topCentresDy[i];
    concatenatedDz[Pstream::myProcNo()][i + bottomPointsDz.size()] = topCentresDz[i];
    concatenatedGlobalIndex[Pstream::myProcNo()][i + bottomPointsDz.size()] = globalIdxTop[i];
}

// Gather values from other processors
Pstream::gatherList<scalarField>(concatenatedPointsX);
Pstream::gatherList<scalarField>(concatenatedPointsY);
Pstream::gatherList<scalarField>(concatenatedPointsZ);
Pstream::gatherList<scalarField>(concatenatedAreas);
Pstream::gatherList<scalarField>(concatenatedDz);
Pstream::gatherList<scalarField>(concatenatedDx);
Pstream::gatherList<scalarField>(concatenatedDy);
Pstream::gatherList<labelField>(concatenatedGlobalIndex);

// Scatter values from other processors
Pstream::scatterList<scalarField>(concatenatedPointsX);
Pstream::scatterList<scalarField>(concatenatedPointsY);
Pstream::scatterList<scalarField>(concatenatedPointsZ);
Pstream::scatterList<scalarField>(concatenatedAreas);
Pstream::scatterList<scalarField>(concatenatedDz);
Pstream::scatterList<scalarField>(concatenatedDx);
Pstream::scatterList<scalarField>(concatenatedDy);
Pstream::scatterList<labelField>(concatenatedGlobalIndex);

// List of interpolation points merged from all the processors
// containing the z=0 mesh points and the z=zMax face centres
scalarField globalPointsX(nGlobalPoints);
scalarField globalPointsY(nGlobalPoints);
scalarField globalPointsZ(nGlobalPoints);
scalarField globalDz(nGlobalPoints);
scalarField globalDx(nGlobalPoints);
scalarField globalDy(nGlobalPoints);
scalarField globalAreas(nGlobalPoints);
    
// Bool for points alreay added    
boolList addedPoint(globalNumPoints, false);
    
label totPoints(0);
    
// Loop over processors to create global list without duplicated points
for (label i = 0; i < Pstream::nProcs(); ++i)
{
    Info << "Merging proc " << i << endl;
    forAll(concatenatedDz[i],pi)
    {
        // The points belonging to more than one processor
        // should not be added twice
        // "accept" becomes false when the point already exists
        label globalI = concatenatedGlobalIndex[i][pi];            

        bool accept(true);
            
        if ( globalI > 0)
        {
            if ( addedPoint[globalI] )
            {
                accept = false;
            }
            else
            {
                addedPoint[globalI] = true;
            }
        }
            
        if (accept)
        {                
            globalPointsX[totPoints] = concatenatedPointsX[i][pi];
            globalPointsY[totPoints] = concatenatedPointsY[i][pi];
            globalPointsZ[totPoints] = concatenatedPointsZ[i][pi];
            globalDz[totPoints] = concatenatedDz[i][pi];
            globalDx[totPoints] = concatenatedDx[i][pi];
            globalDy[totPoints] = concatenatedDy[i][pi];
            globalAreas[totPoints] = concatenatedAreas[i][pi];
            totPoints++;
        }            
    }
}
    
// Reset the size of the field 
globalPointsX.setSize(totPoints);
globalPointsY.setSize(totPoints);
globalPointsZ.setSize(totPoints);
globalDz.setSize(totPoints);
globalDx.setSize(totPoints);
globalDy.setSize(totPoints);
globalAreas.setSize(totPoints);
    
Info << "Global points for deformation " << totPoints << endl;

// Compute the deformation parameter alpha 
// as in Eq.5 from Luke et al. 2012
scalar gamma = 5.0;
scalarField a_n(globalAreas / sum(globalAreas));
scalar dzMean(sum(a_n*globalDz));
scalar alphaAll = gamma / Ldef * max( mag( globalDz - dzMean ) );

Info << "alpha " << alphaAll << endl;

scalar interpDz(0.0);
scalar interpDx(0.0);
scalar interpDy(0.0);

point DeltaInterp;

const label totalPoints = mesh.points().size();
label maxTotalPoints = totalPoints;
reduce(maxTotalPoints, maxOp<label>());

label localCount = 0;  // Count of processed points by this processor

scalar nextPctg(1.0);
scalar percentage(0.0);
  
scalar dxMin_rel;
scalar dxMax_rel;
    
scalar dyMin_rel;
scalar dyMax_rel;

scalar xCoeff;
scalar yCoeff;  

scalar distC;   
scalar distCoeff;   

scalar coeffHor;

Tuple2<scalar, scalar> result;
      
// Loop over all points in the mesh to interpolate vertical deformation
forAll(pDeform,pointi)
{                
    localCount++;      
    percentage = 100.0 * static_cast<scalar>(localCount) / totalPoints;
        
    // The percentage is computed with respect to the maximum number 
    // of points among all the processors
    scalar GlobalPercentage = 100.0 * static_cast<scalar>(localCount) / maxTotalPoints;

    if ( percentage >= 100.0 )
    {
        Sout << "Proc" << Pstream::myProcNo() << " deformation completed" << endl; 
    }

    if ( GlobalPercentage >= nextPctg )
    {
        Info << "Progress: " << nextPctg << "% completed." << endl;
        nextPctg += 1.0;            
    }
  
    pEval = mesh.points()[pointi];

    if ( mag(pEval.z()-zMax) < 1.e-3 )
    {
        interpDz = maxTopo;  
        interpDx = 0.0;
        interpDy = 0.0;
    } 
    else
    {
        if ( pEval.z() < 1.e-3 )
        {
            // New: Compute horizontal deformation using zNeg, dxNeg, dyNeg
            Tuple2<scalar, scalar> negDeform = interpolateNegDeformation(
                pEval.z(), useNegDeformation, zNeg, dxNeg, dyNeg);
            interpDx = negDeform.first();
            interpDy = negDeform.second();    
            // Info << pEval.z() << " " << interpDx << " " << interpDy << endl;                  
            
            // For points on or below the topography consider only (x,y)
            pEval.z() = 0.0;
            
            result = inverseDistanceInterpolationDzBottom(pEval, globalBottomCentresX, 
                 globalBottomCentresY, globalBottomCentresDz, 
                 globalBottomCentresAreas, interpRelRadius);

            interpDz = result.first();
        }    
        else
        {
            // Interpolation based on full 3D weighted inverse distance 
        
            DeltaInterp = inverseDistanceInterpolationDz(Ldef, alphaAll, 
                        coeffVertDeformation, pEval, 
                        globalPointsX, globalPointsY, globalPointsZ, 
                        globalDz, globalDx, globalDy, globalAreas); 


            if ( pEval.z() > noDeformLevel)
            {
                coeffHor = ( zMax - pEval.z() ) / ( zMax - noDeformLevel );
            }
            else 
            {
                coeffHor = 1.0;
            }
            
            interpDx = DeltaInterp.x()*coeffHor;                       
            interpDy = DeltaInterp.y()*coeffHor;                       
            interpDz = coeffHor * DeltaInterp.z() + (1.0 - coeffHor ) * maxTopo;  
        }                                      
    }
  
    // New elevation of the deformed point, used to compute the enlargement
    zNew = pEval.z() + interpDz;

    // Compute coefficient for horizontal enlargement                
    if (dzVert > 0)
    {
        // Enlarge from a fixed height above the maximum
        // topography and the top, thus from an horizontal
        // plane to the top
        z2Rel = max(0, (zNew - zVert) / (zMax + maxTopo - zVert));
    }
    else
    {
        // Enlarge from the topography to the top
        z2Rel = (zNew - interpDz) / (zMax + maxTopo - interpDz);
    }
    z2Rel = std::pow(z2Rel,exp_shape);

    if ( pEval.z() > 1.e-3 )
    {    
        dxMin_rel = (pEval.x() - xMin)/(xMax-xMin);
        dxMax_rel = (xMax - pEval.x())/(xMax-xMin);
    
        dyMin_rel = (pEval.y() - yMin)/(yMax-yMin);
        dyMax_rel = (yMax - pEval.y())/(yMax-yMin);

        distC = Foam::sqrt(pow(pEval.x(),2) + pow(pEval.y(),2));

        distCoeff = max(0.0,min(1.0,(distC-distC1)/(distC2-distC1)));

        xCoeff = min(distCoeff,max(0.0,min(1.0,(min(dxMin_rel,dxMax_rel)-dist_rel1)/(dist_rel2-dist_rel1))));
        yCoeff = min(distCoeff,max(0.0,min(1.0,(min(dyMin_rel,dyMax_rel)-dist_rel1)/(dist_rel2-dist_rel1))));

        pDeform[pointi].x() = xCoeff * interpDx * pEval.z();
        pDeform[pointi].y() = yCoeff * interpDy * pEval.z();     

        // pDeform[pointi].x() = interpDx * pEval.z();
        // pDeform[pointi].y() = interpDy * pEval.z();     


        pDeform[pointi].x() += z2Rel*(expFactor-1.0)*(pEval.x()+pDeform[pointi].x());
        pDeform[pointi].y() += z2Rel*(expFactor-1.0)*(pEval.y()+pDeform[pointi].y());

    }
    else
    {
        pDeform[pointi].x() = interpDx;
        pDeform[pointi].y() = interpDy;    
    }

    pDeform[pointi].z() = interpDz;
} 
      
mesh.setPoints(mesh.points()+pDeform);

Sout << "Proc" << Pstream::myProcNo() << " mesh updated" << endl; 

if ( checkMesh )
{
    const faceList& pFaces = mesh.faces();
    const pointField& pPts = mesh.points();
    const vectorField& pC = mesh.cellCentres();
    const labelList& pOwner = mesh.faceOwner();

    scalar minQ(1.0);

    forAll(z0FaceIndices, facei)
    {
        const face& f = pFaces[z0FaceIndices[facei]];
        scalar flatness = calculateFlatness(f, pPts);
 
        if ( flatness < 0.98 )
        {   
            Sout << "Proc" << Pstream::myProcNo() << " face " << facei 
                 << " flatness " << flatness << endl; 
        }               

        label oCI = pOwner[z0FaceIndices[facei]];

        point oCc = pC[oCI];

        minQ = 1.0;

        forAll(f, faceBasePtI)
        {
            minQ = minQuality(mesh, oCc, z0FaceIndices[facei], true, faceBasePtI);
        }
        
        if (minQ < 1e-15)
        {
            Sout << "Proc" << Pstream::myProcNo() << " face " 
                 << facei << " minQ " << minQ << endl; 
            forAll(f,ip)
            {
                Sout << ip << " coord " << pPts[f[ip]] << endl;
            }
            Sout << " oCc " << oCc << endl;
        }
    }
}
   
Info << "Writing new mesh" << endl;
mesh.setInstance("constant");
mesh.write();

Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;

Info<< "End\n" << endl;
return 0;
} 


// ************************************************************************* //
