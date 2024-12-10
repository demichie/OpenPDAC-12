# Mesh Deformation Procedure in **topoGridNew**

The **topoGridNew** utility deforms the computational mesh to conform to a given topography. The process uses a combination of bilinear interpolation, inverse distance weighting (IDW), and blending techniques to compute vertical displacements and areas. The deformation accounts for contributions from the base ($z = 0$) and the top boundary to compute the deformation at all the points of the mesh.

---

## 1. Bilinear Interpolation for Face Centers with $z = 0$

- **Objective**: Compute vertical displacements ($\Delta z$) and areas for the centers of mesh faces located at $z = 0$, based on the input topography.

- **Procedure**:
  1. The $x, y$ coordinates of each face center are mapped to the corresponding cell in the input topography grid.
  2. **Bilinear interpolation** is performed using the four surrounding grid nodes ($v_{00}, v_{01}, v_{10}, v_{11}$) to estimate the interpolated elevation ($z_{\text{interpolated}}$):
  
     $$z_{\text{interpolated}} = v_{00}(1 - x_{\text{lerp}})(1 - y_{\text{lerp}})
                             + v_{01}x_{\text{lerp}}(1 - y_{\text{lerp}})
                             + v_{10}(1 - x_{\text{lerp}})y_{\text{lerp}}
                             + v_{11}x_{\text{lerp}}y_{\text{lerp}}
     $$
     
  3. The vertical displacement ($\Delta z$) is calculated as:

     $$\Delta z_{\text{face}} = z_{\text{interpolated}} - z_{\text{original}}
     $$

  4. The area of each face ($\text{Area}_{\text{face}}$) is computed based on its geometry and stored for later weighting.

- **Output**:
  - Each face center at $z = 0$ is associated with:
    - Vertical displacement ($\Delta z_{\text{face}}$).
    - Area ($\text{Area}_{\text{face}}$).

---

## 2. First Inverse Distance Weighting (IDW) for Mesh Points at $z = 0$

- **Objective**: Interpolate the vertical displacements ($\Delta z_{\text{face}}$) and areas ($\text{Area}_{\text{face}}$) from the face centers to the mesh points at $z = 0$.

- **Procedure**:
  1. For each mesh point at $z = 0$, calculate its distances ($d_i$) to all $z = 0$ face centers.
  2. Compute the **minimum distance** ($d_{\text{min}}$) to the face centers and define a threshold distance:

     $$d_{\text{threshold}} = \text{interpRelRadius} \cdot d_{\text{min}}
     $$

     - $\text{interpRelRadius}$: A user-defined multiplier (default = 4).
  3. Include only face centers within $d_{\text{threshold}}$ in the summation.
  4. Compute weights ($w_i$) for the included face centers:

     $$w_i = \frac{1}{d_i}, \quad d_i \leq d_{\text{threshold}}
     $$

  5. Interpolate the vertical displacement ($\Delta z$) as:

     $$\Delta z_{\text{p}} = \frac{\sum_i w_i \Delta z_{\text{face},i}}{\sum_i w_i}
     $$

  6. Similarly, interpolate the area ($\text{Area}_{\text{p}}$):

     $$Area_p = \frac{\sum_i w_i \text{Area}_{\text{face},i}}{\sum_i w_i}$$

- **Output**:
  - Each mesh point at $z = 0$ is associated with:
    - Vertical displacement ($\Delta z_{p}$).
    - Interpolated area ($\text{Area}_{p}$).

---

## 3. Merging Displacements and Areas

- **Objective**: Combine the interpolated values from mesh points at $z = 0$ with the fixed displacements and areas associated with the centers of faces of the top boundary.

- **Procedure**:
  1. Prescribed vertical displacements ($\Delta z_{\text{top}}$) are assigned to the top-face centers, typically based on the maximum topography height. The areas are computed from the code. 
  2. Combine these top-face values with the $z = 0$ point values to form a single global list:
     - $\Delta z_{\text{global}}$: Merged list of vertical displacements.
     - $\text{Area}_{\text{global}}$: Merged list of areas.

- **Output**:
  - A single global list of:
    - Vertical displacements ($\Delta z_{\text{global}}$).
    - Areas ($\text{Area}_{\text{global}}$).

---

## 4. Second Inverse Distance Weighting (IDW) for all Mesh Points

- **Objective**: Compute the vertical deformation for all mesh points using the global list of displacements and areas.

- **Procedure**:
1. For each mesh point, compute the weights ($w_i$) for all global points using:

   $$w_i = \text{Area}_{i} \cdot \left[\left(\frac{L}{d_i}\right)^3 + \left(\alpha \cdot \frac{L}{d_i}\right)^5\right]
   $$

   - $L$: Estimated length of the deformation region.
   - $\alpha$: Fraction of $L$, representing the near-body influence region.
   - $d_i$: Euclidean distance to the global point. For the mesh points with $z<0$, we set $z=0$ when computing the distance from the global points.
2. Interpolate the vertical deformation ($\Delta z_{\text{3D}}$):

   $$\Delta z_{\text{3D}} = \frac{\sum_i w_i \Delta z_i}{\sum_i w_i}
   $$


- **Output**:
  - Vertical deformation ($\Delta z_{\text{mesh}}$) for all mesh points.

---

This procedure ensures smooth deformation, accurately capturing local influences at $z = 0$ and global effects across the mesh.

