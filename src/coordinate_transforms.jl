# src/coordinate_transforms.jl
using LinearAlgebra

"""
Utility function that rearranges ephemeris results obtained from Horizons
(e.g., via EphemerisSources.jl) into a 6×N matrix.

# Expected input
- `h.x, h.y, h.z`  : arrays of position components
- `h.ẋ, h.ẏ, h.ż`  : arrays of velocity components

# Returns
- A `6×N` matrix (`N` is the number of epochs)
"""
function horizons_to_matrix(h)

    N = length(h.x)                 # Number of time samples (data points)
    T = eltype(h.x)                 # Element type of h.x (e.g., Float64)

    A = Matrix{T}(undef, 6, N)      # Allocate an uninitialized 6×N matrix
                                     # (rows store x, y, z, vx, vy, vz in this order)

    @views begin                   # Avoid unnecessary array copies during slice assignment

        A[1, :] .= h.x             # Assign x positions to the 1st row
        A[2, :] .= h.y             # Assign y positions to the 2nd row
        A[3, :] .= h.z             # Assign z positions to the 3rd row

        A[4, :] .= h.ẋ             # Assign x-velocity to the 4th row
        A[5, :] .= h.ẏ             # Assign y-velocity to the 5th row
        A[6, :] .= h.ż             # Assign z-velocity to the 6th row
    end

    return A                        # Return the 6×N state matrix
end

"""
Compute the column-wise normalized cross products of two 3×N matrices.

# Arguments
- `a::AbstractMatrix` : A 3×N matrix whose columns represent 3-D vectors.
- `b::AbstractMatrix` : A 3×N matrix whose columns represent 3-D vectors.

# Returns
- A 3×N matrix whose columns are the normalized cross products of the
  corresponding columns of `a` and `b`.
"""
function calculate_normalized_cross_product(a::AbstractMatrix, b::AbstractMatrix)

    @assert size(a,1) == 3 && size(b,1) == 3 "a and b must be 3×N matrices"   # Ensure both inputs have exactly 3 rows
    @assert size(a,2) == size(b,2)           "a and b must have the same number of columns"  # Ensure both inputs have the same number of vectors

    N = size(a,2)                            # Number of column vectors

    # Allocate a 3×N matrix to store the normalized cross products
    # The element type is promoted from the input types and Float64
    result = zeros(promote_type(eltype(a), eltype(b), Float64), 3, N)

    @views for i in 1:N                      # Loop over all columns (use views for slicing)
        c = cross(a[:, i], b[:, i])         # Compute the cross product of the i-th column vectors
        result[:, i] .= c / norm(c)         # Normalize the cross product and store it in the result
    end

    return result                            # Return the 3×N matrix of normalized cross products
end

"""
Compute the column-wise angles between two 3×N matrices of 3-D vectors.

# Arguments
- `a::AbstractMatrix` : A 3×N matrix whose columns represent 3-D vectors.
- `b::AbstractMatrix` : A 3×N matrix whose columns represent 3-D vectors.

# Keyword arguments
- `eps::Real` : Threshold for the product of vector norms.  
  If `norm(a[:, i]) * norm(b[:, i]) <= eps`, the angle is treated as
  undefined or zero according to `zero_behavior`.
- `zero_behavior::Symbol` : Behavior when the norm product is smaller
  than or equal to `eps`.
    - `:nan`  → return `NaN`
    - `:zero` → return `0`

# Returns
- A vector of length `N` containing the angles (in radians) between the
  corresponding columns of `a` and `b`.

# Notes
- The computation uses an `atan2`-style formulation via `atan(y, x)`
  to improve numerical robustness compared to `acos`.
"""
function calculate_angle(a::AbstractMatrix, b::AbstractMatrix;
                         eps::Real = 0.0,
                         zero_behavior::Symbol = :nan)

    @assert size(a,1) == 3 && size(b,1) == 3 "a and b must be 3×N matrices"   # Ensure that both inputs have exactly 3 rows
    @assert size(a,2) == size(b,2)           "a and b must have the same number of columns"  # Ensure that both inputs contain the same number of vectors

    N = size(a,2)                            # Number of column vectors
    T = promote_type(eltype(a), eltype(b), Float64)  # Common element type for the output
    alpha = Vector{T}(undef, N)              # Allocate the output vector of angles

    @inbounds for i in 1:N                   # Loop over all column vectors without bounds checking
        ax, ay, az = a[1,i], a[2,i], a[3,i]  # Components of the i-th vector in a
        bx, by, bz = b[1,i], b[2,i], b[3,i]  # Components of the i-th vector in b

        na = sqrt(ax*ax + ay*ay + az*az)     # Euclidean norm of a[:, i]
        nb = sqrt(bx*bx + by*by + bz*bz)     # Euclidean norm of b[:, i]

        if na * nb <= eps                    # Check for near-zero vectors using the norm product
            alpha[i] = (zero_behavior === :zero) ? zero(T) : T(NaN)  # Handle the degenerate case
        else
            d  = ax*bx + ay*by + az*bz       # Dot product a[:, i] ⋅ b[:, i]

            cx = ay*bz - az*by               # x-component of the cross product
            cy = az*bx - ax*bz               # y-component of the cross product
            cz = ax*by - ay*bx               # z-component of the cross product

            cn = sqrt(cx*cx + cy*cy + cz*cz) # Norm of the cross product

            alpha[i] = atan(cn, d)           # Angle between the vectors (atan2-style)
        end
    end

    return alpha                             # Return the vector of angles
end

"""
Rotate a set of 3-D vectors about column-wise rotation axes using
precomputed sine and cosine values.
The rotation is performed using the Rodrigues rotation formula.

# Arguments
- `coordinates::AbstractMatrix` : A 3×N matrix whose columns are the input
  vectors to be rotated.
- `cosθ::AbstractVector` : A length-N vector containing the cosine of the
  rotation angle for each column.
- `sinθ::AbstractVector` : A length-N vector containing the sine of the
  rotation angle for each column.
- `n::AbstractMatrix` : A 3×N matrix whose columns are the rotation axes.

# Returns
- A 3×N matrix whose columns are the rotated vectors.

# Notes
- Each column is rotated independently.
- This function assumes that each column of `n` is a unit vector.
"""
function rotate_coordinates(coordinates::AbstractMatrix,
                            cosθ::AbstractVector,
                            sinθ::AbstractVector,
                            n::AbstractMatrix)

    @assert size(coordinates,1) == 3 "coordinates must be 3×N"              # Ensure that the input coordinates have three rows
    @assert size(n,1) == 3           "n must be 3×N"                        # Ensure that the rotation axes have three rows
    @assert size(coordinates,2) == size(n,2) "coordinates and n must have same N"  # Ensure that both matrices have the same number of columns
    @assert length(cosθ) == size(coordinates,2)                            # Ensure that cosθ has one value per column
    @assert length(sinθ) == size(coordinates,2)                            # Ensure that sinθ has one value per column

    N = size(coordinates, 2)                                               # Number of vectors (columns)
    T = promote_type(eltype(coordinates), eltype(n), eltype(cosθ), eltype(sinθ), Float64)  # Common element type for the output
    out = Matrix{T}(undef, 3, N)                                           # Allocate the output matrix

    @inbounds for i in 1:N                                                  # Loop over all columns without bounds checking
        x, y, z = coordinates[1,i], coordinates[2,i], coordinates[3,i]    # Components of the input vector to be rotated
        nx, ny, nz = n[1,i], n[2,i], n[3,i]                                 # Components of the rotation axis
        c = cosθ[i]; s = sinθ[i]; one_c = one(T) - c                       # cosθ, sinθ, and (1 − cosθ)

        r11 = nx*nx*one_c + c                                               # (1,1) element of the rotation matrix
        r12 = nx*ny*one_c - nz*s                                            # (1,2) element of the rotation matrix
        r13 = nx*nz*one_c + ny*s                                            # (1,3) element of the rotation matrix

        r21 = ny*nx*one_c + nz*s                                            # (2,1) element of the rotation matrix
        r22 = ny*ny*one_c + c                                               # (2,2) element of the rotation matrix
        r23 = ny*nz*one_c - nx*s                                            # (2,3) element of the rotation matrix

        r31 = nz*nx*one_c - ny*s                                            # (3,1) element of the rotation matrix
        r32 = nz*ny*one_c + nx*s                                            # (3,2) element of the rotation matrix
        r33 = nz*nz*one_c + c                                               # (3,3) element of the rotation matrix

        out[1,i] = r11*x + r12*y + r13*z                                    # x-component of the rotated vector
        out[2,i] = r21*x + r22*y + r23*z                                    # y-component of the rotated vector
        out[3,i] = r31*x + r32*y + r33*z                                    # z-component of the rotated vector
    end

    return out                                                             # Return the rotated vectors as a 3×N matrix
end

"""
Transform spacecraft position vectors into a line-fixed rotating frame
defined by the relative motion of a secondary body with respect to a
primary body.

This transformation is performed in two successive rotation steps
for each epoch (column):

1. **Orbital-plane alignment (plane-aligned frame)**  
   The instantaneous orbital plane of the secondary body, defined by the
   relative position and velocity vectors
   `(secondary − primary)`, is rotated so that its plane coincides with
   the reference *xy*-plane.
   The rotation axis is chosen so that the plane normal is aligned with
   the reference axis `ẑ`.

2. **Line alignment inside the plane (line-fixed frame)**  
   After the orbital plane has been aligned with the *xy*-plane, an
   additional in-plane rotation about the `ẑ` axis is applied so that
   the instantaneous direction from the primary body to the secondary
   body is fixed to the +x axis.

The same two rotations are applied to the spacecraft position vectors,
so that the output represents the spacecraft coordinates in a
secondary-body line-fixed frame centered at the primary body.

# Arguments
- `primary_body_wrt_sun::AbstractMatrix` :
  A 6×N state matrix of the primary body expressed with respect to the Sun
  (or a common inertial origin). The first three rows are position and the
  last three rows are velocity.
- `secondary_body_wrt_sun::AbstractMatrix` :
  A 6×N state matrix of the secondary body expressed with respect to the Sun.
- `sc_wrt_sun::AbstractMatrix` :
  A 6×N state matrix of the spacecraft expressed with respect to the Sun.

# Keyword arguments
- `ẑ::AbstractVector` :
  A constant 3-element reference axis used as the target normal direction
  of the orbital plane. The default is `[0, 0, 1]`.

# Returns
- A 3×N matrix containing the spacecraft position vectors expressed in the
  secondary-body line-fixed frame with origin at the primary body.
"""
function make_line_fixed_frame(primary_body_wrt_sun::AbstractMatrix,
                               secondary_body_wrt_sun::AbstractMatrix,
                               sc_wrt_sun::AbstractMatrix;
                               ẑ::AbstractVector = [0.0, 0.0, 1.0])

    secondary_body_wrt_primary = secondary_body_wrt_sun - primary_body_wrt_sun   # Secondary state with respect to the primary body
    sc_wrt_primary             = sc_wrt_sun - primary_body_wrt_sun              # Spacecraft state with respect to the primary body

    @assert size(secondary_body_wrt_primary,1) == 6 "states must be 6×N"         # Ensure the relative state matrix has 6 rows (position + velocity)
    N = size(secondary_body_wrt_primary,2)                                      # Number of epochs (columns)

    r = secondary_body_wrt_primary[1:3, :]                                      # Relative position vectors of the secondary body (3×N)
    v = secondary_body_wrt_primary[4:6, :]                                      # Relative velocity vectors of the secondary body (3×N)

    @assert length(ẑ) == 3 "ẑ must be a 3-element vector"                     # Ensure that the reference axis has length 3
    ẑ_mat = repeat(ẑ, 1, N)                                                    # Replicate the constant reference axis for all epochs (3×N)

    normal   = calculate_normalized_cross_product(r, v)                         # Unit normal vectors of the instantaneous orbital plane (r × v)
    rot_axis = calculate_normalized_cross_product(normal, ẑ_mat)               # Rotation axes that align the plane normal with ẑ
    α        = calculate_angle(normal, ẑ_mat)                                   # Rotation angles between the plane normal and ẑ

    secondary_body_plane_aligned =
        rotate_coordinates(secondary_body_wrt_primary[1:3, :],                  # Secondary-body positions after orbital-plane alignment
                           cos.(α), sin.(α), rot_axis)

    sc_in_secondary_body_plane =
        rotate_coordinates(sc_wrt_primary[1:3, :],                               # Spacecraft positions after the same plane-alignment rotation
                           cos.(α), sin.(α), rot_axis)

    dist = sqrt.(r[1,:].^2 .+ r[2,:].^2 .+ r[3,:].^2)                            # Distance from the primary body to the secondary body

    cosθ = secondary_body_plane_aligned[1,:] ./ dist                             # Cosine of the in-plane rotation angle for line alignment
    sinθ = -secondary_body_plane_aligned[2,:] ./ dist                            # Sine of the in-plane rotation angle for line alignment

    sc_line_fixed =
        rotate_coordinates(sc_in_secondary_body_plane,                            # Apply in-plane rotation about ẑ
                           cosθ, sinθ, ẑ_mat)                                   # to align the secondary-body direction with +x

    return sc_line_fixed                                                         # Spacecraft positions in the line-fixed frame
end


"""
Rotate spacecraft position vectors into a plane-fixed frame whose orientation
is defined at the initial epoch (t₀) by the relative motion of a secondary body
with respect to a primary body.

For each epoch, the instantaneous orbital plane of the secondary body is
defined by the relative position and velocity vectors
`secondary_body_wrt_sun - primary_body_wrt_sun`.  
However, the rotation used in this function is *frozen at the first epoch*:
the rotation axis and the rotation angle computed at t₀ are applied unchanged
to all epochs.

This produces a coordinate frame that is aligned with the secondary-body
orbital plane at t₀ and remains fixed thereafter.

# Arguments
- `primary_body_wrt_sun::AbstractMatrix` :
  A 6×N state matrix of the primary body expressed with respect to the Sun
  (or a common inertial origin). The first three rows are position and the
  last three rows are velocity.
- `secondary_body_wrt_sun::AbstractMatrix` :
  A 6×N state matrix of the secondary body expressed with respect to the Sun.
- `sc_wrt_sun::AbstractMatrix` :
  A 6×N state matrix of the spacecraft expressed with respect to the Sun.

# Keyword arguments
- `ẑ::AbstractVector` :
  A constant 3-element reference axis defining the target normal direction
  of the orbital plane. The default is `[0.0, 0.0, 1.0]`.

# Returns
- A 3×N matrix containing the spacecraft position vectors expressed in the
  plane-fixed frame defined at the first epoch (t₀), with the origin at the
  primary body.
"""
function rotate_to_plane_fixed_t0(primary_body_wrt_sun::AbstractMatrix,
                                  secondary_body_wrt_sun::AbstractMatrix,
                                  sc_wrt_sun::AbstractMatrix;
                                  ẑ::AbstractVector = [0.0, 0.0, 1.0])

    secondary_body_wrt_primary = secondary_body_wrt_sun - primary_body_wrt_sun   # Secondary body state with respect to the primary body
    sc_wrt_primary             = sc_wrt_sun - primary_body_wrt_sun               # Spacecraft state with respect to the primary body

    @assert size(secondary_body_wrt_primary,1) == 6 "states must be 6×N"         # Ensure the relative state matrix has 6 rows (position + velocity)
    N = size(secondary_body_wrt_primary,2)                                      # Number of epochs (columns)

    r = secondary_body_wrt_primary[1:3, :]                                      # Relative position vectors of the secondary body (3×N)
    v = secondary_body_wrt_primary[4:6, :]                                      # Relative velocity vectors of the secondary body (3×N)

    @assert length(ẑ) == 3 "ẑ must be a 3-element vector"                      # Ensure that the reference axis has three elements
    ẑ_mat = repeat(ẑ, 1, N)                                                    # Replicate the constant reference axis for all epochs (3×N)

    normal = calculate_normalized_cross_product(r, v)                           # Unit normal vectors of the instantaneous orbital plane (r × v)
    axis   = calculate_normalized_cross_product(normal, ẑ_mat)                 # Rotation axes to align the plane normal with the reference axis ẑ
    α      = calculate_angle(normal, ẑ_mat)                                    # Rotation angles between the plane normal and the reference axis ẑ

    α_fixed    = fill(α[1], N)                                                   # Fix the rotation angle to the value at the first epoch (t₀)
    axis_fixed = repeat(axis[:, 1], 1, N)                                       # Fix the rotation axis to the value at the first epoch (t₀)

    sc_in_secondary_body_plane_fixed_t0 =
        rotate_coordinates(sc_wrt_primary[1:3, :],                               # Rotate spacecraft position vectors
                           cos.(α_fixed),                                        # using the rotation angle fixed at t₀
                           sin.(α_fixed),                                        # and
                           axis_fixed)                                           # the rotation axis fixed at t₀

    return sc_in_secondary_body_plane_fixed_t0                                   # Return spacecraft positions in the plane-fixed frame
end

# ----------------------------
# Ecliptic J2000 -> EME2000
# ----------------------------

const EPS_J2000 = deg2rad(23.439291111)          # Obliquity of the ecliptic at J2000 [rad]
const COS_EPS_J2000 = cos(EPS_J2000)             # cos(ε), where ε is the J2000 obliquity
const SIN_EPS_J2000 = sin(EPS_J2000)             # sin(ε), where ε is the J2000 obliquity

"""
Convert state vectors from the Ecliptic J2000 frame to the EME2000 frame
(Earth Mean Equator and Equinox of J2000).

The transformation is a rotation about the x-axis by the J2000 obliquity
of the ecliptic. The same rotation is applied to both position and
velocity components.

# Arguments
- `states::AbstractMatrix{Float64}` :
  A 6×N matrix of state vectors expressed in the Ecliptic J2000 frame.
  The first three rows are position components and the last three rows
  are velocity components.

# Returns
- A 6×N matrix containing the transformed state vectors in the EME2000 frame.

# Notes
- This function assumes a constant obliquity value and does not use SPICE.
- The rotation is applied independently to each column.
"""
function eclJ2000_to_EME2000(states::AbstractMatrix{Float64})

    @assert size(states, 1) == 6 "states must be 6×N"     # Ensure that the input matrix has 6 rows (position + velocity)

    N = size(states, 2)                                  # Number of state vectors (columns)

    out = Matrix{Float64}(undef, 6, N)                   # Allocate output matrix for rotated states

    c = COS_EPS_J2000                                    # Cosine of the obliquity angle
    s = SIN_EPS_J2000                                    # Sine of the obliquity angle

    @inbounds for i in 1:N                               # Loop over all state vectors (disable bounds checking)

        x = states[1,i]; y = states[2,i]; z = states[3,i]   # Extract position components in ECLIPJ2000

        out[1,i] = x                                     # x component is unchanged (rotation about x-axis)
        out[2,i] =  c*y - s*z                            # Rotate y component into EME2000
        out[3,i] =  s*y + c*z                            # Rotate z component into EME2000

        xd = states[4,i]; yd = states[5,i]; zd = states[6,i] # Extract velocity components in ECLIPJ2000

        out[4,i] = xd                                    # x velocity is unchanged
        out[5,i] =  c*yd - s*zd                          # Rotate y velocity component
        out[6,i] =  s*yd + c*zd                          # Rotate z velocity component
    end

    return out                                           # Return the rotated 6×N state matrix
end

#=

# NOTE:
# The same ECLIPJ2000 → J2000 rotation can be performed either by the
# analytical constant-obliquity rotation used in `eclJ2000_to_EME2000`
# or by using SPICE frame transformations (pxform).
#
# We tested both implementations and found that the resulting position
# differences are on the order of ~0.5 mm.
#
# Therefore, the SPICE-based implementation below is currently kept
# deactivated and only retained for verification and reference purposes.

using SPICE

# ECLIPJ2000 -> J2000 rotation matrix (constant for these inertial frames)
const R_ECLIPJ2000_TO_J2000 = SPICE.pxform("ECLIPJ2000", "J2000", 0.0)

function eclJ2000_to_J2000_spice(states::AbstractMatrix{Float64})
    @assert size(states, 1) == 6 "states must be 6×N"
    N = size(states, 2)
    out = Matrix{Float64}(undef, 6, N)
    R = R_ECLIPJ2000_TO_J2000

    @inbounds for i in 1:N
        x  = states[1,i]; y  = states[2,i]; z  = states[3,i]
        xd = states[4,i]; yd = states[5,i]; zd = states[6,i]

        out[1,i] = R[1,1]*x + R[1,2]*y + R[1,3]*z
        out[2,i] = R[2,1]*x + R[2,2]*y + R[2,3]*z
        out[3,i] = R[3,1]*x + R[3,2]*y + R[3,3]*z

        out[4,i] = R[1,1]*xd + R[1,2]*yd + R[1,3]*zd
        out[5,i] = R[2,1]*xd + R[2,2]*yd + R[2,3]*zd
        out[6,i] = R[3,1]*xd + R[3,2]*yd + R[3,3]*zd
    end

    return out
end
=#