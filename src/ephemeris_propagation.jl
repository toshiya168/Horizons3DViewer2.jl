# src/ephemeris_propagation.jl
using LinearAlgebra
using DifferentialEquations
using ProgressMeter

# --- GM [km^3/s^2] ---
const GM_sun     = 1.3271244004193938e11
const GM_earth   = 398600.435436
const GM_moon    = 4902.800066
const GM_jupiter = 126686531.900
const GM_mars    = 37931206.234

# --- Radii [km] ---
const earth_r   = 6378.137
const moon_r    = 1737.4
const mars_r    = 3396.2
const jupiter_r = 69911.0

"""
Container of ephemeris data and time discretization parameters used for
spacecraft propagation driven by tabulated planetary states.

# Fields
- `earth::AbstractMatrix{T}` :
  A 6×N ephemeris state matrix of the Earth.
  The first three rows represent position and the last three rows represent
  velocity, sampled at uniform time steps.
- `moon::AbstractMatrix{T}` :
  A 6×N ephemeris state matrix of the Moon.
- `mars::AbstractMatrix{T}` :
  A 6×N ephemeris state matrix of Mars.
- `jupiter::AbstractMatrix{T}` :
  A 6×N ephemeris state matrix of Jupiter.
- `step_size_sec::S` :
  Sampling interval of the ephemeris data in seconds.
- `sol_from::Int` :
  Column index corresponding to the reference epoch (initial column)
  used to map simulation time to ephemeris indices.

# Notes
- All ephemeris matrices are assumed to share the same time grid and
  indexing convention.
- Linear interpolation between neighboring samples is performed using
  these parameters.
"""
Base.@kwdef struct EphemParams{T, A<:AbstractMatrix{T}, S<:Real}

    earth::A            # Ephemeris state matrix of the Earth (6×N, position and velocity)
    moon::A             # Ephemeris state matrix of the Moon (6×N, position and velocity)
    mars::A             # Ephemeris state matrix of Mars (6×N, position and velocity)
    jupiter::A          # Ephemeris state matrix of Jupiter (6×N, position and velocity)

    step_size_sec::S   # Uniform sampling interval of the ephemeris data [s]
    sol_from::Int      # Starting column index corresponding to t = 0 in the ephemeris tables

end

"""
Linearly interpolate a 3-D position vector from a tabulated 6×N ephemeris
state matrix at a given continuous time.

# Arguments
- `ephem::AbstractMatrix` :
  A 6×N ephemeris state matrix. The first three rows correspond to the
  position components.
- `t::Real` :
  Continuous time (in seconds) measured from the reference epoch
  corresponding to column index `sol_from`.
- `step_size_sec::Real` :
  Sampling interval of the ephemeris data in seconds.
- `sol_from::Int` :
  Column index in `ephem` corresponding to the reference epoch (t = 0).

# Returns
- A length-3 vector containing the linearly interpolated position
  at time `t`.

# Notes
- The interpolation is performed between the two nearest samples
  surrounding the fractional column index.
- Values are clamped to the valid range of the ephemeris table.
"""
@inline function interp_pos_6xN(ephem::AbstractMatrix,   # Ephemeris state matrix (6×N)
                                t::Real,                 # Continuous time since the reference epoch [s]
                                step_size_sec::Real,     # Ephemeris sampling interval [s]
                                sol_from::Int)           # Reference column index corresponding to t = 0

    col = t / step_size_sec + sol_from                   # Fractional column index corresponding to time t
    n = size(ephem, 2)                                   # Total number of ephemeris samples (columns)

    c0 = clamp(floor(Int, col), 1, n - 1)                # Index of the left interpolation sample
    c1 = c0 + 1                                          # Index of the right interpolation sample

    α = col - c0                                         # Interpolation weight in [0, 1]

    @views p0 = ephem[1:3, c0]                           # Position vector at column c0
    @views p1 = ephem[1:3, c1]                           # Position vector at column c1

    return p0 .+ (p1 .- p0) .* α                         # Linearly interpolated 3-D position vector
end

"""
Compute the time derivative of a spacecraft state under point-mass
gravitational attraction from the Sun, Earth, Moon, Mars, and Jupiter,
using tabulated ephemerides for the planetary positions.

The planetary positions are obtained by linear interpolation of the
pre-sampled ephemeris tables stored in `EphemParams`.

# Arguments
- `du` :
  Output vector to be filled with the time derivatives of the state.
- `u` :
  Current spacecraft state vector `[x, y, z, vx, vy, vz]`.
- `p::EphemParams` :
  Structure containing ephemeris tables and interpolation parameters
  for the perturbing bodies.
- `t` :
  Current time (seconds) measured from the reference epoch corresponding
  to `p.sol_from`.

# Returns
- `nothing` (the result is written in-place to `du`).
"""
function spacecraft_motion!(du, u, p::EphemParams, t)

    # Interpolate planetary positions at the current time
    earth_pos   = interp_pos_6xN(p.earth,   t, p.step_size_sec, p.sol_from)   # Earth position [km]
    moon_pos    = interp_pos_6xN(p.moon,    t, p.step_size_sec, p.sol_from)   # Moon position [km]
    mars_pos    = interp_pos_6xN(p.mars,    t, p.step_size_sec, p.sol_from)   # Mars position [km]
    jupiter_pos = interp_pos_6xN(p.jupiter, t, p.step_size_sec, p.sol_from)   # Jupiter position [km]

    @views r_sun = u[1:3]                                                     # Spacecraft position with respect to the Sun
    r_earth   = r_sun .- earth_pos                                           # Spacecraft position with respect to the Earth
    r_moon    = r_sun .- moon_pos                                            # Spacecraft position with respect to the Moon
    r_mars    = r_sun .- mars_pos                                            # Spacecraft position with respect to Mars
    r_jupiter = r_sun .- jupiter_pos                                         # Spacecraft position with respect to Jupiter

    du[1] = u[4]                                                             # dx/dt = vx
    du[2] = u[5]                                                             # dy/dt = vy
    du[3] = u[6]                                                             # dz/dt = vz

    inv_sun3     = inv(norm(r_sun)^3)                                        # 1 / |r_sun|^3
    inv_earth3   = inv(norm(r_earth)^3)                                      # 1 / |r_earth|^3
    inv_moon3    = inv(norm(r_moon)^3)                                       # 1 / |r_moon|^3
    inv_mars3    = inv(norm(r_mars)^3)                                       # 1 / |r_mars|^3
    inv_jupiter3 = inv(norm(r_jupiter)^3)                                    # 1 / |r_jupiter|^3

    # x-component of the gravitational acceleration
    du[4] = -(GM_sun*inv_sun3*r_sun[1] + GM_earth*inv_earth3*r_earth[1] +
              GM_moon*inv_moon3*r_moon[1] + GM_mars*inv_mars3*r_mars[1] +
              GM_jupiter*inv_jupiter3*r_jupiter[1])

    # y-component of the gravitational acceleration
    du[5] = -(GM_sun*inv_sun3*r_sun[2] + GM_earth*inv_earth3*r_earth[2] +
              GM_moon*inv_moon3*r_moon[2] + GM_mars*inv_mars3*r_mars[2] +
              GM_jupiter*inv_jupiter3*r_jupiter[2])

    # z-component of the gravitational acceleration
    du[6] = -(GM_sun*inv_sun3*r_sun[3] + GM_earth*inv_earth3*r_earth[3] +
              GM_moon*inv_moon3*r_moon[3] + GM_mars*inv_mars3*r_mars[3] +
              GM_jupiter*inv_jupiter3*r_jupiter[3])

    return nothing                                                           # In-place update of du
end

"""
Propagate a spacecraft state under multi-body point-mass gravity using
tabulated ephemerides and terminate the propagation when a collision with
any major body is detected.

The spacecraft is propagated in a Sun-centered inertial frame using the
right-hand side defined by `spacecraft_motion!`. The positions of the
Earth, Moon, Mars, and Jupiter are obtained by linear interpolation of
precomputed ephemeris tables.

A termination event is triggered when the spacecraft distance to the
surface of any of the bodies (Earth, Moon, Mars, or Jupiter) becomes zero
or negative.

# Arguments
- `sc1::AbstractMatrix` :
  A 6×N matrix of spacecraft states sampled on the same time grid as the
  ephemerides. The column indexed by `sol_from` is used as the initial
  condition.
- `earth::AbstractMatrix` :
  A 6×N ephemeris state matrix of the Earth.
- `moon::AbstractMatrix` :
  A 6×N ephemeris state matrix of the Moon.
- `mars::AbstractMatrix` :
  A 6×N ephemeris state matrix of Mars.
- `jupiter::AbstractMatrix` :
  A 6×N ephemeris state matrix of Jupiter.

# Keyword arguments
- `sol_from::Int` :
  Column index used as the initial epoch of the propagation.
- `step_size_sec::Real` :
  Sampling interval of the ephemeris data in seconds.
- `reltol` :
  Relative tolerance passed to the ODE solver.
- `abstol` :
  Absolute tolerance passed to the ODE solver.

# Returns
- An `ODESolution` object returned by DifferentialEquations.jl containing
  the propagated spacecraft state history.
"""
function solve_satellite(sc1::AbstractMatrix,                          # Spacecraft ephemeris/state table (6×N)
                         earth::AbstractMatrix,                        # Earth ephemeris (6×N)
                         moon::AbstractMatrix,                         # Moon ephemeris (6×N)
                         mars::AbstractMatrix,                         # Mars ephemeris (6×N)
                         jupiter::AbstractMatrix;                      # Jupiter ephemeris (6×N)
                         sol_from::Int,                                # Initial column index
                         step_size_sec::Real,                          # Ephemeris sampling interval [s]
                         reltol = 1e-10,                               # Relative solver tolerance
                         abstol = 1e-10)                               # Absolute solver tolerance

    @views u0 = vec(sc1[1:6, sol_from])                                 # Initial spacecraft state vector at t = 0

    Nmin = min(size(earth,2), size(moon,2), size(sc1,2),               # Smallest number of samples among all tables
               size(mars,2), size(jupiter,2))

    tmax = (Nmin - sol_from - 1) * step_size_sec                       # Final propagation time [s]
    tspan = (0.0, tmax)                                                 # Time span for the ODE solver

    p = EphemParams(earth = earth,                                     # Pack ephemerides and interpolation parameters
                    moon = moon,
                    mars = mars,
                    jupiter = jupiter,
                    step_size_sec = step_size_sec,
                    sol_from = sol_from)

    # ------------------------------------------------------------------
    # Collision event: terminate when the spacecraft intersects the
    # surface of any body (Earth, Moon, Mars, or Jupiter)
    # ------------------------------------------------------------------
    function collision_g(u, t, integrator)                              # Event function evaluated by the solver

        pp = integrator.p                                               # Retrieve ephemeris parameters

        earth_pos   = interp_pos_6xN(pp.earth,   t,                    # Interpolated Earth position
                                      pp.step_size_sec, pp.sol_from)
        moon_pos    = interp_pos_6xN(pp.moon,    t,                    # Interpolated Moon position
                                      pp.step_size_sec, pp.sol_from)
        mars_pos    = interp_pos_6xN(pp.mars,    t,                    # Interpolated Mars position
                                      pp.step_size_sec, pp.sol_from)
        jupiter_pos = interp_pos_6xN(pp.jupiter, t,                    # Interpolated Jupiter position
                                      pp.step_size_sec, pp.sol_from)

        r = @view u[1:3]                                                # Current spacecraft position vector

        dE  = norm(r .- earth_pos)   - earth_r                         # Signed distance to Earth's surface
        dM  = norm(r .- moon_pos)    - moon_r                          # Signed distance to Moon's surface
        dMa = norm(r .- mars_pos)    - mars_r                          # Signed distance to Mars' surface
        dJu = norm(r .- jupiter_pos) - jupiter_r                       # Signed distance to Jupiter's surface

        return min(dE, dM, dMa, dJu)                                    # Event root when any distance becomes zero
    end

    stop!(integrator) = terminate!(integrator)                          # Terminate integration at the event

    cb = ContinuousCallback(collision_g, stop!; rootfind = true)        # Continuous event callback

    prob = ODEProblem(spacecraft_motion!, u0, tspan, p)                # Define the ODE problem

    ts  = 0.0:step_size_sec:tmax                                        # Output time grid

    sol = solve(prob, Vern7();                                          # Integrate using a 7th-order Runge–Kutta method
                callback = cb,                                         # Collision termination callback
                reltol = reltol,                                       # Relative tolerance
                abstol = abstol,                                       # Absolute tolerance
                saveat = ts)                                           # Save solution at the specified times

    return sol                                                          # Return the ODE solution
end




"""
Propagate "future" spacecraft trajectories for animation, starting from multiple epochs.

This function scans the time index `sol_from` over the input ephemeris/state histories
with a stride `sampling_interval`. For each `sol_from`, it numerically propagates the
spacecraft trajectory forward (via `solve_satellite`) and stores:
- the raw propagated solution (`sol`)
- the spacecraft trajectory expressed in an Earth-centered Earth–Moon line-fixed frame (`sol_em_fixed_e`)
- the spacecraft trajectory expressed in a Moon-centered Earth–Moon line-fixed frame (`sol_em_fixed_m`)
- the spacecraft trajectory expressed in a plane-fixed frame defined at t0 (`sol_lop_e`)

Returns
-------
results_dict::Dict{Int, NamedTuple}
    key: `sol_from` (start index of propagation)
    value: `(sol, sol_em_fixed_e, sol_em_fixed_m, sol_lop_e)`
"""
function propagate_future_trajectories_at_each_epoch_for_animation(  # Propagate future trajectories from many start epochs (for animation)
    sc1_states,       # Spacecraft state history (input to the propagator)
    earth_states,     # Earth state history (used for propagation and frame transforms)
    moon_states,      # Moon state history (used for propagation and frame transforms)
    mars_states,      # Mars state history (optional perturbation / propagator input)
    jupiter_states;   # Jupiter state history (optional perturbation / propagator input)
    sampling_interval::Int,     # Stride (in indices) for choosing propagation start epochs
    step_size_sec,              # Propagation step size [s]
    reltol::Real = 1e-10,       # ODE solver relative tolerance
    abstol::Real = 1e-10,       # ODE solver absolute tolerance
    show_progress::Bool = true, # Whether to show a progress bar
)
    @views begin  # Use views for slices to avoid unnecessary allocations/copies

        N = size(earth_states, 2)  # Total number of time samples in the ephemeris/state histories

        sol_froms = 1:sampling_interval:(N-1)  # List/range of start indices for future propagation

        results_dict = Dict{Int, NamedTuple}()        # Dictionary to store results keyed by start index
        sizehint!(results_dict, length(sol_froms))    # Pre-allocate hash table capacity to reduce rehashing

        p = show_progress ?                               # Create a progress bar only if requested
            Progress(length(sol_froms); desc="Propagating trajectories") :
            nothing                                       # Otherwise keep `p` as `nothing`

        for sol_from in sol_froms  # Loop over each propagation start epoch (index)

            # --- Propagate (in the inertial/propagation frame) ---
            sol = solve_satellite(  # Numerically integrate the spacecraft dynamics from `sol_from` forward
                sc1_states, earth_states, moon_states, mars_states, jupiter_states;
                sol_from = sol_from,           # Start index for propagation
                step_size_sec = step_size_sec, # Time step for integration
                reltol = reltol,               # Relative tolerance
                abstol = abstol                # Absolute tolerance
            )

            sol_N = size(sol, 2)  # Number of samples in the propagated solution `sol` (as stored/returned)

            # --- Safe slicing (avoid indexing past the end of the original histories) ---
            j_end = min(sol_from + sol_N - 1, N)  # Last usable index within the original ephemeris arrays
            sol_N_eff = j_end - sol_from + 1      # Effective number of propagated samples that align with ephemerides

            # --- Slice once and reuse (avoid repeating the same slices) ---
            e = earth_states[1:6, sol_from:j_end]  # 6D Earth state segment aligned with the valid propagation window
            m = moon_states[1:6,  sol_from:j_end]  # 6D Moon state segment aligned with the valid propagation window
            s = sol[1:6, 1:sol_N_eff]              # 6D spacecraft state segment truncated to the valid aligned window

            # --- Frame transforms ---
            s_em_fixed_e = make_line_fixed_frame(e, m, s)  # Spacecraft in Earth-centered Earth–Moon line-fixed frame

            m_em_fixed_e = make_line_fixed_frame(e, m, m)  # Moon itself expressed in the same line-fixed frame (reused)

            s_em_fixed_m = s_em_fixed_e - m_em_fixed_e     # Convert to Moon-centered line-fixed representation

            s_lop_e = rotate_to_plane_fixed_t0(e, m, s)    # Spacecraft in a plane-fixed frame defined at t0

            results_dict[sol_from] = (  # Store results keyed by the propagation start index
                sol            = sol,           # Raw propagated solution (in the propagation/inertial frame representation)
                sol_em_fixed_e = s_em_fixed_e,  # Earth-centered EM line-fixed trajectory
                sol_em_fixed_m = s_em_fixed_m,  # Moon-centered EM line-fixed trajectory
                sol_lop_e      = s_lop_e,       # Plane-fixed trajectory
            )

            show_progress && next!(p)  # Advance the progress bar by one step (only if enabled)
        end

        return results_dict  # Return all propagated trajectories indexed by their start epochs
    end
end

"""
Compute the heliocentric specific mechanical energy time series
for a spacecraft.

The specific mechanical energy is defined as:

    ε = 1/2 |v|² − μ / |r|

where:
- r is the position vector relative to the Sun
- v is the velocity vector relative to the Sun
- μ is the gravitational parameter of the Sun

Arguments
---------
sc1_states::AbstractMatrix
    A 6×N state matrix in a Sun-centered inertial frame.
    Rows 1:3 = position components (x, y, z)
    Rows 4:6 = velocity components (vx, vy, vz)

Keyword Arguments
-----------------
mu::Real (default = GM_sun)
    Gravitational parameter of the Sun in consistent units.

Returns
-------
Vector{T}
    A length-N vector containing the specific mechanical energy
    at each time step.
"""
function calc_heliocentric_specific_energy(
    sc1_states::AbstractMatrix;        # 6×N Sun-centered state matrix
    mu::Real = GM_sun                 # Gravitational parameter of the Sun
)
    @views begin                      # Avoid copying when slicing arrays

        # Ensure the input matrix has at least 6 rows (position + velocity)
        @assert size(sc1_states, 1) ≥ 6 "sc1_states must be at least 6×N"

        # Number of time samples
        N = size(sc1_states, 2)

        # Determine a common numeric type for safe promotion
        T = promote_type(eltype(sc1_states), typeof(mu))

        # Preallocate output energy vector
        ε = Vector{T}(undef, N)

        @inbounds for k in 1:N        # Disable bounds checking for performance

            # Extract heliocentric position components at time step k
            rx = sc1_states[1, k]
            ry = sc1_states[2, k]
            rz = sc1_states[3, k]

            # Extract heliocentric velocity components at time step k
            vx = sc1_states[4, k]
            vy = sc1_states[5, k]
            vz = sc1_states[6, k]

            # Compute squared radial distance |r|^2
            r2 = rx*rx + ry*ry + rz*rz

            # Compute squared speed |v|^2
            v2 = vx*vx + vy*vy + vz*vz

            # Compute specific mechanical energy
            ε[k] = 0.5*v2 - mu / sqrt(r2)
        end

        # Return the time series of specific energy
        return ε
    end
end