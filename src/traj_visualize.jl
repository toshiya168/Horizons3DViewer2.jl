# src/traj_visualize.jl
using GLMakie
using Colors  # to_color を使うため（GLMakieだけで足りれば削ってOK）

"""
Plot a fading trajectory segment and a moving head marker that follow
a frame index observable.

Both the trail and the head marker are automatically updated whenever
`frame_0` changes.

# Arguments
- `ax` :
  A Makie axis (typically `Axis3`) on which the trajectory is drawn.
- `pos_xyz` :
  A 3×M matrix of trajectory positions, where each column is a 3-D point.
- `frame_0::Observable{Int}` :
  Observable index of the current frame (1-based).

# Keyword arguments
- `trail_len::Int` :
  Maximum number of points shown in the trailing segment.
- `linewidth` :
  Line width of the trail.
- `color` :
  Base color of the trajectory (alpha is applied internally).
- `markersize` :
  Marker size of the head point.

# Returns
- A named tuple `(tail, head)` containing the Makie plot objects
  corresponding to the trail line and the head marker.
"""
function traj_fade!(
    ax,                                         # Target axis
    pos_xyz,                                    # 3×M trajectory matrix
    frame_0::Observable{Int};                   # Observable frame index
    trail_len::Int = 800,                       # Length of the fading trail
    linewidth = 1,                              # Line width of the trail
    color = :lime,                              # Base color of the trajectory
    markersize = 10                             # Marker size of the head point
)

    # Observable holding the trailing position segment (3×K)
    pos_obs = @lift begin
        i1 = max(1, $frame_0 - trail_len)       # Starting index of the trail
        @views pos_xyz[1:3, i1:$frame_0]        # View of the trajectory segment
    end

    # Observable holding the per-vertex RGBA colors for the trail (K elements)
    col_obs = @lift begin
        i1 = max(1, $frame_0 - trail_len)       # Starting index of the trail
        n  = $frame_0 - i1 + 1                  # Number of points in the trail

        # Alpha values that increase toward the head of the trail
        alphas = n <= 1 ? Float32[1] :          # Single point → fully opaque
            Float32.(min.(range(0, 2; length=n), 1))

        # Base RGB color converted from the input color
        rgb = RGBf(Makie.to_color(color))

        # Construct RGBA colors with the computed alpha values
        RGBAf.(rgb.r, rgb.g, rgb.b, alphas)
    end

    # Draw the fading trail line
    tail = lines!(ax,                            # Axis on which to draw
                  pos_obs;                       # Observable positions
                  color = col_obs,               # Observable per-vertex colors
                  linewidth,                     # Line width
                  transparency = true)           # Enable alpha blending

    # Draw the head marker at the current frame
    head = scatter!(ax,
        @lift(@views pos_xyz[1:3, $frame_0]'),  # Current position (1×3 for scatter)
        markersize = markersize,                # Marker size
        color = Makie.to_color(color)           # Solid base color for the head
    )

    # Return both plot objects
    return (tail = tail, head = head)
end

"""
Create and configure a 3-D Makie `Axis3` with a consistent default style.

# Arguments
- `parent` :
  The parent layout or figure cell in which the `Axis3` is created.

# Keyword arguments
- `title::AbstractString` :
  Title of the axis.
- `viewmode` :
  Makie camera view mode (default: `:fit`).
- `clip` :
  Whether to clip objects at the axis boundaries.
- `titlegap` :
  Vertical gap between the axis and its title.
- `hide_decorations::Bool` :
  If `true`, hide ticks, spines, and decorations.
- `azimuth` :
  Initial azimuth angle of the camera (in radians).
  If `nothing`, the default value is kept.
- `elevation` :
  Initial elevation angle of the camera (in radians).
  If `nothing`, the default value is kept.
- `scale` :
  If not `nothing`, a symmetric half-range used to set axis limits
  (i.e., `[-scale, scale]` in each direction, before applying `offset`).
- `offset` :
  A 3-tuple `(ox, oy, oz)` added to the center of the axis limits when
  `scale` is specified. This is mainly intended to shift the visible region
  (e.g., along the x direction).

# Returns
- A configured `Axis3` object.
"""
function default_axis3(parent;
    title::AbstractString,                 # Title string for the axis
    viewmode = :fit,                       # Camera view mode
    clip = false,                          # Enable/disable clipping
    titlegap = -20,                        # Gap between title and axis
    hide_decorations::Bool = true,         # Whether to hide ticks and spines
    azimuth = -0.5pi,                      # Initial azimuth angle [rad]
    elevation = 0.5pi,                     # Initial elevation angle [rad]
    scale = nothing,                       # Optional half-range for fixed limits
    offset = (0.0, 0.0, 0.0)                # Optional center offset for limits
    )

    # Create a 3-D axis with common default settings
    ax = Axis3(parent;
        #aspect = :data,                    # Use data aspect ratio
        xlabel = "x", ylabel = "y",        # Axis labels
        viewmode = viewmode,               # Camera view mode
        clip = clip,                       # Clipping option
        title = title,                     # Axis title
        titlegap = titlegap,               # Title gap
        xspinesvisible = false,            # Hide x spines
        yspinesvisible = false,            # Hide y spines
        zspinesvisible = false,            # Hide z spines
        xticksvisible = false,             # Hide x ticks
        yticksvisible = false,             # Hide y ticks
        zticksvisible = false              # Hide z ticks
    )

    # -------------------------
    # decorations on/off
    # -------------------------
    if hide_decorations                      # Check whether decorations should be hidden
        hidedecorations!(ax)                # Hide all decorations of the axis
    end

    # -------------------------
    # view (can be overwritten later)
    # -------------------------
    if azimuth !== nothing                  # If an azimuth angle is specified
        ax.azimuth[] = azimuth              # Set the camera azimuth
    end
    if elevation !== nothing                # If an elevation angle is specified
        ax.elevation[] = elevation          # Set the camera elevation
    end

    # -------------------------
    # scale / limits
    # -------------------------
    if scale !== nothing                    # If a fixed scale is requested
        ox, oy, oz = offset                 # Unpack the center offset
        limits!(ax,                          # Set axis limits
            (-scale + ox, scale + ox),      # x limits
            (-scale + oy, scale + oy),      # y limits
            (-scale + oz, scale + oz)       # z limits
        )
    end

    return ax                               # Return the configured Axis3
end


"""
Set a directional light in a 3-D Makie scene so that the illumination
corresponds to light coming from the Sun direction.

# Arguments
- `ax` :
  A Makie 3-D axis (`Axis3`) whose scene lighting is updated.
- `sun_vec::AbstractVector{<:Real}` :
  Position vector pointing from the origin to the Sun.

# Returns
- `nothing`
"""
function set_sun_light!(ax, sun_vec::AbstractVector{<:Real})

    # Normalize the vector pointing from the origin to the Sun
    dir_to_sun = sun_vec / norm(sun_vec)

    # Light propagation direction is from the Sun toward the origin,
    # hence the negative of the normalized Sun direction
    light_dir  = -Vec3f(dir_to_sun...)

    # Set a single white directional light for the scene
    set_lights!(ax.scene, [DirectionalLight(RGBf(1, 1, 1), light_dir)])

    return nothing
end