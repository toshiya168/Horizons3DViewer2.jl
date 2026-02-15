module Horizons3DViewer2

# 依存は各ファイル内で必要最低限 import してOK
include("coordinate_transforms.jl")
include("ephemeris_propagation.jl")
include("traj_visualize.jl")

# ---- Public API（外に出す関数だけ）----
export horizons_to_matrix,
       make_line_fixed_frame,
       rotate_to_plane_fixed_t0,
       eclJ2000_to_EME2000,
       EphemParams,
       solve_satellite,
       calc_heliocentric_specific_energy,
       propagate_future_trajectories_at_each_epoch_for_animation,
       traj_fade!,
       default_axis3,
       set_sun_light!
end # module