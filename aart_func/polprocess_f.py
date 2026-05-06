from aart_func import *
from params import *
from paramsQU import *
from scipy.spatial import cKDTree

# Calculates velocities

def calculate_kerr_parameters(a):
    """
    Calculates the Horizon radius (rH) and ISCO radius (rms) 
    based on Eqs. B16a-c from the provided PDF.
    """
    # Prevent division by zero for a=1
    if a >= 1.0: a = 0.999999
    
    # Horizon Radius
    rH = 1 + np.sqrt(1 - a**2)
    
    # ISCO Radius (rms) calculation 
    Z1 = 1 + (1 - a**2)**(1/3) * ((1 + a)**(1/3) + (1 - a)**(1/3))
    Z2 = np.sqrt(3 * a**2 + Z1**2)
    rms = 3 + Z2 - np.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2))
    
    return rH, rms

def kerr_disk_velocity(r, a, beta_r, beta_phi, subkep):
    """
    Calculates the 4-velocity components (u^r/u^t, u^phi/u^t) for an accretion flow.
    It automatically detects if 'r' is inside the ISCO and switches to plunging physics.

    Parameters:
    -----------
    r        : float, Radial coordinate
    a        : float, Black hole spin (0 <= a < 1)
    beta_r   : float, Radial velocity mixing (0 = pure infall, 1 = no drift)
    beta_phi : float, Azimuthal velocity mixing
    subkep   : float, Sub-Keplerian factor
    """

    # --- 1. Get Critical Radii ---
    r_H, r_isco = calculate_kerr_parameters(a)

    # --- 2. Metric Potentials (Replaces Delta and PIF) ---
    delta = r**2 - 2*r + a**2
    # The 'A' potential (Eq. B6)
    A_pot = (r**2 + a**2)**2 - a**2 * delta 

    # --- 3. Define Helper: Specific Angular Momentum (lhat) ---
    def get_l_kepler(rad):
        # Calculates angular momentum for a stable circular orbit at radius 'rad'
        return subkep * (rad**2 + a**2 - 2*a*np.sqrt(rad)) / (np.sqrt(rad)*(rad - 2) + a)

    # --- 4. Define "Free-fall" Frame (urbar, Omegabar) ---
    # Radial velocity (pure infall)
    u_r_freefall = -np.sqrt(2*r*(r**2 + a**2)) / (r**2)
    # Angular velocity of the free-fall frame
    omega_freefall = (2*a*r) / A_pot

    # --- 5. Calculate Kinematics ---
    
    if r < r_isco:
        # === PLUNGING REGION (Inside ISCO) ===
        # Physics: Energy and Momentum are conserved from the ISCO values.
        
        # Calculate constants at the ISCO boundary
        L_isco = get_l_kepler(r_isco)
        
        # Calculate Energy at ISCO (Ehat)
        delta_isco = r_isco**2 - 2*r_isco + a**2
        A_pot_isco = (r_isco**2 + a**2)**2 - a**2 * delta_isco
        denom_E = A_pot_isco/(r_isco**2) - (4*a*L_isco)/r_isco - (1 - 2/r_isco)*L_isco**2
        E_isco = np.sqrt(delta_isco / denom_E)

        # 1. Angular Velocity (Omegahat)
        # Driven by the conserved angular momentum L_isco
        omega_kepler = (a + (1 - 2/r)*(L_isco - a)) / (A_pot/(r**2) - (2*a*L_isco)/r)

        # 2. Radial Plunge Velocity (nuhat)
        # Derived from the energy equation with ISCO constants
        potential_term = A_pot/(r**2) - (4*a*L_isco)/r - (1 - 2/r)*L_isco**2 - delta/(E_isco**2)
        # Ensure we don't take sqrt of negative number due to float precision at r=isco
        v_radial_physical = (r / delta) * np.sqrt(np.abs(potential_term))
        # Convert physical radial velocity to 4-velocity component u^r
        u_r_hat = -(delta / r**2) * v_radial_physical * E_isco

        # 3. Apply Mixing
        u_r_mixed = u_r_hat + (1 - beta_r) * (u_r_freefall - u_r_hat)
        omega_mixed = omega_kepler + (1 - beta_phi) * (omega_freefall - omega_kepler)

    else:
        # === STABLE DISK REGION (Outside ISCO) ===
        # Physics: Angular momentum is determined locally by the Keplerian orbit.
        
        L_local = get_l_kepler(r)
        
        # 1. Angular Velocity (Omegahat)
        omega_kepler = (a + (1 - 2/r)*(L_local - a)) / (A_pot/(r**2) - (2*a*L_local)/r)

        # 2. Apply Mixing
        # In the stable disk, the base radial velocity is just the freefall drift
        u_r_mixed = (1 - beta_r) * u_r_freefall
        omega_mixed = omega_kepler + (1 - beta_phi) * (omega_freefall - omega_kepler)

    # --- 6. Normalize 4-Velocity (uttilde) ---
    # We find u^t by enforcing the normalization condition u.u = -1
    norm_numerator = 1 + (u_r_mixed**2 * r**2) / delta
    norm_denominator = 1 - (r**2 + a**2) * omega_mixed**2 - (2/r) * (1 - a*omega_mixed)**2
    
    # Safety check for the denominator (light surface crossing)
    if norm_denominator <= 0:
        return 0.0, 0.0 # Return zero or handle error if velocity is unphysical

    u_t = np.sqrt(norm_numerator / norm_denominator)
    u_phi = u_t * omega_mixed

    # Return normalized 3-velocities (dr/dt, dphi/dt)
    return u_r_mixed / u_t, u_phi / u_t

# TRAJECTORY COMPUTATION

def get_derivatives(t, state, a, beta_r, beta_phi, subkep):
    """
    Computes dr/dt and dphi/dt.
    Corresponding to 'iota' and 'omega'.
    """
    r, phi = state

    # --- Radial Velocity (dr/dt) ---
    iota = kerr_disk_velocity(r,a,beta_r,beta_phi,subkep)[0]
    
    # --- Angular Velocity (dphi/dt) ---
    omega = kerr_disk_velocity(r,a,beta_r,beta_phi,subkep)[1]

    return [iota, omega]

def generate_trajectory(a, beta_r, beta_phi, subkep, r_init, phi_init, n_turns=3, T_max=2000):
    """
    Solves the trajectory and saves to CSV.
    """
    # 1. Get Physics Constants
    rH, rms = calculate_kerr_parameters(a)

    # 2. Define Initial Conditions and Bounds
    r_start = r_init
    phi_start = phi_init
    # Stop just before the Horizon
    r_end = rH * 1.005

    target_phi = phi_init + n_turns * 2 * np.pi

    # EVENT 1: Hit N Turns
    # We want (phi - target) to go from Negative -> Zero.
    # That means the value is INCREASING. Direction = 1.
    def hit_N_turns(t, state, *args):
        return state[1] - target_phi
    hit_N_turns.terminal = True
    hit_N_turns.direction = 1

    # EVENT 2: Hit Horizon
    # We want (r - rH) to go from Positive -> Zero.
    # That means the value is DECREASING. Direction = -1.
    rH = 1 + np.sqrt(1 - a**2) if a < 1.0 else 1.0
    def hit_horizon(t, state, *args):
        return state[0] - (rH * 1.005)
    hit_horizon.terminal = True
    hit_horizon.direction = -1

    # 3. Solver
    y0 = [r_start, phi_start]
    t_span = (0, T_max) # Arbitrary large time, will stop at event

    sol = solve_ivp(
        get_derivatives,
        t_span,
        y0,
        args=(a, beta_r, beta_phi, subkep),
        events=[hit_N_turns, hit_horizon],
        rtol=1e-6,
        atol=1e-6,
        method='RK45'
    )

    # 4. Interpolate for smoother data (optional, similar to PlotPoints->2000)
    # Create 2000 points evenly spaced in time
    t_eval = np.linspace(0, sol.t[-1], 2000)

    # Evaluate solution at these points
    # sol.sol is not available by default unless dense_output=True is set in solve_ivp
    # Rerunning quickly with dense_output or just interpolating manually.
    # To be precise, let's just use the solver's output points if they are dense enough,
    # but re-solving with dense_output is better for smooth curves.

    sol_dense = solve_ivp(
        get_derivatives,
        t_span,
        y0,
        args=(a, beta_r, beta_phi, subkep),
        events=[hit_N_turns, hit_horizon],
        rtol=1e-6,
        atol=1e-6,
        method='RK45',
        dense_output=True
    )

    states = sol_dense.sol(t_eval)
    r_vals = states[0]
    phi_vals = states[1]

    return t_eval, r_vals, phi_vals

def process_polarization(thetao, spin_case, Br, Bth, Bphi, betar, betaphi, sub_kep, gfactor, traj, EVPA_x, EVPA_y, gTotal, n_band=0):
    """
    Iterates through a list of cases, loads H5/CSV data, computes Q&U parameters,
    matches trajectories, and saves the output to text files.
    """

    #print(f"\n--- Processing Case: i={thetao}, a={spin_case}, B=({Br}, {Bth}, {Bphi}), beta_r={betar}, beta_phi={betaphi}, sub_kep={sub_kep} ---")

    # -----------------------------
    # 1. Load HDF5 Files
    # -----------------------------

    # --- Polarization ---
    # polfname=path+"Polarization_a_%f_i_%f_Br_%f_Bth_%f_Bphi_%f_betar_%f_betaphi_%f_sub_kep_%f.h5"%(np.round(spin_case,6),np.round(thetao,6),np.round(Br,6),np.round(Bth,6),np.round(Bphi,6),np.round(betar,6),np.round(betaphi,6),np.round(sub_kep,6))

    # with h5py.File(polfname, 'r') as h5f:
    #     EVPA_x = h5f['EVPA_x'][:]
    #     EVPA_y = h5f['EVPA_y'][:]
    #     gTotal = h5f['gTotal'][:]
    #     isco = h5f['isco'][()]

    # -----------------------------
    # 2. Initialize Variables and Calculate Q&U
    # -----------------------------

    fnrt=path+"Rays_a_%f_i_%f.h5"%(spin_case,thetao)

    h5f = h5py.File(fnrt,'r')
    rs = h5f[f'rs{n_band}'][:]
    phi = h5f[f'phi{n_band}'][:]
    t_ray = h5f[f't{n_band}'][:]
    h5f.close()

    # Calculate Q and U for the whole grid
    ealpha = gTotal**gfactor * EVPA_x
    ebeta = gTotal**gfactor * EVPA_y
    Q_vals = -(ealpha**2 - ebeta**2).real
    U_vals = (-2 * ealpha * ebeta).real #Adding an overall minus sign rotates the image (as seen from below or for an observer with i_case > pi/2)

    # Prepare Grid Points
    r_coords = np.array(rs)
    phi_coords = np.array(phi + np.pi/2) % (2 * np.pi) #Here we add a pi/2 shift to get back to usual convention (raytracing_f has a -pi/2 shift for some reason)
    t_ray_coords = np.array(t_ray)

    # Only include grid points that are not NaN
    mask = ~np.isnan(r_coords) & ~np.isnan(phi_coords) & ~np.isnan(t_ray_coords)
    clean_coords = np.column_stack((r_coords[mask], phi_coords[mask]))
    clean_Q = Q_vals[mask] + Q_shadow
    clean_U = U_vals[mask] + U_shadow
    clean_t_ray = t_ray_coords[mask]

    # -----------------------------
    # 3. Process Trajectory
    # -----------------------------

    tree = cKDTree(clean_coords)

    # Load Data
    t_values = traj[0]
    r_values = traj[1]
    phi_values = (traj[2]) % (2 * np.pi)
    traj_points = np.column_stack((r_values, phi_values))

    # Find the nearest grid point for every trajectory point at once
    # k=1 means "closest neighbor", distance_upper_bound limits the search radius
    dist, indices = tree.query(traj_points, k=1, distance_upper_bound=0.1)

    # Map the results
    # We use a mask for 'inf' distances (where no neighbor was found within the bound)
    valid_hits = dist < np.inf
    Q_out = np.zeros(len(r_values))
    U_out = np.zeros(len(r_values))
    t_ray_out = np.zeros(len(r_values))

    Q_out[valid_hits] = clean_Q[indices[valid_hits]]
    U_out[valid_hits] = clean_U[indices[valid_hits]]
    t_ray_out[valid_hits] = clean_t_ray[indices[valid_hits]]

    return t_values, Q_out, U_out, t_ray_out
