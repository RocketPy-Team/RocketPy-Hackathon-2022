# -*- coding: utf-8 -*-
__author__ = "Franz Masatoshi Yuri, Lucas Kierulff Balabram, Guilherme Fernandes Alves"
__copyright__ = "Copyright 20XX, RocketPy Team"
__license__ = "MIT"

import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

from .Environment import Environment
from .Function import Function


# TODO: Needs tests
def compute_CdS_from_drop_test(
    terminal_velocity, rocket_mass, air_density=1.225, g=9.80665
):
    """Returns the parachute's CdS calculated through its final speed, air
    density in the landing point, the rocket's mass and the force of gravity
    in the landing point.

    Parameters
    ----------
    terminal_velocity : float
        Rocket's speed in m/s when landing.
    rocket_mass : float
        Rocket's dry mass in kg.
    air_density : float, optional
        Air density, in kg/m^3, right before the rocket lands. Default value is 1.225.
    g : float, optional
        Gravitational acceleration experienced by the rocket and parachute during
        descent in m/s^2. Default value is the standard gravity, 9.80665.

    Returns
    -------
    CdS : float
        Number equal to drag coefficient times reference area for parachute.

    """

    return 2 * rocket_mass * g / ((terminal_velocity**2) * air_density)


# TODO: Needs tests


def calculateEquilibriumAltitude(
    rocket_mass,
    CdS,
    z0,
    v0=0,
    env=None,
    eps=1e-3,
    max_step=0.1,
    seeGraphs=True,
    g=9.80665,
    estimated_final_time=10,
):
    """Returns a dictionary containing the time, altitude and velocity of the
    system rocket-parachute in which the terminal velocity is reached.


    Parameters
    ----------
    rocket_mass : float
        Rocket's mass in kg.
    CdS : float
        Number equal to drag coefficient times reference area for parachute.
    z0 : float
        Initial altitude of the rocket in meters.
    v0 : float, optional
        Rocket's initial speed in m/s. Must be negative
    env : Environment, optional
        Environmental conditions at the time of the launch.
    eps : float, optional
        acceptable error in meters.
    max_step: float, optional
        maximum allowed time step size to solve the integration
    seeGraphs : boolean, optional
        True if you want to see time vs altitude and time vs speed graphs,
        False otherwise.
    g : float, optional
        Gravitational acceleration experienced by the rocket and parachute during
        descent in m/s^2. Default value is the standard gravity, 9.80665.
    estimated_final_time: float, optional
        Estimative of how much time (in seconds) will spend until vertical terminal
        velocity is reached. Must be positive. Default is 10. It can affect the final
        result if the value is not high enough. Increase the estimative in case the
        final solution is not founded.


    Returns
    -------
    altitudeFunction: Function
        Altitude as a function of time. Always a Function object.
    velocityFunction:
        Vertical velocity as a function of time. Always a Function object.
    final_sol : dictionary
        Dictionary containing the values for time, altitude and speed of
        the rocket when it reaches terminal velocity.
    """
    final_sol = {}

    if not v0 < 0:
        print("Please set a valid negative value for v0")
        return None

    # TODO: Improve docs
    def check_constant(f, eps):
        """_summary_

        Parameters
        ----------
        f : array, list
            _description_
        eps : float
            _description_

        Returns
        -------
        int, None
            _description_
        """
        for i in range(len(f) - 2):
            if abs(f[i + 2] - f[i + 1]) < eps and abs(f[i + 1] - f[i]) < eps:
                return i
        return None

    if env == None:
        environment = Environment(
            railLength=5.0,
            latitude=0,
            longitude=0,
            elevation=1000,
            date=(2020, 3, 4, 12),
        )
    else:
        environment = env

    # TODO: Improve docs
    def du(z, u):
        """_summary_

        Parameters
        ----------
        z : float
            _description_
        u : float
            velocity, in m/s, at a given z altitude

        Returns
        -------
        float
            _description_
        """
        return (
            u[1],
            -g + environment.density(z) * ((u[1]) ** 2) * CdS / (2 * rocket_mass),
        )

    u0 = [z0, v0]

    us = solve_ivp(
        fun=du,
        t_span=(0, estimated_final_time),
        y0=u0,
        vectorized=True,
        method="LSODA",
        max_step=max_step,
    )

    constant_index = check_constant(us.y[1], eps)

    # TODO: Improve docs by explaining what is happening below with constant_index
    if constant_index is not None:
        final_sol = {
            "time": us.t[constant_index],
            "altitude": us.y[0][constant_index],
            "velocity": us.y[1][constant_index],
        }

    altitudeFunction = Function(
        source=np.array(list(zip(us.t, us.y[0])), dtype=np.float64),
        inputs="Time (s)",
        outputs="Altitude (m)",
        interpolation="linear",
    )

    velocityFunction = Function(
        source=np.array(list(zip(us.t, us.y[1])), dtype=np.float64),
        inputs="Time (s)",
        outputs="Vertical Velocity (m/s)",
        interpolation="linear",
    )

    if seeGraphs:
        altitudeFunction()
        velocityFunction()

    return altitudeFunction, velocityFunction, final_sol


# TODO: Needs tests

def traj(flight, postProcess=True, n=1000):
    '''
    Convert a flight into a trajectory

    flight          Flight
    postProcess     if True, make sure the flight is postProcessed
    n               Number of points for reinterpolation

    Returns
        trajectory  ndarray([(s, x, y, z), ...])
    '''
    if postProcess and not flight.postProcessed:
        flight.postProcess()
    s = np.linspace(0, 1, n)
    tq = flight.x.source[:, 0]
    t0 = tq.min()
    tf = tq.max()
    dt = tf - t0
    t = t0 + dt * s
    x = flight.x(t)
    y = flight.y(t)
    z = flight.z(t)
    return np.stack([s, x, y, z], 1)

def flight_stats(flights, postProcess=True, n=1000):
    '''
    Prepare statistics for a set of flights

    flights         iterable of Flight
    postProcess     if True, make sure each Flight has been postProcessed before proceeding
    n               number of timepoints to sample during reinterpolation

    Returns
        mean        trajectory [(s, x, y, z), ...] representing the mean non-dimensionalized flight
        lower       trajectory [(s, x, y, z), ...] representing a flight one standard deviation lower than the mean
        upper       trajectory [(s, x, y, z), ...] representing a flight one standard deviation higher than the mean
    '''
    # Convert flights into ndarrays
    # Non-dimensionalize time by their respective flight durations
    trajs = np.zeros([len(flights), 1000, 4])
    for i, flight in enumerate(flights):
        trajs[i, :, :] = traj(flight, postProcess=postProcess, n=n)
    
    # Calculate mean trajectory
    mean = sum(trajs) / len(trajs)

    # Calculate stddev at each time point
    stddev = np.sqrt(np.sum((trajs[:, :, 1 :] - mean[:, 1 :])**2, axis=0) / len(trajs))[:, 2]

    # Calculate upper trajectory
    upper = mean.copy()
    upper[:, 3] += stddev
    lower = mean.copy()
    lower[:, 3] -= stddev
    
    return mean, lower, upper

def axes3d():
    '''
    Get a set of axes with 3d projection

    Returns
        ax          axis with 3d projection
    '''
    plt.figure()
    return plt.subplot(projection='3d')

def plot_traj(traj, color, linestyle="-", elevation=None, projections=True, ax=None):
    '''
    Plot trajectory in 3d

    traj            ndarray([(t, x, y, z), ...])
    color           matplotlib color
    linestyle       style of line for plotting, default "-"
    elevation (m)   environment elevation
    projections     (unimplemented) if True, plot projections of trajectory onto principal planes
    ax              axis to plot on

    Returns
        ax          axis which was plotted to
    '''
    # Get axis if one is not provided
    if ax is None:
        plt.figure()
        ax = plt.subplot(projection='3d')
    
    # If environment elevation is not supplied, use minimum z value
    if elevation is None:
        elevation = traj[:, 3].min()
    
    # Does not currently work correctly
    # If plane projections are enabled, plot plane projections first
    # if projections:
    #     ax.plot(traj[:, 1], traj[:, 2], zs=0, zdir="z", linestyle="--", color=color)
    #     ax.plot(traj[:, 1], traj[:, 3] - elevation, zs=elevation, zdir="y", linestyle="--", color=color)
    #     ax.plot(traj[:, 2], traj[:, 2] - elevation, zs=elevation, zdir="x", linestyle="--", color=color)
    
    # Plot 3d trajectory curve
    ax.plot(traj[:, 1], traj[:, 2], traj[:, 3] - elevation, linestyle=linestyle, linewidth="2", color=color)

    # Plot origin
    ax.scatter(0, 0, 0, color="k")

    # Set labels
    ax.set_xlabel("X - East (m)")
    ax.set_ylabel("Y - North (m)")
    ax.set_zlabel("Z - Altitude Above Ground Level (m)")
    ax.set_title("Flight Trajectory")

    # Set view orientation
    ax.view_init(15, 45)
    
    return ax

def plot_flight(flight, color, postProcess=True, n=1000, projections=True, ax=None):
    '''
    Plot flight in 3d

    flight          Flight
    color           matplotlib color
    postProcess     if True, make sure the flight is postProcessed before proceeding
    n               number of points for reinterpolation
    projections     (unimplemented) if True, show projections of trajectory onto principal planes
    ax              axis onto which to plot, if None, an axis is initialized automatically

    Returns
        ax          the axis onto which the flight was plotted
    '''
    return plot_traj(traj(flight, postProcess, n), color, elevation=flight.env.elevation, projections=projections, ax=ax)

def plot_flights(flights, colors=None, elevations=None, postProcess=True, n=1000, projections=True, ax=None):
    '''
    Plot flights or trajectories automatically in 3d

    flights         iterable of Flights or trajectories
    colors          iterable of matplotlib colors, if None, the rainbow is used
    elevations      elevations for trajectories, if None, elevations are inferred
    postProcess     if True, make sure Flights are postProcessed before proceeding
    n               number of points for reinterpolation
    projections     (unimplemented) if True, plot projections of trajectories onto principal planes
    ax              axis onto which to plot, if None axis is intitialized automatically

    Returns
        ax          axis which was plotted to
    '''
    # If colors is not given, select colors automatically
    if colors is None:
        colors = cm.rainbow(np.linspace(0, 1, len(flights)))
    
    # If elevations are not given, do not give elevations
    if elevations is None:
        elevations = [None for _ in flights]
    
    # Plot each flight or trajectory
    for flight, color, elevation in zip(flights, colors, elevations):
        if isinstance(flight, np.ndarray):
            ax = plot_traj(flight, color, elevation=elevation, projections=projections, ax=ax)
        else:
            ax = plot_flight(flight, color, postProcess=postProcess, n=n, projections=projections, ax=ax)
    
    return ax

def plot_dispersion(flights, postProcess=True, n=1000, projections=True, ax=None):
    '''
    Plot the results of a dispersion analysis with mean trajectory

    flights         Flights for which trajectories will be plotted, len(flights) >= 1
    postProcess     if True, make sure flights have been postProcessed before proceeding
    n               number of reinterpolation points
    projections     (unimplemented) if True, plot projections of trajectories onto principal axes
    ax              axis onto which to plot

    Returns
        ax          axis which was plotted to
    '''
    mean, lower, upper = flight_stats(flights, postProcess=postProcess, n=n)
    ax = plot_flights(flights, ['k' for _ in range(len(flights))], postProcess=postProcess, n=n, projections=projections, ax=ax)
    ax = plot_traj(mean, 'r', linestyle="-", elevation=flights[0].env.elevation, projections=projections, ax=ax)
    ax = plot_traj(lower, 'r', linestyle="--", elevation=flights[0].env.elevation, projections=projections, ax=ax)
    ax = plot_traj(upper, 'r', linestyle="--", elevation=flights[0].env.elevation, projections=projections, ax=ax)
    return ax

def set_parabolic_trajectory(flight, t, x, y, zmax, n):
    '''
    Create a parabolic trajectory embedded in 3d space

    flight      Flight object to modify
    t           (t0, tf)
    x           (x0, xf)
    y           (y0, yf)
    zmax        apogee
    n           number of interpolation points

    Returns
        trajectory  [(t, x, y, z), ...]
    '''
    ts = np.linspace(*t, n)
    xs = np.linspace(*x, n)
    ys = np.linspace(*y, n)
    zs = -4 * zmax * (ts - t[0]) * (ts - t[1]) / (t[1] - t[0])**2

    flight.x = Function(np.stack([ts, xs], 1))
    flight.y = Function(np.stack([ts, ys], 1))
    flight.z = Function(np.stack([ts, zs], 1))

    return flight
