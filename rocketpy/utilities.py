# -*- coding: utf-8 -*-
__author__ = "Franz Masatoshi Yuri, Lucas Kierulff Balabram, Guilherme Fernandes Alves"
__copyright__ = "Copyright 20XX, RocketPy Team"
__license__ = "MIT"

from cgi import test
import numpy as np
from scipy.integrate import solve_ivp
import zipfile
import xml.etree.ElementTree as ET
from math import pi

from .Environment import Environment
from .Function import Function


def import_openrocket(
    ork_file_path
):
    """Decompresses an OpenRocket .ork file and saves the uncompressed file for future use.
    
    Parameters
    ----------
    ork_file_path : string
        The directory location of the .ork file, including the file extension.
    
    Returns
    -------
    ork_xml_obj : xml.etree.ElementTree.Element
        An cml Element Tree object of the ork data. 

    """


    # remove .ork extension and add -uncompressed to create new directory name
    ork_root_path = ork_file_path[0:-4] + '-uncompressed' 

    with zipfile.ZipFile(ork_file_path, "r") as compressed_ork_file:
        compressed_ork_file.extractall(ork_root_path)

    with open(ork_root_path+"/rocket.ork", "r") as uncompressed_ork_file:
        ork_xml_data = uncompressed_ork_file.read()

    ork_xml_obj = ET.fromstring(ork_xml_data)

    return ork_xml_obj


def ork_create_environment(
    ork_xml_obj, simulation_number
):
    """Creates a RocketPy Environment class object based on parameters from a given OpenRocket file.
    
    Parameters
    ----------
    ork_xml_obj : xml.etree.ElementTree.Element
        An cml Element Tree object of the ork data. 
    simulation_number : int
        The number simulation in the OpenRocket file that the users wishes to extract environment conditions from. Can be 1-N.  
        
    Returns
    -------
    ork_environment : Environment
        An Environment object containing all openRocket environment parameters
    
    """


    # parse environment parameters
    simulations_list = ork_xml_obj[1].findall('simulation')
    if simulation_number <= 0 or simulation_number > len(simulations_list):
        er_msg = 'simulation_number of {} is not valid. Simulation #{} is not present in the OpenRocket file.'.format(simulation_number)
        raise ValueError(er_msg)

    parameter_simulation = simulations_list[simulation_number]
    ork_railLength =    float(parameter_simulation.find("conditions").find("launchrodlength").text)
    ork_latitude =      float(parameter_simulation.find("conditions").find("launchlatitude").text)
    ork_longitude =     float(parameter_simulation.find("conditions").find("launchlongitude").text)
    ork_elevation=      float(parameter_simulation.find("conditions").find("launchaltitude").text)

    ork_environment = Environment(
        railLength=ork_railLength, latitude=ork_latitude, longitude=ork_longitude, elevation=ork_elevation
    )

    return ork_environment


def ork_add_aero_surfaces(
    rocket, ork_xml_obj, nose_distanceToCM, fin_distanceToCM
):
    """Returns a Rocket object with the added aerodynamic surfaces nosecone and finset. Note: This does not 
    currently support adding tail aero surfaces because OpenRocket does not have a standard tailcone parameter. 
    Additionally, this function requires the distanceToCM values for the nosecone and fins because this is also
    not provided by OpenRocket. 
    
    Parameters
    ----------
    rocket : Rocket
        A RocketPy Rocket object.
    ork_xml_obj : xml.etree.ElementTree.Element
        An cml Element Tree object of the ork data. 
    nose_distanceToCM : int
        Distance from nosecone base to center of mass.
    fin_distanceToCM : int
        Distance from fins to center of mass.

    Returns
    -------
    rocket : Rocket
        A RocketPy Rocket object with nosecone and finset surfaces added.

    """

    ork_nose_length = 0
    ork_nose_kind = ""
    ork_num_fins = 0
    ork_fin_span = 0
    ork_rootChord = 0
    ork_tipChord = 0

    for subcomponent in ork_xml_obj[0].findall('subcomponents'):
        for stage in subcomponent.findall('stage'):
            for stage_subcomponent in stage.findall('subcomponents'):
                for geometry in stage_subcomponent:
                    # Add nosecone parameters
                    if geometry.tag == "nosecone":
                        ork_nose_length = float(geometry.find("length").text)
                        ork_nose_kind = geometry.find("shape").text
                        if ork_nose_kind == "haack":
                            ork_nose_kind = "lvhaack"
                            # TODO, double check from OpenRocket that conical and ogive don't need to be renamed as well
                        else:
                            raise ValueError('Nosecone type of "{}" is not supported'.format(ork_nose_kind))
                    
                    # Add fin geometries
                    # Note: does not support eliptical or freeform fins
                    elif geometry.tag == "bodytube":
                        for geom_subcomponent in geometry.findall('subcomponents'):
                            for finset in geom_subcomponent.findall('trapezoidfinset'):
                                ork_num_fins = float(finset.find("fincount").text)
                                ork_fin_span = float(finset.find("height").text)
                                ork_rootChord = float(finset.find("rootchord").text)
                                ork_tipChord = float(finset.find("tipchord").text)

    NoseCone = rocket.addNose(length=ork_nose_length, kind=ork_nose_kind, distanceToCM=nose_distanceToCM)

    FinSet = rocket.addFins(
        ork_num_fins, span=ork_fin_span, rootChord=ork_rootChord, tipChord=ork_tipChord, distanceToCM=fin_distanceToCM
    )

    return rocket


def ork_add_parachutes(
    rocket, ork_xml_obj, environment
):
    """Returns a Rocket object with the added parachutes.
    
    Parameters
    ----------
    rocket : Rocket
        A RocketPy Rocket object.
    ork_xml_obj : xml.etree.ElementTree.Element
        An cml Element Tree object of the ork data. 
    environment : Environment
        A RocketPy Environment object. 

    Returns
    -------
    rocket : Rocket
        A RocketPy Rocket object with parachutes.

    """

    # TODO find workaround for OpenRocket "auto" value for cd
    for subcomponent in ork_xml_obj[0].findall('subcomponents'):
        for stage in subcomponent.findall('stage'):
            for stage_subcomponent in stage.findall('subcomponents'):
                for bodytube in stage_subcomponent.findall('bodytube'):
                    for body_subcomponent in bodytube.findall('subcomponents'):
                        for parachute in body_subcomponent.findall('parachute'):
                            deployevent_type = parachute.find("deployevent").text
                            if deployevent_type == "apogee":
                                ork_drogue_cd = parachute.find("cd").text
                                if ork_drogue_cd == "auto":
                                    ork_drogue_cd = 1.2 
                                else:
                                    ork_drogue_cd = float(ork_drogue_cd)
                                ork_drogue_diameter = float(parachute.find("diameter").text)
                                ork_drogue_cds = ork_drogue_cd * pi * (ork_drogue_diameter/2)**2 
                            elif deployevent_type == "altitude":
                                ork_main_cd = parachute.find("cd").text
                                if ork_main_cd == "auto":
                                    ork_main_cd = 1.2 
                                else:
                                    ork_main_cd = float(ork_main_cd)
                                ork_main_diameter = float(parachute.find("diameter").text)
                                ork_main_cds = ork_main_cd * pi * (ork_main_diameter/2)**2 
                                ork_main_altitude = float(parachute.find("deployaltitude").text)
                            else:
                                raise ValueError('Unknown parachute deployment type of {}'.format(deployevent_type))


    def drogueTrigger(p, y):
        # p = pressure
        # y = [x, y, z, vx, vy, vz, e0, e1, e2, e3, w1, w2, w3]
        # activate drogue when vz < 0 m/s.
        return True if y[5] < 0 else False      

    # TODO find a way for the main trigger to be set off by ork_main_altitude + environment.elevation instead of constants
    # this may require rewriting how the trigger is used during flight simulation
    def mainTrigger(p, y):
        # p = pressure
        # y = [x, y, z, vx, vy, vz, e0, e1, e2, e3, w1, w2, w3]
        # activate main when vz < 0 m/s and z < 800 + 1400 m (+1400 due to surface elevation).
        return True if y[5] < 0 and y[2] < 800 + 1400 else False         


    Main = rocket.addParachute(
        "Main",
        CdS=ork_main_cds,
        trigger=mainTrigger,
    )

    Drogue = rocket.addParachute(
        "Drogue",
        CdS=ork_drogue_cds,
        trigger=drogueTrigger,
    )

    return rocket


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
