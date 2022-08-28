from rocketpy import Environment, SolidMotor, Rocket, utilities, Flight
import datetime

# Import ork
ork_file_path = 'data/lazarus/lazarus.ork'
ork_xml_obj = utilities.import_openrocket(ork_file_path)

# Create environment
simulation_number = 5
Env = utilities.ork_create_environment(ork_xml_obj, simulation_number)

tomorrow = datetime.date.today() + datetime.timedelta(days=1)
Env.setDate((tomorrow.year, tomorrow.month, tomorrow.day, 12))  # Hour given in UTC time

# Create Motor
# I didn't have time to add motor functionality since the OpenRocket
# motor files are in a serialized Java ArrayList. I attempted to deserialize
# the motor data and failed by the time PRs were due. I have included a deserializer.py
# file but it does not work. 
Pro75M1670 = SolidMotor(
    thrustSource="data/motors/Cesaroni_M1670.eng",
    burnOut=3.9,
    grainNumber=5,
    grainSeparation=5 / 1000,
    grainDensity=1815,
    grainOuterRadius=33 / 1000,
    grainInitialInnerRadius=15 / 1000,
    grainInitialHeight=120 / 1000,
    nozzleRadius=33 / 1000,
    throatRadius=11 / 1000,
    interpolationMethod="linear",
)

# Create rocket
# since intertia and drag curves cannot be extracted from OpenRocket,
# it is best for users to manually define the Rocket Obj. 
# Drag curves can be made from OpenRocket but not extracted from
# an .ork file
ork_rocket = Rocket(
    motor=Pro75M1670,
    radius=127 / 2000,
    mass=19.197 - 2.956,
    inertiaI=6.60,
    inertiaZ=0.0351,
    distanceRocketNozzle=-1.255,
    distanceRocketPropellant=-0.85704,
    powerOffDrag="data/calisto/powerOffDragCurve.csv",
    powerOnDrag="data/calisto/powerOnDragCurve.csv",
)
ork_rocket.setRailButtons([0.2, -0.5])

# Add Aero Surfaces
# The aero surfaces' distance to CM cannot be extracted from .ork files
# and needs to be manually set by the user
nose_distanceToCM = 0.71971
fin_distanceToCM = -1.04956
ork_rocket = utilities.ork_add_aero_surfaces(ork_rocket, ork_xml_obj, nose_distanceToCM, fin_distanceToCM)

# openRocket does not have a standard method for defining tail cones
# so users should define their tailcone separately
Tail = ork_rocket.addTail(
    topRadius=0.0635, bottomRadius=0.0435, length=0.060, distanceToCM=-1.194656
)

# Add parachutes
ork_rocket = utilities.ork_add_parachutes(ork_rocket, ork_xml_obj, Env)

#Test Flight
TestFlight = Flight(rocket=ork_rocket, environment=Env, inclination=85, heading=0)
TestFlight.allInfo()
