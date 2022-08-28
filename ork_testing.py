from rocketpy import Environment, SolidMotor, Rocket
import zipfile
import xml.etree.ElementTree as ET


ork_file_path = "data/lazarus/lazarus.ork"
ork_root_path = ork_file_path[0:-4] + '-uncompressed'

with zipfile.ZipFile(ork_file_path,"r") as compressed_ork_file:
    compressed_ork_file.extractall(ork_root_path)

with open(ork_root_path+"/rocket.ork", "r") as uncompressed_ork_file:
    ork_xml_data = uncompressed_ork_file.read()

openrocket = ET.fromstring(ork_xml_data)

print(openrocket.tag)


# environment conditions:
simulations_list =  openrocket[1].findall('simulation')
last_simulation = simulations_list[-1]
ork_railLength = float(last_simulation.find("conditions").find("launchrodlength").text)
ork_latitude = float(last_simulation.find("conditions").find("launchlatitude").text)
ork_longitude = float(last_simulation.find("conditions").find("launchlongitude").text)
ork_elevation= float(last_simulation.find("conditions").find("launchaltitude").text)


Env = Environment(
    railLength=ork_railLength, latitude=ork_latitude, longitude=ork_longitude, elevation=ork_elevation
)

# motor configuration
# the OpenRocket motor data is a serialized Java ArrayList. I attempted to write code to deserialize
# the arrayList but could not get it to work in time, so for now the motoro must remain hardcoded. For easier use,
# I am using a default RocketPy example engine. 
motor = SolidMotor(
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

Calisto = Rocket(
    motor=motor,
    radius=127 / 2000,
    mass=19.197 - 2.956,
    inertiaI=6.60,
    inertiaZ=0.0351,
    distanceRocketNozzle=-1.255,
    distanceRocketPropellant=-0.85704,
    powerOffDrag="data/calisto/powerOffDragCurve.csv",
    powerOnDrag="data/calisto/powerOnDragCurve.csv",
)

for subcomponent in openrocket[0].findall('subcomponents'):
    print("1" + subcomponent.tag)
    for stage in subcomponent.findall('stage'):
        print("2" + stage.tag)
        for stage_subcomponent in stage.findall('subcomponents'):
            print("3" + stage_subcomponent.tag)
            for geometry in stage_subcomponent:
                print("4" + geometry.tag)
                geometry_name = geometry.find("name").text
                print("5" + geometry_name)

