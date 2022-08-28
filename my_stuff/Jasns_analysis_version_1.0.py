from .Flight import Flight
from .Function import Function
from .Rocket import Rocket

class Analysis:
    def __init__(self, flight: Flight):
        self.flight = flight
        self.rocket = self.flight.Rocket
        self.motor = self.rocket.Motor
        self.env = self.flight.Environment
        
    def apogee_by_mass(self, flight: Flight):
        """
        Returns a RocketPy Function that, for a given Flight configuration,
        estimates the apogee as a function of the dry mass of the rocket.
            Parameters:
                flight (Flight):    Flight object to provide mass-apogee analysis
            Returns:
                RocketPy Function that provides predicted apogee as a function of dry mass
        """
        # Create version of flight that has ambigious mass
        def apogee(mass):
            variable_rocket = Rocket(
                motor = self.motor,
                radius = self.rocket.radius,
                mass = mass,
                inertiaI = self.rocket.inertiaI,
                inertiaZ = self.rocket.inertiaZ,
                distanceRocketNozzle = self.rocket.distanceRocketNozzle,
                distanceRocketPropellant = self.rocket.distanceRocketPropellant,
                powerOffDrag = 'idk', # FIX !!!!!!!!!!
                powerOnDrag = 'idk' # FIX !!!!!!!!!!
            )

            test_flight = Flight(
                rocket=variable_rocket,
                environment=self.env,
                inclination=flight.inclination,
                heading=flight.heading,
                terminateOnApogee=True,
            )

            return test_flight.apogee
        
        
        def chute_radius_finder():
            from numpy import pi
            desiredterminal = float(input('Enter desired landing velocity in m/s '))
            mass = self.rocket.mass # find a way to make sure this calls correctly
            self.calculateDensityProfile()
            d = self.density(1)
        
            parachute_type = input("Enter parachute type (eg: 'toroidal' or 'custom') ")
        
            if parachute_type == 'flat':
                dragcoeff = 0.8
            elif parachute_type == 'toroidal':
                dragcoeff = 2.2
            elif parachute_type == 'spherical':
                dragcoeff = 1.5
            elif parachute_type == 'custom':
                dragcoeff = float(input('Enter custom drag coefficient '))
            else:
                print('Not a valid chute type, enter custom drag coefficient?')
                dragcoeff = float(input())
        
        
            desiredCdS = (2*mass*9.81)/(d*(desiredterminal)**2)
        
            area = desiredCdS / dragcoeff
        
            radius = (area / pi)**(1/2)
        
            print("Estimated required radius: {:.6f} meters".format(radius))
            return(desiredterminal, desiredCdS)
        
        def snatchforce_calculator() :
            density = self.env.density
            CdS = 10.0
            if CdS == True:
                snatchforce = (1/2)* density * Flight.impactVelocity * CdS
            else:
                print ("you do not have a parachute defined!")
            return(snatchforce)
        
        #
        # CdS_finder that uses parameters
        #
        def CdS_finder_paramver(chutetype,ventedchute,oradius,iradius,customCd):
            from numpy import pi
            

            chutetypedict = {
                'flat' : [0.75,0.80],
                'conical' : [0.75,0.90],
                'biconical' : [0.75,0.92],
                'triconical' : [0.80,0.96],
                'polyconical' : [0.80,0.96],
                'extended skirt' : [0.78,0.87],
                'skirt' : [0.78,0.87],
                'hemisphere' : [0.62,0.77],
                'guide surface' : [0.28,0.42],
                'annular' : [0.85,0.95],
            }





            if (chutetype not in chutetypedict.keys()) and ventedchute.lower() == 'n':
                CdShigh = (oradius**2 * pi) * customCd
                CdSlow = (oradius**2 * pi) * customCd
            elif chutetype not in chutetypedict.keys() and ventedchute.lower() == 'y':
                CdShigh = ((oradius**2 * pi)-(iradius**2 * pi)) * customCd
                CdSlow = ((oradius**2 * pi)-(iradius**2 * pi)) * customCd
            elif ventedchute.lower() == 'n':
                CdShigh = (oradius**2 * pi) * chutetypedict[chutetype.lower()][1]
                CdSlow = (oradius**2 * pi) * chutetypedict[chutetype.lower()][0]
            elif ventedchute.lower() == 'y':
                CdShigh = ((oradius**2 * pi)-(iradius**2 * pi)) * chutetypedict[chutetype.lower()][1]
                CdSlow = ((oradius**2 * pi)-(iradius**2 * pi)) * chutetypedict[chutetype.lower()][0]
            
            meanCdS = (CdShigh+CdSlow)/2    
            
            print('With your chute configuradtion, the CdS will be about {:.6f} maximum'.format(CdShigh))
            print('and {:.6f} minimum'.format(CdSlow))
            print('Mean CdS = {:.6}'.format(meanCdS))
            return(CdShigh,CdSlow,meanCdS)
        
        return Function(apogee, inputs="Mass (kg)", outputs="Estimated Apogee (m)")