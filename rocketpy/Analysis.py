from .Flight import Flight
from .Function import Function
from .Rocket import Rocket

class Analysis:
    def __init__(self, flight: Flight):
        self.flight = flight
        self.rocket = self.flight.rocket
        self.motor = self.rocket.motor
        self.env = self.flight.env
        
    def apogee_by_mass(self):
        """
        Returns a RocketPy Function that, for a given Flight configuration,
        estimates the apogee as a function of the dry mass of the rocket.

            Returns:
                RocketPy Function that provides predicted apogee as a function of dry mass
        """
        og_mass = self.rocket.mass # immutable value creates different object unaffected by subsequent code
        # Create version of flight that has variable mass
        def apogee(mass):
            self.rocket.mass = mass

            test_flight = Flight(
                rocket=self.rocket,
                environment=self.env,
                inclination=self.flight.inclination,
                heading=self.flight.heading,
                terminateOnApogee=True,
            )

            return test_flight.apogee
        
        # restore original mass of rocket
        self.rocket.mass = og_mass
        return Function(apogee, inputs="Mass (kg)", outputs="Estimated Apogee (m)")
        
        
    def rail_exit_velocity_by_mass(self, wind_v=-5):
        self.env.setAtmosphericModel(type="CustomAtmosphere", wind_v=wind_v)
        def speed(mass):
            # self.rocket.mass = mass
            
            # ^^^ In theory this could and should replace the variable_rocket definition
            # but due to a bug could not be used; calisto drag curves are used below
            # for demonstration purposes.
            
            # Most likely, some internal process is causing the rail exit velocity to
            # be the same even for different masses. 
            variable_rocket = Rocket(
                motor=self.motor,
                radius=self.rocket.radius,
                mass=mass,
                inertiaI=self.rocket.inertiaI,
                inertiaZ=self.rocket.inertiaZ,
                distanceRocketNozzle=self.rocket.distanceRocketNozzle,
                distanceRocketPropellant=self.rocket.distanceRocketPropellant,
                powerOffDrag="../data/calisto/powerOffDragCurve.csv",
                powerOnDrag="../data/calisto/powerOnDragCurve.csv",
            )
            
            test_flight = Flight(
                rocket=self.rocket,
                environment=self.env,
                inclination=self.flight.inclination,
                heading=self.flight.heading,
                terminateOnApogee=True,
            )

            return test_flight.outOfRailVelocity
                
        return Function(speed, inputs="Mass (kg)", outputs="Out of Rail Speed (m/s)")
            
        
    def chute_radius_finder(self):
        from numpy import pi
        desiredterminal = float(input('Enter desired landing velocity in m/s '))
        mass = self.rocket.mass
        d = self.env.density(1)
    
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
    
        radius = (2 * area) / pi # using formula from fruitychutes
    
        print("Estimated required radius: {:.6f} meters".format(radius))
        
        return radius
    
    def snatchforce_calculator(self, main_CdS):
        """
        Calculates the shock force of parachute deployment based on the CdS
        and velocity while the parachute is deployed. Works only for main parachute,
        does not account for drogue.

        Args:
            main_CdS (number): Drag coefficient times reference area of main chute.

        Returns:
            float: Shock force of main chute
        """
        density = self.env.density
        return (1/2)* density * (self.flight.impactVelocity**2) * main_CdS
        
    
    def CdS_finder(self, chutetype,ventedchute,oradius,iradius,customCd):
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
