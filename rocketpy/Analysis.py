from .Flight import Flight
from .Function import Function
from .Rocket import Rocket

class Analysis:
    def __init__(self, flight: Flight):
        self.flight = flight
        self.rocket = self.flight.Rocket
        self.motor = self.rocket.Motor
        self.env = self.flight.Environment
        
    def apogee_by_mass(self):
        """
        Returns a RocketPy Function that, for a given Flight configuration,
        estimates the apogee as a function of the dry mass of the rocket.

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
                powerOffDrag = 0.5,
                powerOnDrag = 0.5
            )

            test_flight = Flight(
                rocket=variable_rocket,
                environment=self.env,
                inclination=self.flight.inclination,
                heading=self.flight.heading,
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
        
        return Function(apogee, inputs="Mass (kg)", outputs="Estimated Apogee (m)")