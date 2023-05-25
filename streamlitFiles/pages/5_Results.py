import streamlit as st
from rocketpy import Flight
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import shutil
import os

st.set_page_config(
    page_title="Results",
    page_icon="🖥️",
    layout="wide",
    initial_sidebar_state="expanded",
    # menu_items={
    #     "Get Help": "https://www.extremelycoolapp.com/help",
    #     "Report a bug": "https://www.extremelycoolapp.com/bug",
    #     "About": "# This is a header. This is an *extremely* cool app!",
    # },
)
if "rocket" not in st.session_state:
    st.session_state.rocket = None

if "Env" not in st.session_state:
    st.warning("Go back to environment page and initialize environment")

if "motor" not in st.session_state:
    st.warning("Go back to motor page and initialize motor")

if "rocket" not in st.session_state:
    st.warning("Go back to rocket page and initialize rocket")

if "flight" not in st.session_state:
    st.session_state.flight = {
        "inclination": 80.0,
        "heading": 90.0,
        "simulateFight": st.empty,
        "exportKML": False,
    }
col1, col2 = st.columns(2)

with col1:
    st.subheader("Current Configuration")
    st.write(
        "Motor Loaded: **{}**".format(st.session_state.motorInfo["thrustSource"][12:-4])
    )
    st.write("**{}** added".format(st.session_state.rocketInfo["rocketName"]))
    st.write("Radius: **{:.2f} m**".format(st.session_state.rocketInfo["radius"]))
    st.write("Dry Mass: **{:.2f} kg**".format(st.session_state.rocketInfo["mass"]))
    st.write(
        "Parachute: **{}** deployed at: **Apogee** | Delay: **{:.2f} seconds**".format(
            st.session_state.rocket.parachutes[0].name,
            st.session_state.rocketInfo["drogueLag"],
        )
    )
    # st.session_state.rocketInfo[""])
    # st.write(
    #     "Parachute: ",
    #     st.session_state.rocket.parachutes[1].name,
    #     " deployed at: ",
    #     st.session_state.rocketInfo["deploymentHeight"],
    #     " m AGL ",
    #     " Delay: ",
    #     st.session_state.rocketInfo["drogueLag"],
    #     " seconds",
    # )
    st.write(
        "Parachute: **{}** deployed at: **{:.1f} m** AGL | Delay: **{:.2f} seconds**".format(
            st.session_state.rocket.parachutes[1].name,
            st.session_state.rocketInfo["deploymentHeight"],
            st.session_state.rocketInfo["drogueLag"],
        )
    )
    # st.session_state.rocketInfo[""])


with col2:
    with st.form("flightSimulation"):
        st.subheader("Flight Simulation")
        st.session_state.flight["inclination"] = st.number_input(
            "Rail inclination angle relative to ground (°)",
            value=st.session_state.flight["inclination"],
            format="%f",
            help="Rail inclination angle relative to ground, given in degrees. Default is 80.",
        )
        st.session_state.flight["heading"] = st.number_input(
            "Heading angle relative to north (°)",
            value=st.session_state.flight["heading"],
            format="%f",
            help=" Heading angle relative to north given in degrees. Default is 90, which points in the x direction.",
        )
        # add min and max time step, note max time step should be bigger than min
        st.session_state.flight["exportKML"] = st.checkbox(
            "Export KML?", value=st.session_state.flight["exportKML"]
        )

        st.session_state.flight["simulateFight"] = st.form_submit_button("Simulate")
        if st.session_state.flight["simulateFight"]:
            testFlight = Flight(
                rocket=st.session_state.rocket,
                environment=st.session_state.Env,
                inclination=st.session_state.flight["inclination"],
                heading=st.session_state.flight["heading"],
            )
            testFlight.postProcess()

if st.session_state.flight["simulateFight"]:
    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("Initial Conditions")
        st.write(
            "Position - Vx: **{:.2f} m/s** | Vy: **{:.2f} m/s** | Vz: **{:.2f} m/s**".format(
                testFlight.initialSolution[1],
                testFlight.initialSolution[2],
                testFlight.initialSolution[3],
            )
        )
        st.write(
            "Velocity - Vx: **{:.2f} m/s** | Vy: **{:.2f} m/s** | Vz: **{:.2f} m/s**".format(
                testFlight.initialSolution[4],
                testFlight.initialSolution[5],
                testFlight.initialSolution[6],
            )
        )
        st.write(
            "Attitude - e0: **{:.3f}** | e1: **{:.3f}** | e2: **{:.3f}** | e3: **{:.3f}**".format(
                testFlight.initialSolution[7],
                testFlight.initialSolution[8],
                testFlight.initialSolution[9],
                testFlight.initialSolution[10],
            )
        )
        st.write(
            "Euler Angles - Spin φ : **{:.2f}°** | Nutation θ: **{:.2f}°** | Precession ψ: **{:.2f}°**".format(
                testFlight.phi(0),
                testFlight.theta(0),
                testFlight.psi(0),
            )
        )
        st.write(
            "Angular Velocity - ω1: **{:.2f} rad/s**  | ω2: **{:.2f} rad/s** | ω3: **{:.2f} rad/s** ".format(
                testFlight.initialSolution[11],
                testFlight.initialSolution[12],
                testFlight.initialSolution[13],
            )
        )

    with col2:
        st.subheader("Launch rail Orientation")
        st.write("Launch Rail Inclination: **{:.2f}°**".format(testFlight.inclination))
        st.write("Launch Rail Heading: **{:.2f}°**".format(testFlight.heading))

    with col3:
        st.subheader("Surface Wind Conditions")
        st.write(
            "Frontal Surface Wind Speed: **{:.2f} m/s**".format(
                testFlight.frontalSurfaceWind
            )
        )
        st.write(
            "Lateral Surface Wind Speed: **{:.2f} m/s**".format(
                testFlight.lateralSurfaceWind
            )
        )

    col1, col2, col3 = st.columns(3)
    with col1:
        st.subheader("Rail Departure State")
        st.write("Rail Departure Time: **{:.3f} s**".format(testFlight.outOfRailTime))
        st.write(
            "Rail Departure Velocity: **{:.3f} m/s**".format(
                testFlight.outOfRailVelocity
            )
        )
        st.write(
            "Rail Departure Static Margin: **{:.3f} c**".format(
                testFlight.staticMargin(testFlight.outOfRailTime)
            )
        )
        st.write(
            "Rail Departure Angle of Attack: **{:.3f}°**".format(
                testFlight.angleOfAttack(testFlight.outOfRailTime)
            )
        )
        st.write(
            "Rail Departure Thrust-Weight Ratio: **{:.3f}**".format(
                st.session_state.rocket.thrustToWeight(testFlight.outOfRailTime)
            )
        )
        st.write(
            "Rail Departure Reynolds Number: **{:.3e}**".format(
                testFlight.ReynoldsNumber(testFlight.outOfRailTime)
            )
        )

    with col2:
        st.subheader("BurnOut State")
        st.write(
            "BurnOut time: **{:.3f} s**".format(
                st.session_state.rocket.motor.burnOutTime,
            )
        )
        st.write(
            "Altitude at burnout: **{:.3f} m** AGL".format(
                testFlight.z(st.session_state.rocket.motor.burnOutTime),
            )
        )
        st.write(
            "Rocket velocity at burnOut: **{:.3f} m/s**".format(
                testFlight.speed(st.session_state.rocket.motor.burnOutTime),
            )
        )
        st.write(
            "Freestream velocity at burnOut: **{:.3f} m/s**".format(
                (
                    testFlight.streamVelocityX(
                        st.session_state.rocket.motor.burnOutTime
                    )
                    ** 2
                    + testFlight.streamVelocityY(
                        st.session_state.rocket.motor.burnOutTime
                    )
                    ** 2
                    + testFlight.streamVelocityZ(
                        st.session_state.rocket.motor.burnOutTime
                    )
                    ** 2
                )
                ** 0.5,
            )
        )
        st.write(
            "Mach Number at burnOut: **{:.3f}**".format(
                testFlight.MachNumber(st.session_state.rocket.motor.burnOutTime),
            )
        )
        st.write(
            "Kinetic energy at burnOut: **{:.3e} J**".format(
                testFlight.kineticEnergy(st.session_state.rocket.motor.burnOutTime),
            )
        )

    with col3:
        st.subheader("Apogee")
        st.write(
            "Apogee Altitude: **{:.3f} m** (ASL) | **{:.3f} m** (AGL)".format(
                testFlight.apogee,
                testFlight.apogee - st.session_state.Env.elevation,
            )
        )
        st.write(
            "Apogee Time: **{:.3f} s**".format(
                testFlight.apogeeTime,
            )
        )
        st.write(
            "Apogee Freestream Speed: **{:.3f} m/s**".format(
                testFlight.apogeeFreestreamSpeed,
            )
        )

    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("Events")
        if len(testFlight.parachuteEvents) == 0:
            st.write("No Parachute Events Were Triggered.")
        for event in testFlight.parachuteEvents:
            triggerTime = event[0]
            parachute = event[1]
            openTime = triggerTime + parachute.lag
            velocity = testFlight.freestreamSpeed(openTime)
            altitude = testFlight.z(openTime)
            name = parachute.name.title()
            st.write(name + " Ejection Triggered at: **{:.3f} s**".format(triggerTime))
            st.write(name + " Parachute Inflated at: **{:.3f} s**".format(openTime))
            st.write(
                name
                + " Parachute Inflated with Freestream Speed of: **{:.3f} m/s**".format(
                    velocity
                )
            )
            st.write(
                name
                + " Parachute Inflated at Height of: **{:.3f} m** (AGL)".format(
                    altitude - testFlight.env.elevation
                )
            )

    with col2:
        st.subheader("Impact")
        if len(testFlight.impactState) != 0:
            st.write("X Impact: **{:.3f} m**".format(testFlight.xImpact))
            st.write("Y Impact: **{:.3f} m**".format(testFlight.yImpact))
            st.write("Time of Impact: **{:.3f} s**".format(testFlight.tFinal))
            st.write(
                "Velocity at Impact: **{:.3f} m/s**".format(testFlight.impactVelocity)
            )
        elif testFlight.terminateOnApogee is False:
            st.write("End of Simulation")
            st.write("Time: **{:.3f} s**".format(testFlight.solution[-1][0]))
            st.write("Altitude: **{:.3f} m**".format(testFlight.solution[-1][3]))

    with col3:
        st.subheader("Maximum Values")
        st.write(
            "Maximum Speed: **{:.3f} m/s** at ****{:.2f} s****".format(
                testFlight.maxSpeed, testFlight.maxSpeedTime
            )
        )
        st.write(
            "Maximum Mach Number: **{:.3f} Mach** at **{:.2f} s**".format(
                testFlight.maxMachNumber, testFlight.maxMachNumberTime
            )
        )
        st.write(
            "Maximum Reynolds Number: **{:.3e}** at **{:.2f} s**".format(
                testFlight.maxReynoldsNumber, testFlight.maxReynoldsNumberTime
            )
        )
        st.write(
            "Maximum Dynamic Pressure: **{:.3e} Pa** at **{:.2f} s**".format(
                testFlight.maxDynamicPressure, testFlight.maxDynamicPressureTime
            )
        )
        st.write(
            "Maximum Acceleration: **{:.3f} m/s²** at **{:.2f} s**".format(
                testFlight.maxAcceleration, testFlight.maxAccelerationTime
            )
        )
        st.write(
            "Maximum Gs: **{:.3f} g** at **{:.2f} s**".format(
                testFlight.maxAcceleration / testFlight.env.g,
                testFlight.maxAccelerationTime,
            )
        )
        st.write(
            "Maximum Upper Rail Button Normal Force: **{:.3f} N**".format(
                testFlight.maxRailButton1NormalForce
            )
        )
        st.write(
            "Maximum Upper Rail Button Shear Force: **{:.3f} N**".format(
                testFlight.maxRailButton1ShearForce
            )
        )
        st.write(
            "Maximum Lower Rail Button Normal Force: **{:.3f} N**".format(
                testFlight.maxRailButton2NormalForce
            )
        )
        st.write(
            "Maximum Lower Rail Button Shear Force: **{:.3f} N**".format(
                testFlight.maxRailButton2ShearForce
            )
        )

    st.subheader("Numerical Integration Information")
    col1, col2 = st.columns(2)
    with col1:
        st.write("Maximum Allowed Flight Time: **{:f} s**".format(testFlight.maxTime))
        st.write("Maximum Allowed Time Step: **{:f} s**".format(testFlight.maxTimeStep))
        st.write("Minimum Allowed Time Step: **{:e} s**".format(testFlight.minTimeStep))
        st.write("Relative Error Tolerance: ", testFlight.rtol)
        # st.write("Absolute Error Tolerance: ", testFlight.atol)
        st.write("Allow Event Overshoot: ", testFlight.timeOvershoot)

    with col2:
        st.write("Terminate Simulation on Apogee: ", testFlight.terminateOnApogee)
        st.write("Number of Time Steps Used: ", len(testFlight.timeSteps))
        st.write(
            "Number of Derivative Functions Evaluation: ",
            sum(testFlight.functionEvaluationsPerTimeStep),
        )
        st.write(
            "Average Function Evaluations per Time Step: **{:3f}**".format(
                sum(testFlight.functionEvaluationsPerTimeStep)
                / len(testFlight.timeSteps)
            )
        )
    st.write(
        """
        ---
    """
    )
    st.subheader("Trajectory 3d Plot")
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        # Get max and min x and y
        maxZ = max(testFlight.z[:, 1] - testFlight.env.elevation)
        maxX = max(testFlight.x[:, 1])
        minX = min(testFlight.x[:, 1])
        maxY = max(testFlight.y[:, 1])
        minY = min(testFlight.y[:, 1])
        maxXY = max(maxX, maxY)
        minXY = min(minX, minY)

        # Create figure
        fig1 = plt.figure(figsize=(9, 9))
        ax1 = plt.subplot(111, projection="3d")
        ax1.plot(testFlight.x[:, 1], testFlight.y[:, 1], zs=0, zdir="z", linestyle="--")
        ax1.plot(
            testFlight.x[:, 1],
            testFlight.z[:, 1] - testFlight.env.elevation,
            zs=minXY,
            zdir="y",
            linestyle="--",
        )
        ax1.plot(
            testFlight.y[:, 1],
            testFlight.z[:, 1] - testFlight.env.elevation,
            zs=minXY,
            zdir="x",
            linestyle="--",
        )
        ax1.plot(
            testFlight.x[:, 1],
            testFlight.y[:, 1],
            testFlight.z[:, 1] - testFlight.env.elevation,
            linewidth="2",
        )
        ax1.scatter(0, 0, 0)
        ax1.set_xlabel("X - East (m)")
        ax1.set_ylabel("Y - North (m)")
        ax1.set_zlabel("Z - Altitude Above Ground Level (m)")
        ax1.set_title("Flight Trajectory")
        ax1.set_zlim3d([0, maxZ])
        ax1.set_ylim3d([minXY, maxXY])
        ax1.set_xlim3d([minXY, maxXY])
        ax1.view_init(15, 45)
        st.pyplot(fig1)

    st.write(
        """
        ---
    """
    )
    st.subheader("Trajectory Kinematic Plots")
    col1, col2 = st.columns(2)
    with col1:

        st.write("Velocity Magnitude | Acceleration Magnitude  vs Time (s)")
        vaMagDf = pd.DataFrame(
            data=zip(testFlight.speed[:, 1], testFlight.acceleration[:, 1]),
            columns=["Velocity (m/s)", "Acceleration (m/s²)"],
            index=testFlight.acceleration[:, 0],
        )
        st.line_chart(vaMagDf)

        st.write("Velocity Y | Acceleration Y  vs Time (s)")
        vaYDf = pd.DataFrame(
            data=zip(testFlight.vy[:, 1], testFlight.ay[:, 1]),
            columns=["Velocity Y (m/s)", "Acceleration Y (m/s²)"],
            index=testFlight.ay[:, 0],
        )
        st.line_chart(vaYDf)

    with col2:
        st.write("Velocity Z | Acceleration Z  vs Time (s)")
        vaZDf = pd.DataFrame(
            data=zip(testFlight.vz[:, 1], testFlight.az[:, 1]),
            columns=["Velocity Z (m/s)", "Acceleration Z (m/s²)"],
            index=testFlight.az[:, 0],
        )
        st.line_chart(vaZDf)

        st.write("Velocity X | Acceleration X  vs Time (s)")
        vaXDf = pd.DataFrame(
            data=zip(testFlight.vx[:, 1], testFlight.ax[:, 1]),
            columns=["Velocity X (m/s)", "Acceleration X (m/s²)"],
            index=testFlight.ax[:, 0],
        )
        st.line_chart(vaXDf)
    st.write(
        """
        ---
    """
    )
    st.subheader("Angular Position Plots")
    col1, col2 = st.columns(2)
    with col1:
        st.write("Flight Path and Attitude Angle vs Time (s)")
        fpraDf = pd.DataFrame(
            data=zip(testFlight.pathAngle[:, 1], testFlight.attitudeAngle[:, 1]),
            columns=["Flight Path Angle (°)", "Rocket Attitude Angle (°)"],
            index=testFlight.attitudeAngle[:, 0],
        )
        st.line_chart(fpraDf)

    with col2:
        st.write("Lateral Attitude Angle vs Time (s)")
        laaDf = pd.DataFrame(
            data=zip(testFlight.lateralAttitudeAngle[:, 1]),
            columns=["Lateral Attitude Angle (°)"],
            index=testFlight.lateralAttitudeAngle[:, 0],
        )
        st.line_chart(laaDf)

    st.write(
        """
        ---
    """
    )
    st.subheader("Path, Attitude and Lateral Attitude Angle plots")
    col1, col2 = st.columns(2)
    with col1:
        eulerDf = pd.DataFrame(
            data=zip(
                testFlight.e0[:, 1],
                testFlight.e1[:, 1],
                testFlight.e2[:, 1],
                testFlight.e3[:, 1],
            ),
            columns=["e₀", "e₁", "e₂", "e₃"],
            index=testFlight.e0[:, 0],
        )
        st.write("Euler Parameters vs Time (s)")
        st.line_chart(eulerDf)

        eulerPDf = pd.DataFrame(
            data=zip(testFlight.psi[:, 1]),
            columns=["ψ (°)"],
            index=testFlight.psi[:, 0],
        )
        st.write("Euler Precession Angle (°) vs Time (s)")
        st.line_chart(eulerPDf)

    with col2:
        eulerNDf = pd.DataFrame(
            data=zip(testFlight.theta[:, 1]),
            columns=["θ (°)"],
            index=testFlight.theta[:, 0],
        )
        st.write("Euler Nutation Angle (°) vs Time (s)")
        st.line_chart(eulerNDf)

        eulerSDf = pd.DataFrame(
            data=zip(testFlight.phi[:, 1]),
            columns=["φ (°)"],
            index=testFlight.phi[:, 0],
        )
        st.write("Euler Spin Angle (°) vs Time (s)")
        st.line_chart(eulerSDf)

    st.write(
        """
        ---
    """
    )
    st.subheader("Trajectory Angular Velocity and Acceleration Plots")
    col1, col2, col3 = st.columns(3)
    with col1:
        angularVADf = pd.DataFrame(
            data=zip(
                testFlight.w1[:, 1],
                testFlight.alpha1[:, 1],
            ),
            columns=[
                "Angular Velocity - ω₁ (rad/s)",
                "Angular Acceleration - α₁ (rad/s²)",
            ],
            index=testFlight.alpha1[:, 0],
        )
        st.write("Angular Velocity - ω₁ | Angular Acceleration - α₁ vs Time (s)")
        st.line_chart(angularVADf)

    with col2:
        angularVA2Df = pd.DataFrame(
            data=zip(
                testFlight.w2[:, 1],
                testFlight.alpha2[:, 1],
            ),
            columns=[
                "Angular Velocity - ω₂ (rad/s)",
                "Angular Acceleration - α₂ (rad/s²)",
            ],
            index=testFlight.alpha2[:, 0],
        )
        st.write("Angular Velocity - ω₂ | Angular Acceleration - α₂ vs Time (s)")
        st.line_chart(angularVA2Df)

    with col3:
        angularVA3Df = pd.DataFrame(
            data=zip(
                testFlight.w3[:, 1],
                testFlight.alpha3[:, 1],
            ),
            columns=[
                "Angular Velocity - ω₃ (rad/s)",
                "Angular Acceleration - α₃ (rad/s²)",
            ],
            index=testFlight.alpha3[:, 0],
        )
        st.write("Angular Velocity - ω₃ | Angular Acceleration - α₃ vs Time (s)")
        st.line_chart(angularVA3Df)

    # st.write(
    #     """
    #     ---
    # """
    # )
    # st.subheader("Trajectory Force Plots")
    # col1, col2 = st.columns(2)
    # with col1:
    #     railButtonNFDf = pd.DataFrame(
    #         data=zip(
    #             testFlight.railButton1NormalForce[:outOfRailTimeIndex, 1],
    #             testFlight.railButton2NormalForce[:outOfRailTimeIndex, 1],
    #         ),
    #         columns=[
    #             "Upper Rail Button",
    #             "Lower Rail Button",
    #         ],
    #         index=testFlight.railButton1NormalForce[:outOfRailTimeIndex, 0],
    #     )
    #     st.write("Rail Buttons Normal Force vs Time (s)")
    #     st.line_chart(railButtonNFDf)

    # with col2:
    #     railButtonSFDf = pd.DataFrame(
    #         data=zip(
    #             testFlight.railButton1ShearForce[: testFlight.outOfRailTimeIndex, 1],
    #             testFlight.railButton2ShearForce[: testFlight.outOfRailTimeIndex, 1],
    #         ),
    #         columns=[
    #             "Upper Rail Button",
    #             "Lower Rail Button",
    #         ],
    #         index=testFlight.railButton1ShearForce[: testFlight.outOfRailTimeIndex, 0],
    #     )
    #     st.write("Rail Buttons Normal Force vs Time (s)")
    #     st.line_chart(railButtonSFDf)
    st.write(
        """
        ---
    """
    )
    st.subheader("Aerodynamic force and moment plots")
    if len(testFlight.parachuteEvents) > 0:
        eventTime = (
            testFlight.parachuteEvents[0][0] + testFlight.parachuteEvents[0][1].lag
        )
        eventTimeIndex = np.nonzero(testFlight.x[:, 0] == eventTime)[0][0]
    else:
        eventTime = testFlight.tFinal
        eventTimeIndex = -1
    col1, col2 = st.columns(2)
    with col1:

        aeroLiftResDf = pd.DataFrame(
            data=zip(
                testFlight.aerodynamicLift[:eventTimeIndex, 1],
                testFlight.R1[:eventTimeIndex, 1],
                testFlight.R2[:eventTimeIndex, 1],
            ),
            columns=[
                "Resultant",
                "R1",
                "R2",
            ],
            index=testFlight.aerodynamicLift[:eventTimeIndex, 0],
        )
        st.write("Aerodynamic Lift Resultant Force (N) vs Time (s)")
        st.line_chart(aeroLiftResDf)

        aeroDragResDf = pd.DataFrame(
            data=zip(
                testFlight.aerodynamicDrag[:eventTimeIndex, 1],
            ),
            columns=[
                "Drag Force (N)",
            ],
            index=testFlight.aerodynamicDrag[:eventTimeIndex, 0],
        )
        st.write("Aerodynamic Drag Force (N) vs Time (s)")
        st.line_chart(aeroDragResDf)

    with col2:
        aeroBendResDf = pd.DataFrame(
            data=zip(
                testFlight.aerodynamicBendingMoment[:eventTimeIndex, 1],
                testFlight.M1[:eventTimeIndex, 1],
                testFlight.M2[:eventTimeIndex, 1],
            ),
            columns=[
                "Resultant",
                "M1",
                "M2",
            ],
            index=testFlight.aerodynamicBendingMoment[:eventTimeIndex, 0],
        )
        st.write("Aerodynamic Bending Resultant Moment (Nm) vs Time (s)")
        st.line_chart(aeroBendResDf)

        aeroSpinResDf = pd.DataFrame(
            data=zip(
                testFlight.aerodynamicSpinMoment[:eventTimeIndex, 1],
            ),
            columns=[
                "Spin Moment (Nm)",
            ],
            index=testFlight.aerodynamicSpinMoment[:eventTimeIndex, 0],
        )
        st.write("Aerodynamic Spin Moment (Nm) vs Time (s)")
        st.line_chart(aeroSpinResDf)

    st.write(
        """
        ---
    """
    )
    st.subheader("Trajectory Energy Plots")
    # if len(testFlight.parachuteEvents) > 0:
    #     eventTime = (
    #         testFlight.parachuteEvents[0][0] + testFlight.parachuteEvents[0][1].lag
    #     )
    #     eventTimeIndex = np.nonzero(testFlight.x[:, 0] == eventTime)[0][0]
    # else:
    #     eventTime = testFlight.tFinal
    #     eventTimeIndex = -1
    col1, col2 = st.columns(2)
    with col1:

        kecDf = pd.DataFrame(
            data=zip(
                testFlight.kineticEnergy[:, 1],
                testFlight.rotationalEnergy[:, 1],
                testFlight.translationalEnergy[:, 1],
            ),
            columns=[
                "Kinetic Energy",
                "Rotational Energy",
                "Translational Energy",
            ],
            index=testFlight.kineticEnergy[:, 0],
        )
        st.write("Kinetic Energy Components (J) vs Time (s)")
        st.line_chart(kecDf)

        thrustAbsDf = pd.DataFrame(
            data=zip(
                testFlight.thrustPower[:, 1],
            ),
            columns=[
                "Power (W)",
            ],
            index=testFlight.thrustPower[:, 0],
        )
        st.write("Thrust Absolute Power (W) vs Time (s)")
        st.line_chart(thrustAbsDf)

    with col2:
        tmeDf = pd.DataFrame(
            data=zip(
                testFlight.totalEnergy[:, 1],
                testFlight.kineticEnergy[:, 1],
                testFlight.potentialEnergy[:, 1],
            ),
            columns=[
                "Total Energy",
                "Kinetic Energy",
                "Potential Energy",
            ],
            index=testFlight.totalEnergy[:, 0],
        )
        st.write("Total Mechanical Energy Components (J) vs Time (s)")
        st.line_chart(tmeDf)

        dragAbsDf = pd.DataFrame(
            data=zip(
                testFlight.dragPower[:, 1],
            ),
            columns=[
                "Power (W)",
            ],
            index=testFlight.dragPower[:, 0],
        )
        st.write("Drag Absolute Power (W) vs Time (s)")
        st.line_chart(dragAbsDf)

    st.write(
        """
        ---
    """
    )
    st.subheader("Trajectory Fluid Mechanics Plots")
    col1, col2 = st.columns(2)
    with col1:

        machDf = pd.DataFrame(
            data=zip(
                testFlight.MachNumber[:, 1],
            ),
            columns=[
                "Mach Number",
            ],
            index=testFlight.MachNumber[:, 0],
        )
        st.write("Mach Number vs Time (s)")
        st.line_chart(machDf)

        reynoldsDf = pd.DataFrame(
            data=zip(
                testFlight.ReynoldsNumber[:, 1],
            ),
            columns=[
                "Reynolds Number",
            ],
            index=testFlight.ReynoldsNumber[:, 0],
        )
        st.write("Reynolds Number vs Time (s)")
        st.line_chart(reynoldsDf)

    with col2:
        pressureDf = pd.DataFrame(
            data=zip(
                testFlight.dynamicPressure[:, 1],
                testFlight.totalPressure[:, 1],
                testFlight.pressure[:, 1],
            ),
            columns=[
                "Dynamic Pressure",
                "Total Pressure",
                "Static Pressure",
            ],
            index=testFlight.dynamicPressure[:, 0],
        )
        st.write("Total and Dynamic Pressure (Pa) vs Time (s)")
        st.line_chart(pressureDf)

        aoaDf = pd.DataFrame(
            data=zip(
                testFlight.angleOfAttack[:, 1],
            ),
            columns=[
                "Angle of Attack (°)",
            ],
            index=testFlight.angleOfAttack[:, 0],
        )
        st.write("Angle of Attack (°) vs Time (s)")
        st.line_chart(aoaDf)
    with st.container():
        if st.session_state.flight["exportKML"]:
            testFlight.exportKML(
                fileName="{}.kml".format(st.session_state.rocketInfo["rocketName"]),
                extrude=True,
                altitudeMode="relativetoground",
            )
            movePath = "streamlitFiles\__pycache__\{}.kml".format(
                st.session_state.rocketInfo["rocketName"]
            )
            shutil.move(
                "{}.kml".format(st.session_state.rocketInfo["rocketName"]), movePath
            )
            st.write(movePath)
            with open(
                movePath,
                "rb",
            ) as file:
                downloadKML = st.download_button(
                    "Download KML",
                    data=file,
                    mime="application/octet-stream",
                    file_name="{}.kml".format(
                        st.session_state.rocketInfo["rocketName"]
                    ),
                )
                if downloadKML:
                    st.success("File ready")
