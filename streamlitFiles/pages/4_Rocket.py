import streamlit as st
import os
from rocketpy import Rocket, Flight

st.set_page_config(
    page_title="Rocket Configuration",
    page_icon="🚀",
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


def mainTrigger(p, y):
    # p = pressure
    # y = [x, y, z, vx, vy, vz, e0, e1, e2, e3, w1, w2, w3]
    # activate drogue when vz < 0 m/s.
    return (
        True
        if y[5] < 0
        and y[2]
        < st.session_state.rocketInfo["deploymentHeight"]
        + st.session_state.envInfo["elevation"]
        else False
    )


def drogueTrigger(p, y):
    # p = pressure
    # y = [x, y, z, vx, vy, vz, e0, e1, e2, e3, w1, w2, w3]
    # activate drogue when vz < 0 m/s.
    return True if y[5] < 0 else False


if "rocketInfo" not in st.session_state:
    st.session_state.rocketInfo = {
        "rocketName": "rokit",
        "radius": 127 / 2000,
        "mass": 19.197 - 2.956,
        "inertiaI": 6.60,
        "inertiaZ": 0.0351,
        "distanceRocketNozzle": -1.255,
        "distanceRocketPropellant": -0.85704,
        "powerOffDrag": "",
        "powerOnDrag": "",
        "setRailButtons": (-0.5, 0.2),
        "noseLength": 0.55829,
        "noseKind": "Von Karman",
        "noseDistanceToCM": 0.71971,
        "finCount": 4,
        "finSpan": 0.1,
        "finRC": 0.12,
        "finTC": 0.04,
        "finDistanceToCM": 0.71971,
        "tailTopRadius": 0.0635,
        "tailBottomRadius": 0.0435,
        "tailLength": 0.060,
        "tailDistanceToCM": -1.194656,
        "mainCDS": 10.0,
        "drogueCDS": 1.0,
        "deploymentHeight": 800.0,
        # "mainTrigger": None,
        # "drogueTrigger": drogueTrigger,
        "mainLag": 0.0,
        "drogueLag": 0.0,
        "mainObject": None,
        "drogueObject": None,
    }


rocketTab, aeroTab, parachuteTab = st.tabs(
    ["Create a Rocket", "Aerodynamic Surfaces", "Recovery Configuration"]
)

with rocketTab:
    with st.container():
        st.session_state.rocketInfo["rocketName"] = st.text_input(
            "Rocket Name", value=st.session_state.rocketInfo["rocketName"]
        )
        col1, col2, col3 = st.columns(3)
        with col1:
            st.session_state.rocketInfo["inertiaI"] = st.number_input(
                "Unloaded longitudenal moment of inertia (Kg*m^2)",
                value=st.session_state.rocketInfo["inertiaI"],
                format="%f",
                help="Rocket's moment of inertia, without propellant, with respect to to an axis perpendicular to the rocket's axis of cylindrical symmetry, in kg*m^2.",
            )

            st.session_state.rocketInfo["distanceRocketPropellant"] = st.number_input(
                "(Unloaded) Distance between Rocket CG and and Motor CG",
                value=st.session_state.rocketInfo["distanceRocketPropellant"],
                format="%f",
                help="Distance between rocket's center of mass, without propellant, to the motor reference point, which for solid and hybrid motors is the center of mass of solid propellant, in meters.",
            )

        with col2:
            st.session_state.rocketInfo["radius"] = st.number_input(
                "Outer Radius (m)",
                value=st.session_state.rocketInfo["radius"],
                min_value=0.0,
                format="%f",
            )

            st.session_state.rocketInfo["inertiaZ"] = st.number_input(
                "Unloaded rotational moment of inertia (Kg*m^2)",
                value=st.session_state.rocketInfo["inertiaZ"],
                format="%f",
                help="Rocket's moment of inertia, without propellant, with respect to the rocket's axis of cylindrical symmetry, in kg*m^2.",
            )

        with col3:
            st.session_state.rocketInfo["mass"] = st.number_input(
                "Dry Mass (Kg)",
                value=st.session_state.rocketInfo["mass"],
                min_value=0.0,
                format="%f",
                help="Rocket's mass without propellant in kg.",
            )

            st.session_state.rocketInfo["distanceRocketNozzle"] = st.number_input(
                "(Unloaded) Distance between CG and Nozzle Exit (m)",
                value=st.session_state.rocketInfo["distanceRocketNozzle"],
                format="%f",
                help="Distance between rocket's center of mass, without propellant, to the exit face of the nozzle, in meters.",
            )

        col1, col2 = st.columns(2)
        with col1:
            with st.container():

                powerOff = st.file_uploader(
                    label="Power off Drag Curve Data",
                    help="Rocket's drag coefficient as a function of Mach number when the motor is off. Can be found via rasAero (2 column data: mach number and poweroffDrag, keep mach number upto 2)",
                    type=["csv"],
                )
                # st.write(st.session_state.rocket.powerOffDrag)  # debug
                if powerOff is not None:
                    with open(
                        os.path.join(
                            "streamlitFiles\__pycache__",
                            powerOff.name,
                        ),
                        "wb",
                    ) as f:
                        f.write(powerOff.getbuffer())
                        st.session_state.rocketInfo["powerOffDrag"] = os.path.join(
                            "streamlitFiles\__pycache__",
                            powerOff.name,
                        )
                    st.success("File Uploaded")

        with col2:
            with st.container():
                powerOn = st.file_uploader(
                    label="Power on Drag Curve Data",
                    help="Rocket's drag coefficient as a function of Mach number when the motor is on. Can be found via rasAero (2 column data: mach number and powerOnDrag, keep mach number upto 2)",
                    type=["csv"],
                )
                # st.write(st.session_state.rocket.powerOnDrag)  # debug
                if powerOn is not None:
                    with open(
                        os.path.join(
                            "streamlitFiles\__pycache__",
                            powerOn.name,
                        ),
                        "wb",
                    ) as f:
                        f.write(powerOn.getbuffer())
                        st.session_state.rocketInfo["powerOnDrag"] = os.path.join(
                            "streamlitFiles\__pycache__",
                            powerOn.name,
                        )
                    st.success("File Uploaded")

        if "rocketInitialized" not in st.session_state:
            st.session_state.rocketInitialized = False

        st.write(st.session_state.motorInfo["thrustSource"][12:-4])
        rocketSubmitted = st.button("Submit")
        if rocketSubmitted == True or st.session_state.rocketInitialized == True:

            rokit = Rocket(
                motor=st.session_state.motor,
                radius=st.session_state.rocketInfo["radius"],
                mass=st.session_state.rocketInfo["mass"],
                inertiaI=st.session_state.rocketInfo["inertiaI"],
                inertiaZ=st.session_state.rocketInfo["inertiaZ"],
                distanceRocketNozzle=st.session_state.rocketInfo[
                    "distanceRocketNozzle"
                ],
                distanceRocketPropellant=st.session_state.rocketInfo[
                    "distanceRocketPropellant"
                ],
                powerOffDrag=st.session_state.rocketInfo["powerOffDrag"],
                powerOnDrag=st.session_state.rocketInfo["powerOnDrag"],
            )
            st.session_state.rocketInitialized = True
            st.session_state.rocket = rokit

            st.success("Rocket Added")
            with aeroTab:
                st.session_state.rocketInfo["setRailButtons"] = st.slider(
                    "Set Rail button position relative to CG",
                    -5.0,
                    5.0,
                    st.session_state.rocketInfo["setRailButtons"],
                    step=0.01,
                    help="Two values which represent the distance of each of the two rail buttons to the center of mass of the rocket without propellant. If the rail button is positioned above the center of mass, its distance should be a positive value. If it is below, its distance should be a negative value. The order does not matter. All values should be in meters.",
                )
                st.session_state.rocketInfo["setRailButtons"] = list(
                    st.session_state.rocketInfo["setRailButtons"]
                )
                rokit.setRailButtons(st.session_state.rocketInfo["setRailButtons"])
                st.write("""---""")
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.write("Nosecone Setup")
                    st.session_state.rocketInfo["noseLength"] = st.number_input(
                        "Add Nosecone Length",
                        value=st.session_state.rocketInfo["noseLength"],
                        min_value=0.0,
                        format="%f",
                        help="Nose cone length or height in meters. Must be a positive value.",
                    )
                    st.session_state.rocketInfo["noseKind"] = st.selectbox(
                        "Nosecone Shape Type",
                        options=("Von Karman", "conical", "ogive", " lvhaack"),
                    )
                    st.session_state.rocketInfo["noseDistanceToCM"] = st.number_input(
                        "Distance between unloaded rocket CM and base of nose",
                        min_value=0.0,
                        value=st.session_state.rocketInfo["noseDistanceToCM"],
                        format="%f",
                        help="Nose cone position relative to rocket unloaded center of mass, considering positive direction from center of mass to nose cone. Consider the center point belonging to the nose cone base to calculate distance.",
                    )

                    rokit.addNose(
                        length=st.session_state.rocketInfo["noseLength"],
                        kind=st.session_state.rocketInfo["noseKind"],
                        distanceToCM=st.session_state.rocketInfo["noseDistanceToCM"],
                    )

                with col2:
                    # finShape= st.selectbox("Fin Shape Type:", options=("Trapezoidal", "Elliptical"))
                    st.write("Adding Trapezoidal finset")
                    st.session_state.rocketInfo["finCount"] = st.number_input(
                        "# of Fins",
                        min_value=0,
                        value=st.session_state.rocketInfo["finCount"],
                        format="%d",
                    )
                    st.session_state.rocketInfo["finSpan"] = st.number_input(
                        "Fin Span",
                        min_value=0.0,
                        value=st.session_state.rocketInfo["finSpan"],
                        format="%f",
                        help="Fin Span/Height in meters",
                    )
                    st.session_state.rocketInfo["finRC"] = st.number_input(
                        "Root length",
                        min_value=0.0,
                        value=st.session_state.rocketInfo["finRC"],
                        format="%f",
                        help="Fin Root Length in meters",
                    )
                    st.session_state.rocketInfo["finTC"] = st.number_input(
                        "Tip length",
                        min_value=0.0,
                        value=st.session_state.rocketInfo["finTC"],
                        format="%f",
                        help="Fin Tip Length in meters",
                    )
                    st.session_state.rocketInfo["finDistanceToCM"] = st.number_input(
                        "Distance between unloaded rocket CM and base of nose",
                        min_value=0.0,
                        value=st.session_state.rocketInfo["finDistanceToCM"],
                        format="%f",
                        help="Fin set position relative to rocket unloaded center of mass, considering positive direction from center of mass to nose cone. Consider the center point belonging to the top of the fins to calculate distance.",
                    )
                    rokit.addFins(
                        st.session_state.rocketInfo["finCount"],
                        span=st.session_state.rocketInfo["finSpan"],
                        rootChord=st.session_state.rocketInfo["finRC"],
                        tipChord=st.session_state.rocketInfo["finTC"],
                        distanceToCM=(
                            st.session_state.rocketInfo["finDistanceToCM"] * -1
                        ),
                        radius=0,
                        cantAngle=0,
                    )
                    # if finShape == "Trapezoidal":
                    #     rokit.addFin(
                    #         finCount,
                    #         finSpan,
                    #         finRC,
                    #         finTC,
                    #         finDistanceToCM,
                    #         radius=0,
                    #         cantAngle=0,
                    #     )
                    # else:
                    #     rokit.addEllipticalFins(finCount, finRC, finSpan, finDistanceToCM, cantAngle=)

                with col3:
                    # bug fix: add session state to add tail
                    addTailCheck = st.checkbox("Add Tail?")
                    if addTailCheck:
                        st.write("Tail Setup")
                        st.session_state.rocketInfo["tailTopRadius"] = st.number_input(
                            "Tail Top Radius",
                            value=st.session_state.rocketInfo["tailTopRadius"],
                            min_value=0.0,
                            help="Tail top radius in meters, considering positive direction from center of mass to nose cone.",
                        )
                        st.session_state.rocketInfo[
                            "tailBottomRadius"
                        ] = st.number_input(
                            "Tail Bottom Radius",
                            value=st.session_state.rocketInfo["tailBottomRadius"],
                            min_value=0.0,
                            help="Tail bottom radius in meters, considering positive direction from center of mass to nose cone.",
                        )
                        st.session_state.rocketInfo["tailLength"] = st.number_input(
                            "Tail Length/Height",
                            value=st.session_state.rocketInfo["tailLength"],
                            min_value=0.0,
                            help="Tail length or height in meters. Must be a positive value.",
                        )
                        st.session_state.rocketInfo[
                            "tailDistanceToCM"
                        ] = st.number_input(
                            "Tail Distance To CM",
                            value=st.session_state.rocketInfo["tailDistanceToCM"],
                            help="Tail position relative to rocket unloaded center of mass, considering positive direction from center of mass to nose cone. Consider the point belonging to the tail which is closest to the unloaded center of mass to calculate distance.",
                        )

                        tailSubmit = st.button(label="Add Tail")
                        if tailSubmit:
                            rokit.addTail(
                                topRadius=st.session_state.rocketInfo["tailTopRadius"],
                                bottomRadius=st.session_state.rocketInfo[
                                    "tailBottomRadius"
                                ],
                                length=st.session_state.rocketInfo["tailLength"],
                                distanceToCM=(
                                    st.session_state.rocketInfo["tailDistanceToCM"] * -1
                                ),
                            )
            with parachuteTab:
                st.write("Hello")
                col1, col2 = st.columns(2)
                with col1:
                    with st.form("addDrogue", clear_on_submit=False):
                        st.write("Drogue Params")
                        st.session_state.rocketInfo["drogueCDS"] = st.number_input(
                            "Drag coefficient times reference area for parachute. (CdS)",
                            value=st.session_state.rocketInfo["drogueCDS"],
                            format="%f",
                            help=" Drag coefficient times reference area for parachute. It is used to compute the drag force exerted on the parachute by the equation F = ((1/2)*rho*V^2)*CdS, that is, the drag force is the dynamic pressure computed on the parachute times its CdS coefficient. Has units of area and must be given in squared meters.",
                        )
                        chuteDeployment = st.selectbox(
                            "Chute deployed at?",
                            options=[
                                "Apogee",
                                # "Never",
                            ],
                        )
                        # if chuteDeployment == "Apogee":
                        st.session_state.rocketInfo["drogueLag"] = st.number_input(
                            "Lag between deployment (in seconds)",
                            value=st.session_state.rocketInfo["drogueLag"],
                            min_value=0.0,
                        )
                        st.session_state.rocketInfo["drogueTrigger"] = drogueTrigger
                        # st.session_state.rocketInfo["drogueObject"] = None

                        drogue = rokit.addParachute(
                            "Drogue",
                            CdS=st.session_state.rocketInfo["drogueCDS"],
                            trigger=st.session_state.rocketInfo["drogueTrigger"],
                            samplingRate=105,
                            lag=st.session_state.rocketInfo["drogueLag"],
                            noise=(0, 8.3, 0.5),
                        )
                        addDrogueCheck = st.form_submit_button("Add Drogue")
                        if addDrogueCheck:
                            # if drogue is not None:
                            #     rokit.parachutes.remove(drogue)
                            st.success("Drogue Chute Added", icon="✅")

                    # st.write("Drogue Params")
                    # st.session_state.rocketInfo["drogueCDS"] = st.number_input(
                    #     "Drag coefficient times reference area for parachute. (CdS)",
                    #     value=st.session_state.rocketInfo["drogueCDS"],
                    #     format="%f",
                    #     help=" Drag coefficient times reference area for parachute. It is used to compute the drag force exerted on the parachute by the equation F = ((1/2)*rho*V^2)*CdS, that is, the drag force is the dynamic pressure computed on the parachute times its CdS coefficient. Has units of area and must be given in squared meters.",
                    # )
                    # chuteDeployment = st.selectbox(
                    #     "Chute deployed at?",
                    #     options=[
                    #         "Apogee",
                    #         "Never",
                    #     ],
                    # )
                    # if chuteDeployment == "Apogee":
                    #     st.session_state.rocketInfo["drogueLag"] = st.number_input(
                    #         "Lag between deployment (in seconds)",
                    #         value=st.session_state.rocketInfo["drogueLag"],
                    #         min_value=0.0,
                    #     )
                    #     st.session_state.rocketInfo["drogueTrigger"] = drogueTrigger

                    # # st.session_state.rocketInfo["drogueObject"] = None
                    # st.write(st.session_state.rocketInfo["drogueTrigger"])
                    # st.write(type(st.session_state.rocketInfo["drogueTrigger"]))

                    # addDrogueCheck = st.button("Add Drogue")
                    # if addDrogueCheck:
                    #     # if drogue is not None:
                    #     #     rokit.parachutes.remove(drogue)
                    #     testing = rokit.addParachute(
                    #         "Drogue",
                    #         CdS=5,  # st.session_state.rocketInfo["drogueCDS"],
                    #         trigger=drogueTrigger,
                    #         samplingRate=105,
                    #         lag=2,  # st.session_state.rocketInfo["drogueLag"],
                    #         noise=(0, 8.3, 0.5),
                    #     )
                    #     st.success("Drogue Chute Added", icon="✅")
                    #     # st.write(type(st.session_state.rocketInfo["drogueTrigger"]))

                with col2:
                    with st.form("addMain", clear_on_submit=False):
                        st.write("Main Params")
                        st.session_state.rocketInfo["mainCDS"] = st.number_input(
                            "Drag coefficient times reference area for parachute. (CdS)",
                            value=st.session_state.rocketInfo["mainCDS"],
                            format="%f",
                            help=" Drag coefficient times reference area for parachute. It is used to compute the drag force exerted on the parachute by the equation F = ((1/2)*rho*V^2)*CdS, that is, the drag force is the dynamic pressure computed on the parachute times its CdS coefficient. Has units of area and must be given in squared meters.",
                        )
                        chuteDeployment = st.selectbox(
                            "Chute deployed at?",
                            options=[
                                "Specific Height at descent",
                                # "Launch +N seconds",
                            ],
                        )
                        # if chuteDeployment == "Specific Height at descent":
                        st.session_state.rocketInfo[
                            "deploymentHeight"
                        ] = st.number_input(
                            "Deployment Altitude AGL",
                            value=st.session_state.rocketInfo["deploymentHeight"],
                            min_value=0.0,
                        )
                        st.session_state.rocketInfo["mainLag"] = st.number_input(
                            "Lag between deployment (in seconds)",
                            value=st.session_state.rocketInfo["mainLag"],
                            min_value=0.0,
                            key=float,
                            step=0.1,
                        )
                        st.session_state.rocketInfo["mainTrigger"] = mainTrigger

                        # elif chuteDeployment == "Launch +N seconds":
                        #     st.session_state.rocketInfo["deploymentHeight"] = 0.0
                        #     st.session_state.rocketInfo["mainLag"] = st.number_input(
                        #         "Lag between deployment and launch (in seconds)",
                        #         value=st.session_state.rocketInfo["mainLag"],
                        #         min_value=0.0,
                        #     )
                        #     st.session_state.rocketInfo["mainTrigger"] = mainTrigger

                        main = rokit.addParachute(
                            "Main",
                            CdS=st.session_state.rocketInfo["mainCDS"],
                            trigger=st.session_state.rocketInfo["mainTrigger"],
                            samplingRate=105,
                            lag=st.session_state.rocketInfo["mainLag"],
                            noise=(0, 8.3, 0.5),
                        )

                        addMainCheck = st.form_submit_button("Add Main")
                        if addMainCheck:
                            # if main is not None:
                            #     rokit.parachutes.remove(main)

                            st.success("Main Chute Added", icon="✅")

                #     st.write("Main Params")
                #     # st.session_state.rocketInfo["mainCDS"] = st.number_input(
                #     #     "Drag coefficient times reference area for parachute. (CdS)",
                #     #     value=st.session_state.rocketInfo["mainCDS"],
                #     #     format="%f",
                #     #     help=" Drag coefficient times reference area for parachute. It is used to compute the drag force exerted on the parachute by the equation F = ((1/2)*rho*V^2)*CdS, that is, the drag force is the dynamic pressure computed on the parachute times its CdS coefficient. Has units of area and must be given in squared meters.",
                #     # )
                #     # chuteDeployment = st.selectbox(
                #     #     "Chute deployed at?",
                #     #     options=[
                #     #         "Specific Height at descent",
                #     #         "Launch +N seconds",
                #     #     ],
                #     # )
                #     # if chuteDeployment == "Specific Height at descent":
                #     #     st.session_state.rocketInfo[
                #     #         "deploymentHeight"
                #     #     ] = st.number_input(
                #     #         "Deployment Altitude AGL",
                #     #         value=st.session_state.rocketInfo["deploymentHeight"],
                #     #         min_value=0.0,
                #     #     )
                #     #     st.session_state.rocketInfo["mainLag"] = st.number_input(
                #     #         "Lag between deployment (in seconds)",
                #     #         value=st.session_state.rocketInfo["mainLag"],
                #     #         min_value=0.0,
                #     #         key=float,
                #     #         step=0.1,
                #     #     )
                #     #     st.session_state.rocketInfo["mainTrigger"] = mainTrigger

                #     # if chuteDeployment == "Launch +N seconds":
                #     #     st.session_state.rocketInfo["deploymentHeight"] = 0.0
                #     #     st.session_state.rocketInfo["mainLag"] = st.number_input(
                #     #         "Lag between deployment and launch (in seconds)",
                #     #         value=st.session_state.rocketInfo["mainLag"],
                #     #         min_value=0.0,
                #     #     )
                #     #     st.session_state.rocketInfo["mainTrigger"] = mainTrigger

                #     # main = None

                #     # addMainCheck = st.button("Add Main")
                #     # if addMainCheck:
                #     #     # if main is not None:
                #     #     #     rokit.parachutes.remove(main)
                #     #     main = rokit.addParachute(
                #     #         "Main",
                #     #         CdS=st.session_state.rocketInfo["mainCDS"],
                #     #         trigger=st.session_state.rocketInfo["mainTrigger"],
                #     #         samplingRate=105,
                #     #         lag=st.session_state.rocketInfo["mainLag"],
                #     #         noise=(0, 8.3, 0.5),
                #     #     )
                #     #     st.success("Main Chute Added", icon="✅")


# def drogueTrigger(p, y):
#     # p = pressure
#     # y = [x, y, z, vx, vy, vz, e0, e1, e2, e3, w1, w2, w3]
#     # activate drogue when vz < 0 m/s.
#     return True if y[5] < 0 else False


# def mainTrigger(p, y):
#     # p = pressure
#     # y = [x, y, z, vx, vy, vz, e0, e1, e2, e3, w1, w2, w3]
#     # activate main when vz < 0 m/s and z < 800 + 1400 m (+1400 due to surface elevation).
#     return True if y[5] < 0 and y[2] < 800 + 1400 else False


# Main = rokit.addParachute(
#     "Main",
#     CdS=10.0,
#     trigger=mainTrigger,
#     samplingRate=105,
#     lag=1.5,
#     noise=(0, 8.3, 0.5),
# )

# Drogue = rokit.addParachute(
#     "Drogue",
#     CdS=1.0,
#     trigger=drogueTrigger,
#     samplingRate=105,
#     lag=1.5,
#     noise=(0, 8.3, 0.5),
# )
# rokit = rokit
# st.write(rokit.parachutes[:])
# tf = st.button("simulate")
# if tf:
#     st.write(rokit.parachutes[:])
#     testFlight = Flight(
#         rocket=rokit,
#         environment=st.session_state.Env,
#         inclination=85,
#         heading=0,
#     )
#     testFlight.info()
