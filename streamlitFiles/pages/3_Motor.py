import streamlit as st
from rocketpy import Motor
import os
import pandas as pd

st.set_page_config(
    page_title="Motor Configuration",
    page_icon="🔥",
    layout="wide",
    initial_sidebar_state="expanded",
    # menu_items={
    #     "Get Help": "https://www.extremelycoolapp.com/help",
    #     "Report a bug": "https://www.extremelycoolapp.com/bug",
    #     "About": "# This is a header. This is an *extremely* cool app!",
    # },
)
if "Env" not in st.session_state:
    st.warning("Go back to environment page and initialize environment")


def plotGraph(title, tableDF, yAxis):
    st.write(title)
    tableDF.columns = [yAxis]
    st.line_chart(tableDF, use_container_width=True)


# simCheck = False
# if simCheck == True:
#     st.session_state.motorInfo['thrustSource'] = None
if "motor" not in st.session_state:
    st.session_state.motor = None
if "motorInfo" not in st.session_state:
    st.session_state.motorInfo = {
        "thrustSource": "data\motors\Cesaroni_M1670.eng",
        "burnOut": 3.9,
        "grainNumber": 5,
        "grainSeparation": 5 / 1000,
        "grainDensity": 1815.0,
        "grainOuterRadius": 33 / 1000,
        "grainInitialInnerRadius": 15 / 1000,
        "grainInitialHeight": 120 / 1000,
        "nozzleRadius": 33 / 1000,
        "throatRadius": 11 / 1000,
        "interpolationMethod": "linear",
    }
# st.write(st.session_state.motorInfo['thrustSource'])
# thrustSource = st.session_state.motorInfo['thrustSource']
motorList = {}
for root, dir, files in os.walk("data\motors"):
    for file in files:
        filePath = os.path.join(root, file)
        motorList.update({file[:-4]: filePath})  # -4 for removing .eng from name
# st.session_state.motorInfo['thrustSource'] = ""
col1, col2, col3 = st.columns([3, 0.5, 3])
with col3:
    with st.form("Select Engine Form"):
        thrustSourceList = st.selectbox("Pick an Engine", options=motorList.keys())
        submitted = st.form_submit_button("Submit")
        thrustSource = motorList[thrustSourceList]
        if submitted:
            thrustSource = motorList[thrustSourceList]
            st.write(thrustSource)

with col2:
    st.write("#  OR")

with col1:
    with st.form("Upload Engine Form", clear_on_submit=True):
        motorFileUpload = st.file_uploader(label="Upload Engine File", type=["eng"])
        checkboxFlag = st.checkbox("Save to Library")
        # Every form must have a submit button.
        uploaded = st.form_submit_button("Upload")
        if uploaded:
            st.write("Motor", motorFileUpload.name, "Save to Library: ", checkboxFlag)

        if motorFileUpload is not None and checkboxFlag is True:
            fileDetails = {
                "fileName": motorFileUpload.name,
                "fileType": motorFileUpload.type,
            }
            st.write(fileDetails)

            with open(os.path.join("data\motors", motorFileUpload.name), "wb") as f:
                f.write(motorFileUpload.getbuffer())
            thrustSource = os.path.join("data\motors", motorFileUpload.name)

        elif motorFileUpload is not None and checkboxFlag is False:
            with open(
                os.path.join("streamlitFiles\__pycache__", motorFileUpload.name),
                "wb",
            ) as f:
                f.write(motorFileUpload.getbuffer())
                thrustSource = os.path.join(
                    "streamlitFiles\__pycache__", motorFileUpload.name
                )
if uploaded == True or submitted == True:
    st.session_state.motorInfo["thrustSource"] = thrustSource
st.write("**Motor Loaded:**", st.session_state.motorInfo["thrustSource"][12:-4])

col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    st.session_state.motorInfo["burnOut"] = st.number_input(
        "Motor Burnout (s)",
        value=st.session_state.motorInfo["burnOut"],
        min_value=0.0,
        format="%f",
    )
    st.session_state.motorInfo["grainNumber"] = st.number_input(
        "Number of Grains",
        value=st.session_state.motorInfo["grainNumber"],
        min_value=1,
        format="%d",
    )

with col2:
    st.session_state.motorInfo["grainDensity"] = st.number_input(
        "Grain Density (kg/m^3)",
        value=st.session_state.motorInfo["grainDensity"],
        min_value=0.0,
        format="%f",
        help="Density of each grain in kg/meters cubed.",
    )
    st.session_state.motorInfo["grainOuterRadius"] = st.number_input(
        "Grain Outer Radius (m)",
        value=st.session_state.motorInfo["grainOuterRadius"],
        min_value=0.0,
        format="%f",
        help="Outer radius of each grain in meters.",
    )

with col3:

    st.session_state.motorInfo["grainInitialInnerRadius"] = st.number_input(
        "Grain Initial Inner Radius (m)",
        value=st.session_state.motorInfo["grainInitialInnerRadius"],
        min_value=0.0,
        format="%f",
        help="Initial inner radius of each grain in meters.",
    )
    st.session_state.motorInfo["grainInitialHeight"] = st.number_input(
        "Grain Height (m)",
        value=st.session_state.motorInfo["grainInitialHeight"],
        min_value=0.0,
        format="%f",
        help="Initial height of each grain in meters.",
    )

with col4:

    st.session_state.motorInfo["nozzleRadius"] = st.number_input(
        "Nozzle Radius (m)",
        value=st.session_state.motorInfo["nozzleRadius"],
        min_value=0.0,
        format="%f",
        help="Motor's nozzle outlet radius in meters. Used to calculate Kn curve. Optional if the Kn curve is not interesting. Its value does not impact trajectory simulation.",
    )
    st.session_state.motorInfo["throatRadius"] = st.number_input(
        "Nozzle Throat Radius (m)",
        value=st.session_state.motorInfo["throatRadius"],
        min_value=0.0,
        format="%f",
        help="Motor's nozzle throat radius in meters. Its value has very low impact in trajectory simulation, only useful to analyze dynamic instabilities, therefore it is optional.",
    )

with col5:
    st.session_state.motorInfo["grainSeparation"] = st.number_input(
        "Grain Separation (m)",
        value=st.session_state.motorInfo["grainSeparation"],
        min_value=0.0,
        format="%f",
        help="Distance between two grains in meters.",
    )
    st.session_state.motorInfo["interpolationMethod"] = st.selectbox(
        "Interpolation Method",
        ("linear", "akima", "spline"),
        help="Method of interpolation to be used in case thrust curve is given by data set in .csv or .eng, or as an array.",
    )


# check if object is in session state, then execute this as it us
if st.button("Simulate") or st.session_state.motor is not None:
    motir = Motor.SolidMotor(
        thrustSource=st.session_state.motorInfo["thrustSource"],
        burnOut=st.session_state.motorInfo["burnOut"],
        grainNumber=st.session_state.motorInfo["grainNumber"],
        grainSeparation=st.session_state.motorInfo["grainSeparation"],
        grainDensity=st.session_state.motorInfo["grainDensity"],
        grainOuterRadius=st.session_state.motorInfo["grainOuterRadius"],
        grainInitialInnerRadius=st.session_state.motorInfo["grainInitialInnerRadius"],
        grainInitialHeight=st.session_state.motorInfo["grainInitialHeight"],
        nozzleRadius=st.session_state.motorInfo["nozzleRadius"],
        throatRadius=st.session_state.motorInfo["throatRadius"],
        interpolationMethod=st.session_state.motorInfo["interpolationMethod"],
    )
    # saving object to session state
    st.session_state.motor = motir
    col1, col2, col3 = st.columns(3)

    with st.container():

        with col1:
            st.subheader("Nozzle Details")
            st.write("Nozzle Radius: ", motir.nozzleRadius, " m")
            st.write("Nozzle Throat Radius: ", motir.throatRadius, " m")

        with col2:
            st.subheader("Grain Details")
            st.write("Number of Grains: ", motir.grainNumber)
            st.write("Grain Spacing: ", motir.grainSeparation, " m")
            st.write("Grain Density: ", motir.grainDensity, " kg/m3")
            st.write("Grain Mass: ", motir.grainInitialMass, " kg")
            st.write("Grain Outer Radius: ", motir.grainOuterRadius, " m")
            st.write("Grain Inner Radius: ", motir.grainInitialInnerRadius, " m")
            st.write("Grain Height: ", motir.grainInitialHeight, " m")
            st.write("Grain Volume: ", (motir.grainInitialVolume * 1000000), " cm3")

        with col3:
            st.subheader("Motor Details")
            st.write("Total Burning Time: ", motir.burnOutTime, " s")
            st.write("Total Propellant Mass: ", motir.propellantInitialMass, " Kg")
            st.write("Propellant Exhaust Velocity: ", motir.exhaustVelocity, "m/s")
            st.write("Average Thrust: ", motir.averageThrust, " N")
            st.write(
                "Maximum Thrust: ",
                motir.maxThrust,
                " N at ",
                motir.maxThrustTime,
                " s after ignition",
            )
            st.write("Total Impulse: ", motir.totalImpulse, " Ns")

    st.subheader("Plots")
    col4, col5, col6, col7 = st.columns(4)
    with col4:

        thrustDF = pd.DataFrame(
            motir.thrust.source[:, 1], index=motir.thrust.source[:, 0]
        )
        plotGraph("Thrust(N) x Time(s)", thrustDF, "Thrust (N)")

        grainHeightDF = pd.DataFrame(
            motir.grainHeight.source[:, 1],
            index=motir.grainHeight.source[:, 0],
        )
        plotGraph("Grain Height (M) x Time(s)", grainHeightDF, "Grain Height (M)")

        inertiaIDF = pd.DataFrame(
            motir.inertiaI.source[:, 1],
            index=motir.inertiaI.source[:, 0],
        )
        plotGraph(
            "Propellant Inertia I (Kg*M^2) x Time(s)",
            inertiaIDF,
            "Propellant Inertia I (Kg*M^2)",
        )

    with col5:

        massDF = pd.DataFrame(motir.mass.source[:, 1], index=motir.mass.source[:, 0])
        plotGraph(
            "Propellant Total Mass (Kg) x Time(s)", massDF, "Propellant Total Mass (Kg)"
        )

        burnRateDF = pd.DataFrame(
            motir.burnRate.source[:, 1],
            index=motir.burnRate.source[:, 0],
        )
        plotGraph("Burn Rate (M/s) x Time(s)", burnRateDF, "Burn Rate (M/s)")

        inertiaIDotDF = pd.DataFrame(
            motir.inertiaIDot.source[:, 1],
            index=motir.inertiaIDot.source[:, 0],
        )
        plotGraph(
            "Propellant Inertia I Dot (Kg*M^2/s) x Time(s)",
            inertiaIDotDF,
            "Propellant Inertia I Dot (Kg*M^2/s)",
        )

    with col6:

        massDotDF = pd.DataFrame(
            motir.massDot.source[:, 1], index=motir.massDot.source[:, 0]
        )
        plotGraph("Mass Dot (Kg/S) x Time(s)", massDotDF, "Mass Dot (Kg/S)")

        burnAreaDF = pd.DataFrame(
            motir.burnArea.source[:, 1],
            index=motir.burnArea.source[:, 0],
        )
        plotGraph("Burn Area (M^2) x Time(s)", burnAreaDF, "Burn Area (M^2)")

        inertiaZDF = pd.DataFrame(
            motir.inertiaZ.source[:, 1],
            index=motir.inertiaZ.source[:, 0],
        )
        plotGraph(
            "Propellant Inertia Z (Kg*M^2) x Time(s)",
            inertiaZDF,
            "Propellant Inertia Z (Kg*M^2)",
        )

    with col7:

        grainInnerRadiusDF = pd.DataFrame(
            motir.grainInnerRadius.source[:, 1],
            index=motir.grainInnerRadius.source[:, 0],
        )
        plotGraph(
            "Grain Inner Radius (M) x Time(s)",
            grainInnerRadiusDF,
            "Grain Inner Radius (M)",
        )

        KnDF = pd.DataFrame(
            motir.Kn.source[:, 1],
            index=motir.Kn.source[:, 0],
        )
        plotGraph("Kn (M^2/M^2) x Grain Inner Radius(M)", KnDF, "Kn (M^2/M^2)")

        inertiaZDotDF = pd.DataFrame(
            motir.inertiaZDot.source[:, 1],
            index=motir.inertiaZDot.source[:, 0],
        )
        plotGraph(
            "Propellant Inertia Z Dot (Kg*M^2/s) x Time(s)",
            inertiaZDotDF,
            "Propellant Inertia Z Dot (Kg*M^2/s)",
        )
