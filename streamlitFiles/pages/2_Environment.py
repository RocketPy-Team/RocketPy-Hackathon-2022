from datetime import datetime
import streamlit as st
import pandas as pd
from rocketpy import Environment
import numpy as np
import altair as alt

st.set_page_config(
    page_title="Environment Setup",
    page_icon="🏞️",
    layout="wide",
    initial_sidebar_state="expanded",
    # menu_items={
    #     "Get Help": "https://www.extremelycoolapp.com/help",
    #     "Report a bug": "https://www.extremelycoolapp.com/bug",
    #     "About": "# This is a header. This is an *extremely* cool app!",
    # },
)
# def sidebar():
#     # "with" notation
#     with st.sidebar:
#         # st.markdown("""---""")
#         launchPreset = st.radio(
#             label="Pick a Preset:",
#             options=("Spaceport America VLA", "Cape Canaveal, SP", "Custom"),
#         )


if "Env" not in st.session_state:
    st.session_state.Env = None

if "envInfo" not in st.session_state:
    st.session_state.envInfo = {
        "railLength": 5.2,
        "latitude": 32.990254,
        "longitude": -106.974998,
        "elevation": 1400.0,
    }
# sidebar()
# input container
with st.container():
    col1, col2, col3, col4 = st.columns(4)
    # rail length column
    with col1:
        st.session_state.envInfo["railLength"] = st.number_input(
            "Size of Launch Rail in meters",
            value=st.session_state.envInfo["railLength"],
        )

    # latitude column
    with col2:
        st.session_state.envInfo["latitude"] = st.number_input(
            "Latitude",
            value=st.session_state.envInfo["latitude"],
            min_value=-180.0,
            max_value=180.0,
            step=1e-12,
            format="%.5f",
        )
    # longitude column
    with col3:
        st.session_state.envInfo["longitude"] = st.number_input(
            "Longitude",
            value=st.session_state.envInfo["longitude"],
            min_value=-180.0,
            max_value=180.0,
            step=1e-12,
            format="%.5f",
        )
    # elevation column
    with col4:
        st.session_state.envInfo["elevation"] = st.number_input(
            "Mean sea level elevation in meters",
            value=st.session_state.envInfo["elevation"],
        )

    col5, col6, col7, col8 = st.columns(4)
    with col5:
        wenHopDate = st.date_input("When are you launching?")

    with col6:
        wenHopTime = st.time_input("what Time (in UTC)", key=int)

    Env = Environment(
        railLength=st.session_state.envInfo["railLength"],
        latitude=st.session_state.envInfo["latitude"],
        longitude=st.session_state.envInfo["longitude"],
        elevation=st.session_state.envInfo["elevation"],
    )
    Env.setDate(date=datetime.combine(wenHopDate, wenHopTime))
    with col7:
        atmosModel = st.selectbox(
            "Atmospheric Model for Environment Setup",
            (
                "StandardAtmosphere",  # add option to manualy assign wind velocity
                "WyomingSounding",
                "NOAARucSounding",
                "Forecast",
                "Reanalysis",
                "Ensemble",
            ),
            help="Describes the atmospheric model which is being used. Use Forecast for precise results.",
        )

    with col8:
        if atmosModel == "Forecast":
            file_type = st.selectbox(
                "Specify Location of Data", ("GFS", "FV3", "RAP", "NAM")
            )
        elif atmosModel == "Ensemble":
            file_type = st.selectbox("Specify Location of Data", ("GEFS", "CMC"))
        else:
            file_type = ""

        if file_type == "":
            Env.setAtmosphericModel(type=atmosModel)
        else:
            Env.setAtmosphericModel(type=atmosModel, file=file_type)

    df = pd.DataFrame(
        np.random.randn(1, 2) / [25, 25]
        + [st.session_state.envInfo["latitude"], st.session_state.envInfo["longitude"]],
        columns=["lat", "lon"],
    )  # still figuring out how this works, 25,25 seems to be the area wrt the lat and long, and idk how random pays a role here just picked from the documentation
    st.map(
        df,
        use_container_width=True,
        zoom=8,
    )

    environmentTableDict = dict(Env.allInfoReturned())
    environmentPlotsDict = dict(Env.allPlotInfoReturned())
    if st.button("Run Environment Simulation"):
        # saving object to session state
        st.session_state.Env = Env
        col1, col2, col3 = st.columns(3)

        with col1:
            st.subheader("Launch Site Details")
            st.write(
                "Launch Rail length: \t",
                environmentTableDict["launch_rail_length"],
                " m",
            )
            st.write("Launch Date: \t", environmentTableDict["launch_date"], " UTC")
            st.write("Launch Site Latitude: \t", environmentTableDict["lat"], " °")
            st.write("Launch Site Longitude: \t", environmentTableDict["lon"], " °")
            # df = pd.DataFrame([latitude, longitude], columns=["lat", "lon"])
            # df = {"lat": latitude, "lon": longitude}
            st.write(
                "Launch Site UTM coordinates: ",
                environmentTableDict["initialEast"],
                " ",
                environmentTableDict["initialEW"],
                "\t",
                environmentTableDict["initialNorth"],
                " ",
                environmentTableDict["initialHemisphere"],
            )
            st.write("Reference Datum: \t", environmentTableDict["datum"])
            st.write(
                "Launch Site Surface Elevation: \t",
                environmentTableDict["elevation"],
                " m",
            )

            # st.write(environmentTableDict)

        with col2:
            st.subheader("Atmospheric Model Details")
            st.write(
                "Atmospheric Model Type: \t",
                environmentTableDict["modelType"],
            )
            st.write(
                environmentTableDict["modelType"],
                " Maximum Height: \t",
                environmentTableDict["modelTypeMaxExpectedHeight"] / 1000,
                " Km",
            )
            if environmentTableDict["modelType"] != "StandardAtmosphere":
                st.write(
                    environmentTableDict["modelType"],
                    " Time Period: From \t",
                    environmentTableDict["initDate"],
                    " to ",
                    environmentTableDict["endDate"],
                )
                st.write(
                    environmentTableDict["modelType"],
                    " Hour Interval: \t",
                    environmentTableDict["interval"],
                )
                st.write(
                    environmentTableDict["modelType"],
                    " Latitude Range: From \t",
                    environmentTableDict["initLat"],
                    "° to ",
                    environmentTableDict["endLat"],
                    "°",
                )
                st.write(
                    environmentTableDict["modelType"],
                    " Longitude Range: From \t",
                    environmentTableDict["initLon"],
                    "° to ",
                    environmentTableDict["endLon"],
                    "°",
                )

            with col3:
                st.subheader("Surface Atmospheric Conditions")
                st.write(
                    "Surface Wind Speed", environmentTableDict["windSpeed"], " m/s"
                )
                st.write(
                    "Surface Wind Direction", environmentTableDict["windDirection"], "°"
                )
                st.write(
                    "Surface Wind Heading", environmentTableDict["windHeading"], "°"
                )
                st.write(
                    "Surface Pressure", environmentTableDict["surfacePressure"], " hPa"
                )
                st.write(
                    "Surface Temprature",
                    environmentTableDict["surfaceTemperature"],
                    " K",
                )
                st.write(
                    "Surface Air Density",
                    environmentTableDict["surfaceAirDensity"],
                    " kg/m³",
                )
                st.write(
                    "Surface Speed of Sound",
                    environmentTableDict["surfaceSpeedOfSound"],
                    " m/s",
                )

        st.subheader("Atmospheric Plots")
        # st.write(environmentPlotsDict)
        ax1 = pd.DataFrame(environmentPlotsDict)
        # st.vega_lite_chart(ax1,{
        #     'mark'{'type': 'line', 'tooltip':True},
        #     'encoding':{
        #         'x': {'field':'Height above Sea Level (m)'},
        #         'y': {'field':'Wind Speed (m/s) & Wind Direction (°)'}
        #     }
        # })
        # st.line_chart(
        #     ax1,
        #     # x="Height above Sea Level (m)",
        #     # y="Wind Speed (m/s) & Wind Direction (°)",
        # )
        ampCol1, ampCol2 = st.columns(2)
        with ampCol1:
            gridBase = alt.Chart(ax1).encode(
                alt.Y(
                    "grid",
                    axis=alt.Axis(title="Height above Sea Level (m)"),
                )
            )
            windSpeedLine = (
                gridBase.mark_line(stroke="#FF0000", interpolate="monotone")
                .encode(
                    alt.X(
                        "average(windSpeed)",
                        axis=alt.Axis(
                            title="Wind Speed (m/s)",
                            titleColor="#FF0000",
                        ),
                    ),
                    tooltip=["grid", "windSpeed"],
                )
                .interactive()
            )
            windDirLine = (
                gridBase.mark_line(stroke="#0000FF", interpolate="monotone")
                .encode(
                    alt.X(
                        "average(windDirection)",
                        axis=alt.Axis(
                            title="Wind Direction (°)",
                            titleColor="#0000FF",
                        ),
                    ),
                    tooltip=["grid", "windDirection"],
                )
                .interactive()
            )
            chart1 = alt.layer(windSpeedLine, windDirLine).resolve_scale(
                x="independent"
            )
            st.altair_chart(chart1, use_container_width=True)

            windVelXLine = (
                gridBase.mark_line(stroke="#FF0000", interpolate="monotone")
                .encode(
                    alt.X(
                        "average(windVelX)",
                        axis=alt.Axis(
                            title="Wind U (m/s)",
                            titleColor="#FF0000",
                        ),
                    ),
                    tooltip=["grid", "windVelX"],
                )
                .interactive()
            )

            windVelYLine = (
                gridBase.mark_line(stroke="#0000FF", interpolate="monotone")
                .encode(
                    alt.X(
                        "average(windVelY)",
                        axis=alt.Axis(
                            title="Wind V (m/s)",
                            titleColor="#0000FF",
                        ),
                    ),
                    tooltip=["grid", "windVelY"],
                )
                .interactive()
            )

            chart2 = alt.layer(windVelXLine, windVelYLine)
            st.altair_chart(chart2, use_container_width=True)

        with ampCol2:

            densityLine = (
                gridBase.mark_line(stroke="#FF0000", interpolate="monotone")
                .encode(
                    alt.X(
                        "average(density)",
                        axis=alt.Axis(
                            title="Density (kg/m3)",
                            titleColor="#FF0000",
                        ),
                    ),
                    tooltip=["grid", "density"],
                )
                .interactive()
            )
            speedOfSoundLine = gridBase.mark_line(
                stroke="#0000FF",
                interpolate="monotone",
            ).encode(
                alt.X(
                    "average(speedOfSound)",
                    axis=alt.Axis(
                        title="Speed of Sound (m/s)",
                        titleColor="#0000FF",
                    ),
                ),
                tooltip=["grid", "speedOfSound"],
            )
            # Fix Loess Transform for speed of sound line
            # loess = gridBase.transform_loess("grid", "speedOfSound").mark_line(
            #     color="#0000FF"
            # )

            chart3 = alt.layer(densityLine, speedOfSoundLine).resolve_scale(
                x="independent"
            )
            # chart3 = (chart3 + loess).interactive()
            st.altair_chart(chart3, use_container_width=True)

            pressureLine = (
                gridBase.mark_line(stroke="#FF0000", interpolate="monotone")
                .encode(
                    alt.X(
                        "pressure",
                        axis=alt.Axis(
                            title="Pressure (hPa)",
                            titleColor="#FF0000",
                        ),
                    ),
                    tooltip=["grid", "pressure"],
                )
                .interactive()
            )

            tempratureLine = (
                gridBase.mark_line(stroke="#0000FF", interpolate="monotone")
                .encode(
                    alt.X(
                        "average(temperature)",
                        axis=alt.Axis(
                            title="Temprature (K)",
                            titleColor="#0000FF",
                        ),
                    ),
                    tooltip=["grid", "temperature"],
                )
                .interactive()
            )

            chart4 = alt.layer(pressureLine, tempratureLine).resolve_scale(
                x="independent"
            )
            st.altair_chart(chart4, use_container_width=True)
# {
#     "grav": 9.80665,
#     "launch_rail_length": 6,
#     "elevation": 1.019374966621399,
#     "modelType": "Forecast",
#     "modelTypeMaxExpectedHeight": "78862.84",
#     "windSpeed": 7.686906814575195,
#     "windDirection": 236.54574584960938,
#     "windHeading": 56.545753479003906,
#     "surfacePressure": 1000,
#     "surfaceTemperature": 284.6130676269531,
#     "surfaceAirDensity": 1.2240052212166865,
#     "surfaceSpeedOfSound": 338.1990494702038,
#     "launch_date": "2022-12-09 12:00:00",
#     "lat": 65,
#     "lon": 35,
#     "initDate": "2022-12-09 00:00:00",
#     "endDate": "2022-28-09 00:00:00",
#     "interval": 3,
#     "initLat": -90,
#     "endLat": 90,
#     "initLon": 0,
#     "endLon": 359.75,
# }
# Launch Site Details

# Launch Rail Length: 5.2  m
# Launch Date: 2022-09-27 12:00:00 UTC
# Launch Site Latitude: 32.99025°
# Launch Site Longitude: -106.97500°
# Reference Datum: SIRGAS2000
# Launch Site UTM coordinates: 315468.64 W    3651938.65 N
# Launch Site UTM zone: 13S
# Launch Site Surface Elevation: 1471.5 m


# Atmospheric Model Details

# Atmospheric Model Type: Forecast
# Forecast Maximum Height: 79.731 km
# Forecast Time Period: From  2022-09-26 06:00:00  to  2022-10-12 06:00:00  UTC
# Forecast Hour Interval: 3  hrs
# Forecast Latitude Range: From  -90.0 ° To  90.0 °
# Forecast Longitude Range: From  0.0 ° To  359.75 °


# Surface Atmospheric Conditions

# Surface Wind Speed: 0.97 m/s
# Surface Wind Direction: 158.49°
# Surface Wind Heading: 338.49°
# Surface Pressure: 858.87 hPa
# Surface Temperature: 293.07 K
# Surface Air Density: 1.021 kg/m³
# Surface Speed of Sound: 343.19 m/s
