# Simple Algorithm for Evapotranspiration Retrieving (SAFER)
The **SAFER** a new model created in semi-arid parts of Brazil, has been validated using field
data from four flux stations with irrigated crops and natural vegetation (Teixeira 2010). The
model was developed to address the need for real-time operational ET that did not require
anchor pixels, as opposed to **SEBAL** and **METRIC**, to operationalize real ET modeling.It is a simple biophysically realistic model that works well for spatially
mapping surface ET (Santos et al. 2020). **SAFER** needs information from the satellite images,
including the normalized difference vegetation index (NDVI), albedo (α), and Ts. Using suchinformation, the surface ET fraction is determined, which when combined with the reference
evapotranspiration (ETo), enables the mapping of ET using **SAFER** model (ET<sub>SAFER</sub>). In
addition, the energy and water balance components can be calculated both with and without
using the satellite thermal bands (Silva et al. 2019).
This is among very few ET modeling in India. Other than India, most of the work in the **SAFER** model domain has been done for other parts of the world. This study provides the feasibility of the **SAFER** model for ET studies. Several studies on the temporal variability of ET have been carried out in India. To ascertain the decadal trend of ET for the Himalayan Shivalik Mountain range, this study was conducted.
                                                                                          <img width="400" alt="Screenshot 2024-12-12 170940" src="https://github.com/user-attachments/assets/4a8393f0-58dd-4543-9919-7a023c23fed2">

where,
- **α0** represents the surface albedo (dimensionless)
- **NDVI** represents the Normalized Difference Vegetation Index (dimensionless)
- **Ts** represents the surface temperature (<sup>o</sup>C)
- **a** and **b** are *1.8* and *-0.0008* respectively (De Oliveira et al., 2019; Leivas et al., 2015).

The FAO-56 Penman-Monteith technique is commonly applied in the parameterization of **SAFER**. Fewer meteorological variables be required for estimating ETo. Santos et al., (2020) examined ET<sub>SAFER</sub> responses by calculating the ETo determined using the Hargreaves-Samani (HS) and PM methods in the presence of missing data.

## FAO-56-PM
The reference evapotranspiration (ET<sub>o</sub>) represents the evapotranspiration from a hypothesized reference crop (height of 0.12 m, surface resistance of 70 s/m and albedo of 0.23).
The universal standard FAO-56 Penman-Monteith model is a combination-based model that has been approved to be applicable in different climatic conditions without need to a local calibration. The Food FAO-56 Penman-Monteith is considered the standard method to calculate ET<sub>o</sub>.The FAO 56 PM method requires measurements of temperature, relative humidity, wind speed, and solar radiation.
The FAO-56 Penman-Monteith formula for daily ET<sub>o</sub> (mm/ day<sup>-1</sup>) :

<img width="400" alt="Screenshot 2024-03-07 161356" src="https://github.com/AthiraNG/FAO-56-PM/assets/129937610/3fde2d5d-5f1e-4104-bb65-e63f2edb0734">


  Here,
- **ET<sub>o</sub>** denotes the reference evapotranspiration rate (mm day<sup>-1</sup>)
- **Δ** represents the slope of the vapor pressure curve (kPa <sup>∘</sup>C<sup>−1</sup>)
- **Rn** is the net radiation (MJ m<sup>−2</sup> day<sup>−1</sup>)
- **G** is the soil heat flux density (MJ m<sup>−2</sup> day<sup>−1</sup>)
- **𝛾** is the psychrometric constant (kPa <sup>∘</sup>C<sup>−1</sup>)
- **Tmean** is the mean daily air temperature measured at 2 m (<sup>∘</sup>C)
- **u2** is the wind speed at 2 m (ms<sup>−1</sup>)
- **es** is the saturation vapor pressure (kPa)
- **ea** is the actual vapor pressure (kPa)
