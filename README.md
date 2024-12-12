# Simple Algorithm for Evapotranspiration Retrieving (SAFER)
The **SAFER** a new model created in semi-arid parts of Brazil, has been validated using field
data from four flux stations with irrigated crops and natural vegetation (Teixeira 2010). The
model was developed to address the need for real-time operational ET that did not require
anchor pixels, as opposed to **SEBAL** and **METRIC**, to operationalize real ET modeling.It is a simple biophysically realistic model that works well for spatially
mapping surface ET. **SAFER** needs information from the satellite images,
including the normalized difference vegetation index (NDVI), albedo (α), and Ts. Using such information, the surface ET fraction is determined, which when combined with the reference
evapotranspiration (ETo), enables the mapping of ET using **SAFER** model (ET<sub>SAFER</sub>).

In addition, the energy and water balance components can be calculated both with and without
using the satellite thermal bands (Silva et al. 2019).
This is among very few ET modeling in India. Other than India, most of the work in the **SAFER** model domain has been done for other parts of the world. This study provides the feasibility of the **SAFER** model for ET studies. Several studies on the temporal variability of ET have been carried out in India. To ascertain the decadal trend of ET for the Himalayan Shivalik Mountain range, this study was conducted.
                                                                                          <img width="400" alt="Screenshot 2024-12-12 170940" src="https://github.com/user-attachments/assets/4a8393f0-58dd-4543-9919-7a023c23fed2">

where,
- **α0** represents the surface albedo (dimensionless)
- **NDVI** represents the Normalized Difference Vegetation Index (dimensionless)
- **Ts** represents the surface temperature (<sup>o</sup>C)
- **a** and **b** are *1.8* and *-0.0008* respectively.

The FAO-56 Penman-Monteith technique is commonly applied in the parameterization of **SAFER**. Fewer meteorological variables be required for estimating ETo. Santos et al., (2020) examined ET<sub>SAFER</sub> responses by calculating the ETo determined using the Hargreaves-Samani (HS) and PM methods in the presence of missing data.

 ### Output
 
<img width="4000" alt="safer" src="https://github.com/user-attachments/assets/43f25978-8103-45a1-af20-8300f9855bc6">

  ### Resources
  #### [Step by Step Calculation of the Penman-Monteith Evapotranspiration (FAO-56 Method)](https://www.agraria.unirc.it/documentazione/materiale_didattico/1462_2016_412_24509.pdf)
  #### [Comparison of SAFER and METRIC Based Actual Evapotranspiration Models](https://www.researchgate.net/publication/335924830_Comparison_of_SAFER_and_METRIC-Based_Actual_Evapotranspiration_Models_in_a_Subtropical_Area_of_Brazil)
  #### [Energy balance and irrigation performance assessments in lemon orchards by applying the SAFER algorithm to Landsat 8 images](https://www.sciencedirect.com/science/article/abs/pii/S0378377420322691)
