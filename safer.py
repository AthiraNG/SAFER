#!/usr/bin/env python
# coding: utf-8

# In[2]:


pip install earthengine-api


# In[3]:


pip install geemap


# In[4]:


import ee
import geemap
import math


# In[5]:


ee.Authenticate()
ee.Initialize()
Map = geemap.Map()


# In[6]:


# Constants

radians = 0.50841

gsc = 0.0820 # solar constant MJ/m2

# Regression constant, intercept, expressing the fraction of extra-terrestrial radiation reaching the earth on overcast days
As = ee.Number(0.23)
Bs = ee.Number(0.5)


# In[7]:


# Define the year and month
year = 2010
month = 2

# Define start and end dates
start_date = ee.Date.fromYMD(year, month, 1)  # First day of the month
end_date = start_date.advance(1, 'month')    # Advance by one month


# In[8]:


# Define region of interest (roi)
roi = ee.FeatureCollection('projects/ee-athiraganesan314/assets/Shivalik_range')

# Define roi for sampling 
aoi_minmax = ee.FeatureCollection('users/athiraganesan314/aoi_minmax')

# Tower location for Sampling
pts = ee.FeatureCollection('projects/ee-athiraganesan314/assets/tower_locations')


# In[10]:


# Load and process MODIS LST data
lst = (ee.ImageCollection('MODIS/061/MOD11A1')
       .filterDate(start_date, end_date)
       .filterBounds(roi)
       .select('LST_Day_1km')
       .mean()
       .clip(roi)
       .multiply(0.02))  # Convert to Celsius

# Load and process NDVI data
ndvi = (ee.ImageCollection('MODIS/061/MOD13A1')
        .filterDate(start_date, end_date)
        .filterBounds(roi)
        .select('NDVI')
        .median()
        .clip(roi)
        .multiply(0.0001))  # Scale factor

# Load and process ERA5 wind data
wind = (ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
        .filterDate(start_date, end_date)
        .filterBounds(roi)
        .select(['u_component_of_wind_10m', 'v_component_of_wind_10m'])
        .first()
        .clip(roi))

# Surface reflectance dataset and albedo calculation
dataset = (ee.ImageCollection('MODIS/061/MOD09GA')
           .filterDate(start_date, end_date)
           .filterBounds(roi))
bands_alb = (dataset.select(['sur_refl_b01', 'sur_refl_b02'])
             .mean()
             .multiply(0.0001)
             .clip(roi))
alb = bands_alb.expression(
    '0.08 + (0.41 * RED) + (0.14 * NIR)', {
        'NIR': bands_alb.select('sur_refl_b02'),
        'RED': bands_alb.select('sur_refl_b01')
    })

# Load and process Temp data
temp = (ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
        .filterDate(start_date, end_date)
        .filterBounds(roi)
        .select('temperature_2m')
        .first()
        .clip(roi))


# In[13]:


# Calculate crop coefficient (kc)
kc = (lst.subtract(273.15)  # Convert LST from Kelvin to Celsius
          .divide(ndvi.multiply(alb))  # LST / (NDVI * Albedo)
          .multiply(-0.008)
          .add(1.8)
          .exp())

# Sample Min and Max Temperature 
temp_minmax = temp.reduceRegions(
collection = aoi_minmax,
reducer = ee.Reducer.minMax()
)

# Extract the values of some environmental variable(wind,albedo,temperature)

stack=wind.addBands(alb) 
sampled= stack.sampleRegions(
  collection = pts,
  properties = ['point'],
  scale = 9000,
  geometries = True
  
)
mete_values = ee.FeatureCollection([sampled,temp_minmax]).flatten()

uwind= mete_values.aggregate_mean('u_component_of_wind_10m')
vwind= mete_values.aggregate_mean('v_component_of_wind_10m')
dew= mete_values.aggregate_mean('dewpoint_temperature_2m')
alb= mete_values.aggregate_mean('constant')
tmax= mete_values.aggregate_mean('max')
tmin= mete_values.aggregate_mean('min')

# Calculate other atmospheric and environmental variable for the calculation of ETo and ETa

tmean =(tmin.add(tmax)).divide(2)
uh=((uwind.pow(2)).add(vwind.pow(2))).sqrt()
u2=uh.multiply(0.7479511)
ssvp=((((((tmean.subtract(273.15)).multiply(17.27)).divide(((tmean.subtract(273.15)).add(237.3)))).exp()).multiply(0.6108)).multiply(4098)).divide((((tmean.subtract(273.15)).add(237.3)).pow(2)))
p= (ee.Number((293-(0.0065*450))/(293)).pow(ee.Number(5.26))).multiply(101.3)
r= p.multiply(0.000665)
etmax=((((tmax.subtract(273.15)).multiply(17.27)).divide(((tmax.subtract(273.15)).add(237.3)))).exp()).multiply(0.6108)
etmin=((((tmin.subtract(273.15)).multiply(17.27)).divide(((tmin.subtract(273.15)).add(237.3)))).exp()).multiply(0.6108)  
ea= ((((tmin.subtract(273.15)).multiply(17.27)).divide(((tmin.subtract(273.15)).add(237.3)))).exp()).multiply(0.6108)
es= (etmax.add(etmin)).divide(2)
dr =(((ee.Number(((math.pi*2)/365)*1)).cos()).multiply(0.033)).add(1)
sd= (ee.Number(((math.pi*2)/365)-(1.39)).sin()).multiply(0.409)
ws= (((ee.Number(-radians)).tan()).multiply(ee.Number(sd).tan())).acos()
Ra=((ws.multiply(ee.Number(radians).sin()).multiply(sd.sin())).add((ee.Number(radians).cos()).multiply(sd.cos()).multiply(ws.sin()))).multiply(dr).multiply(gsc).multiply(1440/(math.pi))
Rso=Ra.multiply(0.759)
N=(ee.Number(24).divide(math.pi)).multiply(ws)
n=tmean.subtract(273.15).multiply(0.232).add(4.352)
Rs=(ee.Number(n.divide(N))).multiply(Bs).add(As.multiply(Ra))
Rns=((ee.Number(1)).subtract(alb)).multiply(Rs)
Rnl=((((((tmax.subtract(273.15)).add(273.16)).pow(4)).add(((tmin.subtract(273.15)).add(273.16)).pow(4))).divide(2)).multiply(0.000000004903)).multiply(ee.Number(0.34).subtract((ea.sqrt().multiply(0.14)))).multiply(ee.Number(1.35).multiply(Rs.divide(Rso)).subtract(0.35))
Rn=Rns.subtract(Rnl)


# Reference Evapotranspiration (FAO 56)

ETo = (
    ((ee.Number(0.408).multiply(Rn).multiply(ssvp))
     .divide((((u2.multiply(0.34)).add(1)).multiply(r)).add(ssvp)))
    .add(
        ((ee.Number(900).divide((tmean.subtract(273.15)).add(273)))
         .multiply(u2))
        .multiply((es.subtract(ea)).multiply(r))
        .multiply(
            (((ee.Number(0.34).multiply(u2)).add(1))
             .multiply(r))
            .add(ssvp)
        )
    )
)

# Actual Evapotranspiration

ETa =kc.multiply(ETo)

Map.addLayer(ETa,vis_params,"ETa")
Map.centerObject(roi,12)
Map.addLayer(roi,{},"roi")

display(Map)


# In[12]:


vis_params = {'bands': ['LST_Day_1km'], 'palette': ['#5e4fa2', '#5c51a3', '#5b53a4', '#5956a5', '#5758a6', '#555aa7', '#545ca8', '#525fa9', '#5061aa', '#4e63ac', '#4d65ad', '#4b68ae', '#496aaf', '#486cb0', '#466eb1', '#4471b2', '#4273b3', '#4175b4', '#3f77b5', '#3d79b6', '#3b7cb7', '#3a7eb8', '#3880b9', '#3682ba', '#3585bb', '#3387bc', '#3389bd', '#358bbc', '#378ebb', '#3990ba', '#3b92b9', '#3d95b8', '#3f97b7', '#4199b6', '#439bb5', '#459eb4', '#47a0b3', '#49a2b2', '#4ba4b1', '#4ea7b0', '#50a9af', '#52abae', '#54aead', '#56b0ad', '#58b2ac', '#5ab4ab', '#5cb7aa', '#5eb9a9', '#60bba8', '#62bda7', '#64c0a6', '#66c2a5', '#69c3a5', '#6bc4a5', '#6ec5a5', '#71c6a5', '#74c7a5', '#76c8a5', '#79c9a5', '#7ccaa5', '#7ecca5', '#81cda5', '#84cea5', '#86cfa5', '#89d0a4', '#8cd1a4', '#8fd2a4', '#91d3a4', '#94d4a4', '#97d5a4', '#99d6a4', '#9cd7a4', '#9fd8a4', '#a2d9a4', '#a4daa4', '#a7dba4', '#aadca4', '#acdda4', '#aedea3', '#b1dfa3', '#b3e0a2', '#b5e1a2', '#b8e2a1', '#bae3a1', '#bce4a0', '#bfe5a0', '#c1e6a0', '#c3e79f', '#c6e89f', '#c8e99e', '#caea9e', '#cdeb9d', '#cfec9d', '#d1ed9c', '#d3ed9c', '#d6ee9b', '#d8ef9b', '#daf09a', '#ddf19a', '#dff299', '#e1f399', '#e4f498', '#e6f598', '#e7f59a', '#e8f69b', '#e9f69d', '#eaf79e', '#ebf7a0', '#ecf7a1', '#edf8a3', '#eef8a4', '#eff9a6', '#f0f9a7', '#f1f9a9', '#f2faaa', '#f3faac', '#f4faad', '#f5fbaf', '#f6fbb0', '#f7fcb2', '#f8fcb4', '#f9fcb5', '#fafdb7', '#fbfdb8', '#fcfeba', '#fdfebb', '#fefebd', '#ffffbe', '#fffebe', '#fffdbc', '#fffcba', '#fffbb8', '#fffab6', '#fff8b4', '#fff7b2', '#fff6b0', '#fff5ae', '#fff3ac', '#fff2aa', '#fff1a8', '#fff0a6', '#feefa3', '#feeda1', '#feec9f', '#feeb9d', '#feea9b', '#fee999', '#fee797', '#fee695', '#fee593', '#fee491', '#fee28f', '#fee18d', '#fee08b', '#fede89', '#fedc88', '#feda86', '#fed884', '#fed683', '#fed481', '#fed27f', '#fed07e', '#fece7c', '#fecc7b', '#feca79', '#fec877', '#fdc776', '#fdc574', '#fdc372', '#fdc171', '#fdbf6f', '#fdbd6d', '#fdbb6c', '#fdb96a', '#fdb768', '#fdb567', '#fdb365', '#fdb163', '#fdaf62', '#fdad60', '#fcaa5f', '#fca85e', '#fca55d', '#fba35c', '#fba05b', '#fb9d59', '#fa9b58', '#fa9857', '#fa9656', '#f99355', '#f99153', '#f98e52', '#f88c51', '#f88950', '#f8864f', '#f7844e', '#f7814c', '#f67f4b', '#f67c4a', '#f67a49', '#f57748', '#f57547', '#f57245', '#f47044', '#f46d43', '#f36b43', '#f26944', '#f06744', '#ef6645', '#ee6445', '#ed6246', '#eb6046', '#ea5e47', '#e95c47', '#e85b48', '#e75948', '#e55749', '#e45549', '#e3534a', '#e2514a', '#e1504b', '#df4e4b', '#de4c4b', '#dd4a4c', '#dc484c', '#da464d', '#d9444d', '#d8434e', '#d7414e', '#d63f4f', '#d43d4f', '#d23a4e', '#d0384e', '#cd364d', '#cb334d', '#c9314c', '#c72e4c', '#c52c4b', '#c32a4b', '#c1274a', '#be254a', '#bc2249', '#ba2049', '#b81e48', '#b61b48', '#b41947', '#b11747', '#af1446', '#ad1246', '#ab0f45', '#a90d45', '#a70b44', '#a40844', '#a20643', '#a00343', '#9e0142'], 'min': 0.25358117110408673, 'max': 1.234534275434069}


# In[ ]:




