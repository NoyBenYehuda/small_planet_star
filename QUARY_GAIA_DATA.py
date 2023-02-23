import numpy as np
import astropy.io.votable as votable
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
import multiprocessing as mp
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import time

start_time = time.time()
#first we will upload the table of GAIA ID's we found
file = #PATH HERE
table = votable.parse_single_table(file)

# convert the ID's from VOT to list
source_id = []
for row in table.array:
    id = row['source_id']
    source_id.append(str(id))  

 #we will write a GAIA archieve quary for the data we want
Id_list = [f"'{id}'" for id in source_id]
query = f"""
SELECT 
    source_id, ra, dec, phot_g_mean_mag, bp_rp, radial_velocity, vbroad, teff_gspphot, logg_gspphot, mh_gspphot
FROM
    gaiadr3.gaia_source
WHERE 
    source_id IN ({','.join(Id_list)})
"""

# we will quary GAIA on-line
q = Gaia.launch_job_async(query, output_format='votable')
table = q.get_results()


# Add a new columns for data based on the photometry data
distance = 1.0  # Assume 1 kpc distance
g=table['phot_g_mean_mag']
abs_mag = table['phot_g_mean_mag'] - 5 * (np.log10(distance) - 1)
luminosity = 10 ** (-0.4 * (abs_mag - 4.77))
table['luminosity'] = luminosity
tef=table['teff_gspphot']
# Add a new column for effective temperature based on the color index
table['teff'] = 4600 * (1 / (.92 * table['bp_rp'] + 1.7) + 1 / (.92 * table['bp_rp'] + 0.62))
mh=table['mh_gspphot']
#eExport the table
table.write('gaia_data.fits', overwrite=True)
rv=table['radial_velocity']

import numpy as np
from astropy.table import Table
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# Load the data from the FITS file
gaia_data = Table.read('gaia_data.fits')

# Remove any rows with missing values
gaia_data = gaia_data.filled(np.nan)
gaia_data = gaia_data[~np.isnan(gaia_data['teff'])]
gaia_data = gaia_data[~np.isnan(gaia_data['luminosity'])]

# Split the data into training and testing sets
X = gaia_data['ra', 'dec', 'phot_g_mean_mag', 'bp_rp', 'radial_velocity', 'teff', 'mass']
y = gaia_data['luminosity']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Standardize the data
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)
from sklearn.neural_network import MLPRegressor

# Train a multi-layer perceptron regressor on the training data
mlp = MLPRegressor(hidden_layer_sizes=(100, 100), max_iter=500, random_state=42)
mlp.fit(X_train, y_train)
from sklearn.metrics import mean_squared_error

# Evaluate the performance of the model on the testing data
y_pred = mlp.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
print(f"Mean squared error: {mse:.2f}")



fig, ax = plt.subplots()
im = ax.scatter(tef, g, c=mh, cmap='plasma')

# Set the color bar
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('[m/H]')

# Set the labels and title
ax.set_xlabel('T_eff (Solar Masses)')
ax.set_ylabel('g (Solar Luminosities)')
ax.set_title('G vs. T_eff as dependece of the metallicity')
# add annotations for star masses
#for i in range(len(rv)):
 #   ax.annotate(f'{rv[i]} M$_\odot$', (tef[i], luminosity[i]))
plt.show()




end_time = time.time()
total_time = end_time - start_time
print(f"Execution time: {total_time:.2f} seconds")
