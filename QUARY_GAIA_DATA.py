import numpy as np
import astropy.io.votable as votable
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
import multiprocessing as mp
import matplotlib.pyplot as plt
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

end_time = time.time()
total_time = end_time - start_time
print(f"Time passed: {total_time:.2f} seconds")
