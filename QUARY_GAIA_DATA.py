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
    source_id, ra, dec, phot_g_mean_mag, bp_rp, radial_velocity, vbroad, teff_gspspec, logg_gspspec, mh_gspspec, l, b, alphafe_gspspec
FROM
    gaiadr3.gaia_source
WHERE 
    source_id IN ({','.join(Id_list)})
"""

# we will quary GAIA on-line
q = Gaia.launch_job_async(query, output_format='votable')
table = q.get_results()


end_time = time.time()
total_time = end_time - start_time
print(f"Time passed: {total_time:.2f} seconds")
