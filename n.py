import numpy as np
import astropy.io.votable as votable
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
import concurrent.futures
import multiprocessing as mp
from multiprocessing import Pool
import time

start_time = time.time()


# Load the VOTable file

#filename = r'C:\Users\avib-\Desktop\research\my tables\star_names'
filename = r'C:\Users\avib-\Desktop\research\my tables\small_new.VOT'
table = votable.parse_single_table(filename)

# Extract the data from the VOTable and store it in a list of tuples
data = []
for row in table.array:
    star_name = row['hostname']
    ra = row['rastr']
    dec = row['decstr']
    data.append((star_name, ra, dec))


def coo_query(data):
# Set the search radius and Gaia DR3 catalog name
   good_id=[]
   search_radius = 1.0 * u.arcsec
   gaia_dr3 = 'gaiaedr3.gaia_source'
# Loop over the list of stars and query Gaia DR3 for each one
   for star in data:
    # Get the RA and DEC coordinates from the star list
         name, ra, dec = star
         coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

    # Query Gaia DR3 for the closest match to the coordinates
         job = Gaia.cone_search_async(coords, search_radius, table_name=gaia_dr3)
         result = job.get_results()

    # Extract the Gaia ID from the query results
         if len(result) > 0:
             gaia_id = result['source_id'][0]
        #print(f'{name}: Gaia ID = {gaia_id}')
             good_id.append(gaia_id)
         else:
             print(f'{name}: No match found in Gaia DR3')
   gaia_id_list = list(set(good_id))
   return gaia_id_list


#source_id=[704967037090946688, 6794047652729201024, 58200934326315136]


#def pool_handler():

 #   p = Pool(4)

  #  p.map(coo_query(data), data)


#if __name__ == '__main__':
 #   pool_handler()

# create a pool of workers
#pool = mp.Pool(processes=8)
#results = [pool.apply_async(coo_query(data,search_radius,gaia_dr3))]

# retrieve the results
#output = [p.get() for p in results]

# execute the function in parallel
#results = [pool.apply_async(coo_query(data, search_radius, gaia_dr3)) ]

# retrieve the results
#output = [p.get() for p in results]

# combine the results into a single table
#combined_table = Table.vstack(output)
#print(type(output))
#print(output)
# close the pool of workers
#pool.close()
#pool.join()
source_id=(coo_query(data))
#source_id=pool_handler()
#print(source_id)
#tbl = Table(rows=output, names=('source_id'))
tbl = Table(rows=[(x,) for x in source_id], names=('source_id',))

# Define the VOTable file name
filename = r'C:\Users\avib-\Desktop\research\my tables\table.vot'

# Write the table to the VOTable file
votable.writeto(tbl, filename)


# your code here

end_time = time.time()
total_time = end_time - start_time
print(f"Execution time: {total_time:.2f} seconds")

# Upload source IDs as a table
#table = Table([source_id], names=('source_id'))
#job = Gaia.launch_job_async("SELECT source_id, mass, radius FROM gaiaedr3.gaia_source WHERE source_id IN (SELECT source_id FROM tap_upload.my_table)", upload_resource=table)

# Wait for job to finish
#result = job.get_results()
query = """
SELECT 
    source_id, ra, dec, phot_g_mean_mag, bp_rp, radial_velocity
FROM
    gaiadr3.gaia_source
WHERE 
    source_id IN ({','.join(source_id)}

"""
#WHERE
 #   phot_g_mean_mag < 10
# Launch the async job and retrieve the results as an Astropy table
job = Gaia.launch_job_async(query, output_format='votable')
table = job.get_results()
# Print results
print(table)



# Define the ADQL query to retrieve data from Gaia DR2


# Add a new column for luminosity based on the photometry data
distance = 1.0  # Assume 1 kpc distance
abs_mag = table['phot_g_mean_mag'] - 5 * (np.log10(distance) - 1)
luminosity = 10 ** (-0.4 * (abs_mag - 4.77))
table['luminosity'] = luminosity

# Add a new column for effective temperature based on the color index
table['teff'] = 4600 * (1 / (.92 * table['bp_rp'] + 1.7) + 1 / (.92 * table['bp_rp'] + 0.62))

# Export the table to a FITS file
table.write('gaia_data.fits', overwrite=True)
