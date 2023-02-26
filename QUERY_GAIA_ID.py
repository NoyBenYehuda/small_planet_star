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


# Load table(I download from NASA planetary archieve)
filename = #PUT FILE NAME HERE
table = votable.parse_single_table(filename)

# Extract the data from the table
data = []
#we will extract the star's name, RA and DEC
for row in table.array:
    star_name = row['hostname']
    ra = row['rastr']
    dec = row['decstr']
    data.append((star_name, ra, dec))

#Searching for GAIA ID'S using coordinates
def coo_query(data):
# Set the search radius and Gaia DR3 catalog name
   good_id=[]
    #the radius to search around
   search_radius = 1.0 * u.arcsec
   gaia_dr3 = 'gaiaedr3.gaia_source'
# Loop over the list of stars and query Gaia DR3 for each one
   for star in data:
    # the RA and DEC coordinates from list
         name, ra, dec = star
         coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

    # Query Gaia DR3 for the closest match to the coordinates
         job = Gaia.cone_search_async(coords, search_radius, table_name=gaia_dr3)
         result = job.get_results()

    # Extract the Gaia IDs from the query results
         if len(result) > 0:
             gaia_id = result['source_id'][0]
             good_id.append(gaia_id)
         else:
             print(f'{name}: NULL')
    #In case of repeating IDs we wil use only one
   gaia_id_list = list(set(good_id))
   return gaia_id_list


#I HAVE ABOUT 2500 STARS- I made the run faster using parallel prossecing
if __name__ == '__main__':
    mp.set_start_method('spawn')  # set multiprocessing context here

    #1) create a pool of workers and use all the 8 processes in my computer in parallel
    num_processes = mp.cpu_count()
    with mp.Pool(num_processes) as pool:
        results = pool.apply_async(coo_query, (data,))

        # 2)get the results
        output=results.get()
        # close the pool
        pool.close()
        pool.join()
        # export the results to table
        tbl = Table(rows=[(id,) for id in output], names=('source_id',))

        # make a new table with GAIA ID's
        filename = r'C:\Users\avib-\Desktop\research\my tables\table.vot'

        # Write the table to the VOT format
        votable.writeto(tbl, filename)

        end_time = time.time()
        total_time = end_time - start_time
        print(f"Time passed: {total_time:.2f} seconds")




