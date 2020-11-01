from astroquery.vizier import Vizier
from astropy import coordinates
from astropy import units as u

v = Vizier()

c = coordinates.SkyCoord(0, 0, unit=('deg', 'deg'), frame='icrs')
result = v.query_region(c, radius=2*u.deg, catalog = ["I/337","I/345"])

print(result)


#catalog_list = Vizier.find_catalogs('gaia', max_catalogs = 1000)
#print({k:v.description for k,v in catalog_list.items()})
