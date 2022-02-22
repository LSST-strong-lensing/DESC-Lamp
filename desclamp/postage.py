# A few common packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
# We will use astropy's WCS and ZScaleInterval for plotting
from astropy.wcs import WCS
# Also to convert sky coordinates

import galsim
# We will use several stack functions
import lsst.geom
import lsst.afw.display.rgb as rgb
# Source injection
from lsst.pipe.tasks.insertFakes import _add_fake_sources

# And also DESC packages to get the data path
import GCRCatalogs
from GCRCatalogs import GCRQuery
import desc_dc2_dm_data


def catalog_setup(dc2_data_version):
        """ This function fetches catalogs and sets up a butler.

        Parameters
        ----------
        dc2_data_version: str
            data version of the catalog. This is used to instantiate the butler.
        """
        # Fetch GCR catalogs
        GCRCatalogs.get_available_catalogs(names_only=True, name_contains=dc2_data_version)
        cat = GCRCatalogs.load_catalog("dc2_object_run"+dc2_data_version)

        # Sets up butler
        butler = desc_dc2_dm_data.get_butler(dc2_data_version)

        return cat, butler

class Cutout:
    """A class that describes cutout of lens candidates and all the catalog level
    information necessary to identify thhe object as well as lensing-relatted inforrmation."""
    def __init__(self, exposure, catalog):
        self.exposure = exposure
        self.catalog = catalog

    def inject(self, lensed_source, spectra):
        """ A method to do synthetic injection of a lensed source in the cutout
        Parameters
        ----------
        lensed_source: a galsim object or an array
            An image of a lensed source to inject.
        """
        assert len(spectra)==len(self.exposure)
        radec = lsst.geom.SpherePoint(self.catalog["ra"], self.catalog["dec"], lsst.geom.degrees)
        new_exp = copy.deepcopy(self.exposure)
        lensed_obj = galsim.InterpolatedImage(lensed_source, scale = 0.05)
        for i,e in enumerate(new_exp):
            _add_fake_sources(e, [(radec, lensed_obj.withFlux(spectra[i]))])

        return Cutout(new_exp, self.catalog)



class Candidates:
    """ Class that handles catalog querries. Fetches postage stamps and light curves of samples of images and allows visualization """

    def __init__(self, dc2_data_version, skymap='deepCoadd_skyMap'):

        self.cat, self.butler = catalog_setup(dc2_data_version)
        self.skymap = self.butler.get(skymap)


    def catalog_query(self, query, columns = None, tracts = None):
        """ Submits a query and a selection to the catalog. Extract relevant information from catalogs.

        Parameters
        ----------
        query: str
            a set of selection criteria to select galaxy objects
        tracts: list
            list of tract numbers. Used to restrict the search to a small number of tracts.
        """
        # The minimum set of infomation needed about objects in the catalog
        # This will need to included lensing information att some point.

        columns_to_get = ["objectId", "ra", "dec", "tract", "patch"]
        if columns is not None:
            columns_to_get += columns
            columns_to_get = np.unique(columns_to_get)

        assert self.cat.has_quantities(columns_to_get)

        # Submit the query and get catalog of objects
        if tracts is not None:
            filters = f"(tract == {tracts[0]})"
            for t in tracts[1:]:
                filters +=  f" | (tract == {t})"
        objects = self.cat.get_quantities(columns_to_get, filters=GCRQuery(*query), native_filters=filters)

        # make it a pandas data frame for the ease of manipulation.
        # Objects are nont made attributes of the class in case the user wants postage stamps for a smaller set of objects
        objects = pd.DataFrame(objects)
        return objects

    def make_postage_stamps(self, objects, cutout_size=100, bands = 'irg'):
        """ Extracts a coadd postage stamp of an object from the catalog

        Parameters
        ----------
        objects: list of GCR catalog entries
            object to cut a postamp out.
        cutout_size: int
            size of the postage stamp to extract in pixels
        bands: str
            spectral for which patches have to extracted. Default is 'irg'.

        """

        cutout_extent = lsst.geom.ExtentI(cutout_size, cutout_size)
        cutouts = []
        for (_, object_this) in objects.iterrows():
            radec = lsst.geom.SpherePoint(object_this["ra"], object_this["dec"], lsst.geom.degrees)

            center = self.skymap.findTract(radec).getWcs().skyToPixel(radec)
            bbox = lsst.geom.BoxI(lsst.geom.Point2I((center.x - cutout_size*0.5, center.y - cutout_size*0.5)), cutout_extent)

            exposure = [self.butler.get("deepCoadd_sub",
                              bbox=bbox,
                              tract=object_this["tract"],
                              patch=object_this["patch"],
                              filter=band
                             ) for band in bands]

            new_cutout = Cutout(exposure, object_this)
            cutouts.append(new_cutout)

        return cutouts

    def display_cutouts(self, cutouts, figsize=(10,10), data_range = 2, q = 8):
        """ Displays RGB image of cutouts on a mosaic
        """
        n = len(cutouts)

        fig = plt.figure(figsize=figsize, dpi=100)
        if int(np.sqrt(n))**2 == int(n):
            l = int(np.sqrt(n))
        else:
            l = int(np.sqrt(n)+1)

        gs = plt.GridSpec(l, l, fig)

        for i in range(n):
            cutout = cutouts[i]
            image = cutout.exposure
            image_rgb = rgb.makeRGB(*image, dataRange = data_range, Q=q)

            ax = plt.subplot(gs[i],
                             #projection=WCS(cutout.exposure[0].getWcs().getFitsMetadata()),
                             label=str(cutout.catalog["objectId"])
                            )
            ax.imshow(image_rgb, origin='lower')
            del image_rgb  # let gc save some memory for us

            for c in ax.coords:
                c.set_ticklabel(exclude_overlapping=True, size=10)
                c.set_axislabel('', size=0)

        pass