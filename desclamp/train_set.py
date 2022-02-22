from .postage import Candidates
import numpy as np

class TrainSet(Candidates):
    """ Generates training sets for strong gravitational lens images.
    Lensed images are generated using source injection of lensed sources.
    TODO: lensed sources have to be matched to the characteristic of the
    lens galaxy to obtain realistic images.

    Parameters
    ----------
    dc2_data_version: str
        name of the catalog version. For now only for DC2 catalogs, but to be expanded.
    n_samples: int
        number of samples of images to generate.
    batchsize: int
        number of images generated at once in a given batch.
    lens_fraction: float
        the fraction of images that contain lensed-injected features.
    """

    def __init__(self,
                 dc2_data_version,
                 n_samples,
                 batchsize=8,
                 lens_fraction=0.2
                ):

        self.objects = []
        self.n = n_samples
        self.batchsize = batchsize
        self.index = 0
        self.lens_fraction = lens_fraction
        super().__init__(self, dc2_data_version)

    def catalog_query(self, query, columns = None, tracts = None):
        objects = super().catalog_query(query, columns = columns, tracts = tracts)
        self.objects.append(objects)
        return objects

    def __iter__(self):
        return self

    def __next__(self):
        if self.num >= self.batchsize:
            raise StopIteration()
        new_index = np.max([self.index+self.batch_size, self.n_samples])
        cutouts = self.make_postage_stamps(self.objects[self.index, new_index],
                                      cutout_size=100,
                                      bands = 'irg')
        return cutouts



