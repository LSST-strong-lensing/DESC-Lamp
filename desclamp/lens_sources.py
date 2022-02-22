# standard python imports
import numpy as np
import matplotlib.pyplot as plt

# lenstronomy utility functions
import lenstronomy.Util.util as util
import lenstronomy.Util.image_util as image_util
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.Data.imaging_data import ImageData
import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.psf import PSF

import galsim


class Lensing_frame:
    """ Object to described lensed sources"""
    def __init__(self, shape=(100,100), pix = 0.2, wcs = None, hr_factor=1):
        """
        Source object that carries source and lens information and generates lensed source images for injection.
        Parameters
        ----------
        shape: `tuple`
            Shape of the image patch
        pix: `float`
            pixel size in arcseconds. Default is Rubin's pixel.
        wcs: WCS
            wcs information for the source patch.  
        hr_factor: `float`
            the high resolution factor between source and lens plane
        """
        if pix is None: 
            assert wcs is not None
            try:
                model_affine = wcs.wcs.pc
            except AttributeError:
                model_affine = wcs.cd
            pix = np.sqrt(np.abs(model_affine[0, 0])* np.abs(model_affine[1, 1] - model_affine[0, 1] * model_affine[1, 0]))
            
        self.pix = pix
        self.hr_factor = hr_factor
        self.shape = shape
        self._source = None # An image of the source
        self.source_args = None # source arguments dictionary
    
    def lens_source(self, lens_models, lens_args):
        assert self.source is not None, "Please provide a source image."

        kwargs_image = sim_util.data_configure_simple(int(self.shape[0]),  # Issue with rectangle shapes
                                                              self.pix)
        image = ImageData(**kwargs_image)
        
        lensModel = LensModel(lens_models)
        lightModel = LightModel(light_model_list=['INTERPOL'])
        
        imageModel = ImageModel(data_class=image, 
                                         lens_model_class=lensModel, 
                                         source_model_class=lightModel,
                                         psf_class=PSF(psf_type='NONE'),
                                         kwargs_numerics = {'supersampling_factor': self.hr_factor, 
                                                            'supersampling_convolution': False})
        lensed_image = imageModel.image(lens_args, self.source_args)
        return lensed_image
    
    def draw_source(self):
        """ Draws a soource no a grid specified by the parameters of the __init__
        """
        kwargs_list = [
            {'image': self.source, 'scale': self.pix/self.hr_factor, 'phi_G': 0, 'center_x': 0, 'center_y': 0}
        ]
        
        return kwargs_list
    
    @property
    def source(self):
        return self._source
    
    @source.setter
    def source(self, source):
        self._source = source
        self.source_args = self.draw_source()
    
    def from_gsobject(self, gsobject, smooth=0):
        """ Creates a Lensed source from a galsim object.
        Parameters
        ----------
        gsobject: Galsim Object
            image of a galaxy to use in as a source for lensing.
        smooth: bool
            Value of the sigma for a gaussian smoothing kernel. Useful to apply if the input image is noisy or contains sharp features.
        """
        if smooth > 0:
            gso = galsim.Convolve(gsobject, galsim.Gaussian(sigma=smooth))
        #Draws the galsim object on a grid    
        source = gso.drawImage(nx=self.shape[0]*self.hr_factor,
            ny=self.shape[1]*self.hr_factor,
            use_true_center = True,
            method='real_space',
            scale=self.pix/self.hr_factor,
            dtype=np.float64).array
        self.source = source
        return self
        
    def from_lenstronomy(self, kwargs):
        """TO DO"""
        pass
    
    def from_galsim_parametric(self, *args, profile='Spergel', shift=(0,0), shear=(0,0), **kwargs):
        """ Creates a Lensed source object from a source described as a parametric profile.
        Parameters
        ----------
        profile: string
            The galsim profile to use to generate a source galaxy. Currently, Sersic, DeVeaucouleurs and Exponential are implemented
        kwargs:
            parameters of the profile
            
        Returns
        -------
        A `Lensed_source` object.
        """
        
        assert profile in ['Sersic', 'Exponential', 'DeVeaucouleurs', 'Spergel'], "Not a valid profile. Please use 'Sersic', 'Exponential' or 'DeVeaucouleurs'."
        
        if profile == 'Sersic':
            assert len(args)==1
            gso = galsim.Sersic(n= args[0], **kwargs)
        elif profile == 'Spergel':
            assert len(args)==1
            gso = galsim.Spergel(nu= args[0], **kwargs)
        elif profile == 'Exponential':
            gso = galsim.Exponential(*args, *kwargs)
        else:
            gso = galsim.DeVeaucouleurs(*args, *kwargs)
            
        gso = gso.shear(g1=shear[0],g2=shear[1])
        gso = gso.shift(dx=shift[0],dy=shift[1])
        
        #Draws the galsim object on a grid    
        source = gso.drawImage(nx=self.shape[0]*self.hr_factor,
            ny=self.shape[1]*self.hr_factor,
            use_true_center = True,
            method='real_space',
            scale=self.pix/self.hr_factor,
            dtype=np.float64).array
        self.source = source
        return self
    
    
    
    
    