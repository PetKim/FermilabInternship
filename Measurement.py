#Code is courtesy of Dr. A. Nikolakopoulos


import numpy as np
import sys

class Measurement:
    '''
    A class which stores a dataset obtained for a specific flux
    '''

    def __init__(self,E_step,flux,bins,data,E_min=0.):
        self.E_step = E_step
        self.E_min = E_min
        self.lenflux=len(flux)
        self.flux = np.zeros([self.lenflux])
        self.E_vec = np.zeros([self.lenflux])
        for iflx, flxval in enumerate(flux):
            self.flux[iflx] = flxval
            self.E_vec[iflx] = E_min + iflx*E_step

        binshape = bins.shape
        datshape = data.shape
        if binshape[0] != datshape[0]:
            raise ValueError("Dimension of bins and data are not the same")

        if binshape[1] != 4:
            raise ValueError("Bins should be defined by 4 arguments (cosmin, cosmax, Emin, Emax)")
        self.bins=np.zeros(binshape)
        self.data=np.zeros(datshape)

        self.nbins = len(self.data)

        #Hardcopy the bins and data
        for ibin, binvals in enumerate(bins):
            self.data[ibin] = data[ibin]
            self.bins[ibin,:] = bins[ibin,:] #I think this does a hard copy with the ":"

        



class LinComb:
    '''
    A class which stores the linear combination of measurments
    '''

    def __init__(self, Measurements, coeffs):
       
        self.Estep, self.E_min, self.lenflux, self.nbins = self.checkScalars(Measurements)

        self.nMeas = len(Measurements)
        self.Measurements = []
        
        #We make new objects for the measurements, I think this means they won't go out of scope contrary to cpp (test this!)
        for imeas, meas in enumerate(Measurements):
            flux = meas.flux
            bins = meas.bins
            data = meas.data
            self.Measurements.append(Measurement(self.Estep, flux, bins, data, E_min=self.E_min))
        
        if len(coeffs) != self.nMeas:
            raise ValueError("Length of combination coefficients not equal to number of measurements")
        self.bins=Measurements[0].bins
        self.coeffs = coeffs
    
    def checkScalars(self, Measurements):
        E_step = Measurements[0].E_step
        E_min = Measurements[0].E_min
        lenflux = Measurements[0].lenflux
        nbins = Measurements[0].nbins

        for i in range(1,len(Measurements)):
            if E_step != Measurements[i].E_step:
                raise ValueError("Energy step is not identical in different measurements")
            if E_min != Measurements[i].E_min:
                raise ValueError("Energy minimum is not identical in different measurements")
            if lenflux != Measurements[i].lenflux:
                raise ValueError("Flux length is not identical in different measurements")
            if nbins != Measurements[i].nbins:
                raise ValueError("Number of bins is not identical in different measurements")

        return E_step, E_min, lenflux, nbins

    def eval_flux(self):
        flux = np.zeros([self.lenflux])
        for icof,cof in enumerate(self.coeffs):
            flux += cof*self.Measurements[icof].flux
            
        return flux


    def eval_data(self):
        data = np.zeros([self.nbins])
        for icof,cof in enumerate(self.coeffs):
            data += cof*self.Measurements[icof].data
            
        return data

    def eval_Meas(self):
        flux = np.zeros([self.lenflux])
        data = np.zeros([self.nbins])
        for icof,cof in enumerate(self.coeffs):
            flux += cof*self.Measurements[icof].flux
            data += cof*self.Measurements[icof].data

        combi = Measurement(self.Estep, flux, self.bins, data, E_min=self.E_min)

        return combi
