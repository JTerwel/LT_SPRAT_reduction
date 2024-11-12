'''
Extracts spectra from SPRAT pipeline and applies an airmass dependent correction to the flux
Required folders listed for input, output, and plots.
'''

import os
import glob
import pylab as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import argparse

### command line arguments

parser = argparse.ArgumentParser()
parser.add_argument('-z','--z', help='Galaxy lines to plot at z', type=float)
parser.add_argument('-C','--C', help='Cosmir ray removal. Warning! Will remove host lines', 
                    action='store_true')
args = parser.parse_args()
##################

if args.z == None:
    plotgal='n'
else:
    plotgal='y'
    z_gal = args.z

input_location = '/InputSpectra/'
output_location='/OutputSpectra/'
plot_save_location = '/plots/'
#####################
    
def CRReject(x,y, window=37):
    '''
    Cleans 1D spectra of cosmic rays (but also host galaxy lines).
    
    >>> CRReject(x,y, window=37)
    numpy.array(y)
    '''
    
    for j in range(len(x)-window):
        region = y[j:j+window+1]
        
        magic_index = int(j+np.ceil(window/2))
        
        if y[magic_index] > (np.median(region) + 2.5*np.std(region)):
            
            y[magic_index] = np.median(region)
    
    # deal with the beginning
    region = y[:window+1]   
    for idx in range(len(region)):
        if y[idx] > (np.median(region) + 5*np.std(region)):
            y[idx] = np.median(region)
            
    # deal with the beginning
    region = y[j-window:]   
    for idx in range(j-window, len(region)):
        if y[idx] > (np.median(region) + 5*np.std(region)):
            y[idx] = np.median(region)        
            
    return y 

##################

for dir in [input_location,
            output_location,
            plot_save_location,]:
            if os.path.isdir('.%s'%dir)==False:
                print('Making %s' %dir)
                os.mkdir('.%s'%dir)

ext=5

files = glob.glob('.'+input_location+'*.fits')
for file in files:
    hdul = fits.open(file)
    
    if len(hdul) > 5:
        ext=5
    else:
        ext=4
    header = hdul[ext].header
    spectrum = hdul[ext].data[0]
    SN = header['OBJECT']
    dateobs = header['DATE']
    mjd = header['MJD']
    airmass = header['AIRMASS']
    wcs = WCS(hdul[ext].header)
    #wavelength is CDELT*np.arange(NAXIS)+CRVAL
    NAXIS = header['NAXIS1']
    CDELT = wcs.wcs.cdelt[0]
    CRVAL = wcs.wcs.crval[0]
    wl = CDELT*np.arange(NAXIS)+CRVAL
    fl_raw = spectrum
    hdul.close()
    corrections= glob.glob('./CorrectionFiles/*.txt')
    clist = [ float(os.path.basename(c)[:-4])  for c in corrections]
    correction_index = np.argmin(abs(np.array(clist)-airmass))
    corloc = './CorrectionFiles/'+str(clist[correction_index])+'.txt'
    print('Correcting flux for airmass', airmass, 'using', str(clist[correction_index])+'.txt')
        
    s = np.loadtxt(corloc,unpack=True, usecols=(0,1))
    wl_corloc = s[0]
    f_corloc = s[1]

    fl_cor_temp=[]
    for q in range(len(wl)):
        find = np.argmin(abs(wl_corloc-wl[q]))
        fl_cor_temp.append(fl_raw[q]/f_corloc[find])
    fl_cor=np.array(fl_cor_temp)
    
    # Cosmic ray removal
    if args.C:
        fl_cor = CRReject(wl,fl_cor)

    file_name = SN+'_'+dateobs.replace('-', '')
    master=np.column_stack((wl, fl_cor))
    np.savetxt('.'+output_location+os.path.basename(file_name)+'.txt',master,fmt="%1.1f %1.3e", header='Date: '+str(dateobs)+' MJD: '+ str(mjd))

    m = np.nanmean(fl_cor)
    mr = np.nanmean(fl_raw)

    yrange = [ np.min(fl_cor[100:]/m) - 0.5,np.max(fl_cor[100:]/m) + 1 ]
    xrange = [ np.min(wl) - 50,np.max(wl) + 50 ]

    plt.plot(wl,fl_cor/m,color='k',linewidth=1,label=SN)

    ### Plot the galaxy lines
    if (plotgal == 'y') :
        l = [6562.819,4861.333,4340.471]
        ln = ['$H\\alpha$','$H\\beta$','$H\\gamma$']
        for line, line_name in zip(l,ln):            
            plt.axvline(x=line*(1.+z_gal),color='tab:red',linestyle='dashed',zorder=0,linewidth=0.7,label=line_name+' at '+str(z_gal))

        l = [6548.050,6583.460]
        plt.axvline(x=l[0]*(1.+z_gal),color='tab:blue',linestyle='dotted',zorder=0,linewidth=0.7,label='[NII] at '+str(z_gal))
        plt.axvline(x=l[1]*(1.+z_gal),color='tab:blue',linestyle='dotted',zorder=0,linewidth=0.7)

        l = [6716.440,6730.810]
        plt.axvline(x=l[0]*(1.+z_gal),color='tab:green',linestyle='dashed',zorder=0,linewidth=0.7,label='[SII] at '+str(z_gal))
        plt.axvline(x=l[1]*(1.+z_gal),color='tab:green',linestyle='dashed',zorder=0,linewidth=0.7)

        l = [5889.950,5895.924]
        plt.axvline(x=l[0]*(1.+z_gal),color='tab:orange',linestyle='dotted',zorder=0,linewidth=0.7,label='NaI at '+str(z_gal))
        plt.axvline(x=l[1]*(1.+z_gal),color='tab:orange',linestyle='dotted',zorder=0,linewidth=0.7)

        l = [5008.240,4960.295]
        plt.axvline(x=l[0]*(1.+z_gal),color='tab:cyan',linestyle='dotted',zorder=0,linewidth=0.7,label='[OIII] at '+str(z_gal))
        plt.axvline(x=l[1]*(1.+z_gal),color='tab:cyan',linestyle='dotted',zorder=0,linewidth=0.7)


    plt.legend(loc='upper right',fontsize=7)
    plt.ylim(yrange)
    plt.xlim(xrange)
    plt.xlabel('Obs. Wavelength [$\\AA$]')
    plt.ylabel('Scaled flux')
    plt.savefig('.'+plot_save_location+SN+'_'+dateobs+'.pdf',bbox_inches='tight')
    plt.close()

    ### Compare the two spectra
    plt.plot(wl,fl_raw/mr,color='grey',alpha=0.8,linewidth=1,label='Original')
    plt.plot(wl,fl_cor/m,color='tab:red',alpha=1,linewidth=1,label='Corrected')
    plt.legend()
    plt.xlabel('Obs. Wavelength [$\\AA$]')
    plt.ylabel('Scaled flux')
    plt.ylim(yrange)
    plt.xlim(xrange)
    plt.savefig('.'+plot_save_location+SN+'_'+dateobs+'_compare.pdf',bbox_inches='tight')
    plt.close()
