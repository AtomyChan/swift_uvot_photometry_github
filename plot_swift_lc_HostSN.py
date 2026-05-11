import matplotlib.pylab as plt
from astropy.table import Table,Column
import os
import numpy as np

import matplotlib.gridspec as gridspec
from SWIFT_UVOT_photometry_toolkit import swift_rate2mag

# the function for loading the lc data
def get_lcs_files_and_data(snname, lc_master_dir="/home/ping/photometry/SWIFT/", which_version='sum'):
    '''
    which_version: 'sum' or 'no_sum' or 'hybrid
    '''
    sndir = os.path.join(lc_master_dir, snname)
    swift_flts = ['vv', 'bb', 'uu', 'w1', 'm2', 'w2']
    lcfile_dict = {}
    lcdata_dict = {}
    for flt in swift_flts:
        if which_version == 'sum':
            lcfile = os.path.join(sndir, '%s_SWs-%s_uvotsource_ret.txt' % (snname, flt))
        elif which_version == 'no_sum':
            lcfile = os.path.join(sndir, '%s_SWs-%s_uvotsource_ret_no_imgsum.txt' % (snname, flt))
        elif which_version == 'hybrid':
            lcfile = os.path.join(sndir, '%s_SWs-%s_uvotsource_ret_hybrid.txt' % (snname, flt))
        else:
            raise ValueError('which_version should be sum or no_sum or hybrid') 
        
        if not os.path.exists(lcfile):
            lcfile_dict[flt] = None
            lcdata_dict[flt] = None
        else:
            lcfile_dict[flt] = lcfile
            lcdata_flt = Table.read(lcfile, format='ascii.fixed_width')
            if 'sss_factor' not in lcdata_flt.colnames:
                lcdata_flt.add_column(Column([1.0]*len(lcdata_flt), name='sss_factor'))
            lcdata_dict[flt] = lcdata_flt


    return lcfile_dict, lcdata_dict


# the function for displaying the light curve
def plot_lcs_and_prepare_lcfile(lcdata_dict, flts_plt=None, show_flux_density=False, abmag=True, outfile=None, show_on_screen=True, lc_plot_fname=None):

    abvega_offsets = {'vv': -0.01, 'bb': -0.13, 'uu': 1.02, 'w1': 1.51, 'm2': 1.69, 'w2': 1.73} #ab-vega

    fig = plt.figure(figsize=(15, 10))
    gs = gridspec.GridSpec(nrows=2, ncols=3, left=0.08, right=0.95, wspace=0.3, hspace=0.2)

    colors = {'vv': 'r', 'bb': 'y', 'uu': 'c', 'w1': 'g', 'm2': 'b', 'w2': 'k'}

    if flts_plt is None:
        flts_plt = list(lcdata_dict.keys())

    mjd_list = []
    mag_list = []
    magerr_list = []
    limiting_mag_list = []
    magsys_list = []
    filter_list = []
    metadata_list = []

    for i,flt in enumerate(flts_plt):
	
        if flt in ['m2','w2','w1']:
            fltfritz = 'uv'+flt
        else:
            fltfritz=flt[0]
        fltfritz =  "uvot::"+fltfritz
	
        ax = fig.add_subplot(gs[np.divmod(i,3)[0], np.divmod(i,3)[1]])

        color = colors[flt]
        lcdata = lcdata_dict[flt]

        if lcdata is None:
            continue

        if len(lcdata) == 0:
            continue
       
        if abmag or show_flux_density:
            offset = abvega_offsets[flt]
            magsys = 'ab'
        else:
            offset = 0
            magsys = 'vega'

        jds = lcdata['JD']
        mags = lcdata['mag']+offset
        magerrs = np.sqrt(lcdata['mag_err_stat']**2+lcdata['mag_err_sys']**2)
        snrate_err = lcdata['rate_err_corr'].data
        sss_factor = lcdata['sss_factor'].data
        limmags, temp = swift_rate2mag(snrate_err*3, flt)        
        limmags = limmags+offset
        detections = mags<limmags
        nondetections = mags>=limmags
        sss_good = sss_factor==1.0
        sss_bad = sss_factor!=1.0

        if np.sum(nondetections) > 0:
            has_nondetection = 1
        else:
            has_nondetection = 0
       

        for jd,mag,magerr,limmag,detection in zip(jds, mags, magerrs, limmags, detections):
            #print(jd-2400000.5, np.round(mag,3), np.round(magerr,3), np.round(limmag,3), magsys, flt)
            mjd_list.append(np.round(jd-2400000.5,5))
            if detection:
                mag_list.append(np.round(mag,3))
                magerr_list.append(np.round(magerr,3))
            else:
                mag_list.append(None)
                magerr_list.append(None)
            limiting_mag_list.append(np.round(limmag,3))
            magsys_list.append(magsys)
            filter_list.append(fltfritz)
            metadata_list.append(None)

        if show_flux_density:
            ax.errorbar(jds[detections], 10**(-mags[detections]/2.5)*3.631e9,  yerr=magerrs[detections]/1.086*10**(-mags[detections]/2.5)*3.631e9, fmt=color+'o', label=flt) #assuming small uncertainties
            if has_nondetection:
                ax.plot(jds[nondetections], 10**(-limmags[nondetections]/2.5)*3.631e9, 'v', color=color, alpha=0.5)


            ax.set_ylabel('Flux density ($\mu$Jy)', fontsize=15)
        else:
            ax.errorbar(jds[detections*sss_good], mags[detections*sss_good],  yerr=magerrs[detections*sss_good], fmt=color+'o', label=flt + ' (sss good)')
            if has_nondetection:
                ax.plot(jds[nondetections*sss_good], limmags[nondetections*sss_good], 'v', color=color, alpha=0.5)

            ax.errorbar(jds[detections*sss_bad], mags[detections*sss_bad],  yerr=magerrs[detections*sss_bad], fmt=color+'o', mfc='none', mec=color, label=flt + ' (sss bad)')
            if has_nondetection:
                ax.plot(jds[nondetections*sss_bad], limmags[nondetections*sss_bad], 'v', color=color, mfc='none', mec=color, alpha=0.5)

            ax.set_ylabel('Magnitude', fontsize=15)
            ax.invert_yaxis()

        ax.tick_params(which='both', labelsize=15, direction='in')
        ax.set_xlabel('JD', fontsize=15)
        ax.legend(loc=0)

    if lc_plot_fname is not None:
        plt.savefig(lc_plot_fname)

    if show_on_screen:
        plt.show()


    lctable = Table([mjd_list, mag_list, magerr_list, limiting_mag_list, magsys_list, filter_list, metadata_list], names=['mjd','mag','magerr','limiting_mag','magsys','filter','altdata.meta1'])
    print(lctable)
    lctable.sort('mjd')
    if outfile is not None:
        lctable.write(outfile, format='ascii.csv', overwrite=1)

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser()    
    
    def_target_name = ''
    parser.add_option('-t', '--target_name', dest="target_name", type="string", default=def_target_name,
                      help="target name if this pipeline is run for individual source [%s]" % def_target_name)

    
    def_swift_flts = ''
    parser.add_option('-f', '--flts', dest='swift_flts', type="string", default=def_swift_flts,
                      help="filter names for which you want the photometry, one or multiple from vv,bb,uu,w1,m2,w2. Separated by comma if multiple. The defaut is all above filers if not specified")

    def_version = 'sum' # 'sum' or 'no_sum' or 'hybrid'
    parser.add_option('-v', '--version', dest="version", type="string", default=def_version,
                      help="which version of the photometry reuslt to plot [default: %s]. Valid options: sum, no_sum, hybrid" % def_version)
    
    options, remainder = parser.parse_args()

    target_name = options.target_name
    flts = options.swift_flts
    version = options.version 
    
    if flts == '':
        swift_uvot_flts = ['w2', 'm2', 'w1', 'uu', 'bb', 'vv']
    else:
        swift_uvot_flts = flts.split(',')

        
    lc_master_dir="/home/ping/photometry/SWIFT/"
    lcfiles, lcdata = get_lcs_files_and_data(target_name, which_version=version)
    print(lcfiles)

    savefile = os.path.join(lc_master_dir, target_name+'_toFritz.csv')

    plot_lcs_and_prepare_lcfile(lcdata, flts_plt=swift_uvot_flts, outfile=savefile)
#    plot_lcs_and_prepare_lcfile(lcdata, flts_plt=swift_uvot_flts, outfile=None, show_flux_density=True)



