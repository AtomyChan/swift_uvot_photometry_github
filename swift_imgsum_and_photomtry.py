#! /home/ping/anaconda3/bin/python

# swift image sum and uvot photometry

import os
import sys
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table, vstack

from SWIFT_UVOT_photometry_toolkit import prepare_filepaths, prepare_source_reg_func, prepare_bkg_reg_func, sum_swift_images


def aspect_correction_swift_images():
    '''
    align swift images 
    refer to the following for details:
    http://www.swift.ac.uk/analysis/uvot/image.php
    '''

    print("under construction...")


def prepare_source_bkg_region_file(ra, dec, flt, app, sky_in, sky_out, source_reg_filename, bkg_reg_filename, verbose=1, renew_region_file=0):
    '''
    if source region file exists in photometry workplace directory for given target then just copy the target-specific region template as source region files otherwise, create source region from the generic template
    INPUTS:
       ra:
       dec:
       app:
       sky_in:
       sky_out:
       
    '''
    if os.path.exists(source_reg_filename) and (not renew_region_file):
        update_or_create_source_reg = False
    else:
        update_or_create_source_reg = True
        sn_source_reg_template = os.path.join(work_dir, 'source_%sarc.reg'%(str(int(app))))
        if os.path.exists(sn_source_reg_template):
            source_reg_template = sn_source_reg_template
            copy_source_reg_template = True
        else:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            source_reg_template = os.path.join(script_dir, 'source_template.reg')
            copy_source_reg_template = False
            if not os.path.exists(source_reg_template):
                raise IOError("the source region file template does not exist")

            temp = prepare_source_reg_func(ra, dec, app, source_reg_template, source_reg_filename, copy_template_source=copy_source_reg_template, verbose=verbose)

    if os.path.exists(bkg_reg_filename) and (not renew_region_file):
        update_or_create_bkg_reg = False
    else:
        update_or_create_bkg_reg = True
        bkg_reg_flt_specific  = os.path.join( os.path.join(work_dir, 'bkg_region_manual'), 'bkg_%s.reg'%flt) 
        if os.path.exists(bkg_reg_flt_specific):
            sn_bkg_reg_template = bkg_reg_flt_specific
        else:
            sn_bkg_reg_template = os.path.join(work_dir, 'bkg.reg')


        if os.path.exists(sn_bkg_reg_template):
            bkg_reg_template = sn_bkg_reg_template
            copy_bkg_reg_template = True
        else:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            bkg_reg_template = os.path.join(script_dir, 'bkg_template.reg')
            copy_bkg_reg_template = False
            if not os.path.exists(bkg_reg_template):
                raise IOError("the background region file template does not exist")
        temp = prepare_bkg_reg_func(ra, dec, sky_in, sky_out, bkg_reg_template, bkg_reg_filename, copy_template_bkg=copy_bkg_reg_template, verbose=verbose)

    return update_or_create_source_reg, update_or_create_bkg_reg



def swift_uvot_photometry_tpls(target_name, repodir_sn, work_dir, sn_ra, sn_dec, app=5, sky_in=15, sky_out=30, start_over=True, renew_region_file=False, apercorr='CURVEOFGROWTH', verbose=1):
    '''
    work on photometry of template images

    INPUTS:
    
    '''

    tplimg_dir = os.path.join(repodir_sn, 'template')


    # prepare result storage
    DATE_OBS = np.array([])  # start time of the observation
    DATE_END = np.array([])
    JD = np.array([])
    flts = np.array([])
    exptimes = np.array([])
    rates_corr = np.array([])
    rates_err_corr = np.array([])
    fluxs_aa = np.array([])
    fluxs_aa_err_stat = np.array([])
    fluxs_aa_err_sys = np.array([])

    MAG = np.array([])
    MAG_ERR_STAT = np.array([])
    MAG_ERR_SYS = np.array([])
    MAG_ERR = np.array([])

    flts_valid = ['w2', 'm2', 'w1', 'uu', 'bb', 'vv']
    for flt in flts_valid:

        img_name = "tplimg_%s.fits" % flt
        img_abs_name = os.path.join(tplimg_dir, img_name)
        print(img_abs_name)
        if not os.path.isfile(img_abs_name):
            print("template image for filter %s not available" % flt)
            continue

        hdu = fits.open(img_abs_name)
        hdr = hdu[0].header
        date_obs = hdr['DATE-OBS']
        date_end = hdr['DATE-END']
        hdu.close()
        DATE_OBS = np.append(DATE_OBS, date_obs.split('T')[0])
        tstart = Time(date_obs, scale='utc', format='isot').jd
        tend = Time(date_end, scale='utc', format='isot').jd
        t = (tstart+tend)/2  # this is simplied

        # ----------------- prepare source region file and background region file --------------------
        source_reg_filename = os.path.join(work_dir,"source_tpl_%s_sum.reg"%flt)
        bkg_reg_filename = os.path.join(work_dir,"bkg_tpl_%s_sum.reg"%flt)
        update_or_create_source_reg, update_or_create_bkg_reg = prepare_source_bkg_region_file(sn_ra, sn_dec, flt, app, sky_in, sky_out, source_reg_filename, bkg_reg_filename, verbose=verbose, renew_region_file=renew_region_file)
        # ----------------- end of creating source and bkg region ---------------------------

        aperture_size = "_%sarc" %str(int(app))
        log_file = os.path.join(work_dir, img_name.split('.')[0] + aperture_size + '_sum.log')
        if os.path.exists(log_file):
            os.remove(log_file)

        res_file = os.path.join(work_dir, img_name.split('.')[0]+'_res' + aperture_size + '_sum.fits')
        if (not os.path.exists(res_file)) or start_over:
            if os.path.exists(res_file):
                os.remove(res_file)
            phot_command = 'uvotsource image=%s+1 srcreg=%s bkgreg=%s sigma=3 outfile=%s apercorr=%s >%s' % (img_abs_name, source_reg_filename, bkg_reg_filename, res_file, apercorr, log_file)

            print(phot_command)
            os.system(phot_command)
        else:
            print("old photometry result file %s exists and will be used"%res_file)

        JD = np.append(JD, t)
        hdu_res = fits.open(res_file)
        data = hdu_res[1].data

        flts = np.append(flts, flt)

        mag = data['MAG'][0]
        MAG = np.append(MAG, mag)
        mag_err_stat = data['MAG_ERR_STAT'][0]
        MAG_ERR_STAT = np.append(MAG_ERR_STAT, mag_err_stat)
        mag_err_sys = data['MAG_ERR_SYS'][0]
        MAG_ERR_SYS = np.append(MAG_ERR_SYS, mag_err_sys)
        mag_err = np.sqrt(mag_err_stat**2 + mag_err_sys**2)
        MAG_ERR = np.append(MAG_ERR, mag_err)

        exptime = data['EXPOSURE'][0]
        exptimes = np.append(exptimes, exptime)

        rate_corr = data['CORR_RATE'][0]
        rates_corr = np.append(rates_corr, rate_corr)

        rate_err_corr = data['CORR_RATE_ERR'][0]
        rates_err_corr = np.append(rates_err_corr, rate_err_corr)

        flux_aa = data['FLUX_AA'][0]
        fluxs_aa = np.append(fluxs_aa, flux_aa)

        flux_aa_err_stat = data['FLUX_AA_ERR_STAT'][0]
        fluxs_aa_err_stat = np.append(fluxs_aa_err_stat, flux_aa_err_stat)

        flux_aa_err_sys = data['FLUX_AA_ERR_SYS'][0]
        fluxs_aa_err_sys = np.append(fluxs_aa_err_sys, flux_aa_err_sys)

    photret_table = os.path.join(work_dir, target_name+'_tpl_uvotsource_ret.txt')

    photdata = [JD, flts, exptimes, rates_corr, rates_err_corr, fluxs_aa, fluxs_aa_err_stat, fluxs_aa_err_sys, MAG, MAG_ERR_STAT, MAG_ERR_SYS]
    photrettable = Table(photdata, names=['JD', 'flt', 'exptime', 'rate_corr', 'rate_err_corr', 'flux_aa', 'flux_aa_err_stat', 'flux_aa_err_sys', 'mag', 'mag_err_stat', 'mag_err_sys'])
    photrettable.write(photret_table, format='ascii.fixed_width', overwrite=1)


def swift_imgsum_and_uvot_photometry_flt(target_name, data_dir, work_dir,  obs_sequences, flt, sn_ra, sn_dec, sniper_dir, app=5, sky_in=15, sky_out=30, remove_exposure_with_sss_issue=True,  exclude='DEFAULT', photret_table_without_imgsum=None, \
    start_over=True, renew_summed_image=False, renew_region_file=False, check_aspcorr_prior_summing=True, apercorr='CURVEOFGROWTH', verbose=1):
    '''
    summing UVTO images and then perform photometry

    INPUTS:
        target_name:
        data_dir: 
        work_dir:
        obs_sequences: iterable object with items providing the string with 11 digits such as 00016848001 which gives the obsid and the seg number 
        
    '''

    # get all images sequence directories and workplace direcotory for given target
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)

    if remove_exposure_with_sss_issue:
        if (photret_table_without_imgsum is None):
            raise ValueError("the photometry result table without image summing is required to remove exposures with sss issue")
        print(photret_table_without_imgsum)
        print(photret_table_without_imgsum['image']) 
        imgs_all = np.array([img[:-3] for img in photret_table_without_imgsum['image']])  

    # prepare result storage
    obssequences = np.array([])
    DATE_OBS = np.array([])
    JD = np.array([])
    exptimes = np.array([])
    rates_corr = np.array([])
    rates_err_corr = np.array([])
    fluxs_aa = np.array([])
    fluxs_aa_err_stat = np.array([])
    fluxs_aa_err_sys = np.array([])

    MAG = np.array([])
    MAG_ERR_STAT = np.array([])
    MAG_ERR_SYS = np.array([])
    MAG_ERR = np.array([])

    log_file = os.path.join(work_dir, 'photometry_summed_images.log')
    fid = open(log_file, 'a')

    nowtime = Time.now().isot
    fid.write(nowtime+'\n')
    fid.write('%s\n' % target_name)
    fid.write("(%s,%s)\n" % (sn_ra, sn_dec))
    fid.write("aperture size: %s arcs, bkg annulus in: %s arcs, bkg annulus out: %s arcs \n" % ( app, sky_in, sky_out))

    for obs_sequence in obs_sequences:
        obs_sequence_dir = os.path.join(data_dir, obs_sequence)
        # prepare absolute image
        img_dir = os.path.join(obs_sequence_dir, 'uvot/image')
        img_prefix = 'sw'
        img_name = img_prefix + obs_sequence + 'u' + flt + '_sk.img.gz'
        img_abs_name = os.path.join(img_dir, img_name)
        if not os.path.exists(img_abs_name):
            img_abs_name = img_abs_name[:-3] # remove the .gz in the case the compressed version not exist then try to use the uncompressed file exists 

        img_out = flt + '_sum.fits'
        img_out_abs = os.path.join(img_dir, img_out)

        if (not os.path.isfile(img_abs_name)) and (not os.path.exists(img_out_abs)): #something wrong here? to check 
            continue

        if remove_exposure_with_sss_issue:
            mod_img_filename = img_abs_name.replace('img', 'img.mod')
            print(mod_img_filename) 
            
            hduindex_to_remove =[]
            temp = photret_table_without_imgsum[imgs_all == img_abs_name]
            Nhduindex = len(temp) # number of exposures in the image
            if Nhduindex == 0:
                raise ValueError("something wrong in matching image name %s"%img_abs_name)
            for img, sss_factor in temp[['image', 'sss_factor']]:
                index = int(img[-2])
                if sss_factor == -99.9:
                    hduindex_to_remove.append(index)
            if len(hduindex_to_remove) > 0:
                if len(hduindex_to_remove) < Nhduindex:
                    fid.write("obs_sequence %s: exposure(s) %s with sss issue will be removed before summing images \n" % (obs_sequence, ','.join([str(ind) for ind in hduindex_to_remove])))

                    with fits.open(img_abs_name, mode='readonly') as hdulist:
                    # Create a new HDUList excluding the specified extension
                        new_hdulist = fits.HDUList([hdu for i, hdu in enumerate(hdulist) if i not in hduindex_to_remove])
                        new_hdulist.writeto(mod_img_filename, overwrite=True)
                        img_abs_name = mod_img_filename # use the modified image for summing 
                else:
                    fid.write("obs_sequence %s: all exposures have sss issue, no image summing and photometry will be performed \n" % obs_sequence)
                    continue            
            
            
        # summing swift UVOT images
        if (not os.path.exists(img_out_abs)) or renew_summed_image:
            hdu = fits.open(img_abs_name)
            aspcorrs = []
            for i in range(len(hdu)):
                if i == 0:
                    continue
                hdr = hdu[i].header
                aspcorrs.append(hdr['ASPCORR'])
            if 'DIRECT' not in aspcorrs:
                logmessage = ' '.join(aspcorrs)
                fid.write("!!! none of images in %s is aligned to WCS: %s \n" % (img_abs_name, logmessage))

            if 'NONE' in aspcorrs:
                logmessage = ' '.join(aspcorrs)
                fid.write("!!! images without WCS solution exist in %s: %s \n" % (img_abs_name, logmessage))

            sum_swift_images(img_abs_name, img_out_abs, exclude=exclude, verbose=0, ignoreframetime=False, trybestsum=True, check_aspcorr=check_aspcorr_prior_summing, chatter=0)

            hdu.close()

        else:
            fid.write("%s already exists \n" % img_out_abs)

        if not os.path.exists(img_out_abs):
            continue

        hdu = fits.open(img_out_abs)
        hdr = hdu[0].header
        date_obs = hdr['DATE-OBS']
        date_end = hdr['DATE-END']
        hdu.close()

        DATE_OBS = np.append(DATE_OBS, date_obs.split('T')[0])
        tstart = Time(date_obs, scale='utc', format='isot').jd
        tend = Time(date_end, scale='utc', format='isot').jd
        t = (tstart+tend)/2
        # ----------------- prepare source region file and background region file --------------------
        segnum = obs_sequence[-3:]
        source_reg_filename = os.path.join(work_dir, "source_%s%s_sum.reg"%(segnum, flt))
        bkg_reg_filename = os.path.join(work_dir, "bkg_%s%s_sum.reg"%(segnum, flt))
        update_or_create_source_reg, update_or_create_bkg_reg = prepare_source_bkg_region_file(sn_ra, sn_dec, flt, app, sky_in, sky_out, source_reg_filename, bkg_reg_filename, renew_region_file=renew_region_file, verbose=verbose)
        # ----------------- end of creating source and bkg region ---------------------------

        aperture_size = "_%sarc" %str(int(app))
        log_file = os.path.join(work_dir, img_name.split('.')[0] + aperture_size + '_sum.log')
        if os.path.exists(log_file):
            os.remove(log_file)

        res_file = os.path.join(work_dir, img_name.split('.')[0]+'_res' + aperture_size + '_sum.fits')
        if (not os.path.exists(res_file)) or start_over:
            if os.path.exists(res_file):
                os.remove(res_file)
            phot_command = 'uvotsource image=%s+1 srcreg=%s bkgreg=%s sigma=3 outfile=%s apercorr=%s >%s' % (img_out_abs, source_reg_filename, bkg_reg_filename, res_file, apercorr, log_file)

            print(phot_command)
            os.system(phot_command)
        else:
            print("old photometry result file %s exists and will be used"%res_file)

        obssequences = np.append(obssequences, obs_sequence)
        JD = np.append(JD, t)
        hdu_res = fits.open(res_file)
        data = hdu_res[1].data

        mag = data['MAG'][0]
        MAG = np.append(MAG, mag)
        mag_err_stat = data['MAG_ERR_STAT'][0]
        MAG_ERR_STAT = np.append(MAG_ERR_STAT, mag_err_stat)
        mag_err_sys = data['MAG_ERR_SYS'][0]
        MAG_ERR_SYS = np.append(MAG_ERR_SYS, mag_err_sys)
        mag_err = np.sqrt(mag_err_stat**2 + mag_err_sys**2)
        MAG_ERR = np.append(MAG_ERR, mag_err)

        exptime = data['EXPOSURE'][0]
        exptimes = np.append(exptimes, exptime)

        rate_corr = data['CORR_RATE'][0]
        rates_corr = np.append(rates_corr, rate_corr)

        rate_err_corr = data['CORR_RATE_ERR'][0]
        rates_err_corr = np.append(rates_err_corr, rate_err_corr)

        flux_aa = data['FLUX_AA'][0]
        fluxs_aa = np.append(fluxs_aa, flux_aa)

        flux_aa_err_stat = data['FLUX_AA_ERR_STAT'][0]
        fluxs_aa_err_stat = np.append(fluxs_aa_err_stat, flux_aa_err_stat)

        flux_aa_err_sys = data['FLUX_AA_ERR_SYS'][0]
        fluxs_aa_err_sys = np.append(fluxs_aa_err_sys, flux_aa_err_sys)

    photdata = [obssequences, JD, exptimes, rates_corr, rates_err_corr, fluxs_aa, fluxs_aa_err_stat, fluxs_aa_err_sys, MAG, MAG_ERR_STAT, MAG_ERR_SYS]
    photrettable = Table(photdata, names=['obsidseg', 'JD', 'exptime', 'rate_corr', 'rate_err_corr', 'flux_aa', 'flux_aa_err_stat', 'flux_aa_err_sys', 'mag', 'mag_err_stat', 'mag_err_sys'])
    photrettable.sort('JD')
    photrettable['exptime'].info.format = '.3f'
    photrettable['rate_corr'].info.format = '.3f'
    photrettable['rate_err_corr'].info.format = '.3f'
    photrettable['mag'].info.format = '.3f'
    photrettable['mag_err_stat'].info.format = '.3f'
    photrettable['mag_err_sys'].info.format = '.3f'
    photrettable['flux_aa'].info.format = '.3e'
    photrettable['flux_aa_err_stat'].info.format = '.3e'
    photrettable['flux_aa_err_sys'].info.format = '.3e'

    obs_num_real = len(JD)
    lc_data = np.concatenate((JD.reshape((obs_num_real, 1)), MAG.reshape( (obs_num_real, 1)), MAG_ERR.reshape((obs_num_real, 1))), axis=1)
    lc_data = lc_data[np.argsort(lc_data[:,0])]

    fid.close()

    return photrettable, lc_data


if __name__ == '__main__':

    import optparse
    parser = optparse.OptionParser()

    def_target_name = ''
    parser.add_option('-t', '--target_name', dest="target_name", type="string", default=def_target_name,
                      help="target name if this pipeline is run for individual source [%s]" % def_target_name)

    def_sn_work_name = ''
    parser.add_option('-w', '--work_name', dest="sn_work_name", type="string", default=def_sn_work_name,
                      help="workplace and data directory name for the target if this pipeline is run for one given id [%s]" % def_sn_work_name)

    def_swift_id = ''
    parser.add_option('--obsid', '--swift_id', dest='swift_id', type="string", default=def_swift_id,
                      help="observation id assigned to sepcific observation which can be found here www.swift.psu.edu/secure/toop/summary.php")

    def_swift_flts = ''
    parser.add_option('-f', '--flts', dest='swift_flts', type="string", default=def_swift_flts,
                      help="filter names for which you want the photometry, one or multiple from vv,bb,uu,w1,m2,w2. Separated by comma if multiple. The defaut is all above filers if not specified")

    def_RA = '00:00:00'
    parser.add_option('-r', '--ra', '--RA', dest="RA", type="string", default=def_RA,
                      help="Target Right Ascension in hours, sexidecimal [%s]" % def_RA)

    def_DEC = '00:00:00'
    parser.add_option('-d', '--dec', '--DEC', dest="DEC", type="string", default=def_DEC,
                      help="Target Declination in degrees, sexidecimal [%s]" % def_DEC)

    def_app = 5
    parser.add_option('-a', '--app', dest="source_aperture_size", type="float",
                      default=def_app,   help="source aperture size in arcsecond [%s]" % def_app)

    def_sky_in = 100
    parser.add_option('--sky_in', dest="SKY_IN", type="float", default=def_sky_in,
                      help="Inner radius size of background annulus in arcsecond [%s]" % def_sky_in)

    def_sky_out = 150
    parser.add_option('--sky_out', dest="SKY_OUT", type="float", default=def_sky_out,
                      help="Out radius size of background annulus in arcsecond [%s]" % def_sky_out)

    def_startover = False
    parser.add_option('--start_over', dest="start_over", action="store_true", default=def_startover,
                      help="whether or not start over(delete all existed photometry result file)[%s]" % def_startover)

    def_renewsum = False
    parser.add_option('--renew_sum', dest="renew_sum", action="store_true",
                      default=def_renewsum, help="whether or not renew the summed images[%s]" % def_renewsum)

    def_renewreg = False
    parser.add_option('--renew_reg', dest="renew_reg", action="store_true", default=def_renewreg,
                      help="whether or not renew the source/bkg region files[%s]" % def_renewreg)

    def_tplphotonly = False
    parser.add_option('--tplphot_only', dest="tplphot_only", action="store_true",
                      default=def_tplphotonly, help="only do photometry for template images[%s]" % def_tplphotonly)

    def_skip_check_aspcorr = False
    parser.add_option('--skipcheck', dest='skip_check_aspcorr', action="store_true",
                      default=def_skip_check_aspcorr, help="skip checking ASPCORR before summing the images[%s]" % def_skip_check_aspcorr)

    def_remove_sss_issue = False
    parser.add_option('--remove_sss', dest='remove_sss_issue', action="store_true",
                      default=def_remove_sss_issue, help="remove those exposures with Small Scale Sensitivity issue before summing the images[%s]" % def_remove_sss_issue)


    def_exclude = 'DEFAULT' # NONE
    parser.add_option('-e', '--exclude', dest="exclude", type="string", default=def_exclude,
                      help="The option for exclude when performing imgsum [%s]" % def_exclude)

    def_perform_apercorr = True
    parser.add_option('--apercorr', dest="perform_apercorr", action="store_false", default=def_perform_apercorr,
                      help="whether or not perform aperture correction[%s]" % def_perform_apercorr)

    def_verbose = False
    parser.add_option('--verbose', dest="verbose", action="store_true", default=def_verbose,
                      help="whether or not have verbose message[%s]" % def_verbose)

    

    options, remainder = parser.parse_args()

    target_name = options.target_name
    sn_work_name = options.sn_work_name
    obs_id_list = options.swift_id.split(',')
    RA = options.RA
    DEC = options.DEC
    app = options.source_aperture_size
    sky_in = options.SKY_IN
    sky_out = options.SKY_OUT

    flts = options.swift_flts

    startover = options.start_over
    renewsum = options.renew_sum
    renewreg = options.renew_reg
    perform_apercorr=options.perform_apercorr
    if perform_apercorr:
        apercorr = 'CURVEOFGROWTH'
    else:
        print('alert: no aperture correction will be performed...')
        apercorr= 'NONE'

    tplonly = options.tplphot_only
    skipcheck = options.skip_check_aspcorr

    verbose= options.verbose

    if skipcheck:
        check_aspcorr_prior_summing = False
    else:
        check_aspcorr_prior_summing = True

    exclude=options.exclude
    remove_exposure_with_sss_issue = options.remove_sss_issue

    sne_master_table = Table.read('/home/ping/projects/Scripts/swift_uvot_photometry/source_info_SNRADec.csv', format='ascii.csv')

    if RA == '00:00:00':
        if target_name not in sne_master_table['Name']:
            raise IOError("target RA required...")
        else:
            RA = sne_master_table['RA'][sne_master_table['Name'] == target_name].data[0]

    if DEC == '00:00:00':
        if target_name not in sne_master_table['Name']:
            raise IOError("target DEC required...")
        else:
            DEC = sne_master_table['Dec'][sne_master_table['Name'] == target_name].data[0]

    if sn_work_name == '':
        sn_work_name = target_name

    sniper_dir = None  # save the final lc data to a place

    if flts == '':
        swift_uvot_flts = ['w2', 'm2', 'w1', 'uu', 'bb', 'vv']
    else:
        swift_uvot_flts = flts.split(',')

    master_dir = '/home/ping/'  # /mnt/md0/
    data_dir = os.path.join(master_dir, 'image/SWIFT_rawdownload/' + sn_work_name)
    work_dir = os.path.join(master_dir, 'photometry/SWIFT/' + sn_work_name)


    #repodir_swift = os.path.join(master_dir, 'image/SWIFT_rawdownload/')
    #repodir_sn = os.path.join(repodir_swift, target_name)
    #work_dir = os.path.join(master_dir, 'photometry/SWIFT/%s/' % target_name)

    obssequences = np.array([fname for fname in os.listdir(data_dir) if len(fname)==11 and (fname.startswith('000'), fname.startswith('030'))])
    allobsids = np.array([fname[3:8] for fname in obssequences])
    obsids_unique = np.unique(allobsids)
    obsid_seg_dict = {}
    for obsid in obsids_unique:
        obsid_seg_dict[obsid] = obssequences[allobsids==obsid]
    
    if obs_id_list[0]=='':
        print('to work on all available obs ids: ', obsids_unique)
        obs_id_list = obsids_unique



        


    if not tplonly:
        for flt in swift_uvot_flts:
            
            photret_table_file_without_imgsum = os.path.join( work_dir, target_name+'_SWs-' + flt + '_uvotsource_ret_no_imgsum.txt')
            if remove_exposure_with_sss_issue: 
                if not os.path.exists(photret_table_file_without_imgsum):
                    raise IOError("the photometry result file without image summing %s does not exist, please run the swift_photometry_without_imgsum.py to generate this file"%photret_table_file_without_imgsum) 
                photret_table_without_imgsum = Table.read(photret_table_file_without_imgsum, format='ascii.fixed_width') 
                print("read in the photometry result file without image summing %s to remove exposures with sss issue"%photret_table_file_without_imgsum)
                print(photret_table_without_imgsum) 
            else:
                photret_table_without_imgsum = None 
                
            photrettable = None
            lc_data = None
            for obs_id in obs_id_list:
                photrettable_temp, lc_data_temp = swift_imgsum_and_uvot_photometry_flt(target_name, data_dir, work_dir, obsid_seg_dict[obs_id], flt, RA, DEC, \
                    sniper_dir, app=app, sky_in=sky_in, sky_out=sky_out, remove_exposure_with_sss_issue=remove_exposure_with_sss_issue, exclude=exclude, \
                        photret_table_without_imgsum=photret_table_without_imgsum, start_over=startover, renew_summed_image=renewsum, renew_region_file=renewreg, \
                        check_aspcorr_prior_summing=check_aspcorr_prior_summing, apercorr=apercorr, verbose=verbose)

                if photrettable is None:
                    photrettable = photrettable_temp
                else:
                    if len(photrettable_temp)>0:
                        photrettable = vstack([photrettable, photrettable_temp])
                if lc_data is None:
                    lc_data = lc_data_temp
                else:
                    if len(lc_data_temp)>0:
                        lc_data = np.vstack((lc_data, lc_data_temp))

            photret_table_file = os.path.join( work_dir, target_name+'_SWs-' + flt + '_uvotsource_ret.txt')

            if verbose:
                print(photrettable)
            photrettable.write(photret_table_file, format='ascii.fixed_width', overwrite=1)

            lc_file_sndir = os.path.join(work_dir, target_name + '_SWs-' + flt + '.txt')
            fid_sndir = open(lc_file_sndir, 'wt')
            fid_sndir.write('#jd mag(vega) magerr\n') 
            np.savetxt(fid_sndir, lc_data, fmt="%8.4f %6.3f %6.3f")
            fid_sndir.close()

            # don't create empty light curve for SNIPER, which will stop displaying ascii files on SNIPER
            # also remove magnitudes not realistic
            if len(lc_data) > 0 and (sniper_dir is not None):
                lc_file_sniper = os.path.join( sniper_dir, target_name + '_SWs-' + flt + '.txt')
                lc_data_sniper = lc_data[lc_data[:, 1] < 30, :]
                fid_sniper = open(lc_file_sniper, 'wt')
                np.savetxt(fid_sniper, lc_data_sniper, fmt="%8.4f %6.3f %6.3f")

    if tplonly:
        swift_uvot_photometry_tpls(target_name, data_dir, work_dir, RA, DEC, app=app, sky_in=sky_in, sky_out=sky_out, start_over=True, renew_region_file=renewreg, apercorr=apercorr)
