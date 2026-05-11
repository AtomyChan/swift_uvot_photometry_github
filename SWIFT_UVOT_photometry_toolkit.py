#! /usr/bin/python

# first created by Ping Chen on 2016-09-27
# functions for swift photometry
# modified by Ping on 2016-10-23: add start_over argument to swift_photometry_uvot_flt


import os
import sys
import shutil
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table, Row, Column, vstack


def swift_rate2mag(rate, flt, rate_err=None, verbose=0):
    '''
    convert UVOT count rate to magnitudes in UVOT system

    The updated zeropoints presented in https://doi.org/10.1063/1.3621807 will be used	
    '''

    uvot_zpts = {'w2': 17.38, 'm2': 16.85, 'w1': 17.44,
                 'uu': 18.34, 'bb': 19.11, 'vv': 17.89}
    uvot_zpterrs = {'w2': 0.03, 'm2': 0.03,
                    'w1': 0.03, 'uu': 0.02, 'bb': 0.02, 'vv': 0.01}

    mag = -2.5*np.log10(rate) + uvot_zpts[flt]

    if rate_err is not None:
        magerr_stat = 2.5*np.log10(1+rate_err/rate)
    else:
        magerr_stat = 0
        if verbose:
            print("statistical uncertainty on count rate not available; only systematic uncertainty on the output magnitude will be provided")

    magerr = np.sqrt(magerr_stat**2 + uvot_zpterrs[flt]**2)

    return mag, magerr


def swift_photometry_uvot_flt(target_name, obs_id, flt, sn_ra, sn_dec, master_dir='/home/asassn/ASASSN/',  app=5, sky_in=15, sky_out=30, apercorr='CURVEOFGROWTH', start_over=True):
    '''
    perform photometry on swift images with 'uvotsource', see http://www.swift.ac.uk/analysis/uvot/mag.php for how to use 'uvotsource'

    INPUTs:
            target_name: also the directory name storing the images 
            obs_id: SWIFT ID 
            flt: uvot filter for which to extract magnitudes, valid options are 'vv','bb','uu','w1','m2','w2' 
            sn_ra:
            sn_dec:	
            sniper_dir: the directory to store light curves which feed the SNIPER 
            app: source aperture size 
            sky_in: inner annulus radius of background
            sky_out: outer annulus radius of background
            start_over: if True, reproduce the light curve from all images; if False, only perform photometry images from which have not extracted magnitudes 
                        ATTENTIONS: even start_over is true, old source region files will be used
    '''
    if not os.path.exists(master_dir):
        raise IOError("master directory %s not exists..." % master_dir) 
    
    # get all images sequence directories and workplace direcotory for given target
    obs_sequence_dirs, work_dir = prepare_filepaths(target_name, obs_id, master_dir=master_dir)
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)

    # prepare result storage
    DATE_OBS = np.array([])
    EXPOSURE = np.array([])
    JD = np.array([])
    MAG = np.array([])
    MAG_ERR_STAT = np.array([])
    MAG_ERR_SYS = np.array([])
    MAG_ERR = np.array([])
    rates_corr = np.array([])
    rates_err_corr = np.array([])
    fluxs_aa = np.array([])
    fluxs_aa_err_stat = np.array([])
    fluxs_aa_err_sys = np.array([])

    images = np.array([])
    sequences_done = np.array([])
    qualflags = np.array([])
    photflags = np.array([])
    sss_factors = np.array([])

    # check for droped images: this is only relevant if you registered bad images after visual inspections
    drop_flag_file = os.path.join(
        work_dir, target_name + '_SW-' + flt + '_drop.txt')
    if os.path.exists(drop_flag_file):
        droped_imgs = [drop_img.strip()
                       for drop_img in open(drop_flag_file).readlines()]
        drop_specific_obs = True
    else:
        drop_specific_obs = False

    # load old results if exists
    old_lc_file = os.path.join(
        work_dir, target_name + '_SW-' + flt + '_wholeinfo.txt')
    if os.path.exists(old_lc_file):
        print(old_lc_file)
        dtype = {'names': ('JD', 'mag', 'mag_err', 'image', 'sequence'), 'formats': ('f8', 'f8', 'f8', 'S120', 'S20')}
        old_lc = np.loadtxt(old_lc_file, dtype=dtype)
        print(old_lc)
        pass_sequences = np.unique(old_lc['sequence'])
        print(("start over photometey for %s: [%s]" % (target_name, start_over)))
        print(pass_sequences)

    else:
        print(("No old photometry results found for target %s" % target_name))
        if not start_over:
            print(("photometry first try for target %s" % target_name))
            start_over = True

    all_sequences = [seqdir.split('/')[-1] for seqdir in obs_sequence_dirs]
    print(("all observation sequences available for %s:" % target_name))
    print(all_sequences)

    if not start_over:
        todo_sequences_dirs = [seqdir for seqdir in obs_sequence_dirs if seqdir.split('/')[-1] not in pass_sequences]
    else:
        todo_sequences_dirs = obs_sequence_dirs
    print("perform photometry for the following sequences:")
    print(todo_sequences_dirs)

    for obs_sequence_dir in todo_sequences_dirs:

        # prepare absolute image
        img_dir = os.path.join(obs_sequence_dir, 'uvot/image')
        img_prefix = 'sw'
        obs_sequence = obs_sequence_dir.split('/')[-1]
        img_name = img_prefix + obs_sequence + 'u' + flt + '_sk.img.gz'
        img_name2 = img_prefix + obs_sequence + 'u' + flt + '_sk.img'

        img_abs_name = os.path.join(img_dir, img_name)
        print(img_abs_name)
        if not os.path.isfile(img_abs_name):
            img_abs_name =  os.path.join(img_dir, img_name2) #both gz and non-gz files are possible 
            if not os.path.isfile(img_abs_name):
                continue

        # multiple extention images, each snapshoot in one extention
        hdu = fits.open(img_abs_name)
        N = len(hdu)
        print(("N=:", N-1))
        indices = list(range(N))
        for indice in indices[1:]:

            print(indice)
            hdr = hdu[indice].header
            date_obs = hdr['DATE-OBS']
            hdu.close()

            DATE_OBS = np.append(DATE_OBS, date_obs.split('T')[0])
            t = Time(date_obs, scale='utc', format='isot')

            segnum = obs_sequence[-3:]
            source_reg_filename, bkg_reg_filename = prepare_source_bkg_region_filename(obs_sequence, flt, indice)
            source_reg_filename = os.path.join(work_dir, source_reg_filename)
            bkg_reg_filename = os.path.join(work_dir, bkg_reg_filename)

            if os.path.exists(source_reg_filename):
                prepare_source_reg = False
            else:
                prepare_source_reg = True

            if os.path.exists(bkg_reg_filename):
                prepare_bkg_reg = False
            else:
                prepare_bkg_reg = True

            bkg_reg_flt_specific  = os.path.join( os.path.join(work_dir, 'bkg_region_manual'), 'bkg_%s.reg'%flt) 
            if os.path.exists(bkg_reg_flt_specific):
                sn_bkg_reg_template = bkg_reg_flt_specific
            else:
                sn_bkg_reg_template = os.path.join(work_dir, 'bkg.reg')


            sn_source_reg_template = os.path.join(work_dir, 'source_%sarc.reg' % (str(int(app))))

            # if source region file exists in photometry workplace directory for given target
            # then just copy the target-specific region template as source region files
            # otherwise, create source region from the generic template
            if os.path.exists(sn_source_reg_template):
                source_reg_template = sn_source_reg_template
                copy_source_reg_template = True
            else:
                source_reg_template = '/home/ping/projects/Scripts/swift_uvot_photometry/source_template.reg'
                copy_source_reg_template = False

            if os.path.exists(sn_bkg_reg_template):
                bkg_reg_template = sn_bkg_reg_template
                copy_bkg_reg_template = True # use target-specific bkg region template; the sn_bkg_reg_template will be simply copied as bkg region file of the working image  
            else:
                bkg_reg_template = '/home/ping/projects/Scripts/swift_uvot_photometry/bkg_template.reg'
                copy_bkg_reg_template = False

            ra = sn_ra
            dec = sn_dec
            app = app
            sky_in = sky_in
            sky_out = sky_out
             
            if prepare_source_reg:
                prepare_source_reg_func(ra, dec, app, source_reg_template, source_reg_filename, copy_template_source=copy_source_reg_template)
            if prepare_bkg_reg:
                prepare_bkg_reg_func(ra, dec, sky_in, sky_out, bkg_reg_template, bkg_reg_filename, copy_template_bkg=copy_bkg_reg_template)

            aperture_size = "_%sarc" % app
            log_file = os.path.join(work_dir, img_name.split('.')[0] + aperture_size + '_%s.log' % indice)
            if os.path.exists(log_file):
                os.remove(log_file)

            res_file = os.path.join(work_dir, img_name.split('.')[0]+'_res' + aperture_size + '_%s.fits' % indice)
            if os.path.exists(res_file):
                os.remove(res_file)

            phot_command = 'uvotsource image=%s+%s srcreg=%s bkgreg=%s sigma=3 outfile=%s apercorr=%s >%s' % (img_abs_name, indice, source_reg_filename, bkg_reg_filename, res_file, apercorr, log_file)
 
            print(phot_command)
            os.system(phot_command)

            img_label = img_abs_name+"["+str(indice)+"]"

            drop_this = False
            if drop_specific_obs:
                for droped_img in droped_imgs:
                    if droped_img in img_label:
                        print(("!!!! image %s droped !!!!" % droped_img))
                        drop_this = True

            if drop_this:
                continue

            JD = np.append(JD, t.jd)
            hdu_res = fits.open(res_file)
            data = hdu_res[1].data
            EXPOSURE = np.append(EXPOSURE, data['EXPOSURE'][0])
            MAG = np.append(MAG, data['MAG'][0])
            MAG_ERR_STAT = np.append(MAG_ERR_STAT, data['MAG_ERR_STAT'][0])
            MAG_ERR_SYS = np.append(MAG_ERR_SYS, data['MAG_ERR_SYS'][0])
            mag_err = np.sqrt(data['MAG_ERR_STAT'][0]**2 + data['MAG_ERR_SYS'][0]**2)
            MAG_ERR = np.append(MAG_ERR, mag_err)

            rates_corr = np.append(rates_corr,data['CORR_RATE'][0])
            rates_err_corr = np.append(rates_err_corr,data['CORR_RATE_ERR'][0])

            fluxs_aa = np.append(fluxs_aa,data['FLUX_AA'][0])
            fluxs_aa_err_stat = np.append(fluxs_aa_err_stat,data['FLUX_AA_ERR_STAT'][0])
            fluxs_aa_err_sys = np.append(fluxs_aa_err_sys, data['FLUX_AA_ERR_SYS'][0])

            images = np.append(images, img_label)
            sequences_done = np.append(sequences_done, obs_sequence_dir.split('/')[-1])
            qualflags = np.append(qualflags,data['QUALFLAG'][0]) 
            photflags = np.append(photflags,data['PHOTFLAG'][0]) 
            sss_factors = np.append(sss_factors,data['SSS_FACTOR'][0])

    print("============== New Photometry Results ==============")
    print(JD)
    print(MAG)
    print(MAG_ERR)
    print(images)

    if not start_over:
        JD_out           = np.append(old_lc['JD'],         JD)
        EXPOSURE_out     = np.append(old_lc['EXPOSURE'],   EXPOSURE)
        MAG_out          = np.append(old_lc['mag'],        MAG)
        MAG_ERR_out      = np.append(old_lc['mag_err'],    MAG_ERR)
        images_out       = np.append(old_lc['image'],      images)
        sequences_out    = np.append(old_lc['sequence'],   sequences_done)
        qualflags_out    = np.append(old_lc['qualflag'],   qualflags)
        photflags_out    = np.append(old_lc['photflag'],   photflags)
        sss_factors_out  = np.append(old_lc['sss_factor'], sss_factors)
        rates_corr_out       =np.append(old_lc['rate_corr'],        rates_corr)  
        rates_err_corr_out   =np.append(old_lc['rate_err_corr'],    rates_err_corr)
        fluxs_aa_out         =np.append(old_lc['flux_aa'],      fluxs_aa) 
        fluxs_aa_err_stat_out=np.append(old_lc['flux_aa_err_stat'],   fluxs_aa_err_stat)
        fluxs_aa_err_sys_out =np.append(old_lc['flux_aa_err_sys'],   fluxs_aa_err_sys)
        MAG_ERR_STAT_out     =np.append(old_lc['mag_err_stat'],   MAG_ERR_STAT)
        MAG_ERR_SYS_out      =np.append(old_lc['mag_err_sys'], MAG_ERR_SYS)
    else:
        JD_out = JD
        EXPOSURE_out = EXPOSURE
        MAG_out = MAG
        MAG_ERR_out = MAG_ERR
        images_out = images
        sequences_out = sequences_done
        qualflags_out = qualflags
        photflags_out = photflags
        sss_factors_out = sss_factors
        rates_corr_out = rates_corr  
        rates_err_corr_out = rates_err_corr
        fluxs_aa_out = fluxs_aa 
        fluxs_aa_err_stat_out = fluxs_aa_err_stat
        fluxs_aa_err_sys_out = fluxs_aa_err_sys
        MAG_ERR_STAT_out = MAG_ERR_STAT
        MAG_ERR_SYS_out  = MAG_ERR_SYS

    print("============= All Photometry Results ================")
    print(JD_out)
    print(MAG_out)
    print(MAG_ERR_out)
    print(qualflags_out)
    print(photflags_out)
    print(sss_factors_out)

    obs_num_real = len(JD_out)
    lc_data = np.concatenate((JD_out.reshape((obs_num_real, 1)), MAG_out.reshape(
        (obs_num_real, 1)), MAG_ERR_out.reshape((obs_num_real, 1))), axis=1)

    # save photometry result in swift photometry workplace directory
    lc_data_big = np.concatenate((JD_out.reshape((obs_num_real, 1)), EXPOSURE_out.reshape((obs_num_real, 1)), MAG_out.reshape((obs_num_real, 1)), MAG_ERR_out.reshape(
        (obs_num_real, 1)), images_out.reshape((obs_num_real, 1)), sequences_out.reshape((obs_num_real, 1)), qualflags_out.reshape((obs_num_real,1)),
        photflags_out.reshape((obs_num_real, 1)), sss_factors_out.reshape((obs_num_real,1)), rates_corr_out.reshape((obs_num_real,1)), rates_err_corr_out.reshape((obs_num_real,1)), fluxs_aa_out.reshape((obs_num_real,1)), fluxs_aa_err_stat_out.reshape((obs_num_real,1)), fluxs_aa_err_sys_out.reshape((obs_num_real,1)), MAG_ERR_STAT_out.reshape((obs_num_real,1)), MAG_ERR_SYS_out.reshape((obs_num_real,1))), axis=1)

    lc_data_table = Table(data=lc_data_big, names=['JD', 'exposure', 'mag', 'mag_err', 'image', 'sequence', 'qualflag', 'photflag', 'sss_factor', 'rate_corr', 'rate_err_corr', 'flux_aa', 'flux_aa_err_stat', 'flux_aa_err_sys', 'mag_err_stat', 'mag_err_sys'], dtype=['f8', 'f8', 'f8', 'f8', 'S120', 'S20', 'S10', 'S10', 'f8', 'f8', 'f8', 'f8','f8', 'f8', 'f8','f8'])
    lc_data_table['JD'].info.format='.5f'
    lc_data_table['exposure'].info.format='.1f'
    lc_data_table['mag'].info.format='.3f'
    lc_data_table['mag_err'].info.format='.3f'
    lc_data_table['sss_factor'].info.format='.1f'
    lc_data_table['rate_corr'].info.format = '.3f'
    lc_data_table['rate_err_corr'].info.format = '.3f'
    lc_data_table['flux_aa'].info.format = '.3e'
    lc_data_table['flux_aa_err_stat'].info.format = '.3e'
    lc_data_table['flux_aa_err_sys'].info.format = '.3e'
    lc_data_table['mag_err_stat'].info.format = '.3f'
    lc_data_table['mag_err_sys'].info.format = '.3f'


    return lc_data_table, lc_data


def ConfigParser_get_options_given_section(info_file, section):
    '''
    Output:
            info_dict:
    '''

    import configparser
    Config = configparser.ConfigParser()
    Config.read(info_file)

    info_dict = {}
    options = Config.options(section)
    for option in options:
        try:
            info_dict[option] = Config.get(section, option)
            if info_dict[option] == -1:
                print(("Skip: %s" % option))
        except:
            print(("Exception on %s!" % option))
            info_dict[option] = None
    return info_dict


def prepare_filepaths(sn_work_name, obs_id, master_dir='/home/ping'):
    '''
    prepare image directories and workplace directory

    INPUTS:
            sn_work_name: the same as the directory storing image data
            obs_id: the id assigned to specific observation
            master_dir: modified this under different machine	

    OUTPUTS:	
            obs_sequence_dirs: list of image directories for observation sequences 
            work_dir: the swift photometry workplace 	
    '''

    data_dir = os.path.join(master_dir, 'image/SWIFT_rawdownload/' + sn_work_name)
    work_dir = os.path.join(master_dir, 'photometry/SWIFT/' + sn_work_name)

    if len(obs_id) == 8:
        obs_id_format = obs_id
    else:
        obs_id_format = '0'*(8-len(str(obs_id)))+str(obs_id)

    obs_num = 999  # suppose single target maximum observations
    obs_sequences = [obs_id_format+'0' *
                     (3-len(str(x)))+str(x) for x in np.arange(obs_num)+1]

    obs_sequence_dirs = [os.path.join(data_dir, obs_sequence) for obs_sequence in obs_sequences if os.path.exists( os.path.join(data_dir, obs_sequence))]

    return obs_sequence_dirs, work_dir


def prepare_source_bkg_region_filename(obs_sequence, flt, extension):
    '''
    prepare source region file for each image
    INPUTs:
            obs_sequence: observation sequence, for example 00016848001 
            flt: UVOT filters, 'uu','bb','vv','w1','m2','w2'
            extension: extension of fits file in which the working data is stored
    '''
    source_reg_filename = "source_%s%s_%s.reg" % (obs_sequence, flt, extension)
    bkg_reg_filename = "bkg_%s%s_%s.reg" % (obs_sequence, flt, extension)
    return source_reg_filename, bkg_reg_filename

def prepare_source_reg_func(ra, dec, app, source_region_template, source_out, copy_template_source=False, verbose=0):
    '''
    prepare source region file for SWIFT photomery
    source region circle (ra,dec,app)

    INPUTs:
            ra: format 00:00:00:00
            dec: format -00:00:00
            app: aperture size in radius, unit in arcsec for example 5
            source_region_temlate: template file for source region, absolute filename with full path
            source_out: output filename for source region,absolute filename with full path
            copy_template_source: if True, just generate the region file by copying the source region file for specific supernova; otherwise create 
                                  the region file by modifying the general template region file
    '''
    if copy_template_source:
        shutil.copy(source_region_template, source_out)
    else:
        source_reg_template = open(source_region_template, 'rt').read()
        source_reg_ra_modified = source_reg_template.replace('99:99:99.99', ra)
        source_reg_radec_modified = source_reg_ra_modified.replace('-88:88:88.88', dec)
        source_reg_radecapp_modified = source_reg_radec_modified.replace('00.0000', '%s' % app)

        fid_source_region = open(source_out, 'wt')
        fid_source_region.write(source_reg_radecapp_modified)
        fid_source_region.close()

    return 1



def prepare_bkg_reg_func(ra, dec, sky_in, sky_out, bkg_region_template, bkg_out, copy_template_bkg=False, verbose=0):
    '''
    prepare background region file for SWIFT photomery
    bkckground region annulas (ra,dec,sky_in,sky_out)

    INPUTs:
            ra: format 00:00:00:00
            dec: format -00:00:00
            sky_in:	inner radius of annulas in arcsec, for example 10
            sky_out: inner radius of annulas in arcsec, for example 20
            bkg_region_template:  template file for bkg region, absolute filename with full path
            bkg_out: output filename for bkg region,absolute filename with full path
            copy_template_bkg: if True, just generate the region file by copying the bkg region file for specific supernova; otherwise create 
                                  the region file by modifying the general template region file
    '''

    if copy_template_bkg:
        shutil.copy(bkg_region_template, bkg_out)
    else:
        bkg_reg_template = open(bkg_region_template, 'rt').read()
        bkg_reg_ra_modified = bkg_reg_template.replace('99:99:99.99', ra)
        bkg_reg_radec_modified = bkg_reg_ra_modified.replace(
            '-88:88:88.88', dec)
        bkg_reg_radecskyin_modified = bkg_reg_radec_modified.replace(
            '-77.7777', '%s' % sky_in)
        bkg_reg_radecskyinskyout_modified = bkg_reg_radecskyin_modified.replace(
            '-66.6666', '%s' % sky_out)

        fid_bkg_region = open(bkg_out, 'wt')
        fid_bkg_region.write(bkg_reg_radecskyinskyout_modified)
        fid_bkg_region.close()

    return 1


def prepare_source_bkg_reg(ra, dec, app, sky_in, sky_out, source_region_template, bkg_region_template, source_out, bkg_out, copy_template_source=False, copy_template_bkg=False, prepare_source_reg=True, prepare_bkg_reg=True, verbose=0):
    '''
    prepare source region file and background region file for SWIFT photomery

    source region circle (ra,dec,app)
    bkckground region annulas (ra,dec,sky_in,sky_out)

    INPUTs:
            ra: format 00:00:00:00
            dec: format -00:00:00
            app: aperture size in radius, unit in arcsec for example 5
            sky_in:	inner radius of annulas in arcsec, for example 10
            sky_out: inner radius of annulas in arcsec, for example 20
            source_region_temlate: template file for source region, absolute filename with full path
            bkg_region_template:  template file for bkg region, absolute filename with full path
            source_out: output filename for source region,absolute filename with full path
            bkg_out: output filename for bkg region,absolute filename with full path

            copy_template_source: if True, just generate the region file by copying the source region file for specific supernova; otherwise create 
                                  the region file by modifying the general template region file
            copy_template_bkg: same philosophy as 'copy_template_source'
    '''
    if prepare_source_reg:
        if copy_template_source:
            shutil.copy(source_region_template, source_out)
        else:
            source_reg_template = open(source_region_template, 'rt').read()
            source_reg_ra_modified = source_reg_template.replace(
                '99:99:99.99', ra)
            source_reg_radec_modified = source_reg_ra_modified.replace(
                '-88:88:88.88', dec)
            source_reg_radecapp_modified = source_reg_radec_modified.replace(
                '00.0000', '%s' % app)

            fid_source_region = open(source_out, 'wt')
            fid_source_region.write(source_reg_radecapp_modified)
            fid_source_region.close()
    #else:
    #    print(("source region file %s exist" % source_out))

    if prepare_bkg_reg:
        if copy_template_bkg:
            shutil.copy(bkg_region_template, bkg_out)
        else:
            bkg_reg_template = open(bkg_region_template, 'rt').read()
            bkg_reg_ra_modified = bkg_reg_template.replace('99:99:99.99', ra)
            bkg_reg_radec_modified = bkg_reg_ra_modified.replace(
                '-88:88:88.88', dec)
            bkg_reg_radecskyin_modified = bkg_reg_radec_modified.replace(
                '-77.7777', '%s' % sky_in)
            bkg_reg_radecskyinskyout_modified = bkg_reg_radecskyin_modified.replace(
                '-66.6666', '%s' % sky_out)

            fid_bkg_region = open(bkg_out, 'wt')
            fid_bkg_region.write(bkg_reg_radecskyinskyout_modified)
            fid_bkg_region.close()
    #else:
    #    print(("bkg region file %s exists" % bkg_out))

     

def prepare_wget_download_ftp_url(target_id, seg, year, month, wanted_instrument='uvot', wanted_datatype='image', TOO=True):
    '''
    prepare wget command line according to description in http://swift.gsfc.nasa.gov/archive/
    SWIFT ftp: ftp://legacy.gsfc.nasa.gov/swift/data/obs/

    INPUT:
            target_id: target id assigned to specific observation
            seg:	   segment number in given observation with specific target_id
            year:      which year the observation conducted, int for example 2015, or str for example '2015'
            month:     which month the observation conducted, int for example 9 or str for example '9'
            wanted_instrument: instrument level selection, 'bat','xrt','uvot' or 'all'. If 'all', then all data will be downloaded
            wanted_datatype: data type selection for uvot instrument, 'event','hk','image','products' or 'all'. If 'all' then all uvor data will be downloaded 
    '''
    # check input parameters
    if isinstance(year, str):
        year_str = year
    elif isinstance(year, int):
        year_str = str(year)
    else:
        raise IOError(
            'invalid format for input parameter year, only str and int acceptable')

    if len(year_str) != 4:
        raise IOError('please check the input for year')

    month_str = '0'*(2-len(str(month)))+str(month)
    months_acceptable = ['01', '02', '03', '04', '05',
                         '06', '07', '08', '09', '10', '11', '12']

    if month_str not in months_acceptable:
        raise IOError('please check the input for month')

    year_month = year_str + '_' + month_str

    if TOO:
        target_id_str = '0'*(8-len(str(target_id)))+str(target_id)
    else:
        if len(str(target_id)) > 5:
            raise ValueError("please check the input of the target ID")
        # This is not necessarily true but correct for fill-in targets at least which matters when the script is created
        target_id_str = '031' + '0'*(5-len(str(target_id)))+str(target_id)

    seg_str = '0'*(3-len(str(seg)))+str(seg)
    sequence = target_id_str + seg_str

#	wget_command = "wget -r ftp://legacy.gsfc.nasa.gov/swift/data/obs/%s/%s/uvot/image/"%(year_month, sequence)
    if wanted_instrument not in ['bat', 'xrt', 'uvot', 'all']:
        raise ValueError("please check the input for wanted_instrument")

    if wanted_instrument == 'all':
        wget_command = "wget -r ftp://legacy.gsfc.nasa.gov/swift/data/obs/%s/%s/" % (
            year_month, sequence)
    else:
        wget_command = "wget -r ftp://legacy.gsfc.nasa.gov/swift/data/obs/%s/%s/%s/" % (
            year_month, sequence, wanted_instrument)

        if wanted_instrument == 'uvot' and wanted_datatype != 'all':
            wget_command = wget_command + '%s/' % wanted_datatype

    return wget_command


def wget_download_ftp_system_action(wget_command, out_dir, delete_downloaded_filetree=False, to_save_data_level=0):
    '''
    call for system command for wget

    INPUTS:
            wget_command:	command line which can be excuted by system shell 
            out_dir:	directory to store the downloaded data
            delete_downloaded_filetree: delete the downloading file system tree or not
            to_save_data_level: 0 or positive integer smaller than total levels
                                For example, we have ftp://a/b/c/d/e
                                if 0, the lowest level data will be saved, i.e a/b/c/d/e
                                if 3, the third level directory from top will be saved, i.e a/b/c
    '''
    try:
        os.system(wget_command)
    except:
        print("wget command failure...")
        return

    filetree = wget_command.split('//')[1]

    # if the ftp link ending with '/', remove it first
    if filetree[-1] == '/':
        filetree = filetree[:-1]

    dir_segs = filetree.split('/')
    N = len(dir_segs)

    if to_save_data_level == 0:
        to_save_data = filetree
    else:
        if to_save_data_level > N:
            raise IOError(
                "to_save_data_level can not excess total data levels")

        to_save_data = ''
        for indice in range(to_save_data_level):
            to_save_data = os.path.join(to_save_data, dir_segs[indice])

    if os.path.exists(os.path.join(out_dir, to_save_data.split('/')[-1])):
        shutil.rmtree(os.path.join(out_dir, to_save_data.split('/')[-1]))

    shutil.move(to_save_data, out_dir)

    if delete_downloaded_filetree:
        shutil.rmtree(dir_segs[0])


def parse_html_table_new(url_addr, which_one=0):
    from python_tableparser_new import HTMLTableParser
    hp = HTMLTableParser()
    table = hp.parse_url(url_addr)[which_one][1] # Grabbing the table from the tuple
    return table


def parse_html_table(url_addr, which_one=0, use_first_row_as_header=False, header_out=None, dtype_out=None):
    '''
    utilize modified script by Aquil H. Abdullah
    The original code http://devwiki.beloblotskiy.com/index.php5/Generic_HTML_Table_parser_(python)

    INPUTS:
            url_addr: The url address pointing to the html page containing table
            which_one: Here you need to know exactly how many tables are on the page and which one you want. For example, 0 mean first table
            use_first_row_as_header: treat the first row in the parsed table as table header
            header_out: output the table with header? if no, None; if yes, header_out = ['colname1', 'colname2', ...]
    '''
    from urllib.request import Request, urlopen
    from urllib.error import URLError
    from python_tableparser import TableParser
    req = Request(url_addr)
    url = urlopen(req)
    tp = TableParser()
    tp.feed(url.read())

    tables_extracted = tp.get_tables()
    my_table = tables_extracted[which_one]

    try:
        my_table.remove([''])
    except:
        print("no empty str as row")

    try:
        my_table.remove([])
    except:
        print("no empty row")

    Nx1 = len(my_table[0])
    Nx2 = len(my_table[-1])
    if Nx1 != Nx2:
        raise ValueError("not well formatted table...")
    else:
        Nx = Nx1

    Ny = len(my_table)

    if header_out is not None:
        table_colnames = header_out
    elif use_first_row_as_header:
        table_colnames = my_table[0]
    else:
        table_colnames = ['col'+str(i) for i in range(Nx)]

    table_array = np.array(my_table)
    table_array.reshape(Ny, Nx)

    colcharmax = []
    for i in range(Nx):
        colcharnums = [len(x) for x in table_array[:, i]]
        colcharmax.append(np.max(colcharnums))

    if dtype_out is None:
        dtypes = ['S'+str(n) for n in colcharmax]
    else:
        if len(dtype_out) == Nx:
            dtypes = dtype_out
        else:
            raise IOError("%s columns in table, only %s dtypes entries specified..." % (
                Nx, len(dtype_out)))

    out_table = Table(names=table_colnames, dtype=dtypes)

    for row in my_table:
        out_table.add_row(row)

    return out_table


def search_for_new_quicklook_entry(obs_log_table, target_id):
    '''
    search for new observation entry in obs_log_table. 
    the observation entry matched with given target ID will be returned

    INPUTS:
            obs_log_table: parsed table from swift quick-look page 
            target_id:
    OUTPUT:
            exist: True or False
            sequences_selected: list of obs sequence found in quick-look page for given target_id
            redlevel_selected: 
    '''

    target_id = str(target_id)

    info = obs_log_table
    sequences_all = info[info.colnames[0]]
    redlevels_all = info[info.colnames[1]]

    sequence_selected = []
    redlevel_selected = []
    exist = False

    for i, sequence in enumerate(sequences_all):
        redlevel = redlevels_all[i]
        if target_id in sequence:
            if not exist:
                exist = True
            sequence_selected.append(sequence.replace(' ', ''))
            redlevel_selected.append(redlevel.replace(' ', ''))

    return exist, sequence_selected, redlevel_selected


def check_data_archived_status(data_master_dir, sequence):
    '''
    check whether the data for specific sequence already downloaded and archived 

    INPUTS:
            data_dir:
            sequence:	
    '''
    data_check = os.path.join(data_master_dir, sequence)

    if os.path.exists(data_check):
        downloaded = True
    else:
        downloaded = False

    return downloaded


def prepare_wget_download_quicklook_url(sequence, redlevel):
    '''
    the url link to download data from SWIFT quick-look page is well formatte, for example: 
    http://swift.gsfc.nasa.gov/cgi-bin/sdc/download.tar?method=tar&button=Download&top=%2Fdata%2Fswift%2Fsw00034736002.002&dir=sw00034736002.002&dataset=swiftql&%2Fuvot%2Fimage=on

    only '00034736002.002' part changes for different target, and share the format 'sequence(target_id+seg).reduction_level'

    INPUTS:	
                    sequence: 
                    redlevel: swift data reduction level, listed in second column of SWIFT quick-look table
    '''

    urltemplate = 'http://swift.gsfc.nasa.gov/cgi-bin/sdc/download.tar?method=tar&button=Download&top=%2Fdata%2Fswift%2Fsw00000000000.999&dir=sw00000000000.999&dataset=swiftql&%2Fuvot%2Fimage=on'

    # check input for sequence
    if len(sequence) != 11:
        raise IOError(
            "the sequence number for swift is a number with 11 digit, please  check...")
    else:
        sequence_complete = sequence

    redlevel_complete = '0'*(3-len(str(redlevel)))+str(redlevel)

    url_seq = urltemplate.replace('00000000000', sequence_complete)
    url_redlevel = url_seq.replace('999', redlevel_complete)

    url_final = url_redlevel

    return url_final


def sum_swift_images(img_in, img_out, exclude='DEFAULT', ignoreframetime=False, trybestsum=True, check_aspcorr=True, verbose=1, chatter=3):
    '''
    sum swift images with uvotimmsum 	
    refer to the following for the instrctions:
    http://www.swift.ac.uk/analysis/uvot/image.php

    For help on uvotimsum, try
    $ fhelp uvotimsum

    uvotimsum infile=<filename> outfile=<filename> method=<string> pixsize=<float> expmap=<filename> exclude=<string> maskfile=<filename> weightfile=<filename> clobber=<boolean> history=<boolean> cleanup=<boolean> chatter=<integer>

    Notes: Please make sure the images are well aligned before summing
           Pay attention to 'exclude' option:
            (exclude = DEFAULT) [string]
            Comma delimited list of HDU names or numbers (or numbered ranges) to exclude from the sum.
            The special value ASPCORR:NONE can be used to exclude extensions that do not have the ASPCORR keyword. ......
            The special value DEFAULT is treated as 
            exclude=ASPCORR:NONE,DUPEXPID,SETTLING.

    See also: 'uvotskycoor' for aspect correction

    INPUTS:
            img_in:
            img_out:
            exclude:
            ignoreframetime:
            trybestsum: when variation of frametimes are beyond the tolerence and ignoreframetime is no, try to sum frames with the same frametimes and highest total exposure time
    '''

    hdu = fits.open(img_in)

    aspcorrs = []
    frametimes = []
    datesobs = []
    exptimes = []

    for i in range(len(hdu)):
        if i == 0:
            continue

        hdr = hdu[i].header
        if check_aspcorr:
            aspcorrs.append(hdr['ASPCORR'])
        frametimes.append(hdr['FRAMTIME'])
        datesobs.append(hdr['DATE-OBS'])
        exptimes.append(hdr['EXPOSURE'])

    if verbose:
        print(img_in)
        print(("ASPCORR:", aspcorrs))
        print(("FRAMTIME:", frametimes))
        print(("DATE-OBS:", datesobs))
        print(("EXPOSURE:", exptimes))

    if check_aspcorr and ('DIRECT' not in aspcorrs):
        raise IOError("none of the frames has been aligned to WCS")
    else:
        if len(frametimes) > 1:

            frametimes = np.array(frametimes)
            frametimes_diff = np.abs(frametimes - frametimes[0])/frametimes[0]
            if verbose:
                print("the relative difference in frametimes:")
                print(frametimes_diff)

            violence = np.array([timediff < 0.02 for timediff in frametimes_diff])
            if False not in violence:
                uvotimsum_command = "uvotimsum infile=%s outfile=%s exclude=%s chatter=%s clobber=yes < /dev/null > /tmp/uvotimsum_output.log 2>&1" % (img_in, img_out, exclude, chatter)
                os.system(uvotimsum_command)
            else:
                if ignoreframetime:
                    ignoreframetime = 'yes'
                    uvotimsum_command = "uvotimsum infile=%s outfile=%s exclude=%s ignoreframetime=%s chatter=%s clobber=yes < /dev/null > /tmp/uvotimsum_output.log 2>&1" % (img_in, img_out, exclude, ignoreframetime, chatter)
                    os.system(uvotimsum_command)
                else:
                    if trybestsum:
                        exptimes = np.array(exptimes)
                        set1_indice = np.where(violence)[0] + 1
                        set1_totime = np.sum(exptimes[violence])

                        set2_indice = np.where(~violence)[0] + 1
                        set2_totime = np.sum(exptimes[~violence])

                        if set1_totime < set2_totime:
                            set_exclude = set1_indice
                        else:
                            set_exclude = set2_indice

                        exclude_fnums = ','.join([str(idnum) for idnum in set_exclude])
                        if verbose:
                            print(("different frametimes exist and frames %s will be excluded" % exclude_fnums))

                        uvotimsum_command = "uvotimsum infile=%s outfile=%s exclude=%s chatter=%s clobber=yes < /dev/null > /tmp/uvotimsum_output.log 2>&1" % (img_in, img_out, exclude_fnums,  chatter)
                        os.system(uvotimsum_command)

                    else:
                        raise ValueError( "the input images can not be summed under given constraint")

        else:
            uvotimsum_command = "uvotimsum infile=%s outfile=%s exclude=%s chatter=%s clobber=yes < /dev/null > /tmp/uvotimsum_output.log 2>&1" % ( img_in, img_out, exclude, chatter)
            os.system(uvotimsum_command)


def append_swift_images_diff_segs_and_sum(inimage_list, outimage, cleanup=True, chatter=0, clobber=0, check_aspcorr=True):
    '''
    sometimes you want to sum images from different file: first select one image and append all other files to the selected one, then sum all frames together

    the tool will be used for doing the appending is fappend.
    '''

    baseimg = inimage_list[0]
    baseimg_sum = "baseimg_sum_temp.fits"

    if not check_aspcorr:
        exclude='NONE'
    else:
        exclude='DEFAULT'

    sum_swift_images(baseimg, baseimg_sum, exclude=exclude, check_aspcorr=check_aspcorr, chatter=chatter)

    for i, inimg in enumerate(inimage_list):
        if i == 0:
            continue

        inimg_sum = "inimage%s_sum_temp.fits" % i
        sum_swift_images(inimg, inimg_sum, exclude=exclude, check_aspcorr=check_aspcorr, chatter=chatter)
        fappend_command = "fappend %s %s" % (inimg_sum, baseimg_sum)

        os.system(fappend_command)

        if cleanup:
            os.remove(inimg_sum)

    if clobber and os.path.exists(outimage):
        os.remove(outimage)

    sum_swift_images(baseimg_sum, outimage, exclude='NONE', check_aspcorr=False, chatter=chatter)
    if cleanup:
        os.remove(baseimg_sum)


if __name__ == "__main__":
    import sys

    img_in = sys.argv[1]
    img_out = sys.argv[2]

    sum_swift_images(img_in, img_out)
