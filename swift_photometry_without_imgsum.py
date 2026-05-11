import os
import sys
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table, vstack

from SWIFT_UVOT_photometry_toolkit import swift_photometry_uvot_flt


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

    def_sky_in = 15
    parser.add_option('--sky_in', dest="SKY_IN", type="float", default=def_sky_in,
                      help="Inner radius size of background annulus in arcsecond [%s]" % def_sky_in)

    def_sky_out = 30
    parser.add_option('--sky_out', dest="SKY_OUT", type="float", default=def_sky_out,
                      help="Out radius size of background annulus in arcsecond [%s]" % def_sky_out)

    def_startover = False
    parser.add_option('--start_over', dest="start_over", action="store_true", default=def_startover,
                      help="whether or not start over(delete all existed photometry result file)[%s]" % def_startover)


    def_renewreg = False
    parser.add_option('--renew_reg', dest="renew_reg", action="store_true", default=def_renewreg,
                      help="whether or not renew the source/bkg region files[%s]" % def_renewreg)

    def_tplphotonly = False
    parser.add_option('--tplphot_only', dest="tplphot_only", action="store_true",
                      default=def_tplphotonly, help="only do photometry for template images[%s]" % def_tplphotonly)


    def_skip_check_aspcorr = False
    parser.add_option('--skipcheck', dest='skip_check_aspcorr', action="store_true", 
                      default=def_skip_check_aspcorr, help="skip checking ASPCORR before summing the images[%s]"%def_skip_check_aspcorr)


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
    renewreg = options.renew_reg

    tplonly = options.tplphot_only
    skipcheck = options.skip_check_aspcorr
    if skipcheck:
        check_aspcorr_prior_summing = False
    else:
        check_aspcorr_prior_summing = True


    if RA == '00:00:00':
        sne_master_table = Table.read('source_info_SNRADec.csv', format='ascii.csv')
        if target_name not in sne_master_table['Name']:
            raise IOError("target RA required...")
        else:
            RA = sne_master_table['RA'][sne_master_table['Name']== target_name].data[0]

    if DEC == '00:00:00':
        sne_master_table = Table.read('source_info_SNRADec.csv', format='ascii.csv')
        if target_name not in sne_master_table['Name']:
            raise IOError("target DEC required...")
        else:
            DEC = sne_master_table['Dec'][sne_master_table['Name']== target_name].data[0]

    if sn_work_name == '':
        sn_work_name = target_name

    sniper_dir = None  # save the final lc data to a place

    if flts == '':
        swift_uvot_flts = ['w2', 'm2', 'w1', 'uu', 'bb', 'vv']
    else:
        swift_uvot_flts = flts.split(',')




    master_dir = '/home/ping/'  # /mnt/md0/
    data_dir = os.path.join(master_dir, 'image/SWIFT_rawdownload/' + sn_work_name)
    
    if len(obs_id_list) == 1 and obs_id_list[0] == '':
        # use all obsid directories under the data_dir
        obs_id_list = np.unique([f[:8] for f in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, f)) and f.startswith('0') and len(f) == 11])
        
        print('The observation id is not explicitly provided; Use all obsids under %s' % data_dir) 
        if len(obs_id_list) == 0:
            raise ValueError("no valid observation id found...")
        print('Found the following obsids:', obs_id_list)

    work_dir = os.path.join(master_dir, 'photometry/SWIFT/' + sn_work_name)

    if not tplonly:
        for flt in swift_uvot_flts:
            photrettable = None
            lc_data = None
            for obs_id in obs_id_list:
                photrettable_temp, lc_data_temp  = swift_photometry_uvot_flt(target_name, obs_id, flt, RA, DEC, master_dir=master_dir,  app=app, sky_in=sky_in, sky_out=sky_out, start_over=startover)

                if photrettable is None:
                    photrettable = photrettable_temp
                else:
                    photrettable = vstack([photrettable, photrettable_temp])
                if lc_data is None:
                    lc_data = lc_data_temp
                else:
                    lc_data = np.vstack((lc_data, lc_data_temp))

            photret_table_file = os.path.join(work_dir, target_name+'_SWs-' + flt + '_uvotsource_ret_no_imgsum.txt')
            photrettable.write(photret_table_file, format='ascii.fixed_width', overwrite=True)

            lc_file_sndir = os.path.join(work_dir, target_name + '_SWs-' + flt + '_noimgsum.txt')
            fid_sndir = open(lc_file_sndir, 'wt')
            np.savetxt(fid_sndir, lc_data, fmt="%8.4f %6.3f %6.3f")

            # don't create empty light curve for SNIPER, which will stop displaying ascii files on SNIPER
            # also remove magnitudes not realistic
            if len(lc_data) > 0 and sniper_dir is not None:
                lc_file_sniper = os.path.join(sniper_dir, target_name + '_SWs-' + flt + '.txt')
                lc_data_sniper = lc_data[lc_data[:, 1] < 30, :]
                fid_sniper = open(lc_file_sniper, 'wt')
                np.savetxt(fid_sniper, lc_data_sniper, fmt="%8.4f %6.3f %6.3f")


