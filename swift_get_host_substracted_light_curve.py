#! /usr/bin/python

import os
import numpy as np
from astropy.table import Table
from SWIFT_UVOT_photometry_toolkit import swift_rate2mag


def get_hostfree_mags_flt(target_name, flt, master_dir='/home/ping/', which_version='sum', sniper_master_dir=None):

    work_dir = os.path.join(master_dir, 'photometry/SWIFT/%s/' % target_name)

    hostphot_file = os.path.join(work_dir, target_name+'_tpl_uvotsource_ret.txt')
    if which_version == 'sum':
        snphot_file = os.path.join(work_dir, target_name+'_SWs-' + flt + '_uvotsource_ret.txt') 
    elif which_version == 'no_sum':
        snphot_file = os.path.join(work_dir, target_name+'_SWs-' + flt + '_uvotsource_ret_no_imgsum.txt')
    elif which_version == 'hybrid':
        snphot_file = os.path.join(work_dir, target_name+'_SWs-' + flt + '_uvotsource_ret_hybrid.txt')
    else:
         raise ValueError('which_version should be sum or no_sum or hybrid') 

    if not os.path.exists(hostphot_file):
        raise IOError("host photometry not available")

    if not os.path.exists(snphot_file):
        print("no photometry file in filter %s"%flt)
        return

    hostphot = Table.read(hostphot_file, format='ascii.fixed_width')
    snphot = Table.read(snphot_file, format='ascii.fixed_width')

    host_rate = hostphot[hostphot['flt'] == flt]['rate_corr'][0]
    host_rate_err = hostphot[hostphot['flt'] == flt]['rate_err_corr'][0]

    JD = []
    MAG = []
    MAGERR = []
    MAGUPL = [] #upper limit
    for obsphot in snphot:
        jd = obsphot['JD']
        snhostrate = obsphot['rate_corr']
        snhostrate_err = obsphot['rate_err_corr']

        snrate_err = np.sqrt(host_rate_err**2 + snhostrate_err**2)
        magupl, temp = swift_rate2mag(snrate_err*3, flt)        

        if snhostrate < host_rate:
            print("host rate > snrate ...")
            mag = 99.99
            magerr = 99.99
        else:
            snrate = snhostrate-host_rate
            mag, magerr = swift_rate2mag(snrate, flt, rate_err=snrate_err)

        JD.append(jd)
        MAG.append(mag)
        MAGERR.append(magerr)
        MAGUPL.append(magupl)

    lcdata = np.transpose(np.array([JD, MAG, MAGERR, MAGUPL]))

    lc_file_sndir = os.path.join(work_dir, target_name + '_SWs-' + flt + '-HF-'+which_version+'.txt')
    np.savetxt(lc_file_sndir, lcdata, fmt='%12.5f %8.3f %8.3f %8.3f')

    if sniper_master_dir is not None:
        sniper_dir = os.path.join(sniper_master_dir, target_name)
        lc_file_sniper = os.path.join(sniper_dir, target_name + '_SWs-' + flt + '-HF-'+which_version+'.txt')
        np.savetxt(lc_file_sniper, lcdata, fmt='%12.5f %8.3f %8.3f %8.3f')



if __name__ == "__main__":


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

    for flt in swift_uvot_flts:
        get_hostfree_mags_flt(target_name, flt, master_dir='/home/ping/', which_version=version)
