#! /usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd

def get_band_df(path):
	outcar = open(path+'OUTCAR-2', 'r')
	kpoints = open(path+'KPOINTS-2','r')
	    
	kpt=0
	nkpt=0
	band=0
	nbands=0
	ef=0
	start=0

	kx=0
	ky=0
	kz=0

	kdata = []

	for line in outcar:
		sp = line.split()
		if ((len(sp)>13) and (sp[13]=='NBANDS=')):
			nbands = int(float(sp[14]))
		if((len(sp)>13) and (sp[1]=='NKPTS')):
			nkpt = int(float(sp[3]))
		if ((len(sp)>3) and (sp[0]=='E-fermi')):
			ef = float(sp[2])
			start = 1

		if(len(kdata)==(nbands*nkpt)):
			start=0
		
		if ((len(sp)==3) and (start==1) and (sp[0] != 'k-point') and (sp[0]!= 'band') and (sp[0]!='spin')):
			band = int(float(sp[0]))
			energy = float(sp[1])
			data = [kpt, band, energy, energy-ef, kx, ky, kz]
			kdata.append(data)
		
		if ((len(sp)>5) and (sp[0]=='k-point') and (start == 1)):
			kpt = int(float(sp[1]))
			kx = float(sp[3])
			ky = float(sp[4])
			kz = float(sp[5])

	outcar_df = pd.DataFrame(kdata)
	outcar_df.rename(columns = {0:'kpt'
				   ,1:'band'
				   ,2:'raw_energy'
				   ,3:'energy'
				   ,4:'kx'
				   ,5:'ky'
				   ,6:'kz'}, inplace = True)

	kpt_names = []

	for line in kpoints:
		sp = line.split()
		if ((len(sp) > 3) and (sp[3]=='!')):
			kx = float(sp[0])
			ky = float(sp[1])
			kz = float(sp[2])
			kpt_names.append([kx,ky,kz, sp[4]])
		
	#kpt_names_df = pd.DataFrame(kpt_names).drop_duplicates()
	#kpt_names_df.rename(columns = {0:'kx', 1:'ky', 2:'kz', 3:'k_name'},inplace = True)

	#outcar_df0 = outcar_df.copy()
	#outcar_df = pd.merge(outcar_df0, kpt_names_df, how = 'left', on= ['kx','ky','kz'])        
	for i in range(len(kpt_names)):
		outcar_df.loc[(outcar_df.kx == round(kpt_names[i][0],4))&(outcar_df.ky == round(kpt_names[i][1],4))&(outcar_df.kz == round(kpt_names[i][2],4)), 'k_name'] = kpt_names[i][3]
	outcar_df.loc[outcar_df.k_name == '\Gamma', 'k_name'] = '$\Gamma$'

	high_sym_pts = outcar_df.loc[outcar_df.k_name.notnull()][['kpt', 'k_name']].drop_duplicates()

	sym_list = high_sym_pts.kpt.tolist()
	name_list = high_sym_pts.k_name.tolist()
	kpt_to_drop = []
	for i in range(len(sym_list)):
		if ((i > 0) and (name_list[i] == name_list[i-1])):
			kpt_to_drop.append(sym_list[i])
		
	out_df = outcar_df.copy()

	for i in range(len(kpt_to_drop)):
		out_df = out_df.loc[out_df.kpt != kpt_to_drop[i]].copy()

	kpt_df = out_df.loc[(out_df.band == 1)][['kpt']].copy().reset_index()
	kpt_df = kpt_df.drop('index',1)
	kpt_df = kpt_df.reset_index()
	kpt_df.rename(columns = {'index':'new_kpt'},inplace = True)
	kpt_df['new_kpt'] = kpt_df['new_kpt'] + 1

	out_df0 = out_df.copy()
	out_df = pd.merge(out_df0, kpt_df, how = 'left', on = 'kpt')

	out_df_to_plot = out_df.loc[out_df.new_kpt.notnull()].copy()
	out_df_to_plot = get_kcoords(out_df_to_plot)

	return out_df_to_plot

def get_kcoords(df):
        kcoords = np.asarray(df[['new_kpt','kx','ky','kz','k_name']].drop_duplicates())
        name_list = kcoords.T[4]
        kdist = []
        kcum = 0.0
        for i in range(len(kcoords)):
            if ((i == 0) or ((kcoords[i-1][4] in name_list) and (kcoords[i][4] in name_list))):
                kstep = 0.0
            else:
                kstep = ((kcoords[i][1]-kcoords[i-1][1])**2+(kcoords[i][2]-kcoords[i-1][2])**2+(kcoords[i][3]-kcoords[i-1][3])**2)**0.5
            kcum += kstep
            kdist.append([kcoords[i][0], kstep, kcum])

        kdf = pd.DataFrame(kdist)
        kdf.rename(columns = {0:'new_kpt', 1:'kstep', 2:'kcoord'}, inplace = True)
        df1 = pd.merge(df, kdf, on = ['new_kpt'])

        return df1


def get_cbm_vbm(out_df):
	pos_en = out_df.loc[out_df.energy > 0.0].copy()
	neg_en = out_df.loc[out_df.energy <=0.0].copy()
	cbm = min(pos_en.energy)
	vbm = max(neg_en.energy)
	band_gap = cbm - vbm

	out_df['cbm'] = 0
	out_df['vbm'] = 0
	out_df.loc[out_df.energy == cbm, 'cbm'] = 1
	out_df.loc[out_df.energy == vbm, 'vbm'] = 1

	return out_df


def choose_k_pts(df, min_kpt, max_kpt):
	df1 = get_cbm_vbm(df.loc[(df.kpt > (min_kpt - 1))&(df.kpt < (max_kpt +1))]).copy()
	min_new_kpt = min(df1.new_kpt)
	df1['new_kpt_0'] = df1['new_kpt']
	df1['new_kpt'] = df1['new_kpt_0'] - (min_new_kpt - 1)
	return df1


def get_single_band(df, band):
	nparray = np.asarray(df)
	band_list = []
	for i in range(len(nparray)):
		if ((nparray[i,1] == band)):
			band_list.append(nparray[i,3])
	return band_list

def get_high_sym_k_pts(df):
	high_sym_pts1 = df.loc[df.k_name.notnull()][['kcoord', 'k_name']].drop_duplicates()
	high_sym_k_list = high_sym_pts1.kcoord.tolist()
	high_sym_name_list = high_sym_pts1.k_name.tolist()
	return [high_sym_k_list, high_sym_name_list]

def output_bs(df, file_name):
    nbands = max(df.band)
    kcoords = df[['new_kpt','kcoord']].drop_duplicates().kcoord.tolist()
    file = open(file_name, 'w')

    for j in range(nbands):
        band_to_plot = get_single_band(df, (j+1))
        for i in range(len(band_to_plot)):
            file.write(str(kcoords[i]) + ' ' + str(band_to_plot[i]) + '\n')
        file.write('\n')
    file.close()
