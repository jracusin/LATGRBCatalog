#!/usr/bin/env python
"""
------------------------------------------------------------------------

Scripts to collect data for LAT 2nd GRB Catalog

------------------------------------------------------------------------
"""

import numpy as np
import re
from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
from astropy.io import fits
import os
import urllib

def make_cat():

	cat=load_BAtool_results()
	cat=load_GRBOX(cat)
	cat=load_Nicola(cat)

	return cat

def read_catalog():

	### downloaded this manually
	url='http://glast-ground.slac.stanford.edu/LATBA/GRBCatalogReview.jsp'
	filename='LAT_GRB_cat.html'

	f=open(filename,'r')
	lines=f.readlines()

	l=0
	grb=[]
	met=[]
	date=[]
	time=[]
	ra=[]
	dec=[]
	ts=[]
	real=[]
	for line in lines: 
		if "CandidateReview.jsp?NAME=" in line:
			grb=np.append(grb,'GRB'+re.split('<|>|/',line)[4])
			met=np.append(met,float(re.split('<td>|</td>',lines[l+2])[1]))
			date=np.append(date,re.split('<td>|</td>',lines[l+3])[1])
			time=np.append(time,re.split('<td>|</td>',lines[l+4])[1])
			ra=np.append(ra,float(re.split('<td>|</td>',lines[l+5])[1]))
			dec=np.append(dec,float(re.split('<td>|</td>',lines[l+6])[1]))
			ts=np.append(ts,float(re.split('<td>|</td>',lines[l+7])[1]))
			real=np.append(real,re.split('title="|">',lines[l+9])[1])
			# if isreal == 'YES': r=True
			# if isreal == 'NO': r=False
			# real=np.append(real,r)
			
		l=l+1

	d=np.core.defchararray.add(date,' ')
	trigtime=np.core.defchararray.add(d,time)

# 	c1=fits.Column(name='grb',format='10A',array=grb)
# 	c2=fits.Column(name='met',format='D',array=met)
# 	c3=fits.Column(name='trigtime',format='23A',array=trigtime)
# 	c4=fits.Column(name='ra',format='D',array=ra)
# 	c5=fits.Column(name='dec',format='D',array=dec)
# 	c6=fits.Column(name='ts',format='D',array=ts)
# 	c7=fits.Column(name='real',format='L',array=real)
# 	coldefs = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7])
# 	tbhdu = fits.BinTableHDU.from_columns(coldefs)
# 	cat=tbhdu.data
# 	s=np.argsort(cat['grb'])
#	cat=cat[s]

	n=len(grb)
#	catlist=[]
#	for i in range(n):
#		g={'grb':grb[i],'met':met[i],'trigtime':trigtime[i],'ra':ra[i],'dec':dec[i],'ts':ts[i],'real':real[i]}
#		catlist.append(g)

	s=np.argsort(grb)
#	cattable=Table([grb[s],met[s],trigtime[s],ra[s],dec[s],ts[s],real[s]],names=('grb','met','trigtime','ra','dec','ts','real'))

	return grb[s],met[s],trigtime[s],ra[s],dec[s],ts[s],real[s]

def load_BAtool_results():

	from decimal import Decimal
	grb,met,trigtime,ra,dec,ts,real=read_catalog()
	n=len(grb)
	z=np.zeros(n)
	na=np.chararray((n),itemsize=3)
	na[:]='---'
	na10=np.chararray((n),itemsize=10)
	na10[:]='----------'

	cat=Table([grb,na10,met,trigtime,z,z,z,ts,real,z,z,z,\
		z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,na,na],\
		names=['grb','GRBNAME','met','trigtime','ra','dec','err','ts','real','LLE_sig','theta','zenith',\
		'source_events','transient_events','source_events_1GeV','transient_events_1GeV',\
		'max_photon_energy','max_photon_time','photon_flux','photon_flux_err',\
		'energy_flux','energy_flux_err',\
		'LAT_T90','LLE_T90','GBM_T90','LAT_GCN',\
		'z','afterglow','BAT_det'])
		#'PH_index','temp_index',
		
	n=len(cat)
	for i in range(n):
		grb0=re.split('GRB',grb[i])[1]#cat[i]['grb']
		file='/nfs/slac/g/ki/ki08/kocevski/LATBA/DATA/GRBOUT/LTF/'+grb0+'/Prompt/results_'+grb0+'.txt'
#		file='results_'+grb+'.txt'
		if os.path.exists(file):
			print file
			f=open(file,'r')
			lines=f.readlines()

			l=0
			for line in lines:
				tmp=re.split('=|\n',line)
				if l == 0: 
					l=l+1
					continue
				var=tmp[0].replace(' ','')
				val=tmp[1].replace(' ','')
#				print var,val
				if 'FindSrc_DEC1' in line:
					cat[i]['dec']=float(val)
					cat[i]['ra']=float(re.split('=',lines[l+2])[1])
					cat[i]['err']=float(re.split('=',lines[l+1])[1])
				if 'LLE_DetMaxSign' in line:
					cat[i]['LLE_sig']=round(float(val),2)
				# if 'Galacticb' in line:
				# 	cat[i]['Galacticb']=float(val)
				# 	cat[i]['Galacticl']=float(re.split('=',lines[l+1])[1])
				if 'THETA' in line:
					cat[i]['theta']=round(float(val),1)
				if 'ZENITH' in line:
					cat[i]['zenith']=round(float(val),1)
				if var=='NumberOfEvents_S_ROI':
					cat[i]['source_events']=round(float(val))
				if var=='NumberOfEvents_T_ROI':
					cat[i]['transient_events']=round(float(val))
				if var=='NumberOfEvents1GeV_S_ROI':
					cat[i]['source_events_1GeV']=round(float(val))
				if var=='NumberOfEvents1GeV_T_ROI':
					cat[i]['transient_events_1GeV']=round(float(val))
				if 'MaxPhotonEnergy' in line:
					cat[i]['max_photon_energy']=round(float(val),1)
				if 'MaxPhotonEnergyTime' in line:
					cat[i]['max_photon_time']=round(float(val),1)
				if var=='LIKE_MY_FLUX':
					cat[i]['photon_flux']='%.2E' % Decimal(float(val))
					cat[i]['photon_flux_err']='%.2E' % Decimal(float(re.split('=',lines[l+3])[1]))
					cat[i]['energy_flux']='%.2E' % Decimal(float(re.split('=',lines[l+1])[1]))
					cat[i]['energy_flux_err']='%.2E' % Decimal(float(re.split('=',lines[l+2])[1]))



				l=l+1

	ascii.write(cat,'LAT_Cat_Table.dat')
	cat=cat[cat['real']=='YES']
	cat.write('LAT_Cat_Table.html',format='ascii.html')
	return cat

def load_GRBOX(cat):

	from astropy.coordinates import SkyCoord
	from astropy import units as u
	from astropy.time import Time

	#cat=Table.read('LAT_Cat_Table.html',format='ascii.html')
	grbox=ascii.read('grboxtxt.php',format='fixed_width',\
		names=['GRB','UT','T90','RA','DEC','z','det'],data_start=1,\
		col_starts=(0,8,17,23,36,49,55))

	ngrbox=len(grbox['GRB'])

	met0='2001-01-01 00:00:00'
#	gname=[]
#	gfrac=[]
	met=[]
	for i in range(ngrbox):

		if grbox['UT'][i] != '':
			year=grbox['GRB'][i][0:2]
			month=grbox['GRB'][i][2:4]
			day=grbox['GRB'][i][4:6]
			if ':' in grbox['UT'][i][0:2]: q=-1
			else: q=0
			hr=float(grbox['UT'][i][0:2+q])
			mn=float(grbox['UT'][i][3+q:5+q])
			sec=grbox['UT'][i][6+q:8+q]
			if sec != '': sec=float(sec)
			else: sec=0.
			utc='20'+year+'-'+month+'-'+day+' '+grbox[i]['UT']

			times=[met0,utc]
			t=Time(times,format='iso',scale='utc')
			tdiff=t[1]-t[0]
			tdiff.format='sec'
			met=np.append(met,tdiff.value)
		else: met=np.append(met,0.)

	ncat=len(cat)
	cat['afterglow'].astype(str)
	cat[:]['afterglow']='---'

	for i in range(ncat):
		c=SkyCoord(ra=cat[i]['ra']*u.deg,dec=cat[i]['dec']*u.deg)
		d=SkyCoord(ra=grbox['RA']*u.deg,dec=grbox['DEC']*u.deg)
		dist=c.separation(d)
		fsep=abs(cat[i]['met']-met)
		w=np.where((fsep < 60) & (dist<5.*u.deg))
		if (len(w[0])>0):
#			print cat[i]['grb'],grbox[w]['GRB'][0],grbox[w]['z'][0]
			if float(grbox[w]['z'])>0:
				cat[i]['z']=grbox[w]['z'][0]
				cat[i]['GRBNAME']='GRB'+grbox[w[0]]['GRB'][0]
#			print cat[i]['afterglow'],grbox[w[0]]['det']
			# print
			# print grbox[w[0]]['det'][0]
			if grbox[w[0]]['det'].mask == False:
				cat[i]['afterglow']=grbox[w[0]]['det'][0]
				cat[i]['GRBNAME']='GRB'+grbox[w[0]]['GRB'][0]
#			print cat[i]['afterglow']
		# if 'X' in grbox[w]['det']: cat[i]['x_afterglow']='Y'
		# if 'O' in grbox[w]['det']: cat[i]['o_afterglow']='Y'
		# if 'R' in grbox[w]['det']: cat[i]['r_afterglow']='Y'


	cat.write('LAT_Cat_Table.html',format='ascii.html')
	cat.write('LAT_Cat_Table.dat',format='ascii')

	return cat

def load_Nicola(cat):

	import CatalogTools

	results,rtable = CatalogTools.parseCatalogFile_Nicola('CatalogTable_Nicola.txt')

	# match GRB names, fill in table

	ngrb=[]
	for i in range(len(cat)): ngrb=np.append(ngrb,'GRB'+rtable['GCNNAME'][i])

	m1,m2=match(cat['grb'],ngrb)

#	r=Table(results)[m2]

	cat.add_columns(rtable.columns.values())

	cat.write('LAT_Cat_Table.html',format='ascii.html')
	cat.write('LAT_Cat_Table.dat',format='ascii')

	return cat

def load_BAT(cat):

	url='http://swift.gsfc.nasa.gov/results/batgrbcat/summary_cflux/summary_general_info/summary_general.txt'
	file='BAT_summary_general.txt'

	if not os.path.exists(file):
		urllib.urlretrieve(url,file)

	f=open(file,'r')
	lines=f.readlines()
	lines=lines[23:]
	bgrb=[]
	for line in lines:
		d=np.array(re.split('\|| ',line))
		d=d[d!='']
		bgrb=np.append(bgrb,d[0])

	w=np.where(cat['GRBNAME'] != '----------')
	cgrb=np.array(cat['GRBNAME'][w[0]])
	m1,m2=match(cgrb,bgrb)
	m1=w[0][m1]

	cat['BAT_det']='N'
	cat['BAT_det'][m1]='Y'

	cat.write('LAT_Cat_Table.html',format='ascii.html')
	cat.write('LAT_Cat_Table.dat',format='ascii')

	return cat

def match(a, b):

	m1 = []
	m2 = []
	for i in range(len(a)):
		w=np.where(b == a[i])
		if len(w[0]) > 0:
			m2=np.append(m2,w[0])
			m1=np.append(m1,i)

	m1=m1.astype(int)
	m2=m2.astype(int)

	return m1,m2

