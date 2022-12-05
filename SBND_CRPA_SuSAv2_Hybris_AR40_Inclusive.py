#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import sympy as smp
import matplotlib.pyplot as plt
from scipy import integrate
import pandas as pd
from Measurement import Measurement, LinComb


# # Data Imports and Preparation

# ## Importing and Naming Data

# In[2]:


#Reads text files containing data for each flux
#Very unconventional and cumbersome method
#It might be better to use numpy instead of pandas

#Full flux
flux_0 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_0/2d_bins_flux.txt", sep=" ", header=None)
flux_1 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_1/2d_bins_flux.txt", sep=" ", header=None)
flux_2 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_2/2d_bins_flux.txt", sep=" ", header=None)
flux_3 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_3/2d_bins_flux.txt", sep=" ", header=None)
flux_4 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_4/2d_bins_flux.txt", sep=" ", header=None)
flux_5 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_5/2d_bins_flux.txt", sep=" ", header=None)
flux_6 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_6/2d_bins_flux.txt", sep=" ", header=None)
flux_7 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_7/2d_bins_flux.txt", sep=" ", header=None)

el_flux_0 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_0/El_bins_flux.txt", sep=" ", header=None)
el_flux_1 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_1/El_bins_flux.txt", sep=" ", header=None)
el_flux_2 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_2/El_bins_flux.txt", sep=" ", header=None)
el_flux_3 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_3/El_bins_flux.txt", sep=" ", header=None)
el_flux_4 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_4/El_bins_flux.txt", sep=" ", header=None)
el_flux_5 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_5/El_bins_flux.txt", sep=" ", header=None)
el_flux_6 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_6/El_bins_flux.txt", sep=" ", header=None)
el_flux_7 = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid/Flux_7/El_bins_flux.txt", sep=" ", header=None)


#CCQE events only
flux_0_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_0_CCQE/2d_bins_flux_CCQE.txt", sep=" ", header=None)
flux_1_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_1_CCQE/2d_bins_flux_CCQE.txt", sep=" ", header=None)
flux_2_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_2_CCQE/2d_bins_flux_CCQE.txt", sep=" ", header=None)
flux_3_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_3_CCQE/2d_bins_flux_CCQE.txt", sep=" ", header=None)
flux_4_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_4_CCQE/2d_bins_flux_CCQE.txt", sep=" ", header=None)
flux_5_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_5_CCQE/2d_bins_flux_CCQE.txt", sep=" ", header=None)
flux_6_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_6_CCQE/2d_bins_flux_CCQE.txt", sep=" ", header=None)
flux_7_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_7_CCQE/2d_bins_flux_CCQE.txt", sep=" ", header=None)

el_flux_0_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_0_CCQE/El_bins_flux_CCQE.txt", sep=" ", header=None)
el_flux_1_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_1_CCQE/El_bins_flux_CCQE.txt", sep=" ", header=None)
el_flux_2_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_2_CCQE/El_bins_flux_CCQE.txt", sep=" ", header=None)
el_flux_3_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_3_CCQE/El_bins_flux_CCQE.txt", sep=" ", header=None)
el_flux_4_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_4_CCQE/El_bins_flux_CCQE.txt", sep=" ", header=None)
el_flux_5_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_5_CCQE/El_bins_flux_CCQE.txt", sep=" ", header=None)
el_flux_6_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_6_CCQE/El_bins_flux_CCQE.txt", sep=" ", header=None)
el_flux_7_ccqe = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQE/Flux_7_CCQE/El_bins_flux_CCQE.txt", sep=" ", header=None)


#CCQE like
flux_0_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_0_CCQELike/2d_bins_flux_CCQELike.txt", sep=" ", header=None)
flux_1_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_1_CCQELike/2d_bins_flux_CCQELike.txt", sep=" ", header=None)
flux_2_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_2_CCQELike/2d_bins_flux_CCQELike.txt", sep=" ", header=None)
flux_3_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_3_CCQELike/2d_bins_flux_CCQELike.txt", sep=" ", header=None)
flux_4_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_4_CCQELike/2d_bins_flux_CCQELike.txt", sep=" ", header=None)
flux_5_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_5_CCQELike/2d_bins_flux_CCQELike.txt", sep=" ", header=None)
flux_6_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_6_CCQELike/2d_bins_flux_CCQELike.txt", sep=" ", header=None)
flux_7_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_7_CCQELike/2d_bins_flux_CCQELike.txt", sep=" ", header=None)

el_flux_0_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_0_CCQELike/El_bins_flux_CCQELike.txt", sep=" ", header=None)
el_flux_1_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_1_CCQELike/El_bins_flux_CCQELike.txt", sep=" ", header=None)
el_flux_2_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_2_CCQELike/El_bins_flux_CCQELike.txt", sep=" ", header=None)
el_flux_3_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_3_CCQELike/El_bins_flux_CCQELike.txt", sep=" ", header=None)
el_flux_4_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_4_CCQELike/El_bins_flux_CCQELike.txt", sep=" ", header=None)
el_flux_5_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_5_CCQELike/El_bins_flux_CCQELike.txt", sep=" ", header=None)
el_flux_6_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_6_CCQELike/El_bins_flux_CCQELike.txt", sep=" ", header=None)
el_flux_7_ccqe_like = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_7_CCQELike/El_bins_flux_CCQELike.txt", sep=" ", header=None)


#CCQE All
flux_0_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_0/2d_bins_flux_CCQEAll.txt", sep=" ", header=None)
flux_1_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_1/2d_bins_flux_CCQEAll.txt", sep=" ", header=None)
flux_2_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_2/2d_bins_flux_CCQEAll.txt", sep=" ", header=None)
flux_3_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_3/2d_bins_flux_CCQEAll.txt", sep=" ", header=None)
flux_4_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_4/2d_bins_flux_CCQEAll.txt", sep=" ", header=None)
flux_5_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_5/2d_bins_flux_CCQEAll.txt", sep=" ", header=None)
flux_6_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_6/2d_bins_flux_CCQEAll.txt", sep=" ", header=None)
flux_7_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_7/2d_bins_flux_CCQEAll.txt", sep=" ", header=None)

el_flux_0_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_0/El_bins_flux_CCQEAll.txt", sep=" ", header=None)
el_flux_1_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_1/El_bins_flux_CCQEAll.txt", sep=" ", header=None)
el_flux_2_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_2/El_bins_flux_CCQEAll.txt", sep=" ", header=None)
el_flux_3_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_3/El_bins_flux_CCQEAll.txt", sep=" ", header=None)
el_flux_4_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_4/El_bins_flux_CCQEAll.txt", sep=" ", header=None)
el_flux_5_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_5/El_bins_flux_CCQEAll.txt", sep=" ", header=None)
el_flux_6_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_6/El_bins_flux_CCQEAll.txt", sep=" ", header=None)
el_flux_7_ccqe_all = pd.read_csv("./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQEAll/Flux_7/El_bins_flux_CCQEAll.txt", sep=" ", header=None)


#Naming the columns
flux_0.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_1.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_2.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_3.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_4.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_5.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_6.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_7.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']

el_flux_0.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_1.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_2.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_3.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_4.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_5.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_6.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_7.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']


flux_0_ccqe.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_1_ccqe.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_2_ccqe.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_3_ccqe.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_4_ccqe.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_5_ccqe.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_6_ccqe.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_7_ccqe.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']

el_flux_0_ccqe.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_1_ccqe.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_2_ccqe.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_3_ccqe.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_4_ccqe.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_5_ccqe.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_6_ccqe.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_7_ccqe.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']


flux_0_ccqe_like.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_1_ccqe_like.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_2_ccqe_like.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_3_ccqe_like.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_4_ccqe_like.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_5_ccqe_like.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_6_ccqe_like.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_7_ccqe_like.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']

el_flux_0_ccqe_like.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_1_ccqe_like.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_2_ccqe_like.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_3_ccqe_like.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_4_ccqe_like.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_5_ccqe_like.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_6_ccqe_like.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_7_ccqe_like.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']


flux_0_ccqe_all.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_1_ccqe_all.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_2_ccqe_all.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_3_ccqe_all.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_4_ccqe_all.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_5_ccqe_all.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_6_ccqe_all.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']
flux_7_ccqe_all.columns = ['cos_low', ' ', 'cos_up', ' ', 'E_low', ' ', 'E_up', ' ', 'diff_cross']

el_flux_0_ccqe_all.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_1_ccqe_all.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_2_ccqe_all.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_3_ccqe_all.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_4_ccqe_all.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_5_ccqe_all.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_6_ccqe_all.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']
el_flux_7_ccqe_all.columns = ['E_low', ' ', 'E_up', ' ', 'diff_cross']


# ## Flux 0 Cosine Angular Bins

# In[3]:


#Separating out data according to angles for each flux
#It's not efficient since there's so many different conditions we might want to set
#Should create a loop to automate making these
flux_0_1 = flux_0.loc[lambda df: df['cos_low'] == -1, :] #Angle -1.0 < x < -0.9
flux_0_2 = flux_0.loc[lambda df: df['cos_low'] == -.9, :] #Angle -0.9 < x < -0.8
flux_0_3 = flux_0.loc[lambda df: df['cos_low'] == -.8, :] #Angle -0.8 < x < -0.7
flux_0_4 = flux_0.loc[lambda df: df['cos_low'] == -.7, :] #Angle -0.7 < x < -0.6
flux_0_5 = flux_0.loc[lambda df: df['cos_low'] == -.6, :] #Angle -0.6 < x < -0.5
flux_0_6 = flux_0.loc[lambda df: df['cos_low'] == -.5, :] #Angle -0.5 < x < -0.4
flux_0_7 = flux_0.loc[lambda df: df['cos_low'] == -.4, :] #Angle -0.4 < x < -0.3
flux_0_8 = flux_0.loc[lambda df: df['cos_low'] == -.3, :] #Angle -0.3 < x < -0.2
flux_0_9 = flux_0.loc[lambda df: df['cos_low'] == -.2, :] #Angle -0.2 < x < -0.1
flux_0_10 = flux_0.loc[lambda df: df['cos_low'] == -.1, :] #Angle -0.1 < x < 0.0
flux_0_11 = flux_0.loc[lambda df: df['cos_low'] == 0, :] #Angle 0.0 < x < 0.1
flux_0_12 = flux_0.loc[lambda df: df['cos_low'] == .1, :] #Angle 0.1 < x < 0.2 
flux_0_13 = flux_0.loc[lambda df: df['cos_low'] == .2, :] #Angle 0.2 < x < 0.3
flux_0_14 = flux_0.loc[lambda df: df['cos_low'] == .3, :] #Angle 0.3 < x < 0.4
flux_0_15 = flux_0.loc[lambda df: df['cos_low'] == .4, :] #Angle 0.4 < x < 0.5
flux_0_16 = flux_0.loc[lambda df: df['cos_low'] == .5, :] #Angle 0.5 < x < 0.6
flux_0_17 = flux_0.loc[lambda df: df['cos_low'] == .6, :] #Angle 0.6 < x < 0.7
flux_0_18 = flux_0.loc[lambda df: df['cos_low'] == .7, :] #Angle 0.7 < x < 0.8
flux_0_19 = flux_0.loc[lambda df: df['cos_low'] == .8, :] #Angle 0.8 < x < 0.9
flux_0_20 = flux_0.loc[lambda df: df['cos_low'] == .9, :] #Angle 0.9 < x < 1.0

flux_0_1_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == -1, :]
flux_0_2_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == -.9, :]
flux_0_3_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == -.8, :]
flux_0_4_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == -.7, :]
flux_0_5_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == -.6, :]
flux_0_6_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == -.5, :]
flux_0_7_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == -.4, :]
flux_0_8_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == -.3, :]
flux_0_9_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == -.2, :]
flux_0_10_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == -.1, :]
flux_0_11_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == 0, :]
flux_0_12_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == .1, :]
flux_0_13_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == .2, :]
flux_0_14_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == .3, :]
flux_0_15_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == .4, :]
flux_0_16_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == .5, :]
flux_0_17_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == .6, :]
flux_0_18_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == .7, :]
flux_0_19_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == .8, :]
flux_0_20_ccqe = flux_0_ccqe.loc[lambda df: df['cos_low'] == .9, :]

flux_0_1_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == -1, :]
flux_0_2_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == -.9, :]
flux_0_3_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == -.8, :]
flux_0_4_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == -.7, :]
flux_0_5_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == -.6, :]
flux_0_6_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == -.5, :]
flux_0_7_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == -.4, :]
flux_0_8_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == -.3, :]
flux_0_9_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == -.2, :]
flux_0_10_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == -.1, :]
flux_0_11_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == 0, :]
flux_0_12_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == .1, :]
flux_0_13_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == .2, :]
flux_0_14_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == .3, :]
flux_0_15_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == .4, :]
flux_0_16_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == .5, :]
flux_0_17_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == .6, :]
flux_0_18_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == .7, :]
flux_0_19_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == .8, :]
flux_0_20_ccqe_like = flux_0_ccqe_like.loc[lambda df: df['cos_low'] == .9, :]

flux_0_1_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == -1, :]
flux_0_2_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == -.9, :]
flux_0_3_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == -.8, :]
flux_0_4_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == -.7, :]
flux_0_5_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == -.6, :]
flux_0_6_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == -.5, :]
flux_0_7_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == -.4, :]
flux_0_8_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == -.3, :]
flux_0_9_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == -.2, :]
flux_0_10_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == -.1, :]
flux_0_11_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == 0, :]
flux_0_12_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == .1, :]
flux_0_13_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == .2, :]
flux_0_14_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == .3, :]
flux_0_15_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == .4, :]
flux_0_16_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == .5, :]
flux_0_17_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == .6, :]
flux_0_18_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == .7, :]
flux_0_19_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == .8, :]
flux_0_20_ccqe_all = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]


# ## Flux 1 Cosine Angular Bins

# In[4]:


flux_1_1 = flux_1.loc[lambda df: df['cos_low'] == -1, :] #Angle -1.0 < x < -0.9
flux_1_2 = flux_1.loc[lambda df: df['cos_low'] == -.9, :] #Angle -0.9 < x < -0.8
flux_1_3 = flux_1.loc[lambda df: df['cos_low'] == -.8, :] #Angle -0.8 < x < -0.7
flux_1_4 = flux_1.loc[lambda df: df['cos_low'] == -.7, :] #Angle -0.7 < x < -0.6
flux_1_5 = flux_1.loc[lambda df: df['cos_low'] == -.6, :] #Angle -0.6 < x < -0.5
flux_1_6 = flux_1.loc[lambda df: df['cos_low'] == -.5, :] #Angle -0.5 < x < -0.4
flux_1_7 = flux_1.loc[lambda df: df['cos_low'] == -.4, :] #Angle -0.4 < x < -0.3
flux_1_8 = flux_1.loc[lambda df: df['cos_low'] == -.3, :] #Angle -0.3 < x < -0.2
flux_1_9 = flux_1.loc[lambda df: df['cos_low'] == -.2, :] #Angle -0.2 < x < -0.1
flux_1_10 = flux_1.loc[lambda df: df['cos_low'] == -.1, :] #Angle -0.1 < x < 0.0
flux_1_11 = flux_1.loc[lambda df: df['cos_low'] == 0, :] #Angle 0.0 < x < 0.1
flux_1_12 = flux_1.loc[lambda df: df['cos_low'] == .1, :] #Angle 0.1 < x < 0.2 
flux_1_13 = flux_1.loc[lambda df: df['cos_low'] == .2, :] #Angle 0.2 < x < 0.3
flux_1_14 = flux_1.loc[lambda df: df['cos_low'] == .3, :] #Angle 0.3 < x < 0.4
flux_1_15 = flux_1.loc[lambda df: df['cos_low'] == .4, :] #Angle 0.4 < x < 0.5
flux_1_16 = flux_1.loc[lambda df: df['cos_low'] == .5, :] #Angle 0.5 < x < 0.6
flux_1_17 = flux_1.loc[lambda df: df['cos_low'] == .6, :] #Angle 0.6 < x < 0.7
flux_1_18 = flux_1.loc[lambda df: df['cos_low'] == .7, :] #Angle 0.7 < x < 0.8
flux_1_19 = flux_1.loc[lambda df: df['cos_low'] == .8, :] #Angle 0.8 < x < 0.9
flux_1_20 = flux_1.loc[lambda df: df['cos_low'] == .9, :] #Angle 0.9 < x < 1.0

flux_1_1_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == -1, :]
flux_1_2_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == -.9, :]
flux_1_3_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == -.8, :]
flux_1_4_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == -.7, :]
flux_1_5_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == -.6, :]
flux_1_6_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == -.5, :]
flux_1_7_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == -.4, :]
flux_1_8_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == -.3, :]
flux_1_9_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == -.2, :]
flux_1_10_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == -.1, :]
flux_1_11_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == 0, :]
flux_1_12_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == .1, :]
flux_1_13_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == .2, :]
flux_1_14_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == .3, :]
flux_1_15_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == .4, :]
flux_1_16_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == .5, :]
flux_1_17_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == .6, :]
flux_1_18_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == .7, :]
flux_1_19_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == .8, :]
flux_1_20_ccqe = flux_1_ccqe.loc[lambda df: df['cos_low'] == .9, :]

flux_1_1_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == -1, :]
flux_1_2_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == -.9, :]
flux_1_3_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == -.8, :]
flux_1_4_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == -.7, :]
flux_1_5_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == -.6, :]
flux_1_6_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == -.5, :]
flux_1_7_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == -.4, :]
flux_1_8_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == -.3, :]
flux_1_9_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == -.2, :]
flux_1_10_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == -.1, :]
flux_1_11_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == 0, :]
flux_1_12_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == .1, :]
flux_1_13_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == .2, :]
flux_1_14_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == .3, :]
flux_1_15_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == .4, :]
flux_1_16_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == .5, :]
flux_1_17_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == .6, :]
flux_1_18_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == .7, :]
flux_1_19_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == .8, :]
flux_1_20_ccqe_like = flux_1_ccqe_like.loc[lambda df: df['cos_low'] == .9, :]

flux_1_1_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == -1, :]
flux_1_2_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == -.9, :]
flux_1_3_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == -.8, :]
flux_1_4_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == -.7, :]
flux_1_5_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == -.6, :]
flux_1_6_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == -.5, :]
flux_1_7_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == -.4, :]
flux_1_8_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == -.3, :]
flux_1_9_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == -.2, :]
flux_1_10_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == -.1, :]
flux_1_11_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == 0, :]
flux_1_12_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == .1, :]
flux_1_13_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == .2, :]
flux_1_14_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == .3, :]
flux_1_15_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == .4, :]
flux_1_16_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == .5, :]
flux_1_17_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == .6, :]
flux_1_18_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == .7, :]
flux_1_19_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == .8, :]
flux_1_20_ccqe_all = flux_1_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]


# ## Flux 2 Cosine Angular Bins

# In[5]:


flux_2_1 = flux_2.loc[lambda df: df['cos_low'] == -1, :] #Angle -1.0 < x < -0.9
flux_2_2 = flux_2.loc[lambda df: df['cos_low'] == -.9, :] #Angle -0.9 < x < -0.8
flux_2_3 = flux_2.loc[lambda df: df['cos_low'] == -.8, :] #Angle -0.8 < x < -0.7
flux_2_4 = flux_2.loc[lambda df: df['cos_low'] == -.7, :] #Angle -0.7 < x < -0.6
flux_2_5 = flux_2.loc[lambda df: df['cos_low'] == -.6, :] #Angle -0.6 < x < -0.5
flux_2_6 = flux_2.loc[lambda df: df['cos_low'] == -.5, :] #Angle -0.5 < x < -0.4
flux_2_7 = flux_2.loc[lambda df: df['cos_low'] == -.4, :] #Angle -0.4 < x < -0.3
flux_2_8 = flux_2.loc[lambda df: df['cos_low'] == -.3, :] #Angle -0.3 < x < -0.2
flux_2_9 = flux_2.loc[lambda df: df['cos_low'] == -.2, :] #Angle -0.2 < x < -0.1
flux_2_10 = flux_2.loc[lambda df: df['cos_low'] == -.1, :] #Angle -0.1 < x < 0.0
flux_2_11 = flux_2.loc[lambda df: df['cos_low'] == 0, :] #Angle 0.0 < x < 0.1
flux_2_12 = flux_2.loc[lambda df: df['cos_low'] == .1, :] #Angle 0.1 < x < 0.2 
flux_2_13 = flux_2.loc[lambda df: df['cos_low'] == .2, :] #Angle 0.2 < x < 0.3
flux_2_14 = flux_2.loc[lambda df: df['cos_low'] == .3, :] #Angle 0.3 < x < 0.4
flux_2_15 = flux_2.loc[lambda df: df['cos_low'] == .4, :] #Angle 0.4 < x < 0.5
flux_2_16 = flux_2.loc[lambda df: df['cos_low'] == .5, :] #Angle 0.5 < x < 0.6
flux_2_17 = flux_2.loc[lambda df: df['cos_low'] == .6, :] #Angle 0.6 < x < 0.7
flux_2_18 = flux_2.loc[lambda df: df['cos_low'] == .7, :] #Angle 0.7 < x < 0.8
flux_2_19 = flux_2.loc[lambda df: df['cos_low'] == .8, :] #Angle 0.8 < x < 0.9
flux_2_20 = flux_2.loc[lambda df: df['cos_low'] == .9, :] #Angle 0.9 < x < 1.0

flux_2_1_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == -1, :]
flux_2_2_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == -.9, :]
flux_2_3_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == -.8, :]
flux_2_4_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == -.7, :]
flux_2_5_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == -.6, :]
flux_2_6_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == -.5, :]
flux_2_7_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == -.4, :]
flux_2_8_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == -.3, :]
flux_2_9_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == -.2, :]
flux_2_10_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == -.1, :]
flux_2_11_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == 0, :]
flux_2_12_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == .1, :]
flux_2_13_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == .2, :]
flux_2_14_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == .3, :]
flux_2_15_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == .4, :]
flux_2_16_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == .5, :]
flux_2_17_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == .6, :]
flux_2_18_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == .7, :]
flux_2_19_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == .8, :]
flux_2_20_ccqe = flux_2_ccqe.loc[lambda df: df['cos_low'] == .9, :]

flux_2_1_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == -1, :]
flux_2_2_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == -.9, :]
flux_2_3_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == -.8, :]
flux_2_4_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == -.7, :]
flux_2_5_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == -.6, :]
flux_2_6_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == -.5, :]
flux_2_7_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == -.4, :]
flux_2_8_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == -.3, :]
flux_2_9_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == -.2, :]
flux_2_10_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == -.1, :]
flux_2_11_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == 0, :]
flux_2_12_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == .1, :]
flux_2_13_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == .2, :]
flux_2_14_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == .3, :]
flux_2_15_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == .4, :]
flux_2_16_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == .5, :]
flux_2_17_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == .6, :]
flux_2_18_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == .7, :]
flux_2_19_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == .8, :]
flux_2_20_ccqe_like = flux_2_ccqe_like.loc[lambda df: df['cos_low'] == .9, :]

flux_2_1_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == -1, :]
flux_2_2_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == -.9, :]
flux_2_3_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == -.8, :]
flux_2_4_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == -.7, :]
flux_2_5_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == -.6, :]
flux_2_6_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == -.5, :]
flux_2_7_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == -.4, :]
flux_2_8_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == -.3, :]
flux_2_9_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == -.2, :]
flux_2_10_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == -.1, :]
flux_2_11_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == 0, :]
flux_2_12_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == .1, :]
flux_2_13_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == .2, :]
flux_2_14_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == .3, :]
flux_2_15_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == .4, :]
flux_2_16_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == .5, :]
flux_2_17_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == .6, :]
flux_2_18_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == .7, :]
flux_2_19_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == .8, :]
flux_2_20_ccqe_all = flux_2_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]


# ## Flux 3 Cosine Angular Bins

# In[6]:


flux_3_1 = flux_3.loc[lambda df: df['cos_low'] == -1, :] #Angle -1.0 < x < -0.9
flux_3_2 = flux_3.loc[lambda df: df['cos_low'] == -.9, :] #Angle -0.9 < x < -0.8
flux_3_3 = flux_3.loc[lambda df: df['cos_low'] == -.8, :] #Angle -0.8 < x < -0.7
flux_3_4 = flux_3.loc[lambda df: df['cos_low'] == -.7, :] #Angle -0.7 < x < -0.6
flux_3_5 = flux_3.loc[lambda df: df['cos_low'] == -.6, :] #Angle -0.6 < x < -0.5
flux_3_6 = flux_3.loc[lambda df: df['cos_low'] == -.5, :] #Angle -0.5 < x < -0.4
flux_3_7 = flux_3.loc[lambda df: df['cos_low'] == -.4, :] #Angle -0.4 < x < -0.3
flux_3_8 = flux_3.loc[lambda df: df['cos_low'] == -.3, :] #Angle -0.3 < x < -0.2
flux_3_9 = flux_3.loc[lambda df: df['cos_low'] == -.2, :] #Angle -0.2 < x < -0.1
flux_3_10 = flux_3.loc[lambda df: df['cos_low'] == -.1, :] #Angle -0.1 < x < 0.0
flux_3_11 = flux_3.loc[lambda df: df['cos_low'] == 0, :] #Angle 0.0 < x < 0.1
flux_3_12 = flux_3.loc[lambda df: df['cos_low'] == .1, :] #Angle 0.1 < x < 0.2 
flux_3_13 = flux_3.loc[lambda df: df['cos_low'] == .2, :] #Angle 0.2 < x < 0.3
flux_3_14 = flux_3.loc[lambda df: df['cos_low'] == .3, :] #Angle 0.3 < x < 0.4
flux_3_15 = flux_3.loc[lambda df: df['cos_low'] == .4, :] #Angle 0.4 < x < 0.5
flux_3_16 = flux_3.loc[lambda df: df['cos_low'] == .5, :] #Angle 0.5 < x < 0.6
flux_3_17 = flux_3.loc[lambda df: df['cos_low'] == .6, :] #Angle 0.6 < x < 0.7
flux_3_18 = flux_3.loc[lambda df: df['cos_low'] == .7, :] #Angle 0.7 < x < 0.8
flux_3_19 = flux_3.loc[lambda df: df['cos_low'] == .8, :] #Angle 0.8 < x < 0.9
flux_3_20 = flux_3.loc[lambda df: df['cos_low'] == .9, :] #Angle 0.9 < x < 1.0

flux_3_1_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == -1, :]
flux_3_2_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == -.9, :]
flux_3_3_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == -.8, :]
flux_3_4_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == -.7, :]
flux_3_5_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == -.6, :]
flux_3_6_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == -.5, :]
flux_3_7_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == -.4, :]
flux_3_8_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == -.3, :]
flux_3_9_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == -.2, :]
flux_3_10_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == -.1, :]
flux_3_11_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == 0, :]
flux_3_12_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == .1, :]
flux_3_13_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == .2, :]
flux_3_14_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == .3, :]
flux_3_15_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == .4, :]
flux_3_16_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == .5, :]
flux_3_17_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == .6, :]
flux_3_18_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == .7, :]
flux_3_19_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == .8, :]
flux_3_20_ccqe = flux_3_ccqe.loc[lambda df: df['cos_low'] == .9, :]

flux_3_1_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == -1, :]
flux_3_2_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == -.9, :]
flux_3_3_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == -.8, :]
flux_3_4_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == -.7, :]
flux_3_5_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == -.6, :]
flux_3_6_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == -.5, :]
flux_3_7_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == -.4, :]
flux_3_8_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == -.3, :]
flux_3_9_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == -.2, :]
flux_3_10_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == -.1, :]
flux_3_11_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == 0, :]
flux_3_12_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == .1, :]
flux_3_13_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == .2, :]
flux_3_14_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == .3, :]
flux_3_15_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == .4, :]
flux_3_16_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == .5, :]
flux_3_17_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == .6, :]
flux_3_18_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == .7, :]
flux_3_19_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == .8, :]
flux_3_20_ccqe_like = flux_3_ccqe_like.loc[lambda df: df['cos_low'] == .9, :]

flux_3_1_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == -1, :]
flux_3_2_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == -.9, :]
flux_3_3_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == -.8, :]
flux_3_4_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == -.7, :]
flux_3_5_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == -.6, :]
flux_3_6_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == -.5, :]
flux_3_7_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == -.4, :]
flux_3_8_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == -.3, :]
flux_3_9_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == -.2, :]
flux_3_10_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == -.1, :]
flux_3_11_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == 0, :]
flux_3_12_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == .1, :]
flux_3_13_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == .2, :]
flux_3_14_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == .3, :]
flux_3_15_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == .4, :]
flux_3_16_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == .5, :]
flux_3_17_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == .6, :]
flux_3_18_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == .7, :]
flux_3_19_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == .8, :]
flux_3_20_ccqe_all = flux_3_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]


# ## Flux 4 Cosine Angular Bins

# In[7]:


flux_4_1 = flux_4.loc[lambda df: df['cos_low'] == -1, :] #Angle -1.0 < x < -0.9
flux_4_2 = flux_4.loc[lambda df: df['cos_low'] == -.9, :] #Angle -0.9 < x < -0.8
flux_4_3 = flux_4.loc[lambda df: df['cos_low'] == -.8, :] #Angle -0.8 < x < -0.7
flux_4_4 = flux_4.loc[lambda df: df['cos_low'] == -.7, :] #Angle -0.7 < x < -0.6
flux_4_5 = flux_4.loc[lambda df: df['cos_low'] == -.6, :] #Angle -0.6 < x < -0.5
flux_4_6 = flux_4.loc[lambda df: df['cos_low'] == -.5, :] #Angle -0.5 < x < -0.4
flux_4_7 = flux_4.loc[lambda df: df['cos_low'] == -.4, :] #Angle -0.4 < x < -0.3
flux_4_8 = flux_4.loc[lambda df: df['cos_low'] == -.3, :] #Angle -0.3 < x < -0.2
flux_4_9 = flux_4.loc[lambda df: df['cos_low'] == -.2, :] #Angle -0.2 < x < -0.1
flux_4_10 = flux_4.loc[lambda df: df['cos_low'] == -.1, :] #Angle -0.1 < x < 0.0
flux_4_11 = flux_4.loc[lambda df: df['cos_low'] == 0, :] #Angle 0.0 < x < 0.1
flux_4_12 = flux_4.loc[lambda df: df['cos_low'] == .1, :] #Angle 0.1 < x < 0.2 
flux_4_13 = flux_4.loc[lambda df: df['cos_low'] == .2, :] #Angle 0.2 < x < 0.3
flux_4_14 = flux_4.loc[lambda df: df['cos_low'] == .3, :] #Angle 0.3 < x < 0.4
flux_4_15 = flux_4.loc[lambda df: df['cos_low'] == .4, :] #Angle 0.4 < x < 0.5
flux_4_16 = flux_4.loc[lambda df: df['cos_low'] == .5, :] #Angle 0.5 < x < 0.6
flux_4_17 = flux_4.loc[lambda df: df['cos_low'] == .6, :] #Angle 0.6 < x < 0.7
flux_4_18 = flux_4.loc[lambda df: df['cos_low'] == .7, :] #Angle 0.7 < x < 0.8
flux_4_19 = flux_4.loc[lambda df: df['cos_low'] == .8, :] #Angle 0.8 < x < 0.9
flux_4_20 = flux_4.loc[lambda df: df['cos_low'] == .9, :] #Angle 0.9 < x < 1.0

flux_4_1_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == -1, :]
flux_4_2_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == -.9, :]
flux_4_3_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == -.8, :]
flux_4_4_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == -.7, :]
flux_4_5_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == -.6, :]
flux_4_6_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == -.5, :]
flux_4_7_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == -.4, :]
flux_4_8_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == -.3, :]
flux_4_9_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == -.2, :]
flux_4_10_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == -.1, :]
flux_4_11_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == 0, :]
flux_4_12_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == .1, :]
flux_4_13_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == .2, :]
flux_4_14_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == .3, :]
flux_4_15_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == .4, :]
flux_4_16_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == .5, :]
flux_4_17_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == .6, :]
flux_4_18_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == .7, :]
flux_4_19_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == .8, :]
flux_4_20_ccqe = flux_4_ccqe.loc[lambda df: df['cos_low'] == .9, :]

flux_4_1_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == -1, :]
flux_4_2_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == -.9, :]
flux_4_3_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == -.8, :]
flux_4_4_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == -.7, :]
flux_4_5_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == -.6, :]
flux_4_6_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == -.5, :]
flux_4_7_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == -.4, :]
flux_4_8_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == -.3, :]
flux_4_9_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == -.2, :]
flux_4_10_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == -.1, :]
flux_4_11_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == 0, :]
flux_4_12_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == .1, :]
flux_4_13_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == .2, :]
flux_4_14_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == .3, :]
flux_4_15_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == .4, :]
flux_4_16_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == .5, :]
flux_4_17_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == .6, :]
flux_4_18_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == .7, :]
flux_4_19_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == .8, :]
flux_4_20_ccqe_like = flux_4_ccqe_like.loc[lambda df: df['cos_low'] == .9, :]

flux_4_1_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == -1, :]
flux_4_2_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == -.9, :]
flux_4_3_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == -.8, :]
flux_4_4_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == -.7, :]
flux_4_5_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == -.6, :]
flux_4_6_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == -.5, :]
flux_4_7_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == -.4, :]
flux_4_8_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == -.3, :]
flux_4_9_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == -.2, :]
flux_4_10_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == -.1, :]
flux_4_11_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == 0, :]
flux_4_12_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == .1, :]
flux_4_13_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == .2, :]
flux_4_14_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == .3, :]
flux_4_15_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == .4, :]
flux_4_16_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == .5, :]
flux_4_17_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == .6, :]
flux_4_18_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == .7, :]
flux_4_19_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == .8, :]
flux_4_20_ccqe_all = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]


# ## Flux 5 Cosine Angular Bins

# In[8]:


flux_5_1 = flux_5.loc[lambda df: df['cos_low'] == -1, :] #Angle -1.0 < x < -0.9
flux_5_2 = flux_5.loc[lambda df: df['cos_low'] == -.9, :] #Angle -0.9 < x < -0.8
flux_5_3 = flux_5.loc[lambda df: df['cos_low'] == -.8, :] #Angle -0.8 < x < -0.7
flux_5_4 = flux_5.loc[lambda df: df['cos_low'] == -.7, :] #Angle -0.7 < x < -0.6
flux_5_5 = flux_5.loc[lambda df: df['cos_low'] == -.6, :] #Angle -0.6 < x < -0.5
flux_5_6 = flux_5.loc[lambda df: df['cos_low'] == -.5, :] #Angle -0.5 < x < -0.4
flux_5_7 = flux_5.loc[lambda df: df['cos_low'] == -.4, :] #Angle -0.4 < x < -0.3
flux_5_8 = flux_5.loc[lambda df: df['cos_low'] == -.3, :] #Angle -0.3 < x < -0.2
flux_5_9 = flux_5.loc[lambda df: df['cos_low'] == -.2, :] #Angle -0.2 < x < -0.1
flux_5_10 = flux_5.loc[lambda df: df['cos_low'] == -.1, :] #Angle -0.1 < x < 0.0
flux_5_11 = flux_5.loc[lambda df: df['cos_low'] == 0, :] #Angle 0.0 < x < 0.1
flux_5_12 = flux_5.loc[lambda df: df['cos_low'] == .1, :] #Angle 0.1 < x < 0.2 
flux_5_13 = flux_5.loc[lambda df: df['cos_low'] == .2, :] #Angle 0.2 < x < 0.3
flux_5_14 = flux_5.loc[lambda df: df['cos_low'] == .3, :] #Angle 0.3 < x < 0.4
flux_5_15 = flux_5.loc[lambda df: df['cos_low'] == .4, :] #Angle 0.4 < x < 0.5
flux_5_16 = flux_5.loc[lambda df: df['cos_low'] == .5, :] #Angle 0.5 < x < 0.6
flux_5_17 = flux_5.loc[lambda df: df['cos_low'] == .6, :] #Angle 0.6 < x < 0.7
flux_5_18 = flux_5.loc[lambda df: df['cos_low'] == .7, :] #Angle 0.7 < x < 0.8
flux_5_19 = flux_5.loc[lambda df: df['cos_low'] == .8, :] #Angle 0.8 < x < 0.9
flux_5_20 = flux_5.loc[lambda df: df['cos_low'] == .9, :] #Angle 0.9 < x < 1.0

flux_5_1_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == -1, :]
flux_5_2_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == -.9, :]
flux_5_3_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == -.8, :]
flux_5_4_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == -.7, :]
flux_5_5_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == -.6, :]
flux_5_6_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == -.5, :]
flux_5_7_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == -.4, :]
flux_5_8_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == -.3, :]
flux_5_9_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == -.2, :]
flux_5_10_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == -.1, :]
flux_5_11_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == 0, :]
flux_5_12_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == .1, :]
flux_5_13_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == .2, :]
flux_5_14_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == .3, :]
flux_5_15_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == .4, :]
flux_5_16_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == .5, :]
flux_5_17_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == .6, :]
flux_5_18_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == .7, :]
flux_5_19_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == .8, :]
flux_5_20_ccqe = flux_5_ccqe.loc[lambda df: df['cos_low'] == .9, :]

flux_5_1_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == -1, :]
flux_5_2_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == -.9, :]
flux_5_3_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == -.8, :]
flux_5_4_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == -.7, :]
flux_5_5_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == -.6, :]
flux_5_6_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == -.5, :]
flux_5_7_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == -.4, :]
flux_5_8_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == -.3, :]
flux_5_9_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == -.2, :]
flux_5_10_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == -.1, :]
flux_5_11_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == 0, :]
flux_5_12_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == .1, :]
flux_5_13_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == .2, :]
flux_5_14_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == .3, :]
flux_5_15_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == .4, :]
flux_5_16_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == .5, :]
flux_5_17_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == .6, :]
flux_5_18_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == .7, :]
flux_5_19_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == .8, :]
flux_5_20_ccqe_like = flux_5_ccqe_like.loc[lambda df: df['cos_low'] == .9, :]

flux_5_1_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == -1, :]
flux_5_2_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == -.9, :]
flux_5_3_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == -.8, :]
flux_5_4_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == -.7, :]
flux_5_5_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == -.6, :]
flux_5_6_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == -.5, :]
flux_5_7_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == -.4, :]
flux_5_8_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == -.3, :]
flux_5_9_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == -.2, :]
flux_5_10_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == -.1, :]
flux_5_11_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == 0, :]
flux_5_12_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == .1, :]
flux_5_13_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == .2, :]
flux_5_14_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == .3, :]
flux_5_15_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == .4, :]
flux_5_16_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == .5, :]
flux_5_17_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == .6, :]
flux_5_18_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == .7, :]
flux_5_19_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == .8, :]
flux_5_20_ccqe_all = flux_5_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]


# ## Flux 6 Cosine Angular Bins

# In[9]:


flux_6_1 = flux_6.loc[lambda df: df['cos_low'] == -1, :] #Angle -1.0 < x < -0.9
flux_6_2 = flux_6.loc[lambda df: df['cos_low'] == -.9, :] #Angle -0.9 < x < -0.8
flux_6_3 = flux_6.loc[lambda df: df['cos_low'] == -.8, :] #Angle -0.8 < x < -0.7
flux_6_4 = flux_6.loc[lambda df: df['cos_low'] == -.7, :] #Angle -0.7 < x < -0.6
flux_6_5 = flux_6.loc[lambda df: df['cos_low'] == -.6, :] #Angle -0.6 < x < -0.5
flux_6_6 = flux_6.loc[lambda df: df['cos_low'] == -.5, :] #Angle -0.5 < x < -0.4
flux_6_7 = flux_6.loc[lambda df: df['cos_low'] == -.4, :] #Angle -0.4 < x < -0.3
flux_6_8 = flux_6.loc[lambda df: df['cos_low'] == -.3, :] #Angle -0.3 < x < -0.2
flux_6_9 = flux_6.loc[lambda df: df['cos_low'] == -.2, :] #Angle -0.2 < x < -0.1
flux_6_10 = flux_6.loc[lambda df: df['cos_low'] == -.1, :] #Angle -0.1 < x < 0.0
flux_6_11 = flux_6.loc[lambda df: df['cos_low'] == 0, :] #Angle 0.0 < x < 0.1
flux_6_12 = flux_6.loc[lambda df: df['cos_low'] == .1, :] #Angle 0.1 < x < 0.2 
flux_6_13 = flux_6.loc[lambda df: df['cos_low'] == .2, :] #Angle 0.2 < x < 0.3
flux_6_14 = flux_6.loc[lambda df: df['cos_low'] == .3, :] #Angle 0.3 < x < 0.4
flux_6_15 = flux_6.loc[lambda df: df['cos_low'] == .4, :] #Angle 0.4 < x < 0.5
flux_6_16 = flux_6.loc[lambda df: df['cos_low'] == .5, :] #Angle 0.5 < x < 0.6
flux_6_17 = flux_6.loc[lambda df: df['cos_low'] == .6, :] #Angle 0.6 < x < 0.7
flux_6_18 = flux_6.loc[lambda df: df['cos_low'] == .7, :] #Angle 0.7 < x < 0.8
flux_6_19 = flux_6.loc[lambda df: df['cos_low'] == .8, :] #Angle 0.8 < x < 0.9
flux_6_20 = flux_6.loc[lambda df: df['cos_low'] == .9, :] #Angle 0.9 < x < 1.0

flux_6_1_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == -1, :]
flux_6_2_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == -.9, :]
flux_6_3_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == -.8, :]
flux_6_4_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == -.7, :]
flux_6_5_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == -.6, :]
flux_6_6_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == -.5, :]
flux_6_7_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == -.4, :]
flux_6_8_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == -.3, :]
flux_6_9_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == -.2, :]
flux_6_10_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == -.1, :]
flux_6_11_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == 0, :]
flux_6_12_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == .1, :]
flux_6_13_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == .2, :]
flux_6_14_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == .3, :]
flux_6_15_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == .4, :]
flux_6_16_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == .5, :]
flux_6_17_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == .6, :]
flux_6_18_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == .7, :]
flux_6_19_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == .8, :]
flux_6_20_ccqe = flux_6_ccqe.loc[lambda df: df['cos_low'] == .9, :]

flux_6_1_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == -1, :]
flux_6_2_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == -.9, :]
flux_6_3_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == -.8, :]
flux_6_4_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == -.7, :]
flux_6_5_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == -.6, :]
flux_6_6_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == -.5, :]
flux_6_7_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == -.4, :]
flux_6_8_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == -.3, :]
flux_6_9_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == -.2, :]
flux_6_10_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == -.1, :]
flux_6_11_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == 0, :]
flux_6_12_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == .1, :]
flux_6_13_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == .2, :]
flux_6_14_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == .3, :]
flux_6_15_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == .4, :]
flux_6_16_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == .5, :]
flux_6_17_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == .6, :]
flux_6_18_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == .7, :]
flux_6_19_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == .8, :]
flux_6_20_ccqe_like = flux_6_ccqe_like.loc[lambda df: df['cos_low'] == .9, :]

flux_6_1_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == -1, :]
flux_6_2_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == -.9, :]
flux_6_3_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == -.8, :]
flux_6_4_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == -.7, :]
flux_6_5_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == -.6, :]
flux_6_6_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == -.5, :]
flux_6_7_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == -.4, :]
flux_6_8_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == -.3, :]
flux_6_9_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == -.2, :]
flux_6_10_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == -.1, :]
flux_6_11_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == 0, :]
flux_6_12_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == .1, :]
flux_6_13_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == .2, :]
flux_6_14_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == .3, :]
flux_6_15_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == .4, :]
flux_6_16_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == .5, :]
flux_6_17_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == .6, :]
flux_6_18_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == .7, :]
flux_6_19_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == .8, :]
flux_6_20_ccqe_all = flux_6_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]


# ## Flux 7 Cosine Angular Bins

# In[10]:


flux_7_1 = flux_7.loc[lambda df: df['cos_low'] == -1, :] #Angle -1.0 < x < -0.9
flux_7_2 = flux_7.loc[lambda df: df['cos_low'] == -.9, :] #Angle -0.9 < x < -0.8
flux_7_3 = flux_7.loc[lambda df: df['cos_low'] == -.8, :] #Angle -0.8 < x < -0.7
flux_7_4 = flux_7.loc[lambda df: df['cos_low'] == -.7, :] #Angle -0.7 < x < -0.6
flux_7_5 = flux_7.loc[lambda df: df['cos_low'] == -.6, :] #Angle -0.6 < x < -0.5
flux_7_6 = flux_7.loc[lambda df: df['cos_low'] == -.5, :] #Angle -0.5 < x < -0.4
flux_7_7 = flux_7.loc[lambda df: df['cos_low'] == -.4, :] #Angle -0.4 < x < -0.3
flux_7_8 = flux_7.loc[lambda df: df['cos_low'] == -.3, :] #Angle -0.3 < x < -0.2
flux_7_9 = flux_7.loc[lambda df: df['cos_low'] == -.2, :] #Angle -0.2 < x < -0.1
flux_7_10 = flux_7.loc[lambda df: df['cos_low'] == -.1, :] #Angle -0.1 < x < 0.0
flux_7_11 = flux_7.loc[lambda df: df['cos_low'] == 0, :] #Angle 0.0 < x < 0.1
flux_7_12 = flux_7.loc[lambda df: df['cos_low'] == .1, :] #Angle 0.1 < x < 0.2 
flux_7_13 = flux_7.loc[lambda df: df['cos_low'] == .2, :] #Angle 0.2 < x < 0.3
flux_7_14 = flux_7.loc[lambda df: df['cos_low'] == .3, :] #Angle 0.3 < x < 0.4
flux_7_15 = flux_7.loc[lambda df: df['cos_low'] == .4, :] #Angle 0.4 < x < 0.5
flux_7_16 = flux_7.loc[lambda df: df['cos_low'] == .5, :] #Angle 0.5 < x < 0.6
flux_7_17 = flux_7.loc[lambda df: df['cos_low'] == .6, :] #Angle 0.6 < x < 0.7
flux_7_18 = flux_7.loc[lambda df: df['cos_low'] == .7, :] #Angle 0.7 < x < 0.8
flux_7_19 = flux_7.loc[lambda df: df['cos_low'] == .8, :] #Angle 0.8 < x < 0.9
flux_7_20 = flux_7.loc[lambda df: df['cos_low'] == .9, :] #Angle 0.9 < x < 1.0

flux_7_1_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == -1, :]
flux_7_2_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == -.9, :]
flux_7_3_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == -.8, :]
flux_7_4_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == -.7, :]
flux_7_5_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == -.6, :]
flux_7_6_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == -.5, :]
flux_7_7_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == -.4, :]
flux_7_8_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == -.3, :]
flux_7_9_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == -.2, :]
flux_7_10_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == -.1, :]
flux_7_11_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == 0, :]
flux_7_12_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == .1, :]
flux_7_13_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == .2, :]
flux_7_14_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == .3, :]
flux_7_15_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == .4, :]
flux_7_16_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == .5, :]
flux_7_17_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == .6, :]
flux_7_18_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == .7, :]
flux_7_19_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == .8, :]
flux_7_20_ccqe = flux_7_ccqe.loc[lambda df: df['cos_low'] == .9, :]

flux_7_1_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == -1, :]
flux_7_2_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == -.9, :]
flux_7_3_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == -.8, :]
flux_7_4_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == -.7, :]
flux_7_5_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == -.6, :]
flux_7_6_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == -.5, :]
flux_7_7_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == -.4, :]
flux_7_8_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == -.3, :]
flux_7_9_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == -.2, :]
flux_7_10_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == -.1, :]
flux_7_11_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == 0, :]
flux_7_12_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == .1, :]
flux_7_13_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == .2, :]
flux_7_14_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == .3, :]
flux_7_15_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == .4, :]
flux_7_16_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == .5, :]
flux_7_17_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == .6, :]
flux_7_18_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == .7, :]
flux_7_19_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == .8, :]
flux_7_20_ccqe_like = flux_7_ccqe_like.loc[lambda df: df['cos_low'] == .9, :]

flux_7_1_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == -1, :]
flux_7_2_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == -.9, :]
flux_7_3_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == -.8, :]
flux_7_4_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == -.7, :]
flux_7_5_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == -.6, :]
flux_7_6_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == -.5, :]
flux_7_7_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == -.4, :]
flux_7_8_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == -.3, :]
flux_7_9_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == -.2, :]
flux_7_10_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == -.1, :]
flux_7_11_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == 0, :]
flux_7_12_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == .1, :]
flux_7_13_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == .2, :]
flux_7_14_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == .3, :]
flux_7_15_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == .4, :]
flux_7_16_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == .5, :]
flux_7_17_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == .6, :]
flux_7_18_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == .7, :]
flux_7_19_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == .8, :]
flux_7_20_ccqe_all = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]


# # **Flux 0 Plots**

# ## Individual Plots for Flux 0 Cosine Bins

# In[12]:


#It'll be best to figure out a way to make this section look more clean and be efficient of space
dc1 = flux_0_1.diff_cross
edges1 = list(flux_0_1.E_low) + [list(flux_0_1.E_up)[-1]]
plt.figure()
plt.stairs(dc1, edges1)
plt.title("Cross Section Distribution at Angle -1 < x < -0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_1_0.9.pdf")
plt.show()

dc2 = flux_0_2.diff_cross
edges2 = list(flux_0_2.E_low) + [list(flux_0_2.E_up)[-1]]
plt.figure()
plt.stairs(dc2, edges2)
plt.title("Cross Section Distribution at Angle -0.9 < x < -0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.9_0.8.pdf")
plt.show()

dc3 = flux_0_3.diff_cross
edges3 = list(flux_0_3.E_low) + [list(flux_0_3.E_up)[-1]]
plt.figure()
plt.stairs(dc3, edges3)
plt.title("Cross Section Distribution at Angle -0.8 < x < -0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.8_0.7.pdf")
plt.show()

dc4 = flux_0_4.diff_cross
edges4 = list(flux_0_4.E_low) + [list(flux_0_4.E_up)[-1]]
plt.figure()
plt.stairs(dc4, edges4)
plt.title("Cross Section Distribution at Angle -0.7 < x < -0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.7_0.6.pdf")
plt.show()

dc5 = flux_0_5.diff_cross
edges5 = list(flux_0_5.E_low) + [list(flux_0_5.E_up)[-1]]
plt.figure()
plt.stairs(dc5, edges5)
plt.title("Cross Section Distribution at Angle -0.6 < x < -0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.6_0.5.pdf")
plt.show()

dc6 = flux_0_6.diff_cross
edges6 = list(flux_0_6.E_low) + [list(flux_0_6.E_up)[-1]]
plt.figure()
plt.stairs(dc6, edges6)
plt.title("Cross Section Distribution at Angle -0.5 < x < -0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.5_0.4.pdf")
plt.show()

dc7 = flux_0_7.diff_cross
edges7 = list(flux_0_7.E_low) + [list(flux_0_7.E_up)[-1]]
plt.figure()
plt.stairs(dc7, edges7)
plt.title("Cross Section Distribution at Angle -0.4 < x < -0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.4_0.3.pdf")
plt.show()

dc8 = flux_0_8.diff_cross
edges8 = list(flux_0_8.E_low) + [list(flux_0_8.E_up)[-1]]
plt.figure()
plt.stairs(dc8, edges8)
plt.title("Cross Section Distribution at Angle -0.3 < x < -0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.3_0.2.pdf")
plt.show()

dc9 = flux_0_9.diff_cross
edges9 = list(flux_0_9.E_low) + [list(flux_0_9.E_up)[-1]]
plt.figure()
plt.stairs(dc9, edges9)
plt.title("Cross Section Distribution at Angle -0.2 < x < -0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.2_0.1.pdf")
plt.show()

dc10 = flux_0_10.diff_cross
edges10 = list(flux_0_10.E_low) + [list(flux_0_10.E_up)[-1]]
plt.figure()
plt.stairs(dc10, edges10)
plt.title("Cross Section Distribution at Angle -0.1 < x < 0.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.1_0.0.pdf")
plt.show()

dc11 = flux_0_11.diff_cross
edges11 = list(flux_0_11.E_low) + [list(flux_0_11.E_up)[-1]]
plt.figure()
plt.stairs(dc11, edges11)
plt.title("Cross Section Distribution at Angle 0.0 < x < 0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.0_0.1.pdf")
plt.show()

dc12 = flux_0_12.diff_cross
edges12 = list(flux_0_12.E_low) + [list(flux_0_12.E_up)[-1]]
plt.figure()
plt.stairs(dc12, edges12)
plt.title("Cross Section Distribution at Angle 0.1 < x < 0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.1_0.2.pdf")
plt.show()

dc13 = flux_0_13.diff_cross
edges13 = list(flux_0_13.E_low) + [list(flux_0_13.E_up)[-1]]
plt.figure()
plt.stairs(dc13, edges13)
plt.title("Cross Section Distribution at Angle 0.2 < x < 0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.2_0.3.pdf")
plt.show()

dc14 = flux_0_14.diff_cross
edges14 = list(flux_0_14.E_low) + [list(flux_0_14.E_up)[-1]]
plt.figure()
plt.stairs(dc14, edges14)
plt.title("Cross Section Distribution at Angle 0.3 < x < 0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.3_0.4.pdf")
plt.show()

dc15 = flux_0_15.diff_cross
edges15 = list(flux_0_15.E_low) + [list(flux_0_15.E_up)[-1]]
plt.figure()
plt.stairs(dc15, edges15)
plt.title("Cross Section Distribution at Angle 0.4 < x < 0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.4_0.5.pdf")
plt.show()

dc16 = flux_0_16.diff_cross
edges16 = list(flux_0_16.E_low) + [list(flux_0_16.E_up)[-1]]
plt.figure()
plt.stairs(dc16, edges16)
plt.title("Cross Section Distribution at Angle 0.5 < x < 0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1400])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.5_0.6.pdf")
plt.show()

dc17 = flux_0_17.diff_cross
edges17 = list(flux_0_17.E_low) + [list(flux_0_17.E_up)[-1]]
plt.figure()
plt.stairs(dc17, edges17)
plt.title("Cross Section Distribution at Angle 0.6 < x < 0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1500])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.6_0.7.pdf")
plt.show()

dc18 = flux_0_18.diff_cross
edges18 = list(flux_0_18.E_low) + [list(flux_0_18.E_up)[-1]]
plt.figure()
plt.stairs(dc18, edges18)
plt.title("Cross Section Distribution at Angle 0.7 < x < 0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.7_0.8.pdf")
plt.show()

dc19 = flux_0_19.diff_cross
edges19 = list(flux_0_19.E_low) + [list(flux_0_19.E_up)[-1]]
plt.figure()
plt.stairs(dc19, edges19)
plt.title("Cross Section Distribution at Angle 0.8 < x < 0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2750])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.8_0.9.pdf")
plt.show()

dc20 = flux_0_20.diff_cross
edges20 = list(flux_0_20.E_low) + [list(flux_0_20.E_up)[-1]]
plt.figure()
plt.stairs(dc20, edges20)
plt.title("Cross Section Distribution at Angle 0.9 < x < 1.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 3000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_0_plots/SBND_CRPA_SuSAv2_Hybrid_Flux0_0.9_1.0.pdf")
plt.show()


# ## Comparison Plots

# In[ ]:


plt.figure()
plt.stairs(dc20, edges20, label="Angle 0.9 < x < 1.0")
plt.stairs(dc1, edges1, label="Angle -1.0 < x < -0.9")
plt.title("Flux Comparison of Most Off-axis to Most On-axis")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 3000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.legend(loc="upper right")
plt.show()


# # **Flux 1**

# ## Individual Plots for Flux 1 Cosine Bins

# In[ ]:


dc1 = flux_1_1.diff_cross
edges1 = list(flux_1_1.E_low) + [list(flux_1_1.E_up)[-1]]
plt.figure()
plt.stairs(dc1, edges1)
plt.title("Cross Section Distribution at Angle -1 < x < -0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_1_0.9.pdf")
plt.show()

dc2 = flux_1_2.diff_cross
edges2 = list(flux_1_2.E_low) + [list(flux_1_2.E_up)[-1]]
plt.figure()
plt.stairs(dc2, edges2)
plt.title("Cross Section Distribution at Angle -0.9 < x < -0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.9_0.8.pdf")
plt.show()

dc3 = flux_1_3.diff_cross
edges3 = list(flux_1_3.E_low) + [list(flux_1_3.E_up)[-1]]
plt.figure()
plt.stairs(dc3, edges3)
plt.title("Cross Section Distribution at Angle -0.8 < x < -0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.8_0.7.pdf")
plt.show()

dc4 = flux_1_4.diff_cross
edges4 = list(flux_1_4.E_low) + [list(flux_1_4.E_up)[-1]]
plt.figure()
plt.stairs(dc4, edges4)
plt.title("Cross Section Distribution at Angle -0.7 < x < -0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.7_0.6.pdf")
plt.show()

dc5 = flux_1_5.diff_cross
edges5 = list(flux_1_5.E_low) + [list(flux_1_5.E_up)[-1]]
plt.figure()
plt.stairs(dc5, edges5)
plt.title("Cross Section Distribution at Angle -0.6 < x < -0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.6_0.5.pdf")
plt.show()

dc6 = flux_1_6.diff_cross
edges6 = list(flux_1_6.E_low) + [list(flux_1_6.E_up)[-1]]
plt.figure()
plt.stairs(dc6, edges6)
plt.title("Cross Section Distribution at Angle -0.5 < x < -0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.5_0.4.pdf")
plt.show()

dc7 = flux_1_7.diff_cross
edges7 = list(flux_1_7.E_low) + [list(flux_1_7.E_up)[-1]]
plt.figure()
plt.stairs(dc7, edges7)
plt.title("Cross Section Distribution at Angle -0.4 < x < -0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.4_0.3.pdf")
plt.show()

dc8 = flux_1_8.diff_cross
edges8 = list(flux_1_8.E_low) + [list(flux_1_8.E_up)[-1]]
plt.figure()
plt.stairs(dc8, edges8)
plt.title("Cross Section Distribution at Angle -0.3 < x < -0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.3_0.2.pdf")
plt.show()

dc9 = flux_1_9.diff_cross
edges9 = list(flux_1_9.E_low) + [list(flux_1_9.E_up)[-1]]
plt.figure()
plt.stairs(dc9, edges9)
plt.title("Cross Section Distribution at Angle -0.2 < x < -0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.2_0.1.pdf")
plt.show()

dc10 = flux_1_10.diff_cross
edges10 = list(flux_1_10.E_low) + [list(flux_1_10.E_up)[-1]]
plt.figure()
plt.stairs(dc10, edges10)
plt.title("Cross Section Distribution at Angle -0.1 < x < 0.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.1_0.0.pdf")
plt.show()

dc11 = flux_1_11.diff_cross
edges11 = list(flux_1_11.E_low) + [list(flux_1_11.E_up)[-1]]
plt.figure()
plt.stairs(dc11, edges11)
plt.title("Cross Section Distribution at Angle 0.0 < x < 0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.0_0.1.pdf")
plt.show()

dc12 = flux_1_12.diff_cross
edges12 = list(flux_1_12.E_low) + [list(flux_1_12.E_up)[-1]]
plt.figure()
plt.stairs(dc12, edges12)
plt.title("Cross Section Distribution at Angle 0.1 < x < 0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.1_0.2.pdf")
plt.show()

dc13 = flux_1_13.diff_cross
edges13 = list(flux_1_13.E_low) + [list(flux_1_13.E_up)[-1]]
plt.figure()
plt.stairs(dc13, edges13)
plt.title("Cross Section Distribution at Angle 0.2 < x < 0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.2_0.3.pdf")
plt.show()

dc14 = flux_1_14.diff_cross
edges14 = list(flux_1_14.E_low) + [list(flux_1_14.E_up)[-1]]
plt.figure()
plt.stairs(dc14, edges14)
plt.title("Cross Section Distribution at Angle 0.3 < x < 0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.3_0.4.pdf")
plt.show()

dc15 = flux_1_15.diff_cross
edges15 = list(flux_1_15.E_low) + [list(flux_1_15.E_up)[-1]]
plt.figure()
plt.stairs(dc15, edges15)
plt.title("Cross Section Distribution at Angle 0.4 < x < 0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.4_0.5.pdf")
plt.show()

dc16 = flux_1_16.diff_cross
edges16 = list(flux_1_16.E_low) + [list(flux_1_16.E_up)[-1]]
plt.figure()
plt.stairs(dc16, edges16)
plt.title("Cross Section Distribution at Angle 0.5 < x < 0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1400])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.5_0.6.pdf")
plt.show()

dc17 = flux_1_17.diff_cross
edges17 = list(flux_1_17.E_low) + [list(flux_1_17.E_up)[-1]]
plt.figure()
plt.stairs(dc17, edges17)
plt.title("Cross Section Distribution at Angle 0.6 < x < 0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1500])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.6_0.7.pdf")
plt.show()

dc18 = flux_1_18.diff_cross
edges18 = list(flux_1_18.E_low) + [list(flux_1_18.E_up)[-1]]
plt.figure()
plt.stairs(dc18, edges18)
plt.title("Cross Section Distribution at Angle 0.7 < x < 0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.7_0.8.pdf")
plt.show()

dc19 = flux_1_19.diff_cross
edges19 = list(flux_1_19.E_low) + [list(flux_1_19.E_up)[-1]]
plt.figure()
plt.stairs(dc19, edges19)
plt.title("Cross Section Distribution at Angle 0.8 < x < 0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2750])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.8_0.9.pdf")
plt.show()

dc20 = flux_1_20.diff_cross
edges20 = list(flux_1_20.E_low) + [list(flux_1_20.E_up)[-1]]
plt.figure()
plt.stairs(dc20, edges20)
plt.title("Cross Section Distribution at Angle 0.9 < x < 1.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 3000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_1_plots/SBND_CRPA_SuSAv2_Hybrid_Flux1_0.9_1.0.pdf")
plt.show()


# # FLux 2

# ## Individual Plots for Flux 2 Cosine Bins

# In[ ]:


dc1 = flux_2_1.diff_cross
edges1 = list(flux_2_1.E_low) + [list(flux_2_1.E_up)[-1]]
plt.figure()
plt.stairs(dc1, edges1)
plt.title("Cross Section Distribution at Angle -1 < x < -0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_1_0.9.pdf")
plt.show()

dc2 = flux_2_2.diff_cross
edges2 = list(flux_2_2.E_low) + [list(flux_2_2.E_up)[-1]]
plt.figure()
plt.stairs(dc2, edges2)
plt.title("Cross Section Distribution at Angle -0.9 < x < -0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.9_0.8.pdf")
plt.show()

dc3 = flux_2_3.diff_cross
edges3 = list(flux_2_3.E_low) + [list(flux_2_3.E_up)[-1]]
plt.figure()
plt.stairs(dc3, edges3)
plt.title("Cross Section Distribution at Angle -0.8 < x < -0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.8_0.7.pdf")
plt.show()

dc4 = flux_2_4.diff_cross
edges4 = list(flux_2_4.E_low) + [list(flux_2_4.E_up)[-1]]
plt.figure()
plt.stairs(dc4, edges4)
plt.title("Cross Section Distribution at Angle -0.7 < x < -0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.7_0.6.pdf")
plt.show()

dc5 = flux_2_5.diff_cross
edges5 = list(flux_2_5.E_low) + [list(flux_2_5.E_up)[-1]]
plt.figure()
plt.stairs(dc5, edges5)
plt.title("Cross Section Distribution at Angle -0.6 < x < -0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.6_0.5.pdf")
plt.show()

dc6 = flux_2_6.diff_cross
edges6 = list(flux_2_6.E_low) + [list(flux_2_6.E_up)[-1]]
plt.figure()
plt.stairs(dc6, edges6)
plt.title("Cross Section Distribution at Angle -0.5 < x < -0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.5_0.4.pdf")
plt.show()

dc7 = flux_2_7.diff_cross
edges7 = list(flux_2_7.E_low) + [list(flux_2_7.E_up)[-1]]
plt.figure()
plt.stairs(dc7, edges7)
plt.title("Cross Section Distribution at Angle -0.4 < x < -0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.4_0.3.pdf")
plt.show()

dc8 = flux_2_8.diff_cross
edges8 = list(flux_2_8.E_low) + [list(flux_2_8.E_up)[-1]]
plt.figure()
plt.stairs(dc8, edges8)
plt.title("Cross Section Distribution at Angle -0.3 < x < -0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.3_0.2.pdf")
plt.show()

dc9 = flux_2_9.diff_cross
edges9 = list(flux_2_9.E_low) + [list(flux_2_9.E_up)[-1]]
plt.figure()
plt.stairs(dc9, edges9)
plt.title("Cross Section Distribution at Angle -0.2 < x < -0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.2_0.1.pdf")
plt.show()

dc10 = flux_2_10.diff_cross
edges10 = list(flux_2_10.E_low) + [list(flux_2_10.E_up)[-1]]
plt.figure()
plt.stairs(dc10, edges10)
plt.title("Cross Section Distribution at Angle -0.1 < x < 0.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.1_0.0.pdf")
plt.show()

dc11 = flux_2_11.diff_cross
edges11 = list(flux_2_11.E_low) + [list(flux_2_11.E_up)[-1]]
plt.figure()
plt.stairs(dc11, edges11)
plt.title("Cross Section Distribution at Angle 0.0 < x < 0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.0_0.1.pdf")
plt.show()

dc12 = flux_2_12.diff_cross
edges12 = list(flux_2_12.E_low) + [list(flux_2_12.E_up)[-1]]
plt.figure()
plt.stairs(dc12, edges12)
plt.title("Cross Section Distribution at Angle 0.1 < x < 0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.1_0.2.pdf")
plt.show()

dc13 = flux_2_13.diff_cross
edges13 = list(flux_2_13.E_low) + [list(flux_2_13.E_up)[-1]]
plt.figure()
plt.stairs(dc13, edges13)
plt.title("Cross Section Distribution at Angle 0.2 < x < 0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.2_0.3.pdf")
plt.show()

dc14 = flux_2_14.diff_cross
edges14 = list(flux_2_14.E_low) + [list(flux_2_14.E_up)[-1]]
plt.figure()
plt.stairs(dc14, edges14)
plt.title("Cross Section Distribution at Angle 0.3 < x < 0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.3_0.4.pdf")
plt.show()

dc15 = flux_2_15.diff_cross
edges15 = list(flux_2_15.E_low) + [list(flux_2_15.E_up)[-1]]
plt.figure()
plt.stairs(dc15, edges15)
plt.title("Cross Section Distribution at Angle 0.4 < x < 0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.4_0.5.pdf")
plt.show()

dc16 = flux_2_16.diff_cross
edges16 = list(flux_2_16.E_low) + [list(flux_2_16.E_up)[-1]]
plt.figure()
plt.stairs(dc16, edges16)
plt.title("Cross Section Distribution at Angle 0.5 < x < 0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1400])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.5_0.6.pdf")
plt.show()

dc17 = flux_2_17.diff_cross
edges17 = list(flux_2_17.E_low) + [list(flux_2_17.E_up)[-1]]
plt.figure()
plt.stairs(dc17, edges17)
plt.title("Cross Section Distribution at Angle 0.6 < x < 0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1500])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.6_0.7.pdf")
plt.show()

dc18 = flux_2_18.diff_cross
edges18 = list(flux_2_18.E_low) + [list(flux_2_18.E_up)[-1]]
plt.figure()
plt.stairs(dc18, edges18)
plt.title("Cross Section Distribution at Angle 0.7 < x < 0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.7_0.8.pdf")
plt.show()

dc19 = flux_2_19.diff_cross
edges19 = list(flux_2_19.E_low) + [list(flux_2_19.E_up)[-1]]
plt.figure()
plt.stairs(dc19, edges19)
plt.title("Cross Section Distribution at Angle 0.8 < x < 0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2750])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.8_0.9.pdf")
plt.show()

dc20 = flux_2_20.diff_cross
edges20 = list(flux_2_20.E_low) + [list(flux_2_20.E_up)[-1]]
plt.figure()
plt.stairs(dc20, edges20)
plt.title("Cross Section Distribution at Angle 0.9 < x < 1.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 3000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_2_plots/SBND_CRPA_SuSAv2_Hybrid_Flux2_0.9_1.0.pdf")
plt.show()


# # Flux 3

# ## Individual Plots for Flux 3 Cosine Bins

# In[ ]:


dc1 = flux_3_1.diff_cross
edges1 = list(flux_3_1.E_low) + [list(flux_3_1.E_up)[-1]]
plt.figure()
plt.stairs(dc1, edges1)
plt.title("Cross Section Distribution at Angle -1 < x < -0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_1_0.9.pdf")
plt.show()

dc2 = flux_3_2.diff_cross
edges2 = list(flux_3_2.E_low) + [list(flux_3_2.E_up)[-1]]
plt.figure()
plt.stairs(dc2, edges2)
plt.title("Cross Section Distribution at Angle -0.9 < x < -0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.9_0.8.pdf")
plt.show()

dc3 = flux_3_3.diff_cross
edges3 = list(flux_3_3.E_low) + [list(flux_3_3.E_up)[-1]]
plt.figure()
plt.stairs(dc3, edges3)
plt.title("Cross Section Distribution at Angle -0.8 < x < -0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.8_0.7.pdf")
plt.show()

dc4 = flux_3_4.diff_cross
edges4 = list(flux_3_4.E_low) + [list(flux_3_4.E_up)[-1]]
plt.figure()
plt.stairs(dc4, edges4)
plt.title("Cross Section Distribution at Angle -0.7 < x < -0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.7_0.6.pdf")
plt.show()

dc5 = flux_3_5.diff_cross
edges5 = list(flux_3_5.E_low) + [list(flux_3_5.E_up)[-1]]
plt.figure()
plt.stairs(dc5, edges5)
plt.title("Cross Section Distribution at Angle -0.6 < x < -0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.6_0.5.pdf")
plt.show()

dc6 = flux_3_6.diff_cross
edges6 = list(flux_3_6.E_low) + [list(flux_3_6.E_up)[-1]]
plt.figure()
plt.stairs(dc6, edges6)
plt.title("Cross Section Distribution at Angle -0.5 < x < -0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.5_0.4.pdf")
plt.show()

dc7 = flux_3_7.diff_cross
edges7 = list(flux_3_7.E_low) + [list(flux_3_7.E_up)[-1]]
plt.figure()
plt.stairs(dc7, edges7)
plt.title("Cross Section Distribution at Angle -0.4 < x < -0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.4_0.3.pdf")
plt.show()

dc8 = flux_3_8.diff_cross
edges8 = list(flux_3_8.E_low) + [list(flux_3_8.E_up)[-1]]
plt.figure()
plt.stairs(dc8, edges8)
plt.title("Cross Section Distribution at Angle -0.3 < x < -0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.3_0.2.pdf")
plt.show()

dc9 = flux_3_9.diff_cross
edges9 = list(flux_3_9.E_low) + [list(flux_3_9.E_up)[-1]]
plt.figure()
plt.stairs(dc9, edges9)
plt.title("Cross Section Distribution at Angle -0.2 < x < -0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.2_0.1.pdf")
plt.show()

dc10 = flux_3_10.diff_cross
edges10 = list(flux_3_10.E_low) + [list(flux_3_10.E_up)[-1]]
plt.figure()
plt.stairs(dc10, edges10)
plt.title("Cross Section Distribution at Angle -0.1 < x < 0.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.1_0.0.pdf")
plt.show()

dc11 = flux_3_11.diff_cross
edges11 = list(flux_3_11.E_low) + [list(flux_3_11.E_up)[-1]]
plt.figure()
plt.stairs(dc11, edges11)
plt.title("Cross Section Distribution at Angle 0.0 < x < 0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.0_0.1.pdf")
plt.show()

dc12 = flux_3_12.diff_cross
edges12 = list(flux_3_12.E_low) + [list(flux_3_12.E_up)[-1]]
plt.figure()
plt.stairs(dc12, edges12)
plt.title("Cross Section Distribution at Angle 0.1 < x < 0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.1_0.2.pdf")
plt.show()

dc13 = flux_3_13.diff_cross
edges13 = list(flux_3_13.E_low) + [list(flux_3_13.E_up)[-1]]
plt.figure()
plt.stairs(dc13, edges13)
plt.title("Cross Section Distribution at Angle 0.2 < x < 0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.2_0.3.pdf")
plt.show()

dc14 = flux_3_14.diff_cross
edges14 = list(flux_3_14.E_low) + [list(flux_3_14.E_up)[-1]]
plt.figure()
plt.stairs(dc14, edges14)
plt.title("Cross Section Distribution at Angle 0.3 < x < 0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.3_0.4.pdf")
plt.show()

dc15 = flux_3_15.diff_cross
edges15 = list(flux_3_15.E_low) + [list(flux_3_15.E_up)[-1]]
plt.figure()
plt.stairs(dc15, edges15)
plt.title("Cross Section Distribution at Angle 0.4 < x < 0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.4_0.5.pdf")
plt.show()

dc16 = flux_3_16.diff_cross
edges16 = list(flux_3_16.E_low) + [list(flux_3_16.E_up)[-1]]
plt.figure()
plt.stairs(dc16, edges16)
plt.title("Cross Section Distribution at Angle 0.5 < x < 0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1400])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.5_0.6.pdf")
plt.show()

dc17 = flux_3_17.diff_cross
edges17 = list(flux_3_17.E_low) + [list(flux_3_17.E_up)[-1]]
plt.figure()
plt.stairs(dc17, edges17)
plt.title("Cross Section Distribution at Angle 0.6 < x < 0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1500])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.6_0.7.pdf")
plt.show()

dc18 = flux_3_18.diff_cross
edges18 = list(flux_3_18.E_low) + [list(flux_3_18.E_up)[-1]]
plt.figure()
plt.stairs(dc18, edges18)
plt.title("Cross Section Distribution at Angle 0.7 < x < 0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.7_0.8.pdf")
plt.show()

dc19 = flux_3_19.diff_cross
edges19 = list(flux_3_19.E_low) + [list(flux_3_19.E_up)[-1]]
plt.figure()
plt.stairs(dc19, edges19)
plt.title("Cross Section Distribution at Angle 0.8 < x < 0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2750])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.8_0.9.pdf")
plt.show()

dc20 = flux_3_20.diff_cross
edges20 = list(flux_3_20.E_low) + [list(flux_3_20.E_up)[-1]]
plt.figure()
plt.stairs(dc20, edges20)
plt.title("Cross Section Distribution at Angle 0.9 < x < 1.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 3000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_3_plots/SBND_CRPA_SuSAv2_Hybrid_Flux3_0.9_1.0.pdf")
plt.show()


# # Flux 4

# ## Individual Plots for Flux 4 Cosine Bins

# In[ ]:


dc1 = flux_4_1.diff_cross
edges1 = list(flux_4_1.E_low) + [list(flux_4_1.E_up)[-1]]
plt.figure()
plt.stairs(dc1, edges1)
plt.title("Cross Section Distribution at Angle -1 < x < -0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_1_0.9.pdf")
plt.show()

dc2 = flux_4_2.diff_cross
edges2 = list(flux_4_2.E_low) + [list(flux_4_2.E_up)[-1]]
plt.figure()
plt.stairs(dc2, edges2)
plt.title("Cross Section Distribution at Angle -0.9 < x < -0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.9_0.8.pdf")
plt.show()

dc3 = flux_4_3.diff_cross
edges3 = list(flux_4_3.E_low) + [list(flux_4_3.E_up)[-1]]
plt.figure()
plt.stairs(dc3, edges3)
plt.title("Cross Section Distribution at Angle -0.8 < x < -0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.8_0.7.pdf")
plt.show()

dc4 = flux_4_4.diff_cross
edges4 = list(flux_4_4.E_low) + [list(flux_4_4.E_up)[-1]]
plt.figure()
plt.stairs(dc4, edges4)
plt.title("Cross Section Distribution at Angle -0.7 < x < -0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.7_0.6.pdf")
plt.show()

dc5 = flux_4_5.diff_cross
edges5 = list(flux_4_5.E_low) + [list(flux_4_5.E_up)[-1]]
plt.figure()
plt.stairs(dc5, edges5)
plt.title("Cross Section Distribution at Angle -0.6 < x < -0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.6_0.5.pdf")
plt.show()

dc6 = flux_4_6.diff_cross
edges6 = list(flux_4_6.E_low) + [list(flux_4_6.E_up)[-1]]
plt.figure()
plt.stairs(dc6, edges6)
plt.title("Cross Section Distribution at Angle -0.5 < x < -0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.5_0.4.pdf")
plt.show()

dc7 = flux_4_7.diff_cross
edges7 = list(flux_4_7.E_low) + [list(flux_4_7.E_up)[-1]]
plt.figure()
plt.stairs(dc7, edges7)
plt.title("Cross Section Distribution at Angle -0.4 < x < -0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.4_0.3.pdf")
plt.show()

dc8 = flux_4_8.diff_cross
edges8 = list(flux_4_8.E_low) + [list(flux_4_8.E_up)[-1]]
plt.figure()
plt.stairs(dc8, edges8)
plt.title("Cross Section Distribution at Angle -0.3 < x < -0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.3_0.2.pdf")
plt.show()

dc9 = flux_4_9.diff_cross
edges9 = list(flux_4_9.E_low) + [list(flux_4_9.E_up)[-1]]
plt.figure()
plt.stairs(dc9, edges9)
plt.title("Cross Section Distribution at Angle -0.2 < x < -0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.2_0.1.pdf")
plt.show()

dc10 = flux_4_10.diff_cross
edges10 = list(flux_4_10.E_low) + [list(flux_4_10.E_up)[-1]]
plt.figure()
plt.stairs(dc10, edges10)
plt.title("Cross Section Distribution at Angle -0.1 < x < 0.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.1_0.0.pdf")
plt.show()

dc11 = flux_4_11.diff_cross
edges11 = list(flux_4_11.E_low) + [list(flux_4_11.E_up)[-1]]
plt.figure()
plt.stairs(dc11, edges11)
plt.title("Cross Section Distribution at Angle 0.0 < x < 0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.0_0.1.pdf")
plt.show()

dc12 = flux_4_12.diff_cross
edges12 = list(flux_4_12.E_low) + [list(flux_4_12.E_up)[-1]]
plt.figure()
plt.stairs(dc12, edges12)
plt.title("Cross Section Distribution at Angle 0.1 < x < 0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.1_0.2.pdf")
plt.show()

dc13 = flux_4_13.diff_cross
edges13 = list(flux_4_13.E_low) + [list(flux_4_13.E_up)[-1]]
plt.figure()
plt.stairs(dc13, edges13)
plt.title("Cross Section Distribution at Angle 0.2 < x < 0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.2_0.3.pdf")
plt.show()

dc14 = flux_4_14.diff_cross
edges14 = list(flux_4_14.E_low) + [list(flux_4_14.E_up)[-1]]
plt.figure()
plt.stairs(dc14, edges14)
plt.title("Cross Section Distribution at Angle 0.3 < x < 0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.3_0.4.pdf")
plt.show()

dc15 = flux_4_15.diff_cross
edges15 = list(flux_4_15.E_low) + [list(flux_4_15.E_up)[-1]]
plt.figure()
plt.stairs(dc15, edges15)
plt.title("Cross Section Distribution at Angle 0.4 < x < 0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.4_0.5.pdf")
plt.show()

dc16 = flux_4_16.diff_cross
edges16 = list(flux_4_16.E_low) + [list(flux_4_16.E_up)[-1]]
plt.figure()
plt.stairs(dc16, edges16)
plt.title("Cross Section Distribution at Angle 0.5 < x < 0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1400])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.5_0.6.pdf")
plt.show()

dc17 = flux_4_17.diff_cross
edges17 = list(flux_4_17.E_low) + [list(flux_4_17.E_up)[-1]]
plt.figure()
plt.stairs(dc17, edges17)
plt.title("Cross Section Distribution at Angle 0.6 < x < 0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1500])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.6_0.7.pdf")
plt.show()

dc18 = flux_4_18.diff_cross
edges18 = list(flux_4_18.E_low) + [list(flux_4_18.E_up)[-1]]
plt.figure()
plt.stairs(dc18, edges18)
plt.title("Cross Section Distribution at Angle 0.7 < x < 0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.7_0.8.pdf")
plt.show()

dc19 = flux_4_19.diff_cross
edges19 = list(flux_4_19.E_low) + [list(flux_4_19.E_up)[-1]]
plt.figure()
plt.stairs(dc19, edges19)
plt.title("Cross Section Distribution at Angle 0.8 < x < 0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2750])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.8_0.9.pdf")
plt.show()

dc20 = flux_4_20.diff_cross
edges20 = list(flux_4_20.E_low) + [list(flux_4_20.E_up)[-1]]
plt.figure()
plt.stairs(dc20, edges20)
plt.title("Cross Section Distribution at Angle 0.9 < x < 1.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 3000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_4_plots/SBND_CRPA_SuSAv2_Hybrid_Flux4_0.9_1.0.pdf")
plt.show()


# # Flux 5

# ## Individual Plots for Flux 5 Cosine Bins

# In[ ]:


dc1 = flux_5_1.diff_cross
edges1 = list(flux_5_1.E_low) + [list(flux_5_1.E_up)[-1]]
plt.figure()
plt.stairs(dc1, edges1)
plt.title("Cross Section Distribution at Angle -1 < x < -0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_1_0.9.pdf")
plt.show()

dc2 = flux_5_2.diff_cross
edges2 = list(flux_5_2.E_low) + [list(flux_5_2.E_up)[-1]]
plt.figure()
plt.stairs(dc2, edges2)
plt.title("Cross Section Distribution at Angle -0.9 < x < -0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.9_0.8.pdf")
plt.show()

dc3 = flux_5_3.diff_cross
edges3 = list(flux_5_3.E_low) + [list(flux_5_3.E_up)[-1]]
plt.figure()
plt.stairs(dc3, edges3)
plt.title("Cross Section Distribution at Angle -0.8 < x < -0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.8_0.7.pdf")
plt.show()

dc4 = flux_5_4.diff_cross
edges4 = list(flux_5_4.E_low) + [list(flux_5_4.E_up)[-1]]
plt.figure()
plt.stairs(dc4, edges4)
plt.title("Cross Section Distribution at Angle -0.7 < x < -0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.7_0.6.pdf")
plt.show()

dc5 = flux_5_5.diff_cross
edges5 = list(flux_5_5.E_low) + [list(flux_5_5.E_up)[-1]]
plt.figure()
plt.stairs(dc5, edges5)
plt.title("Cross Section Distribution at Angle -0.6 < x < -0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.6_0.5.pdf")
plt.show()

dc6 = flux_5_6.diff_cross
edges6 = list(flux_5_6.E_low) + [list(flux_5_6.E_up)[-1]]
plt.figure()
plt.stairs(dc6, edges6)
plt.title("Cross Section Distribution at Angle -0.5 < x < -0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.5_0.4.pdf")
plt.show()

dc7 = flux_5_7.diff_cross
edges7 = list(flux_5_7.E_low) + [list(flux_5_7.E_up)[-1]]
plt.figure()
plt.stairs(dc7, edges7)
plt.title("Cross Section Distribution at Angle -0.4 < x < -0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.4_0.3.pdf")
plt.show()

dc8 = flux_5_8.diff_cross
edges8 = list(flux_5_8.E_low) + [list(flux_5_8.E_up)[-1]]
plt.figure()
plt.stairs(dc8, edges8)
plt.title("Cross Section Distribution at Angle -0.3 < x < -0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.3_0.2.pdf")
plt.show()

dc9 = flux_5_9.diff_cross
edges9 = list(flux_5_9.E_low) + [list(flux_5_9.E_up)[-1]]
plt.figure()
plt.stairs(dc9, edges9)
plt.title("Cross Section Distribution at Angle -0.2 < x < -0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.2_0.1.pdf")
plt.show()

dc10 = flux_5_10.diff_cross
edges10 = list(flux_5_10.E_low) + [list(flux_5_10.E_up)[-1]]
plt.figure()
plt.stairs(dc10, edges10)
plt.title("Cross Section Distribution at Angle -0.1 < x < 0.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.1_0.0.pdf")
plt.show()

dc11 = flux_5_11.diff_cross
edges11 = list(flux_5_11.E_low) + [list(flux_5_11.E_up)[-1]]
plt.figure()
plt.stairs(dc11, edges11)
plt.title("Cross Section Distribution at Angle 0.0 < x < 0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.0_0.1.pdf")
plt.show()

dc12 = flux_5_12.diff_cross
edges12 = list(flux_5_12.E_low) + [list(flux_5_12.E_up)[-1]]
plt.figure()
plt.stairs(dc12, edges12)
plt.title("Cross Section Distribution at Angle 0.1 < x < 0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.1_0.2.pdf")
plt.show()

dc13 = flux_5_13.diff_cross
edges13 = list(flux_5_13.E_low) + [list(flux_5_13.E_up)[-1]]
plt.figure()
plt.stairs(dc13, edges13)
plt.title("Cross Section Distribution at Angle 0.2 < x < 0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.2_0.3.pdf")
plt.show()

dc14 = flux_5_14.diff_cross
edges14 = list(flux_5_14.E_low) + [list(flux_5_14.E_up)[-1]]
plt.figure()
plt.stairs(dc14, edges14)
plt.title("Cross Section Distribution at Angle 0.3 < x < 0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.3_0.4.pdf")
plt.show()

dc15 = flux_5_15.diff_cross
edges15 = list(flux_5_15.E_low) + [list(flux_5_15.E_up)[-1]]
plt.figure()
plt.stairs(dc15, edges15)
plt.title("Cross Section Distribution at Angle 0.4 < x < 0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.4_0.5.pdf")
plt.show()

dc16 = flux_5_16.diff_cross
edges16 = list(flux_5_16.E_low) + [list(flux_5_16.E_up)[-1]]
plt.figure()
plt.stairs(dc16, edges16)
plt.title("Cross Section Distribution at Angle 0.5 < x < 0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1400])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.5_0.6.pdf")
plt.show()

dc17 = flux_5_17.diff_cross
edges17 = list(flux_5_17.E_low) + [list(flux_5_17.E_up)[-1]]
plt.figure()
plt.stairs(dc17, edges17)
plt.title("Cross Section Distribution at Angle 0.6 < x < 0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1500])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.6_0.7.pdf")
plt.show()

dc18 = flux_5_18.diff_cross
edges18 = list(flux_5_18.E_low) + [list(flux_5_18.E_up)[-1]]
plt.figure()
plt.stairs(dc18, edges18)
plt.title("Cross Section Distribution at Angle 0.7 < x < 0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.7_0.8.pdf")
plt.show()

dc19 = flux_5_19.diff_cross
edges19 = list(flux_5_19.E_low) + [list(flux_5_19.E_up)[-1]]
plt.figure()
plt.stairs(dc19, edges19)
plt.title("Cross Section Distribution at Angle 0.8 < x < 0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2750])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.8_0.9.pdf")
plt.show()

dc20 = flux_5_20.diff_cross
edges20 = list(flux_5_20.E_low) + [list(flux_5_20.E_up)[-1]]
plt.figure()
plt.stairs(dc20, edges20)
plt.title("Cross Section Distribution at Angle 0.9 < x < 1.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 3000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_5_plots/SBND_CRPA_SuSAv2_Hybrid_Flux5_0.9_1.0.pdf")
plt.show()


# # Flux 6

# ## Individual Plots for Flux 6 Cosine Bins

# In[ ]:


dc1 = flux_6_1.diff_cross
edges1 = list(flux_6_1.E_low) + [list(flux_6_1.E_up)[-1]]
plt.figure()
plt.stairs(dc1, edges1)
plt.title("Cross Section Distribution at Angle -1 < x < -0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_1_0.9.pdf")
plt.show()

dc2 = flux_6_2.diff_cross
edges2 = list(flux_6_2.E_low) + [list(flux_6_2.E_up)[-1]]
plt.figure()
plt.stairs(dc2, edges2)
plt.title("Cross Section Distribution at Angle -0.9 < x < -0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.9_0.8.pdf")
plt.show()

dc3 = flux_6_3.diff_cross
edges3 = list(flux_6_3.E_low) + [list(flux_6_3.E_up)[-1]]
plt.figure()
plt.stairs(dc3, edges3)
plt.title("Cross Section Distribution at Angle -0.8 < x < -0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.8_0.7.pdf")
plt.show()

dc4 = flux_6_4.diff_cross
edges4 = list(flux_6_4.E_low) + [list(flux_6_4.E_up)[-1]]
plt.figure()
plt.stairs(dc4, edges4)
plt.title("Cross Section Distribution at Angle -0.7 < x < -0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.7_0.6.pdf")
plt.show()

dc5 = flux_6_5.diff_cross
edges5 = list(flux_6_5.E_low) + [list(flux_6_5.E_up)[-1]]
plt.figure()
plt.stairs(dc5, edges5)
plt.title("Cross Section Distribution at Angle -0.6 < x < -0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.6_0.5.pdf")
plt.show()

dc6 = flux_6_6.diff_cross
edges6 = list(flux_6_6.E_low) + [list(flux_6_6.E_up)[-1]]
plt.figure()
plt.stairs(dc6, edges6)
plt.title("Cross Section Distribution at Angle -0.5 < x < -0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.5_0.4.pdf")
plt.show()

dc7 = flux_6_7.diff_cross
edges7 = list(flux_6_7.E_low) + [list(flux_6_7.E_up)[-1]]
plt.figure()
plt.stairs(dc7, edges7)
plt.title("Cross Section Distribution at Angle -0.4 < x < -0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.4_0.3.pdf")
plt.show()

dc8 = flux_6_8.diff_cross
edges8 = list(flux_6_8.E_low) + [list(flux_6_8.E_up)[-1]]
plt.figure()
plt.stairs(dc8, edges8)
plt.title("Cross Section Distribution at Angle -0.3 < x < -0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.3_0.2.pdf")
plt.show()

dc9 = flux_6_9.diff_cross
edges9 = list(flux_6_9.E_low) + [list(flux_6_9.E_up)[-1]]
plt.figure()
plt.stairs(dc9, edges9)
plt.title("Cross Section Distribution at Angle -0.2 < x < -0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.2_0.1.pdf")
plt.show()

dc10 = flux_6_10.diff_cross
edges10 = list(flux_6_10.E_low) + [list(flux_6_10.E_up)[-1]]
plt.figure()
plt.stairs(dc10, edges10)
plt.title("Cross Section Distribution at Angle -0.1 < x < 0.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.1_0.0.pdf")
plt.show()

dc11 = flux_6_11.diff_cross
edges11 = list(flux_6_11.E_low) + [list(flux_6_11.E_up)[-1]]
plt.figure()
plt.stairs(dc11, edges11)
plt.title("Cross Section Distribution at Angle 0.0 < x < 0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.0_0.1.pdf")
plt.show()

dc12 = flux_6_12.diff_cross
edges12 = list(flux_6_12.E_low) + [list(flux_6_12.E_up)[-1]]
plt.figure()
plt.stairs(dc12, edges12)
plt.title("Cross Section Distribution at Angle 0.1 < x < 0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.1_0.2.pdf")
plt.show()

dc13 = flux_6_13.diff_cross
edges13 = list(flux_6_13.E_low) + [list(flux_6_13.E_up)[-1]]
plt.figure()
plt.stairs(dc13, edges13)
plt.title("Cross Section Distribution at Angle 0.2 < x < 0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.2_0.3.pdf")
plt.show()

dc14 = flux_6_14.diff_cross
edges14 = list(flux_6_14.E_low) + [list(flux_6_14.E_up)[-1]]
plt.figure()
plt.stairs(dc14, edges14)
plt.title("Cross Section Distribution at Angle 0.3 < x < 0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.3_0.4.pdf")
plt.show()

dc15 = flux_6_15.diff_cross
edges15 = list(flux_6_15.E_low) + [list(flux_6_15.E_up)[-1]]
plt.figure()
plt.stairs(dc15, edges15)
plt.title("Cross Section Distribution at Angle 0.4 < x < 0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.4_0.5.pdf")
plt.show()

dc16 = flux_6_16.diff_cross
edges16 = list(flux_6_16.E_low) + [list(flux_6_16.E_up)[-1]]
plt.figure()
plt.stairs(dc16, edges16)
plt.title("Cross Section Distribution at Angle 0.5 < x < 0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1400])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.5_0.6.pdf")
plt.show()

dc17 = flux_6_17.diff_cross
edges17 = list(flux_6_17.E_low) + [list(flux_6_17.E_up)[-1]]
plt.figure()
plt.stairs(dc17, edges17)
plt.title("Cross Section Distribution at Angle 0.6 < x < 0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1500])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.6_0.7.pdf")
plt.show()

dc18 = flux_6_18.diff_cross
edges18 = list(flux_6_18.E_low) + [list(flux_6_18.E_up)[-1]]
plt.figure()
plt.stairs(dc18, edges18)
plt.title("Cross Section Distribution at Angle 0.7 < x < 0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.7_0.8.pdf")
plt.show()

dc19 = flux_6_19.diff_cross
edges19 = list(flux_6_19.E_low) + [list(flux_6_19.E_up)[-1]]
plt.figure()
plt.stairs(dc19, edges19)
plt.title("Cross Section Distribution at Angle 0.8 < x < 0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2750])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.8_0.9.pdf")
plt.show()

dc20 = flux_6_20.diff_cross
edges20 = list(flux_6_20.E_low) + [list(flux_6_20.E_up)[-1]]
plt.figure()
plt.stairs(dc20, edges20)
plt.title("Cross Section Distribution at Angle 0.9 < x < 1.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 3000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_6_plots/SBND_CRPA_SuSAv2_Hybrid_Flux6_0.9_1.0.pdf")
plt.show()


# # Flux 7

# ## Individual Plots for Flux 7 Cosine Bins

# In[ ]:


dc1 = flux_7_1.diff_cross
edges1 = list(flux_7_1.E_low) + [list(flux_7_1.E_up)[-1]]
plt.figure()
plt.stairs(dc1, edges1)
plt.title("Cross Section Distribution at Angle -1 < x < -0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_1_0.9.pdf")
plt.show()

dc2 = flux_7_2.diff_cross
edges2 = list(flux_7_2.E_low) + [list(flux_7_2.E_up)[-1]]
plt.figure()
plt.stairs(dc2, edges2)
plt.title("Cross Section Distribution at Angle -0.9 < x < -0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.9_0.8.pdf")
plt.show()

dc3 = flux_7_3.diff_cross
edges3 = list(flux_7_3.E_low) + [list(flux_7_3.E_up)[-1]]
plt.figure()
plt.stairs(dc3, edges3)
plt.title("Cross Section Distribution at Angle -0.8 < x < -0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.8_0.7.pdf")
plt.show()

dc4 = flux_7_4.diff_cross
edges4 = list(flux_7_4.E_low) + [list(flux_7_4.E_up)[-1]]
plt.figure()
plt.stairs(dc4, edges4)
plt.title("Cross Section Distribution at Angle -0.7 < x < -0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.7_0.6.pdf")
plt.show()

dc5 = flux_7_5.diff_cross
edges5 = list(flux_7_5.E_low) + [list(flux_7_5.E_up)[-1]]
plt.figure()
plt.stairs(dc5, edges5)
plt.title("Cross Section Distribution at Angle -0.6 < x < -0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.6_0.5.pdf")
plt.show()

dc6 = flux_7_6.diff_cross
edges6 = list(flux_7_6.E_low) + [list(flux_7_6.E_up)[-1]]
plt.figure()
plt.stairs(dc6, edges6)
plt.title("Cross Section Distribution at Angle -0.5 < x < -0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.5_0.4.pdf")
plt.show()

dc7 = flux_7_7.diff_cross
edges7 = list(flux_7_7.E_low) + [list(flux_7_7.E_up)[-1]]
plt.figure()
plt.stairs(dc7, edges7)
plt.title("Cross Section Distribution at Angle -0.4 < x < -0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.4_0.3.pdf")
plt.show()

dc8 = flux_7_8.diff_cross
edges8 = list(flux_7_8.E_low) + [list(flux_7_8.E_up)[-1]]
plt.figure()
plt.stairs(dc8, edges8)
plt.title("Cross Section Distribution at Angle -0.3 < x < -0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.3_0.2.pdf")
plt.show()

dc9 = flux_7_9.diff_cross
edges9 = list(flux_7_9.E_low) + [list(flux_7_9.E_up)[-1]]
plt.figure()
plt.stairs(dc9, edges9)
plt.title("Cross Section Distribution at Angle -0.2 < x < -0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.2_0.1.pdf")
plt.show()

dc10 = flux_7_10.diff_cross
edges10 = list(flux_7_10.E_low) + [list(flux_7_10.E_up)[-1]]
plt.figure()
plt.stairs(dc10, edges10)
plt.title("Cross Section Distribution at Angle -0.1 < x < 0.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.1_0.0.pdf")
plt.show()

dc11 = flux_7_11.diff_cross
edges11 = list(flux_7_11.E_low) + [list(flux_7_11.E_up)[-1]]
plt.figure()
plt.stairs(dc11, edges11)
plt.title("Cross Section Distribution at Angle 0.0 < x < 0.1")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.0_0.1.pdf")
plt.show()

dc12 = flux_7_12.diff_cross
edges12 = list(flux_7_12.E_low) + [list(flux_7_12.E_up)[-1]]
plt.figure()
plt.stairs(dc12, edges12)
plt.title("Cross Section Distribution at Angle 0.1 < x < 0.2")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.1_0.2.pdf")
plt.show()

dc13 = flux_7_13.diff_cross
edges13 = list(flux_7_13.E_low) + [list(flux_7_13.E_up)[-1]]
plt.figure()
plt.stairs(dc13, edges13)
plt.title("Cross Section Distribution at Angle 0.2 < x < 0.3")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.2_0.3.pdf")
plt.show()

dc14 = flux_7_14.diff_cross
edges14 = list(flux_7_14.E_low) + [list(flux_7_14.E_up)[-1]]
plt.figure()
plt.stairs(dc14, edges14)
plt.title("Cross Section Distribution at Angle 0.3 < x < 0.4")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.3_0.4.pdf")
plt.show()

dc15 = flux_7_15.diff_cross
edges15 = list(flux_7_15.E_low) + [list(flux_7_15.E_up)[-1]]
plt.figure()
plt.stairs(dc15, edges15)
plt.title("Cross Section Distribution at Angle 0.4 < x < 0.5")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.4_0.5.pdf")
plt.show()

dc16 = flux_7_16.diff_cross
edges16 = list(flux_7_16.E_low) + [list(flux_7_16.E_up)[-1]]
plt.figure()
plt.stairs(dc16, edges16)
plt.title("Cross Section Distribution at Angle 0.5 < x < 0.6")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1400])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.5_0.6.pdf")
plt.show()

dc17 = flux_7_17.diff_cross
edges17 = list(flux_7_17.E_low) + [list(flux_7_17.E_up)[-1]]
plt.figure()
plt.stairs(dc17, edges17)
plt.title("Cross Section Distribution at Angle 0.6 < x < 0.7")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 1500])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.6_0.7.pdf")
plt.show()

dc18 = flux_7_18.diff_cross
edges18 = list(flux_7_18.E_low) + [list(flux_7_18.E_up)[-1]]
plt.figure()
plt.stairs(dc18, edges18)
plt.title("Cross Section Distribution at Angle 0.7 < x < 0.8")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.7_0.8.pdf")
plt.show()

dc19 = flux_7_19.diff_cross
edges19 = list(flux_7_19.E_low) + [list(flux_7_19.E_up)[-1]]
plt.figure()
plt.stairs(dc19, edges19)
plt.title("Cross Section Distribution at Angle 0.8 < x < 0.9")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 2750])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.8_0.9.pdf")
plt.show()

dc20 = flux_7_20.diff_cross
edges20 = list(flux_7_20.E_low) + [list(flux_7_20.E_up)[-1]]
plt.figure()
plt.stairs(dc20, edges20)
plt.title("Cross Section Distribution at Angle 0.9 < x < 1.0")
plt.xlabel("Energy (MeV)")
plt.xlim([0, 3000])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.savefig("./Plots/CRPA_SuSAv2_Hybrid/Flux_7_plots/SBND_CRPA_SuSAv2_Hybrid_Flux7_0.9_1.0.pdf")
plt.show()


# # **Read Files and take Lin. Comb.**

# In[25]:


#Courtesy of Alexis N.
#Small adjustments by myself

def readMeas(folder):
    #Reading the flux and getting the step and E_min from the file flx.data
    flxtxt = np.loadtxt(folder+"El_bins_flux_CCQELike.txt")

    #The minimum is:
    E_min = flxtxt[1,0] 
    
    #The step is
    E_step = flxtxt[2,0] - E_min

    #The flux is:
    flx=np.zeros([len(flxtxt)])
    flx[:] = flxtxt[:,2]


    dattxt = np.loadtxt(folder+"2d_bins_flux_CCQELike.txt")

    #The bins and data are
    bins = np.zeros([len(dattxt),4]) #4 numbers cosmin, max , Emin, Emax for every bin

    dat = np.zeros([len(dattxt)])

    #Explicitly to make structure clear:
    for ibin in range(len(dattxt)):
        bins[ibin,0] = dattxt[ibin,0]
        bins[ibin,1] = dattxt[ibin,1]
        bins[ibin,2] = dattxt[ibin,2]
        bins[ibin,3] = dattxt[ibin,3]

        dat[ibin] = dattxt[ibin,4]

    #Now we can save this info in the class Measurement:

    meas = Measurement(E_step,flx,bins,dat,E_min=E_min)

    return meas
    

def plotflux(combination, fnm="fluxcomb.pdf"):

    #We evaluate the flux for a linear combination and make a plot

    flux_vec = combination.eval_flux()
    E_vec = combination.Measurements[0].E_vec
    flx = el_flux_0.diff_cross

    plt.plot(E_vec, flux_vec, '-b', linewidth='.75', label="Flux 0 CCQE-Like")
    plt.plot(E_vec, flx, '-r', linewidth='.75', label='Flux 0')
    plt.legend()
    plt.grid()
    plt.title("Flux 0")
    plt.xlim([0,2500])
    plt.xlabel("Energy (MeV)")
    plt.ylabel("CS ($10^{-43} cm^{2}/MeV$)")
    plt.savefig(fnm)
    plt.show()

    plt.clf()
    
    
def plotfluxcomb(combination, fnm="fluxdatcomb.pdf"):

    #Add plots that show flux 0 and 7

    flux_vec = combination.eval_flux()
    E_vec = combination.Measurements[0].E_vec
    flx = el_flux_0_ccqe_all.diff_cross

    plt.figure()
    plt.plot(E_vec, flux_vec, '-m', linewidth='.75', label='Flux 7 - 0.3*Flux 0 (CCQE-All)')
    plt.plot(E_vec, flx, '-g', linewidth='.75',  label='Flux 0 (CCQE-All)')
    plt.legend()
    plt.grid()
    plt.title("CCQE Linear Combination")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("CS ($10^{-43} cm^{2}/MeV$)")
    plt.xlim([0, 3000])
    plt.savefig(fnm)
    plt.show()

    plt.clf()
    
#Maybe a second plotfluxcomb function for CCQE specifically and have return values to then plot together

   
def plot_dat(meas,fnm='data.pdf'):
    #We will select only the most forward bin (hard coded, as exampl edit as needed)

    El_vec_1 = []
    CS_El_1 = []


    for ibin, binval in enumerate(meas.bins):
        if binval[0] < 0.85 and binval[1] > 0.85: #In his case we are in [0.9, 1] bin
            El_vec_1.append((binval[2]+binval[3])*0.5) #center of bin
            CS_El_1.append(meas.data[ibin])

    
    plt.figure()
    plt.plot(El_vec_1, CS_El_1, label='Combined Flux 0 and Flux 7')
    #plt.plot(El_vec_1, dc1, label='Flux 0')
    plt.legend()
    plt.xlim([0, 3000])
    plt.grid()
    plt.title(r"0.8 < cos $\theta$ < 0.9")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("CS ($10^{-43} cm^{2}/MeV$)")
    plt.savefig(fnm)
    plt.show()
    plt.clf()
    


    

if __name__=='__main__':
    #We look at following folders:
    folders = ['./Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_0/', './Flux/SBND_CRPA_SuSAv2_Hybrid_CCQELike/Flux_7/']

    #For every folder we create a measurement object and add it to a list
    Measurements = []

    for ifldr,folder in enumerate(folders):
        meas = readMeas(folder)
        
        Measurements.append(meas)

    #We now have a list of measurements, which have info on the flux and a dataset
    #We can make an arbitrary linear combination of them:

    #Start with all coefficients set to 1
    coeffs = np.array([1.,1.])
    combi = LinComb(Measurements, coeffs)
    

    #This function makes a plot of the flux for a certain linear combination
    #Lets just do flux 0 separately by setting different coefficients:
    combi.coeffs[1] = 0. #coefficient for flx7 is zero
    plotflux(combi,fnm='./Python/Plots/CRPA_SuSAv2_Hybrid_CCQE/flux_0_ccqe_like.pdf')
    
    combi.coeffs[0] = 0. #coefficient for flux 0 is 0
    combi.coeffs[1] = 1.
#    plotflux(combi, fnm='./Python/Plots/CRPA_SuSAv2_Hybrid_CCQELike/flux_7_ccqe_like.pdf')

    #Now do an actual combination, flux_7 - 0.3*flx_0
    combi.coeffs[0] = -0.3
    combi.coeffs[1] = 1

    print(f"The coefficients used for the linear combinations are {combi.coeffs}")
    #This result should be similar to Rubens document, the high energy tail is reduced
#    plotfluxcomb(combi,fnm='./Python/Plots/CRPA_SuSAv2_Hybrid_CCQEAll/Linear_combinations/lincom_flx0-7_CCQEAll_w_flx7_CCQE_1.pdf')

    

    #Now lets plot the data in he most forward bin (given that we know how the bins are defined in this cas
    #We can evaluate the combination to give a measurement object:
    
    Meascomb = combi.eval_Meas()
#    plot_dat(Meascomb,fnm='./Python/Plots/CRPA_SuSAv2_Hybrid_CCQEAll/Data_combinations/data_combi_0.8_0.9.pdf')
#    plot_dat2(Meascomb,fnm='data_combi2.pdf')


# # **Comparisons and Ratios**

# In[ ]:


#Comparing full flux to CCQE counterpart
flx0 = flux_0.loc[lambda df: df['cos_low'] == .9, :] 
flx0_e = flux_0_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]

flx7 = flux_7.loc[lambda df: df['cos_low'] == .9, :]
flx7_e = flux_7_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]

flx4 = flux_4.loc[lambda df: df['cos_low'] == .9, :]
flx4_e = flux_4_ccqe_all.loc[lambda df: df['cos_low'] == .9, :]

E_mins = flx0.E_low
E_maxs = flx0.E_up
E_vals = (E_mins + E_maxs)*.5

plt.figure()
plt.plot(E_vals, flx0.diff_cross, '-g', linewidth=.75, label='Flux 0 Full')
plt.plot(E_vals, flx0_e.diff_cross, '--g', linewidth=.75, label='Flux 0 CCQE-All')
plt.plot(E_vals, flx4.diff_cross, '-r', linewidth=.75, label='Flux 4 Full')
plt.plot(E_vals, flx4_e.diff_cross, '--r', linewidth=.75, label='Flux 4 CCQE-All')
plt.plot(E_vals, flx7.diff_cross, '-m', linewidth=.75, label='Flux 7 Full')
plt.plot(E_vals, flx7_e.diff_cross, '--m', linewidth=.75, label='Flux 7 CCQE-All')
plt.xlim([0,3000])
plt.xlabel("Energy (MeV)")
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.title(r'0.9 < cos $\theta$ < 1.0')
plt.grid()
plt.legend()
#plt.savefig('./Plots/CRPA_SuSAv2_Hybrid_CCQEAll/Comparisons/flx_0_4_7_at_0.9_1.0.pdf')
plt.show()


#Comparison of full flux and CCQE flux to corresponding linear combination
flx0_7_lin_comb = flx7.diff_cross - 0.3*flx0.diff_cross
flx0_7_lin_comb_e = flx7_e.diff_cross - 0.3*flx0_e.diff_cross
plt.figure()
plt.plot(E_vals, flx0_7_lin_comb, '-m', label='Flux 7 - 0.3*Flux 0')
plt.plot(E_vals, flx0.diff_cross, '-b', linewidth=.75, label='Flux 0')
plt.plot(E_vals, flx0_7_lin_comb_e, '--m', label='Flux 7 CCQE-All - 0.3*Flux 0 CCQE-All')
plt.plot(E_vals, flx0_e.diff_cross, '--b', linewidth=.75, label='Flux 0 CCQE-All')
plt.xlim([0,3000])
plt.title(r'0.9 < cos $\theta$ < 1.0')
plt.xlabel("Energy (MeV)")
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.legend(prop={'size':7})
#plt.savefig('./Plots/CRPA_SuSAv2_Hybrid_CCQEAll/Full_and_comb/flx_0_7_lincomb_at_0.9_1.0.pdf')
plt.show()


#Comparison of ratios
flx7_ratio = flx7_e.diff_cross/flx7.diff_cross
flx0_ratio = flx0_e.diff_cross/flx0.diff_cross
flx_comb_ratio = flx0_7_lin_comb_e/flx0_7_lin_comb
plt.figure()
plt.plot(E_vals, flx7_ratio, '--r', label='CCQE-All ratio for Flux 7', linewidth=.75)
plt.plot(E_vals, flx0_ratio, '--b', label='CCQE-All ratio for Flux 0', linewidth=.75)
plt.plot(E_vals, flx_comb_ratio, '-m',label='CCQE-All ratio for Flux 7 - 0.3*Flux 0', linewidth=.75)
plt.legend()
plt.title(r'0.9 < cos $\theta$ < 1.0')
plt.xlabel("Energy (MeV)")
plt.legend(prop={'size':8})
plt.grid()
#plt.savefig('./Plots/CRPA_SuSAv2_Hybrid_CCQEAll/Ratios/flx_0_7_ratio_0.9_1.0.pdf')
plt.show()


# # **Normalization of Cross-Sections**

# In[ ]:


#Ratios with respect to flux 0
#Unconventional as well. Making a function or class would be recommended
el_lower = el_flux_0.E_low
el_upper = el_flux_0.E_up
el_values = (el_lower + el_upper)*.5

total_el_flux_0 = integrate.trapezoid(el_flux_0_ccqe.diff_cross, el_values) #Trapezoid-rule integral of the differential cross-section over energy values
total_el_flux_1 = integrate.trapezoid(el_flux_1_ccqe.diff_cross, el_values)
total_el_flux_2 = integrate.trapezoid(el_flux_2_ccqe.diff_cross, el_values)
total_el_flux_3 = integrate.trapezoid(el_flux_3_ccqe.diff_cross, el_values)
total_el_flux_4 = integrate.trapezoid(el_flux_4_ccqe.diff_cross, el_values)
total_el_flux_5 = integrate.trapezoid(el_flux_5_ccqe.diff_cross, el_values)
total_el_flux_6 = integrate.trapezoid(el_flux_6_ccqe.diff_cross, el_values)
total_el_flux_7 = integrate.trapezoid(el_flux_7_ccqe.diff_cross, el_values)

norm_el_flux_0 = el_flux_0_ccqe.diff_cross/total_el_flux_0 #Normalized values
norm_el_flux_1 = el_flux_1_ccqe.diff_cross/total_el_flux_1
norm_el_flux_2 = el_flux_2_ccqe.diff_cross/total_el_flux_2
norm_el_flux_3 = el_flux_3_ccqe.diff_cross/total_el_flux_3
norm_el_flux_4 = el_flux_4_ccqe.diff_cross/total_el_flux_4
norm_el_flux_5 = el_flux_5_ccqe.diff_cross/total_el_flux_5
norm_el_flux_6 = el_flux_6_ccqe.diff_cross/total_el_flux_6
norm_el_flux_7 = el_flux_7_ccqe.diff_cross/total_el_flux_7

non_el_flux_0 = el_flux_0_ccqe.diff_cross #To be used for non-normalized ratios
non_el_flux_1 = el_flux_1_ccqe.diff_cross
non_el_flux_2 = el_flux_2_ccqe.diff_cross
non_el_flux_3 = el_flux_3_ccqe.diff_cross
non_el_flux_4 = el_flux_4_ccqe.diff_cross
non_el_flux_5 = el_flux_5_ccqe.diff_cross
non_el_flux_6 = el_flux_6_ccqe.diff_cross
non_el_flux_7 = el_flux_7_ccqe.diff_cross

norm_ratio_1 = norm_el_flux_1/norm_el_flux_0 #Normalized ratios
norm_ratio_2 = norm_el_flux_2/norm_el_flux_0
norm_ratio_3 = norm_el_flux_3/norm_el_flux_0
norm_ratio_4 = norm_el_flux_4/norm_el_flux_0
norm_ratio_5 = norm_el_flux_5/norm_el_flux_0
norm_ratio_6 = norm_el_flux_6/norm_el_flux_0
norm_ratio_7 = norm_el_flux_7/norm_el_flux_0

non_ratio_1 = non_el_flux_1/non_el_flux_0 #Non-normalized ratio
non_ratio_2 = non_el_flux_2/non_el_flux_0
non_ratio_3 = non_el_flux_3/non_el_flux_0
non_ratio_4 = non_el_flux_4/non_el_flux_0
non_ratio_5 = non_el_flux_5/non_el_flux_0
non_ratio_6 = non_el_flux_6/non_el_flux_0
non_ratio_7 = non_el_flux_7/non_el_flux_0


plt.figure()
plt.plot(el_values, norm_el_flux_0, '--m', label='Norm. diff. cross section for Flux 0', linewidth=.75)
plt.plot(el_values, norm_el_flux_7, '--r', label='Norm. diff. cross section for Flux 7', linewidth=.75)
plt.title("Normalized Differential Cross Section (CCQE)")
plt.legend()
plt.xlabel("Energy (MeV)")
plt.grid()
#plt.savefig('./Plots/CRPA_SuSAv2_Hybrid_Normalized/Norm_diff_cross/flx_7.pdf')
plt.show()

plt.figure()
plt.plot(el_values, norm_ratio_7, 'crimson', label='Flux 7 CCQE/Flux 0 CCQE', linewidth=.75)
plt.plot(el_values, norm_ratio_5, 'deeppink', label='Flux 5 CCQE/Flux 0 CCQE', linewidth=.75)
plt.plot(el_values, norm_ratio_3, 'fuchsia', label='Flux 3 CCQE/Flux 0 CCQE', linewidth=.75)
plt.plot(el_values, norm_ratio_1, 'blueviolet', label='Flux 1 CCQE/Flux 0 CCQE', linewidth=.75)
plt.title("Ratio of Normalized Differential Cross Section (CCQE)")
plt.xlim([0,2000])
plt.xlabel("Energy (MeV)")
plt.legend()
plt.grid()
plt.savefig('./Plots/CRPA_SuSAv2_Hybrid_Normalized/Ratios/odd_fluxes_ccqe.pdf')
plt.show()

plt.figure()
plt.plot(el_values, non_ratio_7, 'crimson', label='Flux 7 CCQE/Flux 0 CCQE', linewidth=.75)
plt.plot(el_values, non_ratio_5, 'deeppink', label='Flux 5 CCQE/Flux 0 CCQE', linewidth=.75)
plt.plot(el_values, non_ratio_3, 'fuchsia', label='Flux 3 CCQE/Flux 0 CCQE', linewidth=.75)
plt.plot(el_values, non_ratio_1, 'blueviolet', label='Flux 1 CCQE/Flux 0 CCQE', linewidth=.75)
plt.title("Ratio of Non-Normalized Differential Cross Section (CCQE)")
plt.xlim([0,2000])
plt.xlabel("Energy (MeV)")
plt.legend()
plt.grid()
plt.savefig('./Plots/CRPA_SuSAv2_Hybrid_Normalized/Ratios/odd_fluxes_non_norm_ccqe.pdf')
plt.show()


# # **Larger Angular Bins**

# ## **Data stuff**

# In[ ]:


#Creating a range of angular bins now
#This is HIGHLY UNCONVENTIONAL

#Full
dcs_1 = flux_0_1.diff_cross #Flux 0
dcs_2 = flux_0_2.diff_cross
dcs_3 = flux_0_3.diff_cross
dcs_4 = flux_0_4.diff_cross
dcs_5 = flux_0_5.diff_cross
dcs_6 = flux_0_6.diff_cross
dcs_7 = flux_0_7.diff_cross
dcs_8 = flux_0_8.diff_cross
dcs_9 = flux_0_9.diff_cross
dcs_10 = flux_0_10.diff_cross
dcs_11 = flux_0_11.diff_cross
dcs_12 = flux_0_12.diff_cross
dcs_13 = flux_0_13.diff_cross
dcs_14 = flux_0_14.diff_cross
dcs_15 = flux_0_15.diff_cross
dcs_16 = flux_0_16.diff_cross
dcs_17 = flux_0_17.diff_cross
dcs_18 = flux_0_18.diff_cross
dcs_19 = flux_0_19.diff_cross
dcs_20 = flux_0_20.diff_cross

dcs7_1 = flux_7_1.diff_cross #Flux 7
dcs7_2 = flux_7_2.diff_cross
dcs7_3 = flux_7_3.diff_cross
dcs7_4 = flux_7_4.diff_cross
dcs7_5 = flux_7_5.diff_cross
dcs7_6 = flux_7_6.diff_cross
dcs7_7 = flux_7_7.diff_cross
dcs7_8 = flux_7_8.diff_cross
dcs7_9 = flux_7_9.diff_cross
dcs7_10 = flux_7_10.diff_cross
dcs7_11 = flux_7_11.diff_cross
dcs7_12 = flux_7_12.diff_cross
dcs7_13 = flux_7_13.diff_cross
dcs7_14 = flux_7_14.diff_cross
dcs7_15 = flux_7_15.diff_cross
dcs7_16 = flux_7_16.diff_cross
dcs7_17 = flux_7_17.diff_cross
dcs7_18 = flux_7_18.diff_cross
dcs7_19 = flux_7_19.diff_cross
dcs7_20 = flux_7_20.diff_cross

np_f_1 = dcs_1.to_numpy() #Flux 0 to numpy
np_f_2 = dcs_2.to_numpy()
np_f_3 = dcs_3.to_numpy()
np_f_4 = dcs_4.to_numpy()
np_f_5 = dcs_5.to_numpy()
np_f_6 = dcs_6.to_numpy()
np_f_7 = dcs_7.to_numpy()
np_f_8 = dcs_8.to_numpy()
np_f_9 = dcs_9.to_numpy()
np_f_10 = dcs_10.to_numpy()
np_f_11 = dcs_11.to_numpy()
np_f_12 = dcs_12.to_numpy()
np_f_13 = dcs_13.to_numpy()
np_f_14 = dcs_14.to_numpy()
np_f_15 = dcs_15.to_numpy()
np_f_16 = dcs_16.to_numpy()
np_f_17 = dcs_17.to_numpy()
np_f_18 = dcs_18.to_numpy()
np_f_19 = dcs_19.to_numpy()
np_f_20 = dcs_20.to_numpy()

np7_f_1 = dcs7_1.to_numpy() #Flux 7 to numpy
np7_f_2 = dcs7_2.to_numpy()
np7_f_3 = dcs7_3.to_numpy()
np7_f_4 = dcs7_4.to_numpy()
np7_f_5 = dcs7_5.to_numpy()
np7_f_6 = dcs7_6.to_numpy()
np7_f_7 = dcs7_7.to_numpy()
np7_f_8 = dcs7_8.to_numpy()
np7_f_9 = dcs7_9.to_numpy()
np7_f_10 = dcs7_10.to_numpy()
np7_f_11 = dcs7_11.to_numpy()
np7_f_12 = dcs7_12.to_numpy()
np7_f_13 = dcs7_13.to_numpy()
np7_f_14 = dcs7_14.to_numpy()
np7_f_15 = dcs7_15.to_numpy()
np7_f_16 = dcs7_16.to_numpy()
np7_f_17 = dcs7_17.to_numpy()
np7_f_18 = dcs7_18.to_numpy()
np7_f_19 = dcs7_19.to_numpy()
np7_f_20 = dcs7_20.to_numpy()


#Different contributions
dcs_e_1 = flux_0_1_ccqe_all.diff_cross #Flux 0
dcs_e_2 = flux_0_2_ccqe_all.diff_cross
dcs_e_3 = flux_0_3_ccqe_all.diff_cross
dcs_e_4 = flux_0_4_ccqe_all.diff_cross
dcs_e_5 = flux_0_5_ccqe_all.diff_cross
dcs_e_6 = flux_0_6_ccqe_all.diff_cross
dcs_e_7 = flux_0_7_ccqe_all.diff_cross
dcs_e_8 = flux_0_8_ccqe_all.diff_cross
dcs_e_9 = flux_0_9_ccqe_all.diff_cross
dcs_e_10 = flux_0_10_ccqe_all.diff_cross
dcs_e_11 = flux_0_11_ccqe_all.diff_cross
dcs_e_12 = flux_0_12_ccqe_all.diff_cross
dcs_e_13 = flux_0_13_ccqe_all.diff_cross
dcs_e_14 = flux_0_14_ccqe_all.diff_cross
dcs_e_15 = flux_0_15_ccqe_all.diff_cross
dcs_e_16 = flux_0_16_ccqe_all.diff_cross
dcs_e_17 = flux_0_17_ccqe_all.diff_cross
dcs_e_18 = flux_0_18_ccqe_all.diff_cross
dcs_e_19 = flux_0_19_ccqe_all.diff_cross
dcs_e_20 = flux_0_20_ccqe_all.diff_cross

dcs7_e_1 = flux_7_1_ccqe_all.diff_cross #Flux 7
dcs7_e_2 = flux_7_2_ccqe_all.diff_cross
dcs7_e_3 = flux_7_3_ccqe_all.diff_cross
dcs7_e_4 = flux_7_4_ccqe_all.diff_cross
dcs7_e_5 = flux_7_5_ccqe_all.diff_cross
dcs7_e_6 = flux_7_6_ccqe_all.diff_cross
dcs7_e_7 = flux_7_7_ccqe_all.diff_cross
dcs7_e_8 = flux_7_8_ccqe_all.diff_cross
dcs7_e_9 = flux_7_9_ccqe_all.diff_cross
dcs7_e_10 = flux_7_10_ccqe_all.diff_cross
dcs7_e_11 = flux_7_11_ccqe_all.diff_cross
dcs7_e_12 = flux_7_12_ccqe_all.diff_cross
dcs7_e_13 = flux_7_13_ccqe_all.diff_cross
dcs7_e_14 = flux_7_14_ccqe_all.diff_cross
dcs7_e_15 = flux_7_15_ccqe_all.diff_cross
dcs7_e_16 = flux_7_16_ccqe_all.diff_cross
dcs7_e_17 = flux_7_17_ccqe_all.diff_cross
dcs7_e_18 = flux_7_18_ccqe_all.diff_cross
dcs7_e_19 = flux_7_19_ccqe_all.diff_cross
dcs7_e_20 = flux_7_20_ccqe_all.diff_cross


np_f_e_1 = dcs_e_1.to_numpy() #Flux 0
np_f_e_2 = dcs_e_2.to_numpy()
np_f_e_3 = dcs_e_3.to_numpy()
np_f_e_4 = dcs_e_4.to_numpy()
np_f_e_5 = dcs_e_5.to_numpy()
np_f_e_6 = dcs_e_6.to_numpy()
np_f_e_7 = dcs_e_7.to_numpy()
np_f_e_8 = dcs_e_8.to_numpy()
np_f_e_9 = dcs_e_9.to_numpy()
np_f_e_10 = dcs_e_10.to_numpy()
np_f_e_11 = dcs_e_11.to_numpy()
np_f_e_12 = dcs_e_12.to_numpy()
np_f_e_13 = dcs_e_13.to_numpy()
np_f_e_14 = dcs_e_14.to_numpy()
np_f_e_15 = dcs_e_15.to_numpy()
np_f_e_16 = dcs_e_16.to_numpy()
np_f_e_17 = dcs_e_17.to_numpy()
np_f_e_18 = dcs_e_18.to_numpy()
np_f_e_19 = dcs_e_19.to_numpy()
np_f_e_20 = dcs_e_20.to_numpy()

np7_f_e_1 = dcs7_e_1.to_numpy() #Flux 7
np7_f_e_2 = dcs7_e_2.to_numpy()
np7_f_e_3 = dcs7_e_3.to_numpy()
np7_f_e_4 = dcs7_e_4.to_numpy()
np7_f_e_5 = dcs7_e_5.to_numpy()
np7_f_e_6 = dcs7_e_6.to_numpy()
np7_f_e_7 = dcs7_e_7.to_numpy()
np7_f_e_8 = dcs7_e_8.to_numpy()
np7_f_e_9 = dcs7_e_9.to_numpy()
np7_f_e_10 = dcs7_e_10.to_numpy()
np7_f_e_11 = dcs7_e_11.to_numpy()
np7_f_e_12 = dcs7_e_12.to_numpy()
np7_f_e_13 = dcs7_e_13.to_numpy()
np7_f_e_14 = dcs7_e_14.to_numpy()
np7_f_e_15 = dcs7_e_15.to_numpy()
np7_f_e_16 = dcs7_e_16.to_numpy()
np7_f_e_17 = dcs7_e_17.to_numpy()
np7_f_e_18 = dcs7_e_18.to_numpy()
np7_f_e_19 = dcs7_e_19.to_numpy()
np7_f_e_20 = dcs7_e_20.to_numpy()


#Additional contributions
dcs_a_1 = flux_0_1_ccqe.diff_cross #Flux 0
dcs_a_2 = flux_0_2_ccqe.diff_cross
dcs_a_3 = flux_0_3_ccqe.diff_cross
dcs_a_4 = flux_0_4_ccqe.diff_cross
dcs_a_5 = flux_0_5_ccqe.diff_cross
dcs_a_6 = flux_0_6_ccqe.diff_cross
dcs_a_7 = flux_0_7_ccqe.diff_cross
dcs_a_8 = flux_0_8_ccqe.diff_cross
dcs_a_9 = flux_0_9_ccqe.diff_cross
dcs_a_10 = flux_0_10_ccqe.diff_cross
dcs_a_11 = flux_0_11_ccqe.diff_cross
dcs_a_12 = flux_0_12_ccqe.diff_cross
dcs_a_13 = flux_0_13_ccqe.diff_cross
dcs_a_14 = flux_0_14_ccqe.diff_cross
dcs_a_15 = flux_0_15_ccqe.diff_cross
dcs_a_16 = flux_0_16_ccqe.diff_cross
dcs_a_17 = flux_0_17_ccqe.diff_cross
dcs_a_18 = flux_0_18_ccqe.diff_cross
dcs_a_19 = flux_0_19_ccqe.diff_cross
dcs_a_20 = flux_0_20_ccqe.diff_cross

dcs7_a_1 = flux_7_1_ccqe.diff_cross #Flux 7
dcs7_a_2 = flux_7_2_ccqe.diff_cross
dcs7_a_3 = flux_7_3_ccqe.diff_cross
dcs7_a_4 = flux_7_4_ccqe.diff_cross
dcs7_a_5 = flux_7_5_ccqe.diff_cross
dcs7_a_6 = flux_7_6_ccqe.diff_cross
dcs7_a_7 = flux_7_7_ccqe.diff_cross
dcs7_a_8 = flux_7_8_ccqe.diff_cross
dcs7_a_9 = flux_7_9_ccqe.diff_cross
dcs7_a_10 = flux_7_10_ccqe.diff_cross
dcs7_a_11 = flux_7_11_ccqe.diff_cross
dcs7_a_12 = flux_7_12_ccqe.diff_cross
dcs7_a_13 = flux_7_13_ccqe.diff_cross
dcs7_a_14 = flux_7_14_ccqe.diff_cross
dcs7_a_15 = flux_7_15_ccqe.diff_cross
dcs7_a_16 = flux_7_16_ccqe.diff_cross
dcs7_a_17 = flux_7_17_ccqe.diff_cross
dcs7_a_18 = flux_7_18_ccqe.diff_cross
dcs7_a_19 = flux_7_19_ccqe.diff_cross
dcs7_a_20 = flux_7_20_ccqe.diff_cross


np_f_a_1 = dcs_a_1.to_numpy() #Flux 0
np_f_a_2 = dcs_a_2.to_numpy()
np_f_a_3 = dcs_a_3.to_numpy()
np_f_a_4 = dcs_a_4.to_numpy()
np_f_a_5 = dcs_a_5.to_numpy()
np_f_a_6 = dcs_a_6.to_numpy()
np_f_a_7 = dcs_a_7.to_numpy()
np_f_a_8 = dcs_a_8.to_numpy()
np_f_a_9 = dcs_a_9.to_numpy()
np_f_a_10 = dcs_a_10.to_numpy()
np_f_a_11 = dcs_a_11.to_numpy()
np_f_a_12 = dcs_a_12.to_numpy()
np_f_a_13 = dcs_a_13.to_numpy()
np_f_a_14 = dcs_a_14.to_numpy()
np_f_a_15 = dcs_a_15.to_numpy()
np_f_a_16 = dcs_a_16.to_numpy()
np_f_a_17 = dcs_a_17.to_numpy()
np_f_a_18 = dcs_a_18.to_numpy()
np_f_a_19 = dcs_a_19.to_numpy()
np_f_a_20 = dcs_a_20.to_numpy()

np7_f_a_1 = dcs7_a_1.to_numpy() #Flux 7
np7_f_a_2 = dcs7_a_2.to_numpy()
np7_f_a_3 = dcs7_a_3.to_numpy()
np7_f_a_4 = dcs7_a_4.to_numpy()
np7_f_a_5 = dcs7_a_5.to_numpy()
np7_f_a_6 = dcs7_a_6.to_numpy()
np7_f_a_7 = dcs7_a_7.to_numpy()
np7_f_a_8 = dcs7_a_8.to_numpy()
np7_f_a_9 = dcs7_a_9.to_numpy()
np7_f_a_10 = dcs7_a_10.to_numpy()
np7_f_a_11 = dcs7_a_11.to_numpy()
np7_f_a_12 = dcs7_a_12.to_numpy()
np7_f_a_13 = dcs7_a_13.to_numpy()
np7_f_a_14 = dcs7_a_14.to_numpy()
np7_f_a_15 = dcs7_a_15.to_numpy()
np7_f_a_16 = dcs7_a_16.to_numpy()
np7_f_a_17 = dcs7_a_17.to_numpy()
np7_f_a_18 = dcs7_a_18.to_numpy()
np7_f_a_19 = dcs7_a_19.to_numpy()
np7_f_a_20 = dcs7_a_20.to_numpy()




range1 = (np_f_1+np_f_2+np_f_3+np_f_4+np_f_5+np_f_6+np_f_7+np_f_8+np_f_9+np_f_10)/.1 #cos in [-1, 0]
range2 = (np_f_11+np_f_12+np_f_13+np_f_14+np_f_15+np_f_16)/.1*.6 #cos in [0, 0.6]
range3 = (np_f_17+np_f_18+np_f_19+np_f_20)/.1*.4 #cos in [0.6, 1]

range7_1 = (np7_f_1+np7_f_2+np7_f_3+np7_f_4+np7_f_5+np7_f_6+np7_f_7+np7_f_8+np7_f_9+np7_f_10)/.1 #cos in [-1, 0]
range7_2 = (np7_f_11+np7_f_12+np7_f_13+np7_f_14+np7_f_15+np7_f_16)/.1*.6 #cos in [0, 0.6]
range7_3 = (np7_f_17+np7_f_18+np7_f_19+np7_f_20)/.1*.4 #cos in [0.6, 1]


range1_e = (np_f_e_1+np_f_e_2+np_f_e_3+np_f_e_4+np_f_e_5+np_f_e_6+np_f_e_7+np_f_e_8+np_f_e_9+np_f_e_10)/.1
range2_e = (np_f_e_11+np_f_e_12+np_f_e_13+np_f_e_14+np_f_e_15+np_f_e_16)/.1*.6
range3_e = (np_f_e_17+np_f_e_18+np_f_e_19+np_f_e_20)/.1*.4

range7_1_e = (np7_f_e_1+np7_f_e_2+np7_f_e_3+np7_f_e_4+np7_f_e_5+np7_f_e_6+np7_f_e_7+np7_f_e_8+np7_f_e_9+np7_f_e_10)/.1
range7_2_e = (np7_f_e_11+np7_f_e_12+np7_f_e_13+np7_f_e_14+np7_f_e_15+np7_f_e_16)/.1*.6
range7_3_e = (np7_f_e_17+np7_f_e_18+np7_f_e_19+np7_f_e_20)/.1*.4


range1_a = (np_f_a_1+np_f_a_2+np_f_a_3+np_f_a_4+np_f_a_5+np_f_a_6+np_f_a_7+np_f_a_8+np_f_a_9+np_f_a_10)/.1
range2_a = (np_f_a_11+np_f_a_12+np_f_a_13+np_f_a_14+np_f_a_15+np_f_a_16)/.1*.6
range3_a = (np_f_a_17+np_f_a_18+np_f_a_19+np_f_a_20)/.1*.4

range7_1_a = (np7_f_a_1+np7_f_a_2+np7_f_a_3+np7_f_a_4+np7_f_a_5+np7_f_a_6+np7_f_a_7+np7_f_a_8+np7_f_a_9+np7_f_a_10)/.1
range7_2_a = (np7_f_a_11+np7_f_a_12+np7_f_a_13+np7_f_a_14+np7_f_a_15+np7_f_a_16)/.1*.6
range7_3_a = (np7_f_a_17+np7_f_a_18+np7_f_a_19+np7_f_a_20)/.1*.4


#Linear combinations
lc_1 = range7_1 - (0.3*range1) #[-1, 0]
lc_2 = range7_2 - (0.3*range2) #[0, 0.6]
lc_3 = range7_3 - (0.3*range3) #[0.6, 1]

lc_e_1 = range7_1_e - (0.3*range1_e) #[-1, 0]
lc_e_2 = range7_2_e - (0.3*range2_e) #[0, 0.6]
lc_e_3 = range7_3_e - (0.3*range3_e) #[0, 0.6]

lc_a_1 = range7_1_a - (0.3*range1_a) #[-1, 0]
lc_a_2 = range7_2_a - (0.3*range2_a) #[0, 0.6]
lc_a_3 = range7_3_a - (0.3*range3_a) #[0, 0.6]


# ## **Plots**

# In[ ]:


plt.figure()
plt.plot(E_vals, range1, '-b', label='Flux 0', linewidth=.8)
plt.plot(E_vals, range1_e, '-r', label='Flux 0 Total CCQE', linewidth=.8)
plt.plot(E_vals, range1_a, '-m', label='Flux 0 Real CCQE', linewidth=.8)
plt.title(r'-1.0 < cos $\theta$ < 0.0')
plt.xlabel("Energy (MeV)")
plt.xlim([0,750])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.legend()
#plt.savefig('./Plots/Angular_bins/Flux_0_ccqelike_w_ccqe_n1.0_0.0.pdf')
plt.show()

plt.figure()
plt.plot(E_vals, range2, '-b', label='Flux 0', linewidth=.8)
plt.plot(E_vals, range2_e, '-r', label='Flux 0 Total CCQE', linewidth=.8)
plt.plot(E_vals, range2_a, '-m', label='Flux 0 Real CCQE', linewidth=.8)
plt.title(r'0.0 < cos $\theta$ < 0.6')
plt.xlabel("Energy (MeV)")
plt.xlim([0,1250])
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.legend()
#plt.savefig('./Plots/Angular_bins/Flux_0_ccqelike_w_ccqe_0.0_0.6.pdf')
plt.show()

plt.figure()
plt.plot(E_vals, range3, '-b', label='Flux 0', linewidth=.8)
plt.plot(E_vals, range3_e, '-r', label='Flux 0 Total CCQE', linewidth=.8)
plt.plot(E_vals, range3_a, '-m', label='Flux 0 Real CCQE', linewidth=.8)
plt.title(r'0.6 < cos $\theta$ < 1.0')
plt.xlabel("Energy (MeV)")
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.legend()
#plt.savefig('./Plots/Angular_bins/Flux_0_ccqelike_w_ccqe_0.6_1.0.pdf')
plt.show()


plt.figure()
plt.plot(E_vals, lc_1, '-b', label='LC Flux 0', linewidth=.8)
plt.plot(E_vals, lc_e_1, '-r', label='LC Flux 0 Total CCQE', linewidth=.8)
plt.plot(E_vals, lc_a_1, '-m', label='LC Flux 0 Real CCQE', linewidth=.8)
plt.title(r'-1.0 < cos $\theta$ < 0.0')
plt.xlim([0,750])
plt.xlabel("Energy (MeV)")
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.legend()
#plt.savefig('./Plots/Angular_bins/LC_Flux_0_ccqelike_w_ccqe_n1.0_0.0.pdf')
plt.show()

plt.figure()
plt.plot(E_vals, lc_2, '-b', label='LC Flux 0', linewidth=.8)
plt.plot(E_vals, lc_e_2, '-r', label='LC Flux 0 Total CCQE', linewidth=.8)
plt.plot(E_vals, lc_a_2, '-m', label='LC Flux 0 Real CCQE', linewidth=.8)
plt.title(r'0.0 < cos $\theta$ < 0.6')
plt.xlim([0,1250])
plt.xlabel("Energy (MeV)")
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.legend()
#plt.savefig('./Plots/Angular_bins/LC_Flux_0_ccqelike_w_ccqe_0.0_0.6.pdf')
plt.show()

plt.figure()
plt.plot(E_vals, lc_3, '-b', label='LC Flux 0', linewidth=.8)
plt.plot(E_vals, lc_e_3, '-r', label='LC Flux 0 Total CCQE', linewidth=.8)
plt.plot(E_vals, lc_a_3, '-m', label='LC Flux 0 Real CCQE', linewidth=.8)
plt.title(r'0.6 < cos $\theta$ < 1.0')
plt.xlabel("Energy (MeV)")
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.legend()
#plt.savefig('./Plots/Angular_bins/LC_Flux_0_ccqelike_w_ccqe_0.6_1.0.pdf')
plt.show()

#for i in range(20):
#    bound1 = zeros(20,1)
#    bound2 = zeros(20,1)
#    bound3 = zeros(20,1)
#    if i < 10:
#        bound1 = bound1 + np_f_(i+1)
#    elif i < 16:
#        b
#    else:
#        c


# In[ ]:


fig, axs = plt.subplots(2, 3, figsize=(15,7))

axs[0, 0].plot(E_vals, range1, '-b', label='Flux 0', linewidth=.8)
axs[0, 0].plot(E_vals, range1_e, '-r', label='Flux 0 CCQE-Like', linewidth=.8)
axs[0, 0].plot(E_vals, range1_a, '-m', label='Flux 0 Real CCQE', linewidth=.8)
axs[0, 0].set_title(r'No Linear Comb. -1.0 < cos $\theta$ < 0.0')
axs[0, 0].set_xlim([0,750])
axs[0, 0].set_ylabel('Cross Section ($10^{-43} cm^{2}/MeV$)')
axs[0, 0].grid()
axs[0, 0].legend(prop={'size':8})

axs[0, 1].plot(E_vals, range2, '-b', label='Flux 0', linewidth=.8)
axs[0, 1].plot(E_vals, range2_e, '-r', label='Flux 0 CCQE-Like', linewidth=.8)
axs[0, 1].plot(E_vals, range2_a, '-m', label='Flux 0 Real CCQE', linewidth=.8)
axs[0, 1].set_title(r'No Linear Comb. 0.0 < cos $\theta$ < 0.6')
axs[0, 1].set_xlim([0,1250])
axs[0, 1].grid()
axs[0, 1].legend(prop={'size':8})

axs[0, 2].plot(E_vals, range3, '-b', label='Flux 0', linewidth=.8)
axs[0, 2].plot(E_vals, range3_e, '-r', label='Flux 0 CCQE-Like', linewidth=.8)
axs[0, 2].plot(E_vals, range3_a, '-m', label='Flux 0 Real CCQE', linewidth=.8)
axs[0, 2].set_title(r'No Linear Comb. 0.6 < cos $\theta$ < 1.0')
axs[0, 2].set_xlim([0,3000])
axs[0, 2].grid()
axs[0, 2].legend(prop={'size':8})

axs[1, 0].plot(E_vals, lc_1, '-b', label='Flux 0', linewidth=.8)
axs[1, 0].plot(E_vals, lc_e_1, '-r', label='Flux 0 CCQE-Like', linewidth=.8)
axs[1, 0].plot(E_vals, lc_a_1, '-m', label='Flux 0 Real CCQE', linewidth=.8)
axs[1, 0].set_title(r'Linear Comb. -1.0 < cos $\theta$ < 0.0')
axs[1, 0].set_xlim([0,750])
axs[1, 0].set_xlabel('Energy (MeV)')
axs[1, 0].set_ylabel('Cross Section ($10^{-43} cm^{2}/MeV$)')
axs[1, 0].grid()
axs[1, 0].legend(prop={'size':8})

axs[1, 1].plot(E_vals, lc_2, '-b', label='Flux 0', linewidth=.8)
axs[1, 1].plot(E_vals, lc_e_2, '-r', label='Flux 0 CCQE-Like', linewidth=.8)
axs[1, 1].plot(E_vals, lc_a_2, '-m', label='Flux 0 Real CCQE', linewidth=.8)
axs[1, 1].set_title(r'Linear Comb. 0.0 < cos $\theta$ < 0.6')
axs[1, 1].set_xlim([0,1250])
axs[1, 1].set_xlabel('Energy (MeV)')
axs[1, 1].grid()
axs[1, 1].legend(prop={'size':8})

axs[1, 2].plot(E_vals, lc_3, '-b', label='Flux 0', linewidth=.8)
axs[1, 2].plot(E_vals, lc_e_3, '-r', label='Flux 0 CCQE-Like', linewidth=.8)
axs[1, 2].plot(E_vals, lc_a_3, '-m', label='Flux 0 Real CCQE', linewidth=.8)
axs[1, 2].set_title(r'Linear Comb. 0.6 < cos $\theta$ < 1.0')
axs[1, 2].set_xlim([0,3000])
axs[1, 2].set_xlabel('Energy (MeV)')
axs[1, 2].grid()
axs[1, 2].legend(prop={'size':8})

#fig.savefig('./Plots/Angular_bins/big_graph.pdf')


# # **Extra**

# In[ ]:


e_low = el_flux_0.E_low
e_upper = el_flux_0.E_up
e_values = (e_low + e_upper)*.5

plt.figure()
plt.plot(e_values, el_flux_0_ccqe_all.diff_cross, 'red', label='Flux 0 CCQE-All')
plt.plot(e_values, el_flux_1_ccqe_all.diff_cross, 'orangered', label='Flux 1 CCQE-All')
plt.plot(e_values, el_flux_2_ccqe_all.diff_cross, 'coral', label='Flux 2 CCQE-All')
plt.plot(e_values, el_flux_3_ccqe_all.diff_cross, 'orange', label='Flux 3 CCQE-All')
plt.plot(e_values, el_flux_4_ccqe_all.diff_cross, 'goldenrod', label='Flux 4 CCQE-All')
plt.plot(e_values, el_flux_5_ccqe_all.diff_cross, 'gold', label='Flux 5 CCQE-All')
plt.plot(e_values, el_flux_6_ccqe_all.diff_cross, 'yellow', label='Flux 6 CCQE-All')
plt.plot(e_values, el_flux_7_ccqe_all.diff_cross, 'greenyellow', label='Flux 7 CCQE-All')
plt.xlim([0,5000])
plt.title('Energy Distribution of CCQE')
plt.xlabel("Energy (MeV)")
plt.ylabel("Cross Section ($10^{-43} cm^{2}/MeV$)")
plt.grid()
plt.legend()
plt.savefig('./Plots/Fluxes_ccqe_all.pdf')
plt.show()


# In[ ]:




