import csep
from csep.core import poisson_evaluations as poisson
from csep.utils import datasets, time_utils, plots, stats
from csep.core import regions, catalog_evaluations
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import datetime
import pandas as pd
import seaborn as sns

### Set yup experiment parameters
start_date = time_utils.strptime_to_utc_datetime('2006-01-01 00:00:00.0')
end_date = time_utils.strptime_to_utc_datetime('2011-01-01 00:00:00.0')

min_mw = 4.95
max_mw = 8.95
dmw = 0.1

# Create space and magnitude regions. The forecast is already filtered in space and magnitude
magnitudes = regions.magnitude_bins(min_mw, max_mw, dmw)
region = regions.california_relm_region()

# Bind region information to the forecast (this will be used for binning of the catalogs)
space_magnitude_region = regions.create_space_magnitude_region(region, magnitudes)

forecast = csep.load_catalog_forecast("Forecasts/SRMS_Km_2511.dat",
                                      start_time = start_date, end_time = end_date, filter_spatial=True,
                                      region = space_magnitude_region)

## Download catalogue of observed events from comcat                                      
comcat_catalog = csep.query_comcat(start_date, end_date, min_magnitude=forecast.min_magnitude)

# Filter observed catalogue using the same region as the forecast
comcat_catalog = comcat_catalog.filter_spatial(forecast.region)

SRMS_full = csep.load_catalog_forecast("/Forecasts/Forecasts_1985_2005/SRMS_Km_2303.dat",
                                      start_time = start_date, end_time = end_date, filter_spatial=True,
                                      region = space_magnitude_region, apply_filters=True) 

SRMS_DC = csep.load_catalog_forecast("/Forecasts/Forecasts_1985_2005/SRMSMs_Km_2303.dat",
                                      start_time = start_date, end_time = end_date, filter_spatial=True,
                                      region = space_magnitude_region, apply_filters=True) 
                                      
SRMSNK_full = csep.load_catalog_forecast("/Forecasts/Forecasts_1985_2005/SRMSNK_Km_2303.dat" ,
                                      start_time = start_date, end_time = end_date, filter_spatial=True,
                                      region = space_magnitude_region, apply_filters=True)  
                                      
SRMSNK_DC = csep.load_catalog_forecast("/Forecasts/Forecasts_1985_2005/SRMSNKMs_Km_2303.dat" ,
                                      start_time = start_date, end_time = end_date, filter_spatial=True,
                                      region = space_magnitude_region, apply_filters=True)
                                      
FDSRMS_full = csep.load_catalog_forecast("/Forecasts/Forecasts_1985_2005/FDSRMS_Km_2303.txt" ,
                                      start_time = start_date, end_time = end_date, filter_spatial=True,
                                      region = space_magnitude_region, apply_filters=True)
                                      
FDSRMS_DC = csep.load_catalog_forecast("/Forecasts/Forecasts_1985_2005/FDSRMSMs_Km_2303.txt" ,
                                      start_time = start_date, end_time = end_date, filter_spatial=True,
                                      region = space_magnitude_region, apply_filters=True)              
                                      
### A function to apply 4 standard tests to a list of catalogue-type forecasts
## Performs pseudo-likelihood test, magnitude test, number test and spatial test 
def alphabet_tests_catalog(forecast_list, catalog):
    LTests = []
    MTests = []
    NTests = []
    STests = []
    
    for i in range(len(forecast_list)):
        print("Running L-tests")
        Lresult = catalog_evaluations.pseudolikelihood_test(forecast_list[i], catalog, verbose=False)
        LTests.append(Lresult)
        
        print("Running S-tests")
        Sresult = catalog_evaluations.spatial_test(forecast_list[i], catalog, verbose=False)
        STests.append(Sresult)
    
        print("Running M-tests")
        Mresult = catalog_evaluations.magnitude_test(forecast_list[i], catalog, verbose=False)
        MTests.append(Mresult)
    
        print("Running N-tests")
        Nresult = catalog_evaluations.number_test(forecast_list[i], catalog, verbose=False)
        NTests.append(Nresult)
    
    return LTests, MTests, NTests, STests  

## Make a list of forecasts
forecast_cats = [SRMS_full, SRMS_DC, FDSRMS_full, FDSRMS_DC, SRMSNK_full, SRMSNK_DC]

## Then run all the tests
LTestCat, MTestCat, NTestCat, STestCat = alphabet_tests_catalog(forecast_cats, comcat_catalog)

start_date2 = time_utils.strptime_to_utc_datetime('2011-01-01 00:00:00.0')
end_date2 = time_utils.strptime_to_utc_datetime('2016-01-01 00:00:00.0')

comcat_catalog2 = csep.query_comcat(start_date2, end_date2, min_magnitude=forecast.min_magnitude)

# Filter observed catalog using the same region as the forecast
comcat_catalog2 = comcat_catalog2.filter_spatial(forecast.region)


start_date3 = time_utils.strptime_to_utc_datetime('2016-01-01 00:00:00.0')
end_date3 = time_utils.strptime_to_utc_datetime('2021-01-01 00:00:00.0')

comcat_catalog3 = csep.query_comcat(start_date3, end_date3, min_magnitude=forecast.min_magnitude)

# Filter observed catalog using the same region as the forecast
comcat_catalog3 = comcat_catalog3.filter_spatial(forecast.region)

LTestCat2, MTestCat2, NTestCat2, STestCat2 = alphabet_tests_catalog(forecast_cats, comcat_catalog2)
LTestCat3, MTestCat3, NTestCat3, STestCat3 = alphabet_tests_catalog(forecast_cats, comcat_catalog3)

## Make dataframes for each forecast test distribution
## There is probably a much tidier way to do this but I'm not a python pro yet...

### N-test
df_Num = pd.DataFrame(dict(srms = NTestCat[0].test_distribution, srmsdc = NTestCat[1].test_distribution, srmsnk = NTestCat[4].test_distribution, srmsnkdc = NTestCat[5].test_distribution, fdsrms = NTestCat[2].test_distribution, fdsrmsdc = NTestCat[3].test_distribution))

df_Num2 = pd.DataFrame(dict(srms = NTestCat2[0].test_distribution, srmsdc = NTestCat2[1].test_distribution, srmsnk = NTestCat2[4].test_distribution, srmsnkdc = NTestCat2[5].test_distribution, fdsrms = NTestCat2[2].test_distribution, fdsrmsdc = NTestCat2[3].test_distribution))

df_Num3 = pd.DataFrame(dict(srms = NTestCat3[0].test_distribution, srmsdc = NTestCat3[1].test_distribution, srmsnk = NTestCat3[4].test_distribution, srmsnkdc = NTestCat3[5].test_distribution, fdsrms = NTestCat3[2].test_distribution, fdsrmsdc = NTestCat3[3].test_distribution))

### S-test
df_spat = pd.DataFrame(dict(srms = STestCat[0].test_distribution, srmsdc = STestCat[1].test_distribution, srmsnk = STestCat[4].test_distribution, srmsnkdc = STestCat[5].test_distribution, fdsrms = STestCat[2].test_distribution, fdsrmsdc = STestCat[3].test_distribution))

df_spat2 = pd.DataFrame(dict(srms = STestCat2[0].test_distribution, srmsdc = STestCat2[1].test_distribution, srmsnk = STestCat2[4].test_distribution, srmsnkdc = STestCat2[5].test_distribution, fdsrms = STestCat2[2].test_distribution, fdsrmsdc = STestCat2[3].test_distribution))

df_spat3 = pd.DataFrame(dict(srms = STestCat3[0].test_distribution, srmsdc = STestCat3[1].test_distribution, srmsnk = STestCat3[4].test_distribution, srmsnkdc = STestCat3[5].test_distribution, fdsrms = STestCat3[2].test_distribution, fdsrmsdc = STestCat3[3].test_distribution))

## Dataframes of observations
spat_obs = [STestCat[0].observed_statistic, STestCat[1].observed_statistic, STestCat[4].observed_statistic, STestCat[5].observed_statistic, STestCat[2].observed_statistic, STestCat[3].observed_statistic]
spat_model = ['srms', 'srmsdc', 'srmsnk', 'srmsnkdc', 'fdsrms', 'fdsrmsdc']
spat_obs_df = pd.DataFrame(list(zip(spat_obs, spat_model)),
                            columns =['obs', 'model'])
spat_obs2 = [STestCat2[0].observed_statistic, STestCat2[1].observed_statistic, STestCat2[4].observed_statistic, STestCat2[5].observed_statistic, STestCat2[2].observed_statistic, STestCat2[3].observed_statistic]
spat_obs_df2 = pd.DataFrame(list(zip(spat_obs2, spat_model)),
                            columns =['obs', 'model'])
                            
spat_obs3 = [STestCat3[0].observed_statistic, STestCat3[1].observed_statistic, STestCat3[4].observed_statistic, STestCat3[5].observed_statistic, STestCat3[2].observed_statistic, STestCat3[3].observed_statistic]
spat_obs_df3 = pd.DataFrame(list(zip(spat_obs3, spat_model)),
                            columns =['obs', 'model'])

### Pseudo-likelihood test distributions
df_pl = pd.DataFrame(dict(srms = LTestCat[0].test_distribution, srmsdc = LTestCat[1].test_distribution, srmsnk = LTestCat[4].test_distribution, srmsnkdc = LTestCat[5].test_distribution, fdsrms = LTestCat[2].test_distribution, fdsrmsdc = LTestCat[3].test_distribution))

df_pl2 = pd.DataFrame(dict(srms = LTestCat2[0].test_distribution, srmsdc = LTestCat2[1].test_distribution, srmsnk = LTestCat2[4].test_distribution, srmsnkdc = LTestCat2[5].test_distribution, fdsrms = LTestCat2[2].test_distribution, fdsrmsdc = LTestCat2[3].test_distribution))

df_pl3 = pd.DataFrame(dict(srms = LTestCat3[0].test_distribution, srmsdc = LTestCat3[1].test_distribution, srmsnk = LTestCat3[4].test_distribution, srmsnkdc = LTestCat3[5].test_distribution, fdsrms = LTestCat3[2].test_distribution, fdsrmsdc = LTestCat3[3].test_distribution))

## And observations
pl_obs = [LTestCat[0].observed_statistic, LTestCat[1].observed_statistic, LTestCat[4].observed_statistic, LTestCat[5].observed_statistic, LTestCat[2].observed_statistic, LTestCat[3].observed_statistic]
pl_obs_df = pd.DataFrame(list(zip(pl_obs, spat_model)),
                            columns =['obs', 'model'])

pl_obs2 = [LTestCat2[0].observed_statistic, LTestCat2[1].observed_statistic, LTestCat2[4].observed_statistic, LTestCat2[5].observed_statistic, LTestCat2[2].observed_statistic, LTestCat2[3].observed_statistic]
pl_obs_df2 = pd.DataFrame(list(zip(pl_obs2, spat_model)),
                            columns =['obs', 'model'])

pl_obs3 = [LTestCat3[0].observed_statistic, LTestCat3[1].observed_statistic, LTestCat3[4].observed_statistic, LTestCat3[5].observed_statistic, LTestCat3[2].observed_statistic, LTestCat3[3].observed_statistic]
pl_obs_df3 = pd.DataFrame(list(zip(pl_obs3, spat_model)),
                            columns =['obs', 'model'])
                            
import matplotlib.lines as mlines
### Set up figure
fig, axes = plt.subplots(1, 3, sharey=True, figsize=(11,11))

# plot data as a boxenplot - shows 95% of results (k_depth=proportion) with 5% specified as outliers (outlier_prop)
sns.boxenplot(ax=axes[0], data=df_Num, palette="Set3", linewidth=0.5, orient="h",  k_depth='proportion', outlier_prop=0.05, showfliers=False)

## Draw lines to indicate observed numbers of events for each testing period
axes[0].axvline(NTestCat[0].observed_statistic, ls='--', color='#f2665b', label= '2006-11')
axes[0].axvline(NTestCat2[0].observed_statistic, ls=':', color='#5b7df2', label='2011-16')
axes[0].axvline(NTestCat3[0].observed_statistic, ls='-.', color='#7fc086', label='2016-21')

axes[0].set_title("Forecast n-test")
axes[0].set(xlabel='forecast # events')
axes[0].xaxis.get_label().set_fontsize(8)

## Same as for N-test (above)
sns.boxenplot(ax=axes[1], data=df_spat, palette="Set3", linewidth=0.5, orient="h", k_depth='proportion', outlier_prop=0.05, showfliers=False)

## Plot observed S-test statistics as symbols, label them so you can make a legend later!
sns.stripplot(ax=axes[1], x="obs", y="model", data= spat_obs_df, size=8, color="red", marker="*", linewidth=1, label= 'T1_mark')
sns.stripplot(ax=axes[1], x="obs", y="model", data= spat_obs_df2, size=6, color="blue", marker="D", linewidth=1, label='T2_mark')
sns.stripplot(ax=axes[1], x="obs", y="model", data= spat_obs_df3, size=6, color="green", marker="o", linewidth=1, label='T3_mark')

ax=axes[1].set_title("Spatial likelihood test")
axes[1].set(xlabel='normalised pseudo-likelihood')
axes[1].set(ylabel=' ')
axes[1].xaxis.get_label().set_fontsize(8)

## And again for PL-tests
sns.boxenplot(ax=axes[2], data=df_pl, palette="Set3", linewidth=0.5, orient="h",  k_depth='proportion', outlier_prop=0.05, showfliers=False)

### Observed PL as mark
sns.stripplot(ax=axes[2], x="obs", y="model", data= pl_obs_df, size=8, color="red", linewidth=1, marker="*")
sns.stripplot(ax=axes[2], x="obs", y="model", data= pl_obs_df2, size=6, color="blue", linewidth=1, marker="D")
sns.stripplot(ax=axes[2], x="obs", y="model", data= pl_obs_df3, size=6, color="green", linewidth=1, marker="o")

ax=axes[2].set_title("Pseudolikelihood")
axes[2].set(xlabel='normalised pseudo-likelihood')
axes[2].set(ylabel=' ')
axes[2].xaxis.get_label().set_fontsize(8)

### Make a legend explaining those symbols/coloured lines
T1line = mlines.Line2D([], [], color='#f2665b', ls='--', label='2006-2011')
T1mark = mlines.Line2D([], [], marker='*', ls='None', color = 'red', linewidth=1, markersize=6)
T2line = mlines.Line2D([], [], color='#5b7df2', ls=':', label='2011-2016')
T2mark = mlines.Line2D([], [], marker='D', ls='None', color = 'blue', markersize=6)
T3line = mlines.Line2D([], [], color='#7fc086', ls='-.', label='2016-2021')
T3mark = mlines.Line2D([], [], marker='o', ls='None', color = 'green', markersize=6)

lgd = plt.legend(handles=[T1line, T1mark, T2line, T2mark, T3line, T3mark], bbox_to_anchor=(1.05, 1))

plt.show()

print("Gridded forecasts")

SRMS = csep.load_gridded_forecast("Forecasts/SRMS_gridded_Km_2511.dat",
                                  start_date=start_date,
                                  end_date=end_date,
                                  name='SRMS')
## Plot to check this looks right!                                    
#SRMS.plot()
#plt.show()

SRMSms = csep.load_gridded_forecast("/Forecasts/Forecasts_1985_2005/SRMSms_gridded_Km_2303.dat",
                                  start_date=start_date,
                                  end_date=end_date,
                                  name='SRMS declustered')
                                  
FDSRMS = csep.load_gridded_forecast("/Forecasts/Forecasts_1985_2005/FDSRMS_gridded_Km_2303.dat",
                                  start_date=start_date,
                                  end_date=end_date,
                                  name='FDSRMS')
                                  
FDSRMSms = csep.load_gridded_forecast("/Forecasts/Forecasts_1985_2005/FDSRMSms_gridded_Km_2303.dat",
                                  start_date=start_date,
                                  end_date=end_date,
                                  name='FDSRMS declustered')
                                  
SRMSNK = csep.load_gridded_forecast("/Forecasts/Forecasts_1985_2005/SRMSNK_gridded_Km_2303.dat",
                                  start_date=start_date,
                                  end_date=end_date,
                                  name='SRMSNK')
                                  
SRMSNKms = csep.load_gridded_forecast("/Forecasts/Forecasts_1985_2005/SRMSNKms_gridded_Km_2303.dat",
                                  start_date=start_date,
                                  end_date=end_date,
                                  name='SRMSNK declustered')

## Load Helmstetter catalogues from pycsep
Helmstetter = csep.load_gridded_forecast(datasets.helmstetter_aftershock_fname,
                                         start_date=start_date,
                                         end_date=end_date,
                                         name='helmstetter_aftershock')

Helmstetter_dec = csep.load_gridded_forecast(datasets.helmstetter_mainshock_fname,
                                         start_date=start_date,
                                         end_date=end_date,
                                         name='helmstetter declustered')

                                         
## Check your catalogues line up with Helmstetter or your results will be inconsistent
r = Helmstetter.region
SR = SRMS.region
np.testing.assert_allclose(r.midpoints(), SR.midpoints())

## Define a function to run L, M. N and CL tests.
## I've set number of simulations here, you might want to adjust this if you're in a hurry but you should get more reproducible results with a high num_simulations

def alphabet_tests_gridded(forecast_list, catalog):
    LTests = []
    MTests = []
    NTests = []
    STests = []
    
    for i in range(len(forecast_list)):
        Lresult = poisson.conditional_likelihood_test(forecast_list[i], catalog, num_simulations=100000)
        LTests.append(Lresult)
    
        Mresult = poisson.magnitude_test(forecast_list[i], catalog, num_simulations=100000)
        MTests.append(Mresult)
    
        Nresult = poisson.number_test(forecast_list[i], catalog)
        NTests.append(Nresult)
    
        Sresult = poisson.spatial_test(forecast_list[i], catalog, num_simulations=100000)
        STests.append(Sresult)
    
    return LTests, MTests, NTests, STests

### Now list forecasts and apply the function to them. Takes a wee while with num_simulations so high.
my_forecasts = [SRMS, SRMSms, FDSRMS, FDSRMSms, SRMSNK, SRMSNKms, Helmstetter, Helmstetter_dec]

LtestT1, MtestT1, NtestT1, StestT1 = alphabet_tests_gridded(my_forecasts, comcat_catalog)
LtestT2, MtestT2, NtestT2, StestT2 = alphabet_tests_gridded(my_forecasts, comcat_catalog2)
LtestT3, MtestT3, NtestT3, StestT3 = alphabet_tests_gridded(my_forecasts, comcat_catalog3)

fig = plt.figure(figsize=(25, 2))
#subfigs = fig.subfigures(1, 2, wspace=0.07)
#fig.suptitle('Figure title')


# create 3x1 subfigs
subfigs = fig.subfigures(nrows=3, ncols=1)
subfigs[0].suptitle('2006-2011')

# create 1x3 subplots per subfig
axs = subfigs[0].subplots(nrows=1, ncols=3, sharey=True)
ax = plots.plot_poisson_consistency_test(NtestT1,
                                         plot_args={'xlabel': 'likelihood', 'tight_layout': False}, axes=axs[0])
                                        
ax = plots.plot_poisson_consistency_test(StestT1, one_sided_lower=True,
                                         plot_args={'xlabel': 'Spatial likelihood', 'tight_layout':False}, axes=axs[1])
ax =plots.plot_poisson_consistency_test(LtestT1,one_sided_lower=True,
                                         plot_args={'xlabel': 'forecast likelihood', 'tight_layout':False}, axes=axs[2])  

subfigs[1].suptitle('2011-2016')
# create 1x3 subplots per subfig
axs = subfigs[1].subplots(nrows=1, ncols=3, sharey=True)
ax = plots.plot_poisson_consistency_test(NtestT2,
                                         plot_args={'xlabel': 'likelihood', 'tight_layout':False}, axes=axs[0])
ax = plots.plot_poisson_consistency_test(StestT2, one_sided_lower=True,
                                         plot_args={'xlabel': 'Spatial likelihood', 'tight_layout':False}, axes=axs[1])
ax =plots.plot_poisson_consistency_test(LtestT2,one_sided_lower=True,
                                         plot_args={'xlabel': 'forecast likelihood', 'tight_layout':False}, axes=axs[2])    

subfigs[2].suptitle('2016-2021')
# create 1x3 subplots per subfig
axs = subfigs[2].subplots(nrows=1, ncols=3, sharey=True)
ax = plots.plot_poisson_consistency_test(NtestT3,
                                         plot_args={'xlabel': 'likelihood', 'tight_layout':False}, axes=axs[0])
ax = plots.plot_poisson_consistency_test(StestT3, one_sided_lower=True,
                                         plot_args={'xlabel': 'Spatial likelihood', 'tight_layout':False}, axes=axs[1])
ax =plots.plot_poisson_consistency_test(LtestT3,one_sided_lower=True,
                                         plot_args={'xlabel': 'forecast likelihood', 'tight_layout':False}, axes=axs[2])    

fig.set_size_inches(10, 10)
fig.savefig('Gridded_plots_v1.pdf')

## Set up list of forecasts
my_objects = [SRMS, SRMSms, FDSRMS, FDSRMSms, SRMSNK, SRMSNKms]

## For each time period, compare each model with Helmstetter model
paired_test_result1 = []

for i in range(len(my_objects)):
  result = poisson.paired_t_test(my_objects[i], Helmstetter, comcat_catalog)
  paired_test_result1.append(result)

paired_test_result2 = []

for i in range(len(my_objects)):
  result = poisson.paired_t_test(my_objects[i], Helmstetter, comcat_catalog2)
  paired_test_result2.append(result)


paired_test_result3 = []

for i in range(len(my_objects)):
  result = poisson.paired_t_test(my_objects[i], Helmstetter, comcat_catalog3)
  paired_test_result3.append(result)

## Plot args for each figure 
comp_args1 = {
    'xlabel':'model',
    'title':'2006-2011',
    'ylabel': 'information gain per earthquake',
    'xlabel_fontsize' : '10',
    'ylabel_fontsize' : '10'
}
comp_args2 = {
    'xlabel':'model',
    'title':'2011-2016',
    'ylabel': 'information gain per earthquake',
    'xlabel_fontsize' : '10',
    'ylabel_fontsize' : '10'
}
comp_args3 = {
    'xlabel':'model',
    'title':'2016-2021',
    'ylabel': 'information gain per earthquake',
    'xlabel_fontsize' : '10',
    'ylabel_fontsize' : '10'
}

fig, axes = plt.subplots(1, 3, figsize=(40,50)) 
ax1 = plots.plot_comparison_test(paired_test_result1, plot_args= comp_args1, axes=axes[0])
ax2 = plots.plot_comparison_test(paired_test_result2, plot_args= comp_args2, axes=axes[1])
ax3 = plots.plot_comparison_test(paired_test_result3, plot_args= comp_args3, axes=axes[2])

plt.show()

import cartopy 

f, axes = plt.subplots(3, 1, figsize=(11, 6))

args_dict1 = {'title': '2006-2011',
    'grid_labels': True,
    'borders': True,
    'feature_lw': 0.5,
    'basemap': 'ESRI_topo',
    'projection': cartopy.crs.Mercator(),
    'markercolor': 'black'}

axes[0] = plots.plot_catalog(comcat_catalog, ax=None, show=True, extent=None, plot_args=args_dict1)
plt.savefig('CalCat06-11.pdf')

args_dict1 = {'title': '2011-2016',
    'grid_labels': True,
    'borders': True,
    'feature_lw': 0.5,
    'basemap': 'ESRI_topo',
    'projection': cartopy.crs.Mercator(),
    'markercolor': 'black'}

axes[1] = plots.plot_catalog(comcat_catalog2, ax=None, show=True, extent=None, plot_args=args_dict1)
plt.savefig('CalCat11-16.pdf')

args_dict1 = {'title': '2016-2021',
    'grid_labels': True,
    'borders': True,
    'feature_lw': 0.5,
    'basemap': 'ESRI_topo',
    'projection': cartopy.crs.Mercator(),
    'markercolor': 'black'}

axes[2] = plots.plot_catalog(comcat_catalog3, ax=None, show=True, extent=None, plot_args=args_dict1)
plt.savefig('CalCat16-21.pdf')

plt.show()

start_date4 = time_utils.strptime_to_utc_datetime('1984-01-01 00:00:00.0')
end_date4 = time_utils.strptime_to_utc_datetime('2021-01-01 00:00:00.0')

comcat_catalog4 = csep.query_comcat(start_date4, end_date4, min_magnitude=SRMS.min_magnitude)

# Filter observed catalog using the same region as the forecast
comcat_catalog4 = comcat_catalog4.filter_spatial(SRMS.region)
plots.plot_catalog(comcat_catalog4, ax=None, show=True, extent=None)

cati = pd.DataFrame(comcat_catalog4.catalog)
plt.clf()

cati['dt'] = pd.to_datetime(cati['origin_time'], unit='ms')

plt.axvspan(pd.to_datetime("2006-1-1"), pd.to_datetime("2011-1-1"), facecolor='#f8ab9a', alpha=0.3, zorder=1)
plt.axvline(pd.to_datetime("2006-1-1"), color="red", alpha=0.5, ls="--", zorder=2)
plt.axvspan(pd.to_datetime("2011-1-1"), pd.to_datetime("2016-1-1"), facecolor='#9ab8f8', alpha=0.3, zorder=1)
plt.axvline(pd.to_datetime("2011-1-1"), color="blue", alpha=0.5, ls=":", zorder=2)
plt.axvspan(pd.to_datetime("2016-1-1"), pd.to_datetime("2021-1-1"), facecolor='#abfbb3', alpha=0.3, zorder=1)
plt.axvline(pd.to_datetime("2016-1-1"), color="green", alpha=0.5, ls="-.", zorder=2)
sns.scatterplot(data=cati,  x="dt", y="magnitude", zorder=3, linewidth=2, color="darkblue")

plt.ylabel("Magnitude")
plt.xlabel("Year")

#plt.xlabel('Date')
#plt.ylabel('magnitude')
#plt.show()
#plt.savefig('Cat.pdf')

cati['year'] = pd.DatetimeIndex(cati['dt']).year
yr_counts = cati[['magnitude','year']].groupby('year').count()

plt.clf()
#yrs= cati['year'].unique()
#yr_df = pd.DataFrame(dict(year = yrs, count = yr_counts))

plt.axvspan(2006, 2011, facecolor='#f8ab9a', alpha=0.3, zorder=1)
plt.axvline(2006, color="red", alpha=0.5, ls="--", zorder=2)
plt.axvspan(2011, 2016, facecolor='#9ab8f8', alpha=0.3, zorder=1)
plt.axvline(2011, color="blue", alpha=0.5, ls=":", zorder=2)
plt.axvspan(2016, 2021, facecolor='#abfbb3', alpha=0.3, zorder=1)
plt.axvline(2016, color="green", alpha=0.5, ls="-.", zorder=2)
sns.scatterplot(data=yr_counts,  x="year", y="magnitude", zorder=3, linewidth=2, color="darkblue")

plt.ylabel( "Events per year with M > 4.95")
plt.xlabel("Year")
plt.show()
plt.savefig('events_per_year.pdf')
