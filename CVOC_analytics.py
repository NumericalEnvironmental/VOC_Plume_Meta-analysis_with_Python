###################################################################################
#
# CVOC_analytics.py
#
# by Walt McNab
#
# exploration of chlorinated hydrocarbon plume data from GeoTracker/GAMA database
# GAMA data = California Groundwater Ambient Monitoring and Assessment Program
#
###################################################################################

from numpy import *
from pandas import *
from datetime import *
from sklearn.cluster import DBSCAN
import pyproj
from scipy import stats
from scipy.spatial import distance_matrix, Delaunay
import seaborn as sns


# options.mode.chained_assignment = None


#####################
#
# class definitions
#
#####################


class SiteSummary:

    # holds data frame for sites statistical summary

    def __init__(self, first_site_df, analytes, phi):
        # initialize summary data frame
        self.wgs84=pyproj.Proj("+init=EPSG:4326")           # geographic coordinate system references
        self.utm10N=pyproj.Proj("+init=EPSG:32610")         # note that UTM units = meters
        self.utm11N=pyproj.Proj("+init=EPSG:32611")
        self.phi = phi                                      # porosity (universal); used to calculate mass/depth
        self.summary_df = self.SummaryStats(first_site_df, analytes)        

    def AppendSite(self, site_df, analytes):
        # process new site and append to summary data frame
        new_summary_df = self.SummaryStats(site_df, analytes)
        self.summary_df = self.summary_df.append(new_summary_df)

    def SummaryStats(self, site_df, analytes):
        # summarize site characteristics; return to append to all-sites summary
        site_summary_series = site_df[analytes].max(axis=0)
        site_summary_df = site_summary_series.to_frame()
        site_summary_df = site_summary_df.transpose()
        for name in analytes: site_summary_df.rename(columns={name: name+'_maxC'}, inplace=True)        
        site_summary_df['site_index'] = array(site_df['cluster'])[0]
        site_summary_df['start_date'] = site_df['DATE'].min()
        wells_df, long_mean, lat_mean = self.WellSpatial(site_df)        # extract well coordinates (median concentrations returned)
        site_summary_df['num_wells'] = len(wells_df)
        site_summary_df['longitude'] = long_mean
        site_summary_df['latitude'] = lat_mean        
        loc_set = wells_df[['x', 'y']].values                   
        dm = distance_matrix(loc_set, loc_set)                      # find maximum interwell distance
        site_summary_df['max_extent'] = dm.max()
        tot_area, mass_df = self.Triangles(loc_set, wells_df[analytes], analytes)
        site_summary_df['tot_area'] = tot_area
        site_summary_df = concat([site_summary_df, mass_df], axis=1)
        #re-order columns ...
        header = ['site_index', 'longitude', 'latitude', 'start_date', 'num_wells', 'max_extent', 'tot_area']
        for name in analytes: header.append(name+'_maxC')
        for name in analytes: header.append(name+'_mass_L')		
        site_summary_df = site_summary_df[header]
        return site_summary_df

    def Triangles(self, loc_set, median_analytes_df, analytes):
        # parse site wells into delauney traingles, calculate areas
        median_analytes_df.fillna(value = 0.0, inplace=True)
        tot_area = 0.                       # total area serviced by well network
        tri = Delaunay(loc_set)
        tri_indices = tri.simplices         # index numbers of triangles
        for i, idx in enumerate(tri_indices):
            p = loc_set[idx, :]
            area = self.PolygonArea(p)
            c_vertices = median_analytes_df.iloc[idx, :]               # median historical concs, at vertices
            c_rep_series = c_vertices.median()
            if not i: mass_series = c_rep_series * area * 1e-6 *self.phi                   # initial mass Series (per analyte); convert ug/L to kg/m^3
            else: mass_series = mass_series.add(c_rep_series * area * 1e-6 * self.phi)                     # add to mass, per analyte, for this site
            tot_area += area
        mass_df = mass_series.to_frame()
        mass_df = mass_df.transpose()
        for name in analytes: mass_df.rename(columns={name: name+'_mass_L'}, inplace=True) 
        return tot_area, mass_df
        
    def PolygonArea(self, p):
        # classic cross-product planimeter approach for calculating area of polygon by vertices (p = array of points)
        pc = append(p, [p[0]], axis=0)
        s = (pc[1: -1,0] * (pc[2:, 1] - pc[0:-2, 1])).sum()
        s += pc[-1, 0] * (pc[1, 1] - pc[-2, 1])
        return abs(s/2)        

    def SubZone(self, sub_df, geo_sys):
        # return transformed coordinates for single (UTM) zone
        lon = array(sub_df['APPROXIMATE LONGITUDE'])
        lat = array(sub_df['APPROXIMATE LATITUDE'])
        x, y = pyproj.transform(self.wgs84, geo_sys, lon, lat)
        sub_df['x'] = x
        sub_df['y'] = y
        return sub_df

    def WellSpatial(self, site_df):
        # return well information for site (locations in UTM)
        wells_df = site_df.groupby('WELL NAME').median()
        wells_df.reset_index(inplace=True)
        long_mean = wells_df['APPROXIMATE LONGITUDE'].mean()       # site centroid
        lat_mean = wells_df['APPROXIMATE LATITUDE'].mean()        
        wells_w_df = wells_df[wells_df['APPROXIMATE LONGITUDE']<=-120.]      # for UTM Zone 10N
        wells_w_df = self.SubZone(wells_w_df, self.utm10N)                     
        wells_e_df = wells_df[wells_df['APPROXIMATE LONGITUDE']>-120.]      # for UTM Zone 11N
        wells_e_df = self.SubZone(wells_e_df, self.utm11N)
        wells_df = concat([wells_w_df, wells_e_df], axis=0)                    # recombine data frames
        wells_df.reset_index(inplace=True)
        return wells_df, long_mean, lat_mean
        

########################
#
# supporting functions
#
########################


def Histogram(data_set, bins, x_title, label, file_name, fit='none'):
    # draw histogram(s) of data set(s)
    if fit == 'none':
        for i, data in enumerate(data_set): ax = sns.distplot(data, bins=bins, kde=False, norm_hist=False, hist_kws={'label': label[i]})    
        ax.set(xlabel=x_title, ylabel='Count')
    else:
        for i, data in enumerate(data_set): ax = sns.distplot(data, bins=bins, kde=False, fit=fit, hist_kws={'label': label[i]})    
        ax.set(xlabel=x_title, ylabel='Probability density')
    ax.legend()
    fig = ax.get_figure()
    fig.savefig(file_name)    
    sns.plt.show()     


def ScatterPlot(data, x_col, y_col, x_label, y_label, file_name, hexplot='none'):
    # create scatter plot, with bivariate distributions and KDE
    if hexplot == 'none':
        g = (sns.jointplot(x=x_col, y=y_col, data=data).plot_joint(sns.kdeplot, zorder=0, n_levels=6))
    else:
        g = sns.jointplot(x=x_col, y=y_col, data=data, kind='hex')
    g.set_axis_labels(x_label, y_label)
    g.savefig(file_name) 
    sns.plt.show() 


def ReadGama(analytes, headings, file_name):

    # read entire gama data file
    gama_df = read_csv(file_name, sep='\t')
    mask = gama_df['CHEMICAL'].isin(analytes)                   # limit chemicals in dataframe to those on analytes list
    gama_df = gama_df[mask]
    gama_df['DATE'] = to_datetime(gama_df['DATE'])                # format dates
    gama_df = gama_df[gama_df['DATASET']=='EDF']                  # exclude data NOT from shallow environmental monitor wells

    # handle unit conversions
    ppm_match = array(gama_df['UNITS']=='MG/L')                                                  
    gama_df['CONC'] = 1000.*ppm_match*gama_df['RESULT'] + (1-ppm_match)*gama_df['RESULT']

    # exclude any non-detects
    gama_df = gama_df[(gama_df['CONC']>0.) & (gama_df['QUALIFIER']!='<')]

    # create pivot table
    gama_df = pivot_table(gama_df, values='CONC', index=headings, columns=['CHEMICAL'])
    gama_df.reset_index(inplace=True)

    return gama_df


def ReadList(input_file):
    # read in a simple list of strings from file and return
    item = []
    input_file = open(input_file, 'r')
    for line in input_file:
        name = line.split()
        item.append(name[0])
    input_file.close()
    return item


def CheckMinWells(site_df, min_wells):
    # check to see if minimum number of distinct wells exist at site (return boolean)
    loc_array = site_df[['APPROXIMATE LATITUDE', 'APPROXIMATE LONGITUDE']].values
    well_list = tuple(map(tuple, loc_array))
    num_wells = len(list(set(well_list)))
    return (num_wells>=min_wells)


#####################################
#
# main script
#
#####################################


def Process_gama(read_opt):

    # definitions
    data_files = ReadList('counties.txt')            # list of files containing complete county GAMA data (one file per county)
    analytes = ReadList('analytes.txt')              # list of analytes (CVOCs) in GAMA data to be processed
    headings = ['WELL NAME', 'APPROXIMATE LATITUDE', 'APPROXIMATE LONGITUDE', 'DATE', 'DATASET', 'COUNTY']          # for dataframe processing
    min_wells = 5                                   # minimum number of monitoring wells to qualify a site
    phi = 0.25                                      # aquifer porosity (universal)
    start_date = '01-01-2000'
    end_date = '01-01-2020'
    min_date = to_datetime(start_date)
    max_date = to_datetime(end_date)

    if not read_opt:        # read GAMA data sets
        for i, data_file in enumerate(data_files):              # extract GAMA data from multiple county data files
            print 'Reading', data_file
            county_gama_df = ReadGama(analytes, headings, data_file)
            if not i: gama_df = county_gama_df
            else:
                if len(county_gama_df): gama_df = concat([gama_df, county_gama_df], axis=0)
        gama_df.to_csv('gama_df.csv', index=False)      # write to summary data file
    else:                                           # read summary data file
        print 'Reading combined data file'
        gama_df = read_csv('gama_df.csv', sep=',')
        gama_df['DATE'] = to_datetime(gama_df['DATE'])                # format dates

    # filter date window
    gama_df = gama_df[(gama_df['DATE']>=min_date)]
    gama_df = gama_df[(gama_df['DATE']<=max_date)]

    # parse sites with DBSCAN cluster analysis
    print 'Performing DBSCAN to parse sites ...'
    X = gama_df[['APPROXIMATE LATITUDE', 'APPROXIMATE LONGITUDE']].values
    dbscan_analysis = DBSCAN(eps=0.00075, min_samples=5)                          # eps = 0.002 degree seems to work (min distance between clsuter wells) ...
    cluster = dbscan_analysis.fit_predict(X)                                                
    gama_df['cluster'] = cluster
    cluster_array = arange(0, cluster.max()+1)
    
    # drop data for which (1) no cluster has been assigned, or (2) there are fewer than min_wells distinct wells
    print 'Disqualifying posited sites, as warranted ...'
    gama_df = gama_df[gama_df['cluster']!=-1]
    remove_cluster = []
    for c in cluster_array:
        site_df = gama_df[gama_df['cluster']==c]
        keep_flag = CheckMinWells(site_df, min_wells)
        if not keep_flag: remove_cluster.append(c)            # add to list to be removed
    remove_cluster = array(remove_cluster)
    cluster_array = delete(cluster_array, remove_cluster, None)
    mask = gama_df['cluster'].isin(cluster_array)                                                   # limit chemicals in dataframe to those on analytes list
    gama_df = gama_df[mask]         
    gama_df.to_csv('gama_processed_df.csv')		# write to embellished and edited summary file

    # run summary statistics on lumped site data sets
    for i, c in enumerate(cluster_array):
        site_df = gama_df[gama_df['cluster']==c]
        if not i: site_summary = SiteSummary(site_df, analytes, phi)      # instantantiate all-sites summary object
        else: site_summary.AppendSite(site_df, analytes)
    site_summary.summary_df['tot_area'] = 0.0001*site_summary.summary_df['tot_area']    # convert from m^2 to hectares
    site_summary.summary_df['well_density'] = site_summary.summary_df['num_wells']/site_summary.summary_df['tot_area']         

    # write site summary characteristics to file
    site_summary.summary_df.to_csv('all_sites_summary.csv')


    #####################################################
    #
    # generate example plots
    #
    #####################################################


    sns.set(color_codes=True)
    
    # (1a) Maximum VOC plume extent histogram
    data_set = [log10(site_summary.summary_df['max_extent'])]
    bins = linspace(1., 4., num=21, endpoint=True)   
    Histogram(data_set, bins, 'Plume extent (log m)', ['plume extent'], 'plume_extent.png', fit=stats.norm)

    # (1b) Plume area footprint histogram
    data_set = [log10(site_summary.summary_df['tot_area'])]
    bins = linspace(-2., 3., num=21, endpoint=True)   
    Histogram(data_set, bins, 'Plume area (log hectares)', ['plume area'], 'plume_area.png', fit=stats.norm)

    # (2a) TCE - cis-12-DCE - vinyl chloride (reductive dehalogenation sequence) max. concs histograms
    data_set = []
    rd_halo = ['TCE', 'DCE12C', 'DCE12T', 'VC']
    for voc in rd_halo:
        col_name = voc + '_maxC'
        df = site_summary.summary_df[site_summary.summary_df[col_name]>0.]    
        data_set.append(log10(df[col_name]))
    bins = linspace(-1., 6., num=21, endpoint=True)
    Histogram(data_set, bins, 'Maximum concentration (log ug/L)', rd_halo, 'rd_vocs.png')

    # (2b) TCE and PCE mass/unit-depth histograms
    data_set = []
    mass_depth = ['TCE', 'PCE']
    for voc in mass_depth:
        col_name = voc + '_mass_L'
        df = site_summary.summary_df[site_summary.summary_df[col_name]>0.]    
        data_set.append(log10(df[col_name]))
    bins = linspace(-3., 3., num=21, endpoint=True)
    Histogram(data_set, bins, 'Mass/aquifer thickness (log kg/m)', mass_depth, 'mass_depth.png')
    
   # (3) well density histogram for site with a total area >= 1 hectare
    larger_sites_df = site_summary.summary_df[site_summary.summary_df['tot_area']>=1.0]
    data_set = [larger_sites_df['well_density']]
    bins = linspace(0., 50., num=31, endpoint=True)   
    Histogram(data_set, bins, 'Well density (no. of wells/hectare)', ['well density'], 'well_density.png', fit=stats.gamma)

    # (4a) 1,1-DCE versus 1,1,1-TCA max concs, per site where both are present (example)
    x_col = 'TCA111_maxC'
    y_col = 'DCE11_maxC'
    df = site_summary.summary_df[[x_col, y_col]]
    df = df[(df[x_col]>0.) & (df[y_col]>0.)]    
    ScatterPlot(log10(df), x_col, y_col, '111-TCA max conc. (log ug/L)', '11-DCE max conc. (log ug/L)', 'DCE11_vs_TCA111.png')

    # (4b) 1,1-DCA versus PCE max concs, per site where both are present (example)
    x_col = 'PCE_maxC'
    y_col = 'DCA11_maxC'
    df = site_summary.summary_df[[x_col, y_col]]
    df = df[(df[x_col]>0.) & (df[y_col]>0.)]    
    ScatterPlot(log10(df), x_col, y_col, 'PCE max conc. (log ug/L)', '11-DCA max conc. (log ug/L)', 'DCA11_vs_PCE.png')

    # (5a) no. of wells versus PCE max concs, per site (example)
    x_col = 'PCE_maxC'
    y_col = 'num_wells'
    df = site_summary.summary_df[[x_col, y_col]]
    df = df[df[x_col]>0.]    
    ScatterPlot(log10(df), x_col, y_col, 'PCE max conc. (log ug/L)', 'No. of wells (log)', 'num_wells vs PCE_maxC.png')
    
    # (5b) PCE mass versus no. of wells, per site (example)
    x_col = 'num_wells'
    y_col = 'PCE_mass_L'
    df = site_summary.summary_df[[x_col, y_col]]
    df = df[df[y_col]>0.]    
    ScatterPlot(log10(df), x_col, y_col, 'No. of wells (log)', 'PCE mass (log kg/m)', 'PCE_mass vs wells.png', hexplot='hex')
    
    # (6a) 1,4-dioxane versus 1,1,1-TCA mass per unit depth
    x_col = 'TCA111_mass_L'
    y_col = 'DIOXANE14_mass_L'
    df = site_summary.summary_df[[x_col, y_col]]
    df = df[(df[x_col]>0.) & (df[y_col]>0.)]    
    ScatterPlot(log10(df), x_col, y_col, '1,1,1-TCA mass (log kg/m)', '1,4-dioxane mass (log kg/m)', 'dioxane_vs_TCA.png')
   
    # (6b) 1,4-dioxane versus 1,1-DCE mass per unit depth
    x_col = 'DCE11_mass_L'
    y_col = 'DIOXANE14_mass_L'
    df = site_summary.summary_df[[x_col, y_col]]
    df = df[(df[x_col]>1e-6) & (df[y_col]>1e-6)]    
    ScatterPlot(log10(df), x_col, y_col, '1,1-DCE mass (log kg/m)', '1,4-dioxane mass (log kg/m)', 'dioxane_vs_DCE11.png')
    
    # (6c) 1,4-dioxane versus TCE mass per unit depth
    x_col = 'PCE_mass_L'
    y_col = 'DIOXANE14_mass_L'
    df = site_summary.summary_df[[x_col, y_col]]
    df = df[(df[x_col]>0.) & (df[y_col]>0.)]    
    ScatterPlot(log10(df), x_col, y_col, 'PCE mass (log kg/m)', '1,4-dioxane mass (log kg/m)', 'dioxane_vs_PCE.png')
    
    print 'Done.'


##### run script #####
#Process_gama(read_opt)
read_opt = 0 		# --> read from GAMA database file (by county)
read_opt = 1 		# --> read from CVOC summary buffer file
Process_gama(1)

