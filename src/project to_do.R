

#notes from 2/11 gates meeting:
#-proxy project to support better DALY estimations
#before the positive control analysis, test predictions on VMNIS subnational data

#merge HFID data for all countries
#"C:\Users\andre\OneDrive\Documents\mn-proxies\data\HFID\hfid_hv1.csv"

#merge Cadre Harmonize data for all countries
#"C:\Users\andre\OneDrive\Documents\mn-proxies\data\CadreHarmonise\cadre_harmonise_caf_ipc_dec25.xlsx"


#If I'm going to get subnational predictions in other countries, I really think
#I'm going to need subnational prevalence models. Aggregate the surveys at admin2
#Predict the prevalences directly (maybe use old SL package). New country data
#could be aggregated at that level and then predictions made (bootstrap CIs)?
#Use DHS + Cadre harmonize + can I get GEE averages by region.
#Add midpoint GPS as a covariate

#Do positive control analysis with anemia or with wasting


#merge in Cadre Harmonise data into all 4 countries
#https://data.humdata.org/dataset/harmonized-food-insecurity-dataset-hfid
#get the admin2 DHS estimates
#get GEE soil estimates



#merge in Individual level MICS data into Gambia (it was a subsample)




#Figure out why the cluster_SL code has extreme MSE and errors



# Harmonized Food Insecurity Dataset (HFID)
#
# The HFID dataset consists of files with administrative units geometries and one file with HFID variables. The geometries are found (1) in the GADM geometries with distinct administrative levels 0,1 and 2 (GADM.zip), (2) in the geopackage file with combined geometries at relevant administrative level with exact geometries (hfid_all_geom.gpkg) or simplified geometries (simplified_hfid_geom.gpkg). The GADM.zip file contains the shapefiles with administrative level 0 (gadm_410_L0), level 1 (gadm_410_L1) and administrative level 2 (gadm_410_L2) reference geometries.
#
#
# The data table with the HFID variables (HFID_hv1.csv) is composed of several columns:
#
#   â€¢ year_month: indicates the year and month of the record in the format "%Y-%m".
# â€¢ ADMIN_0: name of the administrative level 0 (country) as defined in the GADM.
# â€¢ ADMIN_1: name of the administrative level 1 (region) as defined in the GADM.
# â€¢ ADMIN_2: name of the administrative level 2 (sub-region) as defined in the GADM.
# â€¢ ipc_phase_fews: FEWS NET-IPC compatible phase classification (phase-FEWS variable) if available, else NaN. Possible values: 1,2,3,4,5.
# â€¢ ha_fews: 1 if humanitarian aid has been distributed (owing to FEWS NET source), else NaN.
# â€¢ ipc_phase_ipcch: IPC/CH Phase classification (phase-IPC/CH variable) if available, else NaN. Possible values: 1,2,3,4,5,6.
# â€¢ ha_ipcch: 1 if humanitarian aid has been distributed (owing to IPC API source), else NaN.
# â€¢ set_ipcch: 1 if HouseHold Group or Internally Displaced People or Urban settlements are found in the area (owing to
#                                                                                                                IPC API source), else NaN.
# â€¢ rfg_ipcch: 1 if refugees are found in the area (owing to IPC API source), else NaN.
# â€¢ fcs_lit: monthly average of population prevalence of insufficient food consumption score from the WFP-LIT source (FCS-LIT variable) if available, else NaN. Possible values: between 0 and 1.
# â€¢ rcsi_lit monthly average of population prevalence of crisis or above reduced Coping Strategy Index from the WFP-LIT source (FCS-LIT variable) if available, else NaN. Possible values: between 0 and 1.
# â€¢ fcs_rt mean: monthly average of population prevalence of insufficient food consumption score from the WFP-RT source (FCS-RT variable) if available, else NaN. Possible values: between 0 and 1.
# â€¢ fcs_rt max: monthly maximum of population prevalence of insufficient food consumption score from the WFP-RT sourceif available, else NaN. Possible values: between 0 and 1.
# â€¢ fcs_rt min: monthly minimum of population prevalence of insufficient food consumption score from the WFP-RT source if available, else NaN. Possible values: between 0 and 1.
# â€¢ rcsi_rt mean: monthly average of population prevalence of crisis or above reduced Coping Strategy Index from the WFP-RT source (FCS-RT variable) if available, else NaN. Possible values: between 0 and 1.
# â€¢ rcsi_rt max: monthly maximum of population prevalence of crisis or above reduced Coping Strategy Index from the WFP-RT source if available, else NaN. Possible values: between 0 and 1.
# â€¢ rcsi_rt min: monthly minimum of population prevalence of crisis or above reduced Coping Strategy Index from the WFP-RT source if available, else NaN. Possible values: between 0 and 1.
# â€¢ iso2: ISO 3166-1 code for country in two letters.
# â€¢ iso3: ISO 3166-1 code for country in three letters.
# â€¢ region: region name from the United Nations Geoscheme.
#
# Cite this article
# Machefer, M., Ronco, M., Thomas, AC. et al. A monthly sub-national Harmonized Food Insecurity Dataset for comprehensive analysis and predictive modeling. Sci Data 12, 741 (2025). https://doi.org/10.1038/s41597-025-05034-4
#
#
# Cite this dataset
# Machefer, M., Ronco, M., Thomas, A.-C., Assouline, M., Rabier, M., Corbane, C., & Rembold, F. (2025). Harmonized Food Insecurity Dataset (HFID): data and code (v1.1.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.15017473
