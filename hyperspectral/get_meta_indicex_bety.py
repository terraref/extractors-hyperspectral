#!/usr/bin/python
import os
import sys
import json
import urllib2


#new_var_List=["EVI","ARVI", "GEMI","GARI","DVI", "GNDVI","GRVI","IPVI","LAI","RDVI","MSR","NLI","MNLI","SAVI","TDVI","VARI","RENDVI","mRESR","mRENDVI","VOG1","VOG2","VOG3","MCARI","MCARI1","MCARI2","MTVI", "MTVI2","GMI1","GMI2","G","Lic1","Lic2","Lic3","PRI","NDNI","NRI510","NRI850","NDL1","CAI","PSRI","CRI1","CRI2","ARI1","ARI2","SRPI","NPQI","NPCI","WBI","NDWI","MSI","NDII","NMDI","HI","CLSI","SBRI","PMI","DWSI","Crt1","Crt2","BIG2","LSI","BRI"]



# var_List=["NDVI","R", "OSAVI","CHL","msR705","TCARI","CarChap","Car1Black","Car2Black","PRI570","IPI","antGamon","antGitelson","CHLDela","CI","PRI586","PRI512","FRI1","FRI2","NDVI1","RERI","ZM","REP","NDRE","TVI"];


rootUrl='https://terraref.ncsa.illinois.edu/bety/api/beta/variables?key=9999999999999999999999999999999999999999&limit=none&name='



queryAllUrl='https://terraref.ncsa.illinois.edu/bety/api/beta/variables?key=9999999999999999999999999999999999999999&limit=none&type=%7EReflectance+index'

field_List=["long_name","standard_name", "description", "notes","units"]



DEBUG=0


# returns the text content of an http request
def myGet(Url):
    req = urllib2.Request(Url);
    f = urllib2.urlopen(req)
    #req = http.request('GET', Url)
    sLines=f.read();
  
    if DEBUG:
        sys.stderr.write("json-response\n");
        sys.stderr.write(sLines);
          
    if len(sLines)<10:
        sys.stderr.write("requested json response too short\n")
        return "";
    else:
        return sLines

def prnVar(myVar):


    nameInd=myVar['name']   

    if 'label' in myVar:
       myVar['long_name']=myVar['label']
  

    for field in field_List:
        try:
            value=myVar[field].encode('utf8')
        except KeyError:
            continue;

        if len(value)>0:
            print '\'{0}@{1}\'="{2}";'.format(nameInd, field, value)



#get ajson data for all indices
sLines=myGet(queryAllUrl)
js=json.loads(sLines)

for myD in js["data"]:
    prnVar (myD["variable"])




##########################  example json request ##########################################################################################################

# {
#   "metadata": {
#     "URI": "https://terraref.ncsa.illinois.edu/bety/api/beta/variables?key=9999999999999999999999999999999999999999&limit=none&name=OSAVI",
#     "timestamp": "2017-03-21T05:43:17-05:00",
#     "count": 1
#   },
#   "data": [
#     {
#       "variable": {
#         "id": 6000000022,
#         "description": "Optimized Soil Adjusted Vegetation Index Formula: (1 + 0.16) * (R800 - R670)/(R800 + R670 + 0.16) Rondeaux et al. (1996)",
#         "units": "ratio",
#         "notes": "Optimized Soil Adjusted Vegetation Index OSAVI",
#         "created_at": "2017-02-08T15:19:59-06:00",
#         "updated_at": "2017-03-15T14:31:42-05:00",
#         "name": "OSAVI",
#         "max": "Infinity",
#         "min": "0",
#         "standard_name": "optimized_soil_adjusted_vegetation_index",
#         "standard_units": null,
#         "label": "Optimized Soil Adjusted Vegetation Index",
#         "type": "Reflectance Index",
#         "number of associated covariates": 0,
#         "number of associated traits": 0,
#         "number of associated formats_variables": 0,
#         "number of associated formats": 0,
#         "number of associated priors": 0,
#         "number of associated likelihoods": 0,
#         "view_url": "https://terraref.ncsa.illinois.edu/bety/variables/6000000022",
#         "edit_url": "https://terraref.ncsa.illinois.edu/bety/variables/6000000022/edit"
#       }
#     }
#   ]
# }


