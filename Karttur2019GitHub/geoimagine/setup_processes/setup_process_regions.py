'''
Created on 21 feb. 2018

@author: thomasgumbricht
'''

from geoimagine.postgresdb import ManageProcess, SelectUser, ManageMODIS, ManageRegion, ManageAncillary, ManageSentinel, ManageLandsat
#from geoimagine.setup_processes import ProcessProcess
from geoimagine.ancillary import ProcessAncillary
from geoimagine.modis import ProcessModis
from geoimagine.sentinel import ProcessSentinel
from geoimagine.kartturmain import runxml, MainProc
from geoimagine.support import ConvertLandsatScenesToStr
import csv 
#import geoimagine.gis.mj_gis_v80 as mj_gis

from os import path, makedirs
from sys import exit

#import geoimagine.gis.gis as mj_gis
import geoimagine.gis.mj_gis_v80 as mj_gis

class ProcessDefaultRegions:
    '''
    '''  
    def __init__(self, process,session,verbose):
        #self.session = session
        self.process = process
        self.verbose = verbose
        #print ('self.process.proc.userProj.system',self.process.proc.system.system)
        print ('self.process.proc.userProj.srcsystem',self.process.system.srcsystem)
        print ('self.process.proc.userProj.dstsystem',self.process.system.dstsystem)
        print ('self.process.proc.userProj.dstsystem',self.process.system.system)

        if (self.process.system.dstsystem) != 'system':
            exit('setup_process_class - setting default regions requires system = system ')

        #direct to subprocess
        if self.process.proc.processid == 'regioncategories':
            self._InsertRegionCat(session)
        elif self.process.proc.processid == 'defaultregion':
            self._DefaultRegion(session)
        elif self.process.proc.processid == 'defaultregionfromvector':
            self._DefaultRegFromVec(session)
        elif self.process.proc.processid == 'linkregionswrs':
            if self.process.proj.projectid == "karttur" and self.process.proj.system == 'system': 
                self.LinkAllRegionsToWRS()
        elif self.process.proc.processid == 'linkregionsmodtiles':
            if self.process.proj.projectid == "karttur" and self.process.proj.system == 'system': 
                self.LinkAllRegionsToMODIS()
            else:
                exit('Only superuser can link wrs to regions')
        else:
            exitstr = 'No process %s under Processregion' %(self.process.processid)
            exit(exitstr)
            
    def _InsertRegionCat(self,session): 
        #query = dict((x,y)  for x,y in self.process.proc.parameters.__dict__.items())
        session._InsertRegionCat(self.process, self.process.proc.paramsD)
        
    def _DefaultRegion(self,session):
        '''
        '''
        for locus in self.process.dstLayerD:
            for datum in self.process.dstLayerD[locus]:
                for comp in self.process.dstLayerD[locus][datum]:
                    query = self.process.proc.paramsD
                    fieldDD = self._SetfieldD( query['regionid'], query['regionname'], query['regioncat'], query['stratum'], query['parentid'], query['parentcat'])
                    layer = self.process.dstLayerD[locus][datum][comp]

                    #The destination region must be forced,this is because the locus to be created did not exists when chekcing for the feault locus
                    if self.verbose:
                        print ('        forcing region from to',locus, self.process.proc.paramsD['regionid'])
                    
                    if layer.locus.locus != self.process.proc.paramsD['regionid']:
                        layer.locus.locus = self.process.proc.paramsD['regionid']
                    
                    if layer.locus.path != self.process.proc.paramsD['regionid']:
                        layer.locus.path = self.process.proc.paramsD['regionid']
                                                                      

                    layer._SetPath()
                    
                    if self.verbose:
                        print ('        filepath', layer.FPN)
   
                    layer.CreateVectorAttributeDef(fieldDD)
                    layer._SetBounds(query['epsg'],query['minlon'],query['minlat'],query['maxlon'],query['maxlat'])

                    projection = mj_gis.MjProj()
                    projection.SetFromEPSG(query['epsg'])

                    if not layer._Exists() or self.process.proc.overwrite:
                        mj_gis.CreateESRIPolygonPtL(layer.FPN, layer.fieldDefL, layer.BoundsPtL, projection.proj_cs, query['regionid'])          
                    boundsD = mj_gis.GetFeatureBounds(layer.FPN,'REGIONID')
                    
                    #Set lonlat projection
                    lonlatproj = mj_gis.MjProj()
                    lonlatproj.SetFromEPSG(4326)
                    
                    #Get the corners in lonlat
                    llD = mj_gis.ReprojectBounds(layer.BoundsPtL,projection.proj_cs,lonlatproj.proj_cs)

                    session._InsertDefRegion(self.process, layer, query, boundsD[query['regionid']], llD, self.process.proc.overwrite, self.process.proc.delete )

    def _DefaultRegFromVec(self,session):
        for locus in self.process.dstLayerD:
            if self.verbose:
                print ('locus',locus)
            for datum in self.process.dstLayerD[locus]:
                if self.verbose:
                    print ('datum',datum)
                for comp in self.process.dstLayerD[locus][datum]:

                    dstLayer = self.process.dstLayerD[locus][datum][comp]
                    #if not dstLayer._Exists() or self.process.proc.overwrite:
                    srcLayer = self.process.srcLayerD[locus][datum][comp]
                    if not path.isfile(srcLayer.FPN):
                        exitstr = 'setup_process_class: No source layer in _DefaultRegFromVec', srcLayer.FPN
                        exit(exitstr)
                    p = self.process.params
                    fieldL = [p.idcol, p.namecol, p.categorycol, p.parentidcol, p.parentcatcol, p.stratumcol,p.titlecol,p.labelcol]
                    fieldD = mj_gis.GetFeatureAttributeList(srcLayer.FPN, fieldL, p.idcol)
                    if not fieldD:
                        exit('setup_process_class: fieldD failed in _DefaultRegFromVec')
                    for key in fieldD:
                        fieldDD = self._SetfieldD( str(fieldD[key][p.idcol]), str(fieldD[key][p.namecol]), str(fieldD[key][p.categorycol]), int(fieldD[key][p.stratumcol]), str(fieldD[key][p.parentidcol]),str(fieldD[key][p.parentcatcol]) )
                        regionid, regioncat, parentid, regionname, parentcat = fieldD[key][p.idcol],fieldD[key][p.categorycol],fieldD[key][p.parentidcol],fieldD[key][p.namecol],fieldD[key][p.parentcatcol]
                        parentid = str(parentid.lower()).replace(' ', '-')
                        #regionname = str(regionname.lower()).replace(' ', '-')
                        regionid = str(regionid.lower()).replace(' ', '-')
                        if self.verbose:
                            print ('        forcing dst region from to',locus, regionid)

                        dstLayer.locus.locus = regionid
                        dstLayer.locus.path = regionid
                        if self.verbose:
                            print ('        forcing comp band and prefix region to "roi"')
                        #print ('dstLayer.comp.band',dstLayer.comp.band)
                        dstLayer.comp.band = dstLayer.comp.prefix = 'roi'
                        dstLayer.comp._SetCompid()
                        #print ('dstLayer.comp.compid',dstLayer.comp.compid)
                        dstLayer._SetPath()

                        dstLayer.CreateVectorAttributeDef(fieldDD)
                        fieldname = p.idcol
                        valueLL = [[fieldD[key][p.idcol]]]
                        if not dstLayer._Exists() or self.process.proc.overwrite:
                            mj_gis.ExtractFeaturesToNewDS(srcLayer.FPN, dstLayer.FPN,fieldname,valueLL, dstLayer.fieldDefL)
                            fieldname = 'REGIONID'
                            #Get the epsg and bounds
                            boundsD = mj_gis.GetFeatureBounds(dstLayer.FPN,fieldname) 
                            if len(boundsD) != 1:
                                exitstr = 'Default regions must consist on only one (1) feature (polygon or multipolygon): %s' %(dstLayer.FPN)
                                exit(exitstr)
                            projection = mj_gis.GetVectorProjection(dstLayer.FPN)
                            k = list(boundsD)[0]
                            bounds = boundsD[k]
    
                            dstLayer._SetBounds(projection.epsg,boundsD[k][0], boundsD[k][1], boundsD[k][2], boundsD[k][3] )  
                  
                            #Set lonlat projection
                            lonlatproj = mj_gis.MjProj()
                            lonlatproj.SetFromEPSG(4326)
                            #Get the corners in lonlat
                            llD = mj_gis.ReprojectBounds(dstLayer.BoundsPtL,projection.proj_cs,lonlatproj.proj_cs)
    
                            title = label = 'default region %s' %(regionid)
                            query = {'regionname':regionname,'regioncat':regioncat, 'parentid':parentid, 'parentcat':parentcat,'regionid':regionid, 'title':title,'label':label,'epsg':projection.epsg}
                            session._InsertDefRegion(self.process, dstLayer, query, bounds, llD,self.process.overwrite,self.process.delete )
                                                        
                        else:
                            boundsD = mj_gis.GetFeatureBounds(dstLayer.FPN,"REGIONID")
                            key = list(boundsD.keys())[0]
                            bounds = boundsD[key]
                            minx,miny,maxx,maxy = bounds
                            llparams = ['ullon','ullat','urlon','urlat','lrlon','lrlat','lllon','lllat']
                            llvalues = [minx,maxy,maxx,maxy,maxx,miny,minx,miny]
                            llD = dict(zip(llparams,llvalues))
                            title = label = 'default region %s' %(regionid)
                            query = {'regionname':regionname,'regioncat':regioncat, 'parentid':parentid, 'parentcat':parentcat,'regionid':regionid, 'title':title,'label':label,'epsg':4326}

                            session._InsertDefRegion(self.process, dstLayer, query, bounds, llD,self.process.overwrite,self.process.delete )

                            if self.verbose:
                                printstr = '    Layer %s already exists, skipping' %(dstLayer.FPN)
                                print (printstr)

    def _SetfieldD(self,regionid,regionname,regioncat,stratum,parentid,parentcat):
        #TGTODO SHOULD BE FROM DB
        fieldDD = {}
        fieldDD['REGIONID'] = {'name':'REGIONID', 'type':'string','width':32,'precision':0,'transfer':'constant','source':regionid }
        fieldDD['NAME'] = {'name':'NAME', 'type':'string','width':64,'precision':0,'transfer':'constant','source':regionname }
        fieldDD['CATEGORY'] = {'name':'CATEGORY', 'type':'string','width':32,'precision':0,'transfer':'constant','source':regioncat }
        fieldDD['STRATUM'] = {'name':'STRATUM', 'type':'integer','width':4,'precision':0,'transfer':'constant','source':stratum }
        fieldDD['PARENTID'] = {'name':'PARENTID', 'type':'string','width':32,'precision':0,'transfer':'constant','source':parentid }
        fieldDD['PARENTCAT'] = {'name':'PARENTCAT', 'type':'string','width':32,'precision':0,'transfer':'constant','source':parentcat }
        return fieldDD
           
def DegreeTiles():
    '''Setup default region tiles covering 1deg by 1deg, this is like tile system for any data
    '''
    session = SelectUser()
    userid = 'karttur'
    userData = session._SelectUserCreds(userid,'')
    session._Close() 
    if userData[2] < 10:
        exitstr = 'The user does not have sufficient rights for this process'
        exit(exitstr)
      
    headL = ['regioncat','regionid','regionname','parentid','title','label']
    headL =['regioncat','regionid','regiontype','epsg','ullat','ullon','urlat','urlon','lrlat',
            'lrlon','lllat','lllon','monx','miny','maxx','maxy']
    csvL = []
    regionL = []
    #the tile "region" naming is based on the lower left corner
    for lat in range(-90,90):
        for lon in range(-180,180): 
            if lat < 0:  
                latnr = abs(lat)
                latsign = 'S'
                if lon < 0:
                    lonnr = abs(lon)
                    lonsign = 'W'
                    #tileid = 'S%(lat)dW%(lon)d' %{'lat':abs(lat),'lon':abs(lon)}
                else:
                    lonnr = lon
                    lonsign = 'E'
                    #tileid = 'S%(lat)dE%(lon)d' %{'lat':lat,'lon':abs(lon)}
            else:
                latnr = lat
                latsign = 'N'
                if lon < 0:
                    lonnr = abs(lon)
                    lonsign = 'W'
                    #tileid = 'S%(lat)dW%(lon)d' %{'lat':abs(lat),'lon':abs(lon)}
                else:
                    lonnr = lon
                    lonsign = 'E'
                    #tileid = 'S%(lat)dE%(lon)d' %{'lat':lat,'lon':abs(lon)}
                '''
                if lon < 0:
                    tileid = 'N%(lat)dW%(lon)d' %{'lat':lat,'lon':abs(lon)}
                else:
                    tileid = 'N%(lat)dE%(lon)d' %{'lat':lat,'lon':lon}
                '''
            if latnr < 10:
                latstr = '0%(l)d' %{'l':latnr}
            else:
                latstr = '%(l)d' %{'l':latnr}
            if lonnr < 10:
                lonstr = '00%(l)d' %{'l':lonnr}
            elif lonnr < 100:
                lonstr = '0%(l)d' %{'l':lonnr}
            else:
                lonstr = '%(l)d' %{'l':lonnr}
            tileid = '%(latsign)s%(latstr)s%(lonsign)s%(lonstr)s' %{'latsign':latsign, 
                        'latstr':latstr, 'lonsign':lonsign, 'lonstr':lonstr}
            #if lonnr ==  180: 
            #    print (lat,lon,tileid)
            #    queryD = {'regioncat':'global', 'regionid':tileid,'regionname':tileid,'parentid':'global',
            #              'title':'1 deg square global','label': '1 deg square global'}
            csvL.append(['global',tileid,tileid,'global','1deqsquare','1 deg square global']) 
            regionL.append(['global', tileid, 'D', 4326, lat+1.0, lon, lat+1.0, lon+1.0, lat,
                       lon+1.0, lat, lon, lon, lat, lon+1.0, lat+1.0])
            ['regioncat','regionid','regiontype','epsg','ullat','ullon','urlat','urlon','lrlat',
            'lrlon','lllat','lllon','minx','miny','maxx','maxy']
            

    home = path.expanduser("~")      
    tmpFPN = path.join(home,'1degglobaltiles.csv')          
    with open(tmpFPN, 'w') as csvfile:
        wr = csv.writer(csvfile)
        wr.writerows(csvL) 
    regFPN = path.join(home,'1degglobalregions.csv')
    with open(regFPN, 'w') as csvfile:
        wr = csv.writer(csvfile)
        wr.writerows(regionL)           
    #Open the db session for MODIS
    session = ManageRegion()

    #session._LoadBulkDefregions(tmpFPN)
    session._LoadBulkRegions(regFPN)
    session._Close()

                    
def ModisTileCoords():
    '''
    '''
    #Check  the karttur password (local .netrc file
    session = SelectUser()
    userid = 'karttur'
    userData = session._SelectUserCreds(userid,'')
    if userData[2] < 10:
        exitstr = 'The user does not have sufficient rights for this process'
        exit(exitstr)
        
    #Open the db session for MODIS
    session = ManageRegion()
    SINproj = mj_gis.MjProj()
    SINproj.SetFromProj4('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs')
    LatLonproj = mj_gis.MjProj()
    LatLonproj.SetFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
    ptL = []
    for lon in range(360):
        ptL.append((lon-180,90))
    for lat in range(180):
        ptL.append((180,-1*(lat-90)))
    for lon in range(360):
        ptL.append((-1*(lon-180),-90))
    for lat in range(180):
        ptL.append((-180,lat-90))
    #worldgeom = mj_gis.ShapelyPoygon('polygon',ptL) 
    #print ('ptL',ptL)
    #worldgeom = mj_gis.Geometry()
    #worldgeom.PointsToPolygonGeom(ptL)
    #print (worldgeom.shapelyGeom)
    #SNULLEKOLLA
    worldgeom = mj_gis.ShapelyPolyGeom(ptL)
    
    worldgeom.ShapelyToOgrGeom()
    worldgeom.GeoTransform(LatLonproj,SINproj)
    worldgeom.OgrGeomToShapely()
    
    home = path.expanduser("~")
    tarShpFP =  path.join(path.dirname(__file__),'data')

    if not path.exists(tarShpFP):
        makedirs(tarShpFP)
    FN = 'modtiles-multi_karttur_global_epsg6842.shp'
    tarShpFPN = path.join(tarShpFP,FN) 
    fieldDefD = {'type':'string','transfer':'constant','source':'globe','width':8}
    fieldDefL = [mj_gis.FieldDef('name',fieldDefD)]
    mj_gis.CreateESRIPolygonGeom(tarShpFPN, fieldDefL, worldgeom, SINproj.proj_cs, 'globe')
    # Create a shape file for all individual tiles in SIN proj
    FN = 'modtiles-single_karttur_global_epsg6842.shp'
    tarShpFPN = path.join(tarShpFP,FN)
    tarDS,tarLayer = mj_gis.ESRICreateDSLayer(tarShpFPN, SINproj.proj_cs, 'polygon', 'tiles', fieldDefL)
    # Create a shape file for all individual tiles in Geographic coordinates
    FN = 'modtiles_karttur_global_0.shp'
    tarShpFPN = path.join(tarShpFP,FN)
    tarDSLonLat,tarLayerLonLat = mj_gis.ESRICreateDSLayer(tarShpFPN, LatLonproj.proj_cs, 'polygon', 'tiles', fieldDefL)
    
    #create a region with all tiles
    tlen = 20015109.3539999984204769
    tlen /= 18
    #regioncat = 'globe'
    #regionid = 'globe'
    for h in range(36):
        minx = tlen*(18-36)+h*tlen
        maxx = minx+tlen
        for v in range(18):
            maxy = tlen*(9-18)+(18-v)*tlen
            miny = maxy-tlen
            ptL = [(minx,maxy),(maxx,maxy),(maxx,miny),(minx,miny)]
            #TGTODO MOVE TO mj_gis
            tilegeom = mj_gis.ShapelyMultiPointGeom(ptL)
            #convert to ogr
            tilegeom.ShapelyToOgrGeom()
            
            #write target feature 

            tilegeom.GeoTransform(SINproj,LatLonproj)
            tilegeom.OgrGeomToShapely()
            coordL = []
            for point in [ptgeom for ptgeom in tilegeom.shapelyGeom]:
                coordL.extend([list(point.coords)[0][0],list(point.coords)[0][1]])
            ullon, ullat, urlon, urlat, lrlon, lrlat, lllon, lllat = coordL    
            tilepoly = mj_gis.ShapelyPolyGeom([(minx, maxy), (maxx, maxy), (maxx, miny), (minx,miny)])
            #Test if this tile is inside the globe

            if tilepoly.shapelyGeom.intersects(worldgeom.shapelyGeom): 
                if h < 10:
                    htile = 'h0%s' %(h)
                else:
                    htile = 'h%s' %(h)
                if v < 10:
                    vtile = 'v0%s' %(v)
                else:
                    vtile = 'v%s' %(v)
                hvtile = '%s%s' %(htile,vtile)
                polytilegeom = mj_gis.ShapelyPolyGeom(ptL)
                polytilegeom.ShapelyToOgrGeom()
                fieldDefD = {'type':'string','transfer':'constant','source':hvtile,'width':8}
                fieldDefL = [mj_gis.FieldDef('name',fieldDefD)]
                #create target feature
                tarFeat = mj_gis.ogrFeature(tarLayer)
                print (polytilegeom.shapelyGeom)
                tarFeat.CreateOgrFeature(polytilegeom.ogrGeom, fieldDefL) 
                if h == 17:
                    pass
                else:
                    #to be correct 5 points are needed and also the lat must be fitted
                    if h < 18 and ullon > 0:
                        ullon = -180
                        #continue
                    if h < 18 and lllon > 0:
                        lllon = -180
                        #continue
                    if h < 18 and urlon > 0:
                        urlon = -180
                        #continue
                    if h < 18 and lrlon > 0:
                        lrlon = -180
                        #continue
                    if h > 18 and urlon < 0:
                        urlon = 180
                        #continue 
                    if h > 18 and lrlon < 0:
                        lrlon = 180
                        #continue
                    if h > 18 and ullon < 0:
                        ullon = 180
                        #continue 
                    if h > 18 and lllon < 0:
                        lllon = 180
                    if hvtile == 'h24v01':
                        urlon = 180
                    if hvtile == 'h24v16':
                        lrlon = 180  
                    if hvtile == 'h11v01':
                        ullon = -180
                    if hvtile == 'h11v16':
                        lllon = -180 
                        
                if ullon > urlon:
                    print ('ERROR','ullon > urlon',hvtile,ullon,urlon)

                if lllon > lrlon:
                    print ('ERROR','lllon > lrlon',hvtile, lllon, lrlon)

                #
                polytilegeom = mj_gis.ShapelyPolyGeom([(ullon, ullat), (urlon, urlat), (lrlon, lrlat), (lllon,lllat)])
                polytilegeom.ShapelyToOgrGeom()
                #polytilegeom.GeoTransform(SINproj,LatLonproj)
                #create target feature
                tarLonLatFeat = mj_gis.ogrFeature(tarLayerLonLat)
                tarLonLatFeat.CreateOgrFeature(polytilegeom.ogrGeom, fieldDefL)
                west,south,east,north = polytilegeom.shapelyGeom.bounds

                session._InsertModisTileCoord(hvtile,h,v,minx,maxy,maxx,miny,west,south,east,north,ullat,ullon,lrlon,lrlat,urlon,urlat,lllon,lllat)
                #print ('session.name',session.name)
                query = {'system':'system','table':'regions','h':h,'v':v,'hvtile':hvtile,'regionid':'global','regioncat':'global','regiontype':'default','delete':False}
                session._InsertModisRegionTile(query)
    tarDS.CloseDS()
    tarDSLonLat.CloseDS()
    session._Close()
    print ('Check the shaoe file',tarShpFPN)
    return (tarShpFPN)

def SentinelMGRSCoordsOld(srcShpFPN,maskShpFPN):
    '''Open sentinel file and loop all tiles to set the coords
    '''
    srcDS, srcVector, fieldDefL = mj_gis.ESRIOpenGetLayer(srcShpFPN,'read')
    #maskDS, maskLayer, maskDefL = mj_gis.ESRIOpenGetLayer(maskShpFPN,'read')
    featureCount = srcVector.layer.GetFeatureCount()
    print ("Number of features in %s: %d" % (srcShpFPN,featureCount))
    for feature in srcVector.layer:
        print (feature.GetField("Name"))
        geom = feature.GetGeometryRef()
        print (geom.Centroid().ExportToWkt())
    
def GetSentinelTilesDict(session):
    '''
    '''
    recs = session._SelectSentinelTileCoords({})
    senTileD ={}
    for rec in recs:
        epsg,mgrs,utmzone,mgrsid,minx,miny,maxx,maxy,ullat,ullon,lrlat,lrlon,urlat,urlon,lllat,lllon = rec
        llptL = ((ullon,ullat),(urlon,urlat),(lrlon,lrlat),(lllon,lllat))
        sentilegeom = mj_gis.Geometry()
        sentilegeom.PointsToPolygonGeom(llptL)
        west, south, east, north = sentilegeom.shapelyGeom.bounds
        senTileD[mgrs] = {'mgrs':mgrs,'mgrsid':mgrsid,'utmzone':utmzone,'geom':sentilegeom,
                              'west':west,'south':south,'east':east,'north':north}
    return senTileD
    
def LinkSentineModisTiles():
    session = SelectUser()
    userid = 'karttur'
    userData = session._SelectUserCreds(userid,'')
    if userData[2] < 10:
        exitstr = 'The user does not have sufficient rights for this process'
        exit(exitstr)    
    #Open the db session for MODIS
    sessionModis = ManageMODIS()
    #Open the session for sentinel
    
    recs = sessionModis._SelectModisTileCoords({})
    sessionModis._Close()
    sessionSentinel = ManageSentinel()
    sentileD = GetSentinelTilesDict(sessionSentinel)
    for rec in recs:
        print ('rec',rec)
        hvtile,h,v,minxsin,minysin,maxxsin,maxysin,ullat,ullon,lrlat,lrlon,urlat,urlon,lllat,lllon = rec
        llptL = ((ullon,ullat),(urlon,urlat),(lrlon,lrlat),(lllon,lllat))
        modtilegeom = mj_gis.Geometry()
        modtilegeom.PointsToPolygonGeom(llptL)
        west, south, east, north = modtilegeom.shapelyGeom.bounds
        #Get all sentinel tiles falling inside the bounding box (nothing is omitted?)
        sentiles = sessionSentinel._SearchTilesFromWSEN(west, south, east, north)
        for sentile in sentiles:
            mgrs,tilewest,tilesouth,tileeast,tilenorth,ullon,ullat,urlon,urlat,lrlon,lrlat,lllon,lllat, minx, miny, maxx, maxy = sentile  
            llptL = ((ullon,ullat),(urlon,urlat),(lrlon,lrlat),(lllon,lllat))
            sentilegeom = mj_gis.Geometry()
            sentilegeom.PointsToPolygonGeom(llptL)

            if modtilegeom.shapelyGeom.intersects(sentilegeom.shapelyGeom):
                #print ('    overlap',productoverlap)
                #query = {'system':'system', 'regionid':regionid,'regiontype':'default', 'overwrite':False, 'delete':False, 'mgrs':mgrs,'utmzone':self.senTileD[mgrs]['utmzone'], 'mgrsid':self.senTileD[mgrs]['mgrsid']}
                #self.session._InsertSentinelRegionTile(query)
                print ()
                queryD = {'system':'system', 'regionid':hvtile, 'regiontype':'sintile', 'mgrs':mgrs, 'utmzone':sentileD[mgrs]['utmzone'], 'mgrsid':sentileD[mgrs]['mgrsid'], 'overwrite':False, 'delete':False}
                print ('in')
                #query = {'system':'system', 'regionid':regionid,'regiontype':'default', 'overwrite':False, 'delete':False, 'mgrs':mgrs,'utmzone':self.senTileD[mgrs]['utmzone'], 'mgrsid':self.senTileD[mgrs]['mgrsid']}
                sessionSentinel._InsertSentinelMODISTile(queryD)
                #sessionSentinel._InsertSingleSentinelRegion(queryD)
                print ('out')
    sessionSentinel._Close()

def LinkLandsatModisTiles():
    session = SelectUser()
    userid = 'karttur'
    userData = session._SelectUserCreds(userid,'')
    if userData[2] < 10:
        exitstr = 'The user does not have sufficient rights for this process'
        exit(exitstr)    
    #Open the db session for MODIS
    sessionModis = ManageMODIS()
    #Open the session for sentinel
    sessionLandsat = ManageLandsat()
    recs = sessionModis._SelectModisTileCoords({})
    for rec in recs:
        hvtile,h,v,minxsin,minysin,maxxsin,maxysin,ullat,ullon,lrlat,lrlon,urlat,urlon,lllat,lllon = rec
        llptL = ((ullon,ullat),(urlon,urlat),(lrlon,lrlat),(lllon,lllat))
        modtilegeom = mj_gis.Geometry()
        modtilegeom.PointsToPolygonGeom(llptL)
        west, south, east, north = modtilegeom.shapelyGeom.bounds
        #Get all sentinel tiles falling inside the bounding box (nothing is omitted?)
        landsattiles = sessionLandsat._SearchTilesFromWSEN(west, south, east, north)
        #loop
        for landsattile in landsattiles:
            wrs,ldir,p,r,tilewest,tilesouth,tileeast,tilenorth,ullon,ullat,urlon,urlat,lrlon,lrlat,lllon,lllat, minx, miny, maxx, maxy = landsattile  
            llptL = ((ullon,ullat),(urlon,urlat),(lrlon,lrlat),(lllon,lllat))
            tilegeom = mj_gis.Geometry()
            tilegeom.PointsToPolygonGeom(llptL)
            prD = ConvertLandsatTilesToStr(p,r)
            overlapGeom = tilegeom.ShapelyIntersection(modtilegeom)
            productoverlap = overlapGeom.area/tilegeom.shapelyGeom.area
            if productoverlap > 0:    
                queryD = {'regionid':hvtile, 'regiontype':'sintile', 'prstr':prD['prstr'], 'dir':ldir, 'wrs':wrs, 'path':p, 'row':r}
                sessionLandsat._InsertSingleLandsatRegion(queryD)
                       
def LinkSentinelLandsatTiles():
    session = SelectUser()
    userid = 'karttur'
    userData = session._SelectUserCreds(userid,'')
    if userData[2] < 10:
        exitstr = 'The user does not have sufficient rights for this process'
        exit(exitstr)    
    #Open the db session for MODIS
    sessionLandsat = ManageLandsat()
    #Open the session for sentinel
    sessionSentinel = ManageSentinel()
    #Select all landsat tiles
    recs = sessionLandsat._SelectLandsatTileCoords({})
    doneL = sessionSentinel._SelectSentineRegions({'regiontype':'lsattile'})

    for rec in recs:
        wrs,dir,wrspath,wrsrow,minx,miny,maxx,maxy,ullat,ullon,lrlat,lrlon,urlat,urlon,lllat,lllon = rec
        prD = ConvertLandsatTilesToStr(wrspath,wrsrow)
        if prD['prstr'] in doneL:
            continue
        print ('adding sentinel region:',prD['prstr'])
        llptL = ((ullon,ullat),(urlon,urlat),(lrlon,lrlat),(lllon,lllat))
        landsattilegeom = mj_gis.Geometry()
        landsattilegeom.PointsToPolygonGeom(llptL)
        west, south, east, north = landsattilegeom.shapelyGeom.bounds
        #Get all sentinel tiles falling inside the bounding box (nothing is omitted?)
        sentiles = sessionSentinel._SearchTilesFromWSEN(west, south, east, north)
        print ('adding sentinel region:',prD['prstr'],len(sentiles))
        for sentile in sentiles:
            mgrs,tilewest,tilesouth,tileeast,tilenorth,ullon,ullat,urlon,urlat,lrlon,lrlat,lllon,lllat, minx, miny, maxx, maxy = sentile  
            llptL = ((ullon,ullat),(urlon,urlat),(lrlon,lrlat),(lllon,lllat))
            tilegeom = mj_gis.Geometry()
            tilegeom.PointsToPolygonGeom(llptL)
            #Get the overlap
            overlapGeom = tilegeom.ShapelyIntersection(landsattilegeom)
            productoverlap = overlapGeom.area/tilegeom.shapelyGeom.area
            if productoverlap > 0:
                #Replace this with a COPY command, faster and assures that a tile is always complete
                queryD = {'regionid':prD['prstr'], 'regiontype':'lsattile', 'mgrs':mgrs, 'utm':int(mgrs[0:2]), 'mgrsid':mgrs[2:5]}
                sessionSentinel._InsertSingleSentinelRegion(queryD)
 
def Setup(relpath,projFN,verbose):
    '''
    Setup processes
    '''
    srcFP = path.join(path.dirname(__file__),relpath)
    projFPN = path.join(srcFP,projFN)
    #Read the textfile with links to the xml files defining schemas and tables
    procLL = runxml.ReadXMLProcesses(projFPN,verbose) 

    for procL in procLL:
        for proc in procL:
            #At this stage the proc (process) is a well defined and cheked dictionary  
            #change to local path
            if hasattr(proc, 'srcpathD'):
                if proc.srcpathD['volume'] == 'None':
                            proc.srcpathD['volume'] = path.dirname(__file__)
            if hasattr(proc, 'dstpathD'):
                if proc.dstpathD['volume'] == 'None':
                            proc.dstpathD['volume'] = path.dirname(__file__)

            if verbose: 
                print ('    PROCESS', proc.processid, proc.rootprocid)
            if proc.rootprocid == 'manageprocess':
                #Mangeprocess does not contain composition definitions and timestep settings and can be run without further settings
                #Connect to the Postgres Server for managing processes
                session = ManageProcess()
                ProcessProcess(session,proc)
            else:

                if proc.rootprocid == 'ManageRegion':
                    session = ManageRegion()
                    #proj, proc, session, verbose
                    process = MainProc(proc, session, verbose)
                    #Connect to the Postgres Server for managing regions
                    ProcessDefaultRegions(process,session,verbose)  
                                  
                elif proc.rootprocid == 'Ancillary':
                    #Connect to the Postgres Server for managing regions
                    session = ManageAncillary()

                    

                    process = MainProc(proc,session,verbose)                  
                    ProcessAncillary(process,session,verbose)
                elif proc.rootprocid == 'MODISProc':
                    #Connect to the Postgres Server for managing regions
                    session = ManageMODIS()

                    process = MainProc(proc,session,verbose)                  
                    ProcessModis(process,session,verbose)
                elif proc.rootprocid == 'SentinelProcess':
                    #Connect to the Postgres Server for managing regions
                    session = ManageSentinel()
                    process = MainProc(proc,session,verbose)                  
                    ProcessSentinel(process,session,verbose)

                else:
                    exitstr = 'Unrecognized root process %(r)s' %{'r':proc.rootprocid}
                    exit(exitstr)
            #Close server connection
            
            session._Close()

if __name__ == "__main__":
    '''
    '''
    DegreeTiles()
    TETST
    DefaultRegions = True
    MODIS = False
    Landsat = False
    Sentinel = False
    Ancillary = False
    Climate = False
    verbose = True
    
    '''Link to project file that sets up default regions, arbitrary regions and special regions. 
    '''
    if DefaultRegions:
        projFN = 'regions_karttur_setup_20181116.txt'
        Setup('regiondoc',projFN,verbose)

    if MODIS:
        '''Stand alone script that defines the MODIS tile coordinates'''
        
        FPN = ModisTileCoords()
        '''
        exitstr = 'The script ModisTileCoords() produced a shape with all MODIS SIN tiles projected to Geographic coordiantes.\n \
            Copy the shape data sourse: %(fpn)s,\n and edit the xml file for importing this layer.\n \
            Then comment out the "exit" command and re-run the module.' %{'fpn':FPN}
        exit(exitstr)
        '''
        projFN = 'modis_karttur_setup_20181116.txt'
        Setup('modisdoc',projFN,verbose)
    
    if Sentinel:
        '''Link to project file that sets up the Sentinel tiling system'''
        projFN = 'sentinel_karttur_setup_2018116.txt'
        Setup('sentineldoc',projFN,verbose)
        
    if Landsat:
        pass
    
    if MODIS and Sentinel:
        '''Stand alone script that links sentinel and modis, requires that all sentinel tiles are in the db'''
        LinkSentineModisTiles()

    if MODIS and Landsat:
        LinkLandsatModisTiles()

    if Sentinel and Landsat:
        LinkSentinelLandsatTiles()

    if Ancillary:
        ''' link to project file that imports default ancillary data'''
        projFN = 'ancillary_karttur_setup_20180221_0.txt'
        Setup('ancildoc',projFN,verbose)
    
    if Climate:
        ''' Climate data'''
        projFN = 'climate_karttur_setup_20181116.txt'
        Setup('climatedoc',projFN,verbose)
    