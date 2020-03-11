'''
Created on 21 feb. 2018

@author: thomasgumbricht
'''
#Import the session connector
#from geoimagine.postgresdb import PGsession

from sys import exit
import geoimagine.gis.gis as mj_gis
from os import path
#from geoimagine.gis.gis import GetVectorProjection
#import subprocess
#from copy import deepcopy


class ProcessDefaultRegions:
    'class for all region management'   
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
        #query = dict((x,y)  for x,y in self.process.proc.parameters.__dict__.items())
        for locus in self.process.dstLayerD:
            for datum in self.process.dstLayerD[locus]:
                for comp in self.process.dstLayerD[locus][datum]:
                    query = self.process.proc.paramsD
                    fieldDD = self._SetfieldD( query['regionid'], query['regionname'], query['regioncat'], query['stratum'], query['parentid'], query['parentcat'])
                    layer = self.process.dstLayerD[locus][datum][comp]

                    #The destination region must be forced
                    if self.verbose:
                        print ('        forcing region from to',locus, self.process.proc.paramsD['regionid'])
                    layer.locus = self.process.proc.paramsD['regionid']
                    layer.locuspath = self.process.proc.paramsD['regionid']
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

                    session._InsertDefRegion(self.process, layer, query, boundsD[query['regionid']], llD )
   
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
                        dstLayer.locus = regionid
                        dstLayer.locuspath = regionid
                        if self.verbose:
                            print ('        forcing comp band and prefix region to "roi"')
                        dstLayer.comp.band = dstLayer.comp.prefix = 'roi'
                        dstLayer.comp._SetCompid()

                        dstLayer._SetPath()

                        dstLayer.CreateVectorAttributeDef(fieldDD)
                        fieldname = p.idcol
                        valueLL = [[fieldD[key][p.idcol]]]
                        if not dstLayer._Exists() or self.process.proc.overwrite: #or overwrite
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
                            session._InsertDefRegion(self.process, dstLayer, query, bounds, llD )
                                                        
                        else:
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

