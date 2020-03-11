'''
Created on 21 feb. 2018

@author: thomasgumbricht
'''

from geoimagine.postgresdb import ManageProcess, ManageLayout, ManageRegion
from geoimagine.layout import ProcessLayout
#from geoimagine.kartturmain import UserProj, SetXMLProcess, MainProc
from geoimagine.kartturmain import runxml, MainProc

from os import path

class ProcessProcess:  
    """"class for processes defining other processes"""  
    def __init__(self, session, process): 
        """"The constroctur requires an instance of the main process, and the xml elements in the tags used for defining the particular process to run"""  
        self.process = process
        self.session = session
        if self.process.processid == 'addrootproc':
            self.session._ManageRootProcess(self.process)
        elif self.process.processid == 'addsubproc':
            self.session._ManageSubProcess(self.process)
        else:
            exitstr = 'subprocess %s not defined in manageprocess' %(self.process.processid)
            exit( exitstr )
                
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
            #At this stage the proc (process) is a well defined and checked dictionary  
            if verbose: 
                print ('    PROCESS', proc.processid, proc.rootprocid)
            if proc.rootprocid == 'manageprocess':
                #Mangeprocess does not contain composition definitions and timestep settings and can be run without further settings
                #Connect to the Postgres Server for managing processes
                session = ManageProcess()
                ProcessProcess(session,proc)
            elif proc.rootprocid == 'LayoutProc':
                #Mangeprocess does not contain composition definitions and timestep settings and can be run without further settings
                #Connect to the Postgres Server for managing processes
                session = ManageLayout()
                process = MainProc(proc,session,verbose)
                
                ProcessLayout(process,session,verbose)
            else:
                exitstr = 'Unrecognized root process %(r)s' %{'r':proc.rootprocid}
                exit(exitstr)
            #Close server connection
            session._Close()

if __name__ == "__main__":
    verbose = True
    ''' Link to project file that sets up all processes'''
    projFN = 'process_karttur_setup_20181116.txt'
    Setup('dbdoc',projFN,verbose)
