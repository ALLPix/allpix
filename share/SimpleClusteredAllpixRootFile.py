from ROOT import *
from math import fsum,fabs
from array import array
from scipy.cluster.hierarchy import fclusterdata
import sys

pitchX = 0.25
pitchY = 0.05
npix_X = 80
npix_Y = 336

halfChip_X = npix_X*pitchX/2.
halfChip_Y = npix_Y*pitchY/2.

# sigma = 0.015
um = 1e-3
mm = 1
cm =10


someData_in_um = 1000*um

scaler =1


class Cluster:

    col = []
    row = []
    tot = []
    energyGC = []
    energyPbPC = []
    sizeX = 0
    sizeY = 0
    size  = 0
    totalTOT =0
    aspectRatio = 0

    # local coordinates
    relX = 0.
    relY = 0.
    relZ = 0.
    relX_energyGC = 0.
    relY_energyGC = 0.
    relZ_energyGC = 0.
    relX_energyPbPC = 0.
    relY_energyPbPC = 0.
    relZ_energyPbPC = 0.

    # telescope coordinates
    absX =-10000.
    absY =-10000.
    absZ =-10000.
    absX_energyGC =-10000.
    absY_energyGC =-10000.
    absZ_energyGC =-10000.
    absX_energyPbPC =-10000.
    absY_energyPbPC =-10000.
    absZ_energyPbPC =-10000.

    resX = -10000.
    resY = -10000.
    resX_energyGC = -10000.
    resY_energyGC = -10000.
    resX_energyPbPC = -10000.
    resY_energyPbPC = -10000.

    id = 0

    # track number
    tracknum = -1

    def __init__(self):
        self.sizeX = 0
        self.sizeY = 0
        self.size  = 0
        self.col = []
        self.row = []
        self.tot = []
        self.energyGC = []
        self.energyPbPC = []

    def addPixel(self,col,row,tot):
        self.col.append(col)
        self.row.append(row)
        self.tot.append(tot)
        #self.energyGC.append(energyGC)
        #self.energyPbPC.append(energyPbPC)

    def Print(self):
        for i in range(len(self.col)):
            print "x:%d y%d tot:%.3f"%(self.col[i],self.row[i],self.tot[i])
        print "Cluster Total size = %d , in X = %d Y = %d , Aspect Ratio = %.3f , Total Energy (keV) = %.1f ID = %i"%(self.size,self.sizeX,self.sizeY,self.aspectRatio,self.totalTOT,self.id)
        print "Position in sensor X = %.3f Y = %.3f"%(self.absX,self.absY)

    def Statistics(self) :
        self.totalTOT=fsum(self.tot)
        self.size=len(self.col)
        self.sizeX=max(self.col)-min(self.col)+1
        self.sizeY=max(self.row)-min(self.row)+1
        self.aspectRatio=float(self.sizeY)/self.sizeX
        #self.totalEnergyGC = fsum(self.energyGC)
        #self.totalEnergyPbPC = fsum(self.energyPbPC)

#
#compute the hit position as the mean of the fired pixels positions weighted by their deposited energy
#
    def GetQWeightedCentroid(self) :
        self.relX=0.
        self.relY=0.
        self.relX_energyGC=0.
        self.relY_energyGC=0.
        self.relX_energyPbPC=0.
        self.relY_energyPbPC=0.

        if(self.totalTOT==0):
            self.GetDigitalCentroid()
        else :
            for index,tot_tmp in enumerate(self.tot) :
                self.relX+=(self.col[index]*pitchX)*tot_tmp
                self.relY+=(self.row[index]*pitchY)*tot_tmp

            self.relX = self.relX/self.totalTOT + pitchX/2.
            self.relY = self.relY/self.totalTOT + pitchY/2.

            self.absX=self.relX - halfChip_X
            self.absY=self.relY - halfChip_Y
            self.absZ=0


#
#compute the hit position as the mean of the fired pixels positions (digital method)
#
    def GetDigitalCentroid(self) :
        self.relX=0.
        self.relY=0.
        for index,col_tmp in enumerate(self.col) :
            self.relX+=(self.col[index]*pitchX)
            self.relY+=(self.row[index]*pitchY)

        self.relX = self.relX/len(self.col) + pitchX/2.
        self.relY = self.relY/len(self.row) + pitchY/2.

        self.absX=self.relX - halfChip_X
        self.absY=self.relY - halfChip_Y
        self.absZ=0

#
#compute the hit position as the center of the fired pixel with the highest energy
#
    def GetMaxTOTCentroid(self) :
        maxTOTindex_tmp=0
        maxTOT_tmp=self.tot[0]
        for index,tot_tmp in enumerate(self.tot) :
            if self.tot[index]>maxTOT_tmp:
                maxTOT_tmp=self.tot[index]
                maxTOTindex_tmp=index

        self.relX=self.col[maxTOTindex_tmp]*pitchX + pitchX/2.
        self.relY=self.row[maxTOTindex_tmp]*pitchY + pitchY/2.

        self.absX=self.relX - halfChip_X
        self.absY=self.relY - halfChip_Y
        self.absZ=0

#
#compute the hit position as the Qweighted method but adding an eta correction du to the charge sharing between the fired pixels
#
    def GetEtaCorrectedQWeightedCentroid(self,sigma=0.003,sigmaGC=0.003,sigmaPbPC=0.003) :

        self.relX = -1000
        self.relY = -1000
        self.absX = -1000
        self.absY = -1000



        if(self.size==2) :
            if(self.sizeX==2 and self.sizeY==1) :
                # cluster size 2x1
                Qrel = self.tot[self.col.index(min(self.col))] / self.totalTOT
                self.relX = max(self.col)*pitchX - shiftLat(sigma,Qrel)
                self.relY = self.row[0]*pitchY + pitchY/2.


            elif(self.sizeX==1 and self.sizeY==2) :
                # cluster size 1x2
                Qrel = self.tot[self.row.index(min(self.row))] / self.totalTOT
                self.relX = self.col[0]*pitchX + pitchX/2.
                self.relY = max(self.row)*pitchY - shiftLat(sigma,Qrel)
                

            elif(self.sizeX==2 and self.sizeY==2) :
                # cluster size 2 with sizeX = 2 and sizeY = 2 i.e. 2 pixels on a diagonal
                self.GetMaxTOTCentroid()


        elif(self.size==4) :
            if(self.sizeX==2 and self.sizeY==2) :
                for i in xrange(len(self.row)):
                    if self.row[i] == min(self.row) and self.col[i] == min(self.col):
                        bottomlefti = i
                    if self.row[i] == max(self.row) and self.col[i] == max(self.col):
                        toprighti = i
                    if self.row[i] == min(self.row) and self.col[i] == max(self.col):
                        toplefti = i
                    if self.row[i] == max(self.row) and self.col[i] == min(self.col):
                        bottomrighti = i
                Qrel1 = self.tot[bottomlefti] / (self.tot[bottomlefti] + self.tot[toprighti])
                Qrel2 = self.tot[bottomrighti] / (self.tot[bottomrighti] + self.tot[toplefti])

                shift1X,shift1Y = shiftDiag(sigma,Qrel1)
                shift2X,shift2Y = shiftDiag(sigma,Qrel2)

                self.relX = max(self.col)*pitchX - shift1X - shift2X
                self.relY = max(self.row)*pitchY - shift1Y + shift2Y

            else : # not 2x2 -> using the simple Qweighted centroid
                self.GetQWeightedCentroid()


        elif(self.size==3) :
            if(self.sizeX==2 and self.sizeY==2) :
                # copy original row, col, tot to keep safe
                orig_row = list(self.row)
                orig_col = list(self.col)
                orig_tot = list(self.tot)
                orig_energyGC = list(self.energyGC)
                orig_energyPbPC = list(self.energyPbPC)

                # calculate and add missing pixel
                for i in xrange(len(self.row)):
                    if self.row.count(self.row[i]) == 1:
                        self.row.append(self.row[i])
                    if self.col.count(self.col[i]) == 1:
                        self.col.append(self.col[i])
                        self.tot.append(1.0)
                        self.energyGC.append(1.0)
                        self.energyPbPC.append(1.0)

                # proceed as 4 cluster case
                for i in xrange(len(self.row)):
                    if self.row[i] == min(self.row) and self.col[i] == min(self.col):
                        bottomlefti = i
                    if self.row[i] == max(self.row) and self.col[i] == max(self.col):
                        toprighti = i
                    if self.row[i] == min(self.row) and self.col[i] == max(self.col):
                        toplefti = i
                    if self.row[i] == max(self.row) and self.col[i] == min(self.col):
                        bottomrighti = i
                Qrel1 = self.tot[bottomlefti] / (self.tot[bottomlefti] + self.tot[toprighti])
                Qrel2 = self.tot[bottomrighti] / (self.tot[bottomrighti] + self.tot[toplefti])

                shift1X,shift1Y = shiftDiag(sigma,Qrel1)
                shift2X,shift2Y = shiftDiag(sigma,Qrel2)

                self.relX = max(self.col)*pitchX - shift1X - shift2X
                self.relY = max(self.row)*pitchY - shift1Y + shift2Y

                # put back original row, col, tot
                self.row = orig_row
                self.col = orig_col
                self.tot = orig_tot

            else : # not 2x2 -> using the simple Qweighted centroid
                self.GetQWeightedCentroid()


        else : # other cluster sizes -> using the simple Qweighted centroid
            self.GetQWeightedCentroid()


        self.absX = self.relX - halfChip_X
        self.absY = self.relY - halfChip_Y
        self.absZ = 0


        if (self.relX == -1000 or self.relY == -1000 or self.absX == -1000 or self.absY == -1000):
            print "WARNING GetEtaCorrectedQWeightedCentroid didn't calculate centroid for some cluster"
            print "WARNING This should never happen - review code"



    def GetResiduals(self,x,y) :
        self.resX = self.absX-(x)
        self.resY = self.absY-(y)



    def GetPixelResiduals(self,trackx,tracky) :
        # compute the x, y distances between a track 
        # and the centre of each pixel in the cluster
        # return the smallest combined

        dr = []
        dx = []
        dy = []
        for i in xrange(self.size):
            resX = self.col[i]*pitchX + pitchX/2. - halfChip_X - trackx
            resY = self.row[i]*pitchY + pitchY/2. - halfChip_Y - tracky

            dr.append(sqrt(resX**2 + resY**2))
            dx.append(resX)
            dy.append(resY)

        return min(dr), dx[dr.index(min(dr))], dy[dr.index(min(dr))]




class AllpixData:
    
    m_RootFile = 0
    m_RawDataTree = 0
    
    def __init__(self,filename):
        
        self.m_RootFile=TFile(filename)
        self.m_RawDataTree = self.m_RootFile.Get("tree")
    
    
    
    
    def GetEvent(self,i):
         
        self.m_RawDataTree.GetEntry(i)
        
#         for data in self.m_RawDataTree.posX : 
#             print data
#         for data in self.m_RawDataTree.posY : 
#             print data           
#             
#         for data in self.m_RawDataTree.energyTotal : 
#             print data
#         for data in self.m_RawDataTree.TOT : 
#             print data        

        
    
    
    
    
    
    def ClusterEvent(self,i,method="QWeighted",sigma=0.003,sigmaGC=0.003,sigmaPbPC=0.003):

        self.GetEvent(i)

        row_tmp = [s for s in self.m_RawDataTree.posX]
        col_tmp = [s for s in self.m_RawDataTree.posY]
        tot_tmp = [s for s in self.m_RawDataTree.TOT]
        
        print row_tmp
        print col_tmp
        print tot_tmp
        
        #energyGC_tmp = [s for s in self.p_energyGC]
        #energyPbPC_tmp = [s for s in self.p_energyPbPC]

#         # ------------------------------------------------------------------------------------#
#         # Temporary solution for pixels hit several times. Include TOA in the future analysis
#         # ------------------------------------------------------------------------------------#
#         if (SensorType=="Timepix3" or SensorType=="CLICpix"):
#             indexPixelsToRemove=[]
#             for index in range(0, len(row_tmp)):
#                 row_temp=row_tmp[index]
#                 col_temp=col_tmp[index]            
#                 for index2 in range(index+1, len(row_tmp)):
#                     if(row_temp==row_tmp[index2] and col_temp==col_tmp[index2]):
#                         indexPixelsToRemove.append(index)
#                         indexPixelsToRemove.append(index2)
# 
#             row_tmp=[ row_tmp[k] for k in range(0, len(row_tmp)) if k not in indexPixelsToRemove ]
#             col_tmp=[ col_tmp[k] for k in range(0, len(col_tmp)) if k not in indexPixelsToRemove ]   
#             tot_tmp=[ tot_tmp[k] for k in range(0, len(tot_tmp)) if k not in indexPixelsToRemove ]
#             energyGC_tmp=[ energyGC_tmp[k] for k in range(0,len(energyGC_tmp)) if k not in indexPixelsToRemove ]
#             energyPbPC_tmp=[ energyPbPC_tmp[k] for k in range(0,len(energyPbPC_tmp)) if k not in indexPixelsToRemove ]
        # ------------------------------------------------------------------------------------#


        # set a maximum number of hit pixels to be clustered (skips large events)
        if len(col_tmp) < 5000:
#             try : 
            print "clustering event %i"%i
            clusters = self.RecursiveClustering(col_tmp,row_tmp,tot_tmp)
            print "%i clusters found"%len(clusters)
#             except : 
#                 print "except"
#                 clusters=[]
        else:
            print "Event", i, "not beng clustered,", len(col_tmp), "hit pixels"
            clusters=[]


        for cluster in clusters :
            cluster.Statistics()    
    
        clusters = [cluster for cluster in clusters if cluster.totalTOT>0]
        clusterid = 0

        for cluster in clusters :
            if (method=="QWeighted"):
                cluster.GetQWeightedCentroid()
            elif (method=="DigitalCentroid"):
                cluster.GetDigitalCentroid()
            elif (method=="maxTOT"):
                cluster.GetMaxTOTCentroid()
            elif (method=="EtaCorrection"):
                cluster.GetEtaCorrectedQWeightedCentroid(sigma,sigmaGC,sigmaPbPC)

            cluster.id=clusterid
            clusterid+=1
            cluster=0

    
        return clusters


    def SciPyClustering(self,col,row,tot):

        print "inside scipy"
        pixels = [[col[i],row[i]] for i,x in enumerate(col)]
        if(len(pixels)>1):
            result=fclusterdata(pixels,sqrt(2.),criterion="distance")
            clusters=[Cluster() for i in range(max(result))]
            [clusters[x-1].addPixel(col[j],row[j],tot[j]) for j,x in enumerate(result)]
        else:
            if(len(pixels)==1):
                c=Cluster()
                c.addPixel(col[0],row[0],tot[0])
                clusters=[c]        
    
        print len(clusters)
        return clusters

    def RecursiveClustering(self,row,col,tot) :
        clusters = []
        while(len(row)!=0) :

            cluster = Cluster()
            cluster.addPixel(col[0], row[0], tot[0])
            print "[DEBUG] adding pixel col=%d row=%d as seed"%(col[0],row[0])
            row.pop(0)
            col.pop(0)
            tot.pop(0)
            while(self.addNeighbor(cluster, col,row, tot)>0):
                pass
            clusters.append(cluster)
        return clusters
    
    
    def addNeighbor(self,cluster,col,row,tot):
        counter =0
        i=0
        j=0
        len_col=len(col)
        len_clu_col=len(cluster.col)
        while(i<len_col):
            j=0
            while(j<len_clu_col):

                if((col[i]-cluster.col[j])**2>1) :
                    j+=1
                    continue

                if((row[i]-cluster.row[j])**2>1) :
                    j+=1
                    continue

                cluster.addPixel(col[i],row[i],tot[i])

            #print "[DEBUG] after adding pixel col=%d row=%d to existing cluster as neighbor to x=%d y=%d "%(col[i],row[i],cluster.col[j],cluster.row[j])

                col.pop(i)
                row.pop(i)
                tot.pop(i)
                counter+=1
                i+=-1
                len_col=len(col)
                len_clu_col=len(cluster.col)
                break
            i+=1
        return counter   
    
    def DumpClusterTree(self,outfilename):

        outfile = TFile(outfilename,"RECREATE")
        outfile.cd()

        clusterTree = TTree('clusters','cluster tree')

        maxn = 500
        col = array( 'i', maxn*[ 0 ] )
        row = array( 'i', maxn*[ 0 ] )
        tot = array( 'f', maxn*[ 0. ] )

        event = array( 'i', [ 0 ] )
        sizeX = array( 'i', [ 0 ] )
        sizeY = array( 'i', [ 0 ] )
        size  = array( 'i', [ 0 ] )
        totalTOT =array( 'f', [ 0. ] )

        aspectRatio = array( 'f', [ 0. ] )
        relX = array( 'f', [ 0. ] )
        relY = array( 'f', [ 0. ] )
        absX = array( 'f', [ 0. ] )
        absY = array( 'f', [ 0. ] )
        resX = array( 'f', [ 0. ] )
        resY = array( 'f', [ 0. ] )
        id = array( 'f', [ 0. ] )

        clusterTree.Branch( 'event', event, 'event/I' )
        clusterTree.Branch( 'size', size, 'size/I' )
        clusterTree.Branch( 'sizeX', sizeX, 'sizeX/I' )
        clusterTree.Branch( 'sizeY', sizeY, 'sizeY/I' )
        clusterTree.Branch( 'totalTOT', totalTOT, 'totalTOT/F' )

        clusterTree.Branch( 'aspectRatio', aspectRatio, 'aspectRatio/F' )
        clusterTree.Branch( 'relX', relX, 'relX/F' )
        clusterTree.Branch( 'relY', relY, 'relY/F' )
        clusterTree.Branch( 'absX', absX, 'absX/F' )
        clusterTree.Branch( 'absY', absY, 'absY/F' )
        clusterTree.Branch( 'resX', resX, 'resX/F' )
        clusterTree.Branch( 'resY', resY, 'resY/F' )

        clusterTree.Branch( 'id', id, 'id/I' )
        clusterTree.Branch( 'col', col, 'col[size]/I' )
        clusterTree.Branch( 'row', row, 'row[size]/I' )
        clusterTree.Branch( 'tot', tot, 'tot[size]/F' )



        for j in range(self.m_RawDataTree.GetEntries()) :
            clusters = self.ClusterEvent(j)
            for cluster in clusters :

                if(len(cluster.col)<maxn) :
        
                    for i in range(len(cluster.col)) :
                        col[i]=cluster.col[i]
                    for i in range(len(cluster.row)) :
                        row[i]=cluster.row[i]
                    for i in range(len(cluster.tot)) :
                        tot[i]=cluster.tot[i]
            
                else : 
                    for i in range(maxn) :
                        col[i]=cluster.col[i]
                    for i in range(maxn) :
                        row[i]=cluster.row[i]
                    for i in range(maxn) :
                        tot[i]=cluster.tot[i]

                    

                sizeX[0]=cluster.sizeX
                sizeY[0]=cluster.sizeY
                size[0]=cluster.size
                totalTOT[0]=cluster.totalTOT

                aspectRatio[0]=cluster.aspectRatio
                relX[0]=cluster.relX
                relY[0]=cluster.relY
                resX[0]=cluster.resX
                resY[0]=cluster.resY
                absX[0]=cluster.absX
                absY[0]=cluster.absY

                id[0]=cluster.id
                event[0]=j
 
                clusterTree.Fill()
        outfile.Write()
    
  
#########################################################################################
def usage():
    print 'Converts allpix generated ROOT files to a clustered ROOT file (pyEudetNTuple format) https://github.com/pyEudetAnalysis/pyEudetAnalysis'
    print 'Usage:\n python %s <input ROOT file> <output ROOT file>' % ( sys.argv[0] )

#########################################################################################  
  
    
if __name__ == "__main__":
    if len( sys.argv ) < 3:
        usage()
        sys.exit( 1 )    
        
    filename,outfile = sys.arv[1],sys.argv[2]    
    data = AllpixData(filename)
    
#     for i in range(10): 
#         clusters = data.ClusterEvent(i)
#         for cluster in clusters : 
#             cluster.Print()

    data.DumpClusterTree(outfile)        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
