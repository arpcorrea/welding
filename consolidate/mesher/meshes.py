# -*- coding: utf-8 -*-
import numpy as np

class MeshTwoPlates():

    def __init__(self, deck,geometry):
        self.deck = deck
        self.geometry=geometry
        self.set_mesh_grid() 
        self.set_dx()
        self.set_dy()
        self.set_temperatures()
        self.set_conductivity()
        self.set_density()
        self.set_specific_heat()
        self.set_diffusivity()
        self.set_viscosity()   
        self.set_dic()
        self.set_heat_density()
        self.dist()
        
        
        
    def set_mesh_grid(self):
        self.nx = int(self.deck.doc["Simulation"]["Number of Elements X"])
        self.ny1 = int(self.deck.doc["Simulation"]["Number of Elements Y1"])
        self.ny2=int(self.deck.doc["Simulation"]["Number of Elements Y2"])
        self.ny=self.ny1+self.ny2+1
    
        self.dx=self.geometry.Lx/self.nx
        self.dy1=self.geometry.Ly1/self.ny1
        self.dy2=self.geometry.Ly2/self.ny2
        self.dhe=self.geometry.Lhe
        
        
        X, Y = np.meshgrid(np.arange(0, self.ny+1), np.arange(0, self.nx))
        X=X[1,:].copy()
        Y=Y[:,1].copy()
        self.X = X
        self.Y = Y

    
        
        
        
    def set_dx(self):
        Mdx=np.zeros((self.ny, self.nx)) 
        Mdx[0:self.ny, 0:self.nx]=self.dx
        self.Mdx=Mdx.copy()
    
    def set_dy(self):
        Mdy=np.zeros((self.ny, self.nx)) 
        Mdy[0:self.ny1, 0:self.nx]=self.dy1
        Mdy[self.ny1:self.ny, 0:self.nx]=self.dy2
        Mdy[self.ny1, 1:-1]=self.dhe
        self.Mdy=Mdy.copy()
        
    def set_temperatures(self):
        T = np.zeros((self.ny, self.nx))        
        T[0:self.ny1, 0:self.nx] = self.deck.doc["Materials"]["Material1"]["Initial Temperature"] # Set array size and set the interior value with Tini
        T[self.ny1:self.ny, 0:self.nx] = self.deck.doc["Materials"]["Material2"]["Initial Temperature"] # Set array size and set the interior value with Tini
        T[self.ny1, 1:-1] = self.deck.doc["Processing Parameters"]["Temperature"]

        self.T = T.copy()
        self.T0=T.copy()
        
      
    def set_conductivity(self):
        KtotalX= np.zeros((self.ny, self.nx)) 
        KtotalX[0:self.ny1, 0:self.nx] = self.deck.doc["Materials"]["Material1"]["Thermal Conductivity X"]
        KtotalX[self.ny1:self.ny, 0:self.nx] = self.deck.doc["Materials"]["Material2"]["Thermal Conductivity X"]
        KtotalX[self.ny1, 1:-1]=3.0
        self.KtotalX=KtotalX
                                                                                         
        KtotalY= np.zeros((self.ny, self.nx)) 
        KtotalY[0:self.ny1, 0:self.nx] = self.deck.doc["Materials"]["Material1"]["Thermal Conductivity Y"]
        KtotalY[self.ny1:self.ny, 0:self.nx] = self.deck.doc["Materials"]["Material2"]["Thermal Conductivity Y"]      
        KtotalY[self.ny1, 1:-1]=3.0                                                                             
        self.KtotalY=KtotalY         

    def set_density(self):
        RhoTotal= np.zeros((self.ny, self.nx)) 
        RhoTotal[0:self.ny1, 0:self.nx] = self.deck.doc["Materials"]["Material1"]["Density"]
        RhoTotal[self.ny1:self.ny, 0:self.nx] = self.deck.doc["Materials"]["Material2"]["Density"]    
        RhoTotal[self.ny1, 1:-1]=4600                                                                            
        self.RhoTotal=RhoTotal  

    def set_specific_heat(self):
        CpTotal= np.zeros((self.ny, self.nx)) 
        CpTotal[0:self.ny1, 0:self.nx] = self.deck.doc["Materials"]["Material1"]["Cp"]
        CpTotal[self.ny1:self.ny, 0:self.nx] = self.deck.doc["Materials"]["Material2"]["Cp"]     
        CpTotal[self.ny1, 1:-1]=700                                                                              
        self.CpTotal=CpTotal  

    def set_diffusivity(self):                                                                                  
        DiffTotalX = np.zeros((self.ny, self.nx)) 
        DiffTotalX[0:,0:]=self.KtotalX[0:,0:]/(self.RhoTotal[0:,0:]*self.CpTotal[0:,0:])
        self.DiffTotalX = DiffTotalX.copy()
        
        DiffTotalY = np.zeros((self.ny, self.nx)) 
        DiffTotalY[0:,0:]=self.KtotalY[0:,0:]/(self.RhoTotal[0:,0:]*self.CpTotal[0:,0:])
        self.DiffTotalY = DiffTotalY.copy()

    def set_viscosity(self):         
        Visc=np.zeros((self.ny, self.nx))
        Visc[0:, 0:]=1.14*10**(-12)*np.exp(26300/self.T[0:, 0:])
        self.Visc=Visc.copy()

    def set_dic(self):
        Dic=np.ones((self.ny, self.nx))
        self.dic=1/(1+0.45)
        Dic[self.ny1,1:-1]=self.dic
        self.Dic0=Dic.copy()
        self.Dic=Dic.copy()   


    def set_heat_density(self):       
        self.q=float(self.deck.doc["Processing Parameters"]["Power Density"])
        Q=np.zeros((self.ny, self.nx))
        # Q[int(self.ny/2), 0:] = self.q
        # Q[int(self.ny/2-1), 0:] = self.q
        self.Q=Q.copy()


    def dist(self):
        
        distY=np.zeros((self.ny+1))
        distY[0:self.ny1+1]=self.X[0:self.ny1+1]*self.dy1
        distY[self.ny1+1]=distY[self.ny1]+self.dhe
        distY[self.ny1+2:self.ny+1]=distY[self.ny1+1]+self.X[1:self.ny2+1]*self.dy2
        # DistY[self.ny1, 1:-1]=self.dhe
      
        self.distY=distY.copy()
        