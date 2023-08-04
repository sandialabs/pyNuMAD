import matplotlib.pyplot as plt
import numpy as np


def plot_airfoil(airfoil):
    """Plot airfoil"""
    fig, ax = plt.subplots()
    # ax[0].plot(self.x,self.y,'.-')
    ax.plot(airfoil.coordinates[:, 0], airfoil.coordinates[:, 1], ".-")
    ax.plot(airfoil.c, airfoil.camber)
    # mtx = self.maxthick * np.array([1,1])
    # kn = find(self.c >= self.maxthick,1)
    # mty = self.camber(kn) + self.thickness(kn) * np.array([0.5,- 0.5])
    # line(mtx,mty,'LineStyle',':','Color','k')
    # else:
    fig.show()
    return ax

def plot_regions(blade):
    fig, ax = plt.subplots()
    
    n = blade.keypoints.shape[0]
    for kn in range(n):
        blade.hgKeypoints[kn] = ax.plot(
            blade.keypoints[kn,2,:],
            blade.keypoints[kn,0,:],
            blade.keypoints[kn,1,:])
    fig.show()
    return ax

def plot_geometry(self):
    fig, ax = plt.subplots()
    n = self.geometry.shape[2]
    for k in range(n):
        self.hgGeometry[k] = ax.plot(
            self.geometry[:,2,k],
            self.geometry[:,0,k],
            self.geometry[:,1,k])
    fig.show()
    return ax

def plot_profile(blade, k):
    fig, ax = plt.subplots()
    ax.plot(
        blade.profiles[:,0,k],
        blade.profiles[:,1,k],
        '.-')
    fig.show()
    return ax

# def plotbom(self,k = None): 
#     if not ('k' is not None)  or len(k)==0:
#         k = 1
#     if k=='cb':
#         k = np.rint(get(self.bomPlot.uisliderHP,'Value'))
#     if len(self.bomPlot.hgLinesHP)==0 or not np.all(ishandle(self.bomPlot.hgLinesHP)) :
#         clf
#         self.bomPlot.axHP = axes('Position',np.array([0.1,0.6,0.8,0.3]))
#         self.bomPlot.hgLinesHP = plt.plot(self.ispan,self.keyarcs(1,:) - self.HParcx0,'k-.',self.ispan,self.keyarcs(7,:) - self.HParcx0,'k-.')
#         self.bomPlot.hTitleHP = plt.title('','Interpreter','none')
#         self.bomPlot.axLP = axes('Position',np.array([0.1,0.2,0.8,0.3]))
#         self.bomPlot.hgLinesLP = plt.plot(self.ispan,self.keyarcs(13,:) - self.LParcx0,'k-.',self.ispan,self.keyarcs(7,:) - self.LParcx0,'k-.')
#         self.bomPlot.hTitleLP = plt.title('','Interpreter','none')
#         n = self.bom['hp'].shape[1-1]
#         self.bomPlot.uisliderHP = uicontrol('Style','slider','Min',1,'Max',n,'Value',1,'SliderStep',np.array([1 / (n - 1),10 / (n - 1)]),'Units','normalized','Position',np.array([0.02,0.02,0.96,0.04]),'Callback',lambda src = None,evt = None: self.plotbom('cb'))
#     k = np.amin(k,self.bom['hp'].shape[1-1])
#     self.bomPlot.kLayer = k
#     set(self.bomPlot.uisliderHP,'Value',k)
#     str = sprintf('HP Layer #d: #s',k,self.bom['hp'][k,3])
#     set(self.bomPlot.hTitleHP,'String',str)
#     str = sprintf('LP Layer #d: #s',k,self.bom['lp'][k,3])
#     set(self.bomPlot.hTitleLP,'String',str)
#     hp = self.bomIndices.hp(k,:)
#     x = self.ispan(np.arange(hp(1),hp(2)+1))
#     y1 = self.keyarcs(hp(3),np.arange(hp(1),hp(2)+1)) - self.HParcx0(np.arange(hp(1),hp(2)+1))
#     y2 = self.keyarcs(hp(4),np.arange(hp(1),hp(2)+1)) - self.HParcx0(np.arange(hp(1),hp(2)+1))
#     if ishandle(self.bomPlot.hgPatchHP):
#         os.delete(self.bomPlot.hgPatchHP)
#     # 
#     axes(self.bomPlot.axHP)
#     self.bomPlot.hgPatchHP = patch(np.array([x,fliplr(x)]),np.array([y1,fliplr(y2)]),'b')
#     lp = self.bomIndices.lp(k,:)
#     x = self.ispan(np.arange(lp(1),lp(2)+1))
#     y1 = self.keyarcs(lp(3),np.arange(lp(1),lp(2)+1)) - self.LParcx0(np.arange(lp(1),lp(2)+1))
#     y2 = self.keyarcs(lp(4),np.arange(lp(1),lp(2)+1)) - self.LParcx0(np.arange(lp(1),lp(2)+1))
#     if ishandle(self.bomPlot.hgPatchLP):
#         os.delete(self.bomPlot.hgPatchLP)
#     # 
#     axes(self.bomPlot.axLP)
#     self.bomPlot.hgPatchLP = patch(np.array([x,fliplr(x)]),np.array([y1,fliplr(y2)]),'b')
#     return

