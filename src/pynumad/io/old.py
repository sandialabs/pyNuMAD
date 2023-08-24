def writeBOMxls(self,file = None): 
        # This method writes the bill-of-materials out to a spreadsheet.
        
        # Example:
        
        #   ``bladeDef.writeBOMxls('bom.xlsx')``
        
        m_to_mm = 1000.0
        if os.path.exist(str('BOM_template.xlsx')):
            copyfile('BOM_template.xlsx',file)
        
        header = np.array([['Layer #','Material ID','Component','Begin Station','End Station','Max width','Mean width','3D area','Layer Thickness','Computed layer weight'],['','','','(m)','(m)','(m)','(m)','(m^2)','(mm)','(g)']])
        # LP skin table
        array_ = np.array([[header],[self.bom['lp']]])
        xlswrite(file,array_,'LP skin')
        # HP skin table
        array_ = np.array([[header],[self.bom['hp']]])
        xlswrite(file,array_,'HP skin')
        # shear web table
        array_ = np.array([np.array([['SW #'],['']]),header])
        
        for k in np.arange(1,np.asarray(self.bom.sw).size+1).reshape(-1):
            nr = self.bom.sw[k].shape[1-1]
            array_ = np.array([[array_],[np.array([np.matlib.repmat(np.array([k]),nr,1),self.bom.sw[k]])]])
        
        xlswrite(file,array_,'shear webs')
        # bond line lengths
        array_ = np.array([['','Length'],['','(mm)'],['Root diameter',self.ichord(1) * m_to_mm],['LE bond',np.round(self.bom.lebond)],['TE bond',np.round(self.bom.tebond)]])
        for r in np.arange(1,2+1).reshape(-1):
            for k in np.arange(1,np.asarray(self.bom.swbonds).size+1).reshape(-1):
                surfs = np.array(['HP','LP'])
                str = sprintf('#s bond, SW #d',surfs[r],k)
                cellrow = np.array([str,np.round(self.bom.swbonds[k](r))])
                array_ = np.array([[array_],[cellrow]])
        
        xlswrite(file,array_,'lengths')
        return
    
    
def writePlot3D(self,file = None,breakpoints = None):
        #NOTE ask team about this
        # ignore for now  

        # Write the current blade geometry in Plot3D format.
        # breakpoints is a list of chord fractions at which the
        # surface geometry is divided into blocks
        
        # Examples:
        
        #   ``BladeDef.writePlot3D(filename,[breakpoints])``
        
        #   ``BladeDef.writePlot3D('file.p3d',[-.3, .3]);``
        
        if not ('breakpoints' is not None) :
            breakpoints = []
        
        indicesOfBreakpoints = np.zeros((1,np.asarray(breakpoints).size))
        # get the chordwise spacing of points, assuming identical
        # spacing for all stations
        chordspacing = self.cpos[:,1]
        for kBreakpoint in np.arange(1,np.asarray(breakpoints).size+1).reshape(-1):
            bp = breakpoints(kBreakpoint)
            __,ind = np.amin(np.sqrt((chordspacing - bp) ** 2))
            indicesOfBreakpoints[kBreakpoint] = ind
        
        N,M = self.cpos.shape
        INCLUDE_REPEATS = False
        if INCLUDE_REPEATS:
            indicesOfBreakpoints = unique(np.array([1,indicesOfBreakpoints,N]))
        else:
            indicesOfBreakpoints = unique(np.array([2,indicesOfBreakpoints,N - 1]))
        
        fid = open(file,'wt')
        
        if (fid == - 1):
            raise Exception('Could not open file "#s"',file)
        
        # output the data in Plot3d format
        #TODO
        # try:
        #     nBlocks = np.asarray(indicesOfBreakpoints).size - 1
        #     fid.write('#d\n' # (nBlocks))
        #     for kblock in np.arange(1,nBlocks+1).reshape(-1):
        #         a = indicesOfBreakpoints(kblock)
        #         b = indicesOfBreakpoints(kblock + 1)
        #         fid.write('#d  #d  #d\n' # (1 + b - a,M,1))
        #     columnsPerLine = 5
        #     for kblock in np.arange(1,nBlocks+1).reshape(-1):
        #         a = indicesOfBreakpoints(kblock)
        #         b = indicesOfBreakpoints(kblock + 1)
        #         self.fprintf_matrix(fid,self.geometry(np.arange(a,b+1),1,:),columnsPerLine)
        #         self.fprintf_matrix(fid,self.geometry(np.arange(a,b+1),2,:),columnsPerLine)
        #         self.fprintf_matrix(fid,self.geometry(np.arange(a,b+1),3,:),columnsPerLine)
        # finally:
        #     pass
        
        # fid.close()
        # return
    
    def fprintf_matrix(self,fid = None,matrixData = None,columnsPerLine = None): 
        kColumn = 1
        for kData in np.arange(1,np.asarray(matrixData).size+1).reshape(-1):
            fid.write('%g ', (matrixData(kData)))
            kColumn = kColumn + 1
            if kColumn > columnsPerLine:
                fid.write('\n')
                kColumn = 1
        
        return