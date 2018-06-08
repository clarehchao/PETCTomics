# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 19:42:54 2015

@author: rharnish
"""

#%%

import numpy as np
import dicom_series

#%%

class ImageGeometry:
    """ ImageGeometry """


    def __init__(self,*args):

        # origin always represented in Patient XYZ coords. like in DICOM
        self._origin = None
        
        # span in XYZ
        self._span = None

        # numpy matrices are indexed Slice=Plane, Row, Col, hence SRC
        self._shapeSRC     = None
        self._samplingSRC  = None
        
        # ITK and MATLAB matrices are indexed Row, Col, SLice, hence RCS
        self._shapeRCS     = None
        self._samplingRCS  = None
        
        # direction cosines give directions along each dimension where
        # each is a vector [X,Y,Z]
        self._dir_as_col_idx_increases   = None # --> along rows as column index increases
        self._dir_as_row_idx_increases   = None # 
        self._dir_as_slice_idx_increases = None # from slice to slice
        
        # axial, coronal, sagittal, oblique
        self._orientation = None
        
        # dir cosine matrix
        self._direction_cosine_matrix = None
        
        # allow for ImageGeometry to be initialized with 
        # a DicomSeries object as arg[0]
        if args:
            if isinstance(args[0],dicom_series.DicomSeries):
                series = args[0]
                self.get_geometry_from_dicom_series(series)
            inputTypeString = type(args[0]).__name__
            if 'itkImage' in inputTypeString:
                itkImage = args[0]
                self.get_geometry_from_itkImage(itkImage)
        
        
    # ---------------------------------------------------------------------    
    # Properties
    # --------------------------------------------------------------------- 
    @property
    def shapeSRC(self):
        """ """
        return self._shapeSRC

    @property
    def samplingSRC(self):
        """ """
        return self._samplingSRC
        
    @property
    def shapeRCS(self):
        """ """
        return self._shapeRCS
    
    @property
    def samplingRCS(self):
        """ """
        return self._samplingRCS
        
    @property
    def origin(self):
        """ The XYZ position of the center of the (0,0,0) voxel """
        return self._origin
        
    @property
    def span(self):
        """ The XYZ exent of the image """
        return self._span
        
    @property
    def dir_as_col_idx_increases(self):
        """ """
        return self._dir_as_col_idx_increases
        
    @property
    def dir_as_row_idx_increases(self):
        """ """
        return self._dir_as_row_idx_increases

    @property
    def dir_as_slice_idx_increases(self):
        """ """
        return self._dir_as_slice_idx_increases
        
    @property 
    def direction_cosine_matrix(self):
        """ """
        return self._direction_cosine_matrix
        
    @property
    def orientation(self):
        """ """
        return self._orientation
        
        
           
    # --------------------------------------------------------------------- 
    # Private methods
    # --------------------------------------------------------------------- 
    def _determine_orientation(self):
        
        orientation = 'undetermined'     
        orvec = np.cross(self._dir_as_col_idx_increases,self._dir_as_row_idx_increases)

        if abs(np.dot(orvec,[0.0, 0.0, 1.0])) > 0.9999:
            orientation = 'axial'
            
        if abs(np.dot(orvec,[0.0, 1.0, 0.0])) > 0.9999:
            orientation = 'coronal'
            
        if abs(np.dot(orvec,[1.0, 0.0, 0.0])) > 0.9999:
            orientation = 'sagittal'
            
        # TODO add oblique
            
        self._orientation = orientation
        
        return orientation
        
        
#    def _generate_direction_cosine_matrix(self):
#        """ Use row, col, and slice dir to generate direcion cosine matrix """
#        irs = self._dir_as_row_idx_increases
#        ics = self._dir_as_col_idx_increases
#        iss = self._dir_as_slice_idx_increases
#        dir_cos_matrix = np.matrix([irs, ics, iss])
#        self._direction_cosine_matrix = dir_cos_matrix
        
    def _generate_direction_cosine_matrix(self):
        """ Use row, col, and slice dir to generate direcion cosine matrix """
        irs = self._dir_as_row_idx_increases
        ics = self._dir_as_col_idx_increases
        iss = self._dir_as_slice_idx_increases
        dir_cos_matrix = np.matrix([irs, ics, iss]).T
        self._direction_cosine_matrix = dir_cos_matrix
        
    def _generate_I2W_affine(self):
        self._generate_direction_cosine_matrix()
        A = np.zeros([4,4])
        A[0:3,0] = [float(i) * self.samplingRCS[0] for i in self._dir_as_row_idx_increases]   
        A[0:3,1] = [float(i) * self.samplingRCS[1] for i in self._dir_as_col_idx_increases]   
        A[0:3,2] = [float(i) * self.samplingRCS[2] for i in self._dir_as_slice_idx_increases]         
        A[0:3,3] = [float(i) for i in self._origin]
        A[3,3]   = 1
        self._I2W = np.matrix(A)
        return self._I2W
        
    def _generate_W2I_affine(self):
        I2W = self._I2W
        W2I = np.linalg.inv(I2W)
        print I2W
        print W2I
        self._W2I = W2I
        return self._W2I
        
    def idx_to_coords(self,row_idx,col_idx,slice_idx):
        '''
        x,y,z = idx_to_coords(r,c,s)
        '''
        row_spacing   = self._samplingRCS[0]
        col_spacing   = self._samplingRCS[1]
        slice_spacing = self._samplingRCS[2]
        
        increasing_rows   = self._dir_as_row_idx_increases
        increasing_cols   = self._dir_as_col_idx_increases
        increasing_slices = self._dir_as_slice_idx_increases
        
        offsetmmrow   = np.multiply(row_idx * row_spacing,     increasing_rows)
        offsetmmcol   = np.multiply(col_idx * col_spacing,     increasing_cols)
        offsetmmslice = np.multiply(slice_idx * slice_spacing, increasing_slices)
        offsetmm      = np.add( offsetmmslice, np.add(offsetmmrow,offsetmmcol) )        
        
        return np.add(self._origin,offsetmm) 
      
    def coords_to_idx(self,x,y,z):
        '''
        r,c,s = coords_to_idx(x,y,z)
        '''      
        
        # reverse translation by origin
        x = x - self.origin[0]  
        y = y - self.origin[1] 
        z = z - self.origin[2]
        
        increasing_rows   = self._dir_as_row_idx_increases
        increasing_cols   = self._dir_as_col_idx_increases
        increasing_slices = self._dir_as_slice_idx_increases
        
        '''
        Get contribution by row to [x,y,z] by dot product
        of row_direction with [x,y,z] divided by spacing
        in the row direction. Do the equivalent for columns
        and slices
        '''
        r = np.dot( np.array([x,y,z]), increasing_rows )   / self.samplingRCS[0]     
        c = np.dot( np.array([x,y,z]), increasing_cols )   / self.samplingRCS[1]
        s = np.dot( np.array([x,y,z]), increasing_slices ) / self.samplingRCS[2]
        
        # round to nearest integer
        # return int(r),int(c),int(s)
        return int(round(r)), int(round(c)), int(round(s))
      
    def _determine_span(self):
        span_near = self.idx_to_coords(0,0,0)
        s = self._shapeRCS
        span_far  = self.idx_to_coords(s[0],s[1],s[2])
#        print 'span near : ', span_near
#        print 'span far  : ', span_far
        span = np.subtract(span_far,span_near)
        span = [abs(s) for s in span]
        self._span = span
        
        
        
    # --------------------------------------------------------------------- 
    # Public methods
    # --------------------------------------------------------------------- 
    def get_geometry_from_dicom_series(self,series):
        # DICOM IPP <==> position of first voxel (0,0,0) <==> "origin"
        self._origin = [float(i) for i in series.info.ImagePositionPatient]
        
        # from DICOM IOP
        iop = [float(i) for i in series.info.ImageOrientationPatient]
        self._dir_as_col_idx_increases = iop[0:3] # direction along rows (with col idx increased)
        self._dir_as_row_idx_increases = iop[3:]  # direction along columns (with row index increased)
        self._dir_as_slice_idx_increases = series.slice_direction # direction from slice to slice (slice idx increased)

        # DicomSeries object returns shape in numpy (SRC) order
        self._shapeSRC    = [int(i) for i in series.shape]
        self._samplingSRC = [float(i) for i in series.sampling]
        
        # create new lists representing shape and sampling in RCS order
        self._shapeRCS    = self._shapeSRC[1:3]    + [self._shapeSRC[0]]
        self._samplingRCS = self._samplingSRC[1:3] + [self._samplingSRC[0]]
                
        self._determine_orientation()
        self._determine_span()
        self._generate_direction_cosine_matrix()

    def get_geometry_from_itkImage(self,itkImage):
        
        # use accessors to get info from itkImage
        origin     = itkImage.GetOrigin()
        spacing    = itkImage.GetSpacing()
        region     = itkImage.GetLargestPossibleRegion()
        size       = region.GetSize()
        vnl_matrix = itkImage.GetDirection().GetVnlMatrix()         
        
        # ITK image gives origin as XYZ
        self._origin = [float(i) for i in origin]
             
        # TODO: really get this straightened out     
        # columns of ITK vnl_matrix represent directions of R, C, then S
        self._dir_as_row_idx_increases   = [vnl_matrix.get(i,0) for i in np.arange(3)]        
        self._dir_as_col_idx_increases   = [vnl_matrix.get(i,1) for i in np.arange(3)]
        self._dir_as_slice_idx_increases = [vnl_matrix.get(i,2) for i in np.arange(3)]
              
        # sampling and size given in RCS order by ITK images      
        self._samplingRCS = [float(i) for i in spacing]
        self._shapeRCS    = [float(i) for i in size]

        # permute RCS --> SRC
        self._shapeSRC    = [self._shapeRCS[2]]    + self._shapeRCS[0:2]
        self._samplingSRC = [self._samplingRCS[2]] + self._samplingRCS[0:2]                
                
        self._determine_orientation()
        self._determine_span()
        self._generate_direction_cosine_matrix()
        
        
    def print_self(self):
        print 'origin          : ', self.origin
        print 'span            : ', self.span
        print
        print 'shapeRCS        : ', self.shapeRCS
        print 'samplingRCS     : ', self.samplingRCS
        print
        print 'shapeSRC        : ', self.shapeSRC
        print 'samplingSRC     : ', self.samplingSRC
        print
        print 'dir_as_row_idx_increases   : ', self._dir_as_row_idx_increases
        print 'dir_as_col_idx_increases   : ', self._dir_as_col_idx_increases
        print 'dir_as_slice_idx_increases : ', self._dir_as_slice_idx_increases
        print 
        print 'orientation     : ', self._orientation
        print
        
        
        