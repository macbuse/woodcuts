#!/usr/bin/python

import inkex, sys, simplestyle, cubicsuperpath
import  numpy as np
from itertools import chain

def flatten(x):
   return list(chain(*x))

def return_copy(func):
   def wrapper(self, *args):
        payload = func(self,*args)
        tt = self.copy()
        tt.flat = payload
        tt.update()
        return tt
    
   return wrapper

class Path(list):
    def __init__(self, pth=None, copy=False):
        if copy: return
        self.obj = pth
        if pth == None:
           return 
        self.base_list = cubicsuperpath.parsePath(pth)
        self.shape = [ len(subp) for subp in self.base_list]
        
        self.flat = []
        for subp in self.base_list:
            for node in subp:
                tt = [np.array(x) for x in node]
                self.flat.extend(tt)
                
        self.center = self.flat[1]
        self.start = self.flat[1]
        self.end = self.flat[-2]
        
    def update(self):
        self.center = self.flat[1]
        self.start = self.flat[1]
        self.end = self.flat[-2]
        
        
    def copy(self):
        xx = Path(None,copy=True)
        xx.shape = self.shape[:]
        xx.flat = self.flat[:]
        return xx
        
    def paths(self):
        u  = self.flat
        nodes = [[ list(x) for x in node]
                           for node in zip(u[::3],u[1::3],u[2::3]) ]
        
        tt = [ [ nodes.pop(0) for j in range(k)]
                              for k in self.shape]
        return cubicsuperpath.formatPath(tt)
    
    @return_copy
    def __add__(self,other):
        return  [ v + other for v in self.flat]

    @return_copy
    def __rmul__(self,other):
        if isinstance(other,(int,float)):
            return [ v*other for v in self.flat]

    @return_copy
    def __div__(self,other):
        if isinstance(other,(int,float)): 
            return [ v/other for v in self.flat]

    @return_copy     
    def __mul__(self,other):
        
        if isinstance(other,(int,float)):
            return [ (v - self.center)*other + self.center   for v in self.flat]
      
        if isinstance(other,np.ndarray): 
            return [ np.dot(other, v - self.center) + self.center for v in self.flat]
                   
    def __str__(self):
        return str(self.shape)
    

def mk_plugs(CV, width=100, N=2):
   
    tx = CV.end - CV.start
    R = np.array([[-1,0],[0,-1]])
    
    motifs =  flatten([CV.flat , reversed( (CV*R + 2*tx).flat) ])
    payload = flatten([ [ y + 2*k*tx for y in motifs] for k in range(N) ] )
    
    scale = width/(2*N*tx[0])
    cp = Path()
    cp.flat = payload
    cp.update()
    cp.shape = [len(payload)/3]
    cp = cp*scale
    return cp
   
   
def slitsZ(xnum = 6, ynum = 5,
           ww = 100, hh = 6,
           inter_cut = 8, 
           width=100,
           diff=3.,
           top_left=np.array([0,0]),
           LRdiff=0,
           hinge_offset=0) :
   
    #make the initial slit
    hh += .5*diff
    inter_row =  hh
    pts = [ [0,-hh], [.5*ww-hh  ,-hh], [.5*ww+hh,hh],  [ww,hh] ]
    pts = [np.array(x) for x in pts]
    pts.extend([x - diff*np.array([0,1]) for x in reversed(pts)])
    pts.append(pts[0])
    
    #slits are interleaved so you don't have to move them so much to center
    offs = .125*np.array([ (inter_cut + ww)*xnum + ww, 0])
    subp = [ x - offs for x in pts]
    hh -= .5*diff
    
    #set up some stuff for doing affine transformations
    tx = .5*np.array([ inter_cut + ww, 0])
    ty1 =  2*np.array([0, inter_row + hh ])
    ty2 =  np.array([0, 4*inter_row + 4*hh + 2*diff ])
    
    SX = np.array([[1,0],[0,-1]])
    
    row = []

    for j in range(xnum):    
        row.extend([x + j*tx for x in subp ])
        #we could calulate these but why bother ?
        if j == 0 : min_x= row[2][0]
        if j == xnum - 1: max_x = row[-2][0]
        row.extend([np.dot(SX,x) + j*tx + ty1 for x in subp ])
        
    rows = []
    for k in range(ynum):
        rows.extend([x + k*ty2 for x in row ])
    
    #trim to fit a rectangle min_x and max_x are calculated above
    for k,x in enumerate(rows):
      #continue
      rows[k][0] = min( max( x[0],min_x), max_x)
    
    scale = width/(max_x - min_x)
    real_top_left = top_left + np.array([0,hinge_offset])
    rows = [scale*(x - rows[0]) +  real_top_left for x in rows]
    #the last 2 points are raised so take the third last
    
    bottom = (rows[-3] - top_left)[1]
    rows = [[rows.pop(0) for x in subp] for j in range(len(rows)/len(subp) -1 )]
    
    for k,subp in enumerate(rows):
      subp.append(subp[-2])
      tt  = [a for a,b in zip(subp,subp[1:]) if sum(abs(a-b)) > .01]
      uu = [a[0] for a in tt]
      rows[k] = None if max(uu) == min(uu) else tt
    rows = [x for x in rows if x is not None ]
    
    edge1 = [ [0,0], [0, bottom] ]
    edge1 = [np.array(x) + top_left for x in edge1]
    edge2 = [x + np.array([width,0]) for x in edge1]
    edge2[0][1] += LRdiff
    rows.extend([edge1,edge2] )

    return rows
   
def pts2curve(pairs):
    '''makes a polyline path element from a list of np array
    '''

    pth = [ '%.2f, %.2f '%tuple(z) for z in pairs]
    return 'M '+ ''.join(pth)


class Project(inkex.Effect):
    def __init__(self):
        inkex.Effect.__init__(self)
        #load up the parameters from the GUI
        
        def add_opt(var, var_type, default):
            self.OptionParser.add_option("", "--" + var,
                                     action="store", type=var_type,
                                     dest=var, default=default,
                                     help="command line help")
       
        # list of parameters defined in the .inx file
        #nothing fancy here so delegate 
        add_opt('nx','int',3)
        add_opt('ny','int',3)
        add_opt('width','int',100)
        add_opt('num_pegs','int',1)
        add_opt('depth','int',30)
        add_opt('spacing','int',8)
        add_opt('stack','float',2.)
        add_opt('diff','float',3.)

        
        # I won't touch this as it gives errors if I cut it out
        self.OptionParser.add_option("", "--active-tab",
                                     action="store", type="string",
                                     dest="active_tab", default='title', # use a legitmate default
                                     help="Active tab.")
   
        
    def effect(self):
      
        if self.selected == {}:
           inkex.errormsg('Select a curve and try again')
           sys.exit(1)
        
        obj = self.selected[self.options.ids[0]]
        
        #should be some checking here
        pth = Path(obj.get('d'))
        
        plugs = mk_plugs(pth,
                      width=self.options.width,
                      N=self.options.num_pegs)
        
        LRdiff = (plugs.flat[-2] - plugs.flat[1])[1]
        
        pegs_attribs = { 'd': plugs.paths()}
        
        payload = slitsZ(xnum=self.options.nx+2,
                         ynum=self.options.ny,
                         width=self.options.width,
                         ww=100 - self.options.spacing,
                         inter_cut=self.options.spacing,
                         diff=self.options.diff,
                         hh=self.options.stack,
                         top_left=pth.center,
                         LRdiff=LRdiff,
                         hinge_offset=self.options.depth)
        
        slitz_attribs = { 'd': ' '.join([pts2curve(x) for x in payload])}

        gp = self.mk_group()
        inkex.etree.SubElement( gp,
                               inkex.addNS('path','svg'),
                               pegs_attribs )
        
        inkex.etree.SubElement( gp,
                               inkex.addNS('path','svg'),
                               slitz_attribs )
        
        
    def mk_group(self,tag='greg_gp'):

            style = { 'stroke': '#FF0000',
                      'fill': 'none' ,
                      'stroke-width': 1.}
            
            style = simplestyle.formatStyle(style) 
 
            # Make a nice useful name
            g_attribs = { inkex.addNS('label','inkscape'): tag,
                          inkex.addNS('transform-center-x','inkscape'): str(0),
                          inkex.addNS('transform-center-y','inkscape'): str(0),
                          #'transform': 'translate(%s,%s)' % self.view_center[:2],
                          'style' : style,
                          'info':'N: '+ tag}
            
            # add the group to the document's current layer
            return inkex.etree.SubElement(self.current_layer, 'g', g_attribs )
            
    
if __name__ == '__main__':
    e = Project()
    e.affect()

