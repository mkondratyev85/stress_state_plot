import operator
from math import cos, sin, acos, asin, atan, tan, radians, degrees, sqrt

class Plane(object):
    ''' The Plane is a class of a surface of a fracture or a tranch in a therms 
    of a structural geology (a Dir angle and a Dip angle). This class could 
    perform relationship between Cos1--Cos2-Cos3 and Dir and Dip angles. 
    Cos1-Cos2-Cos3 are also coordinats of a 3D vector of a unit lenght. The 
    direction of this vector is coincides with the dipping of a plane.
    The direction of a Cos1 axis is coincides with the North. The direction of a
    Cos3 axis is coincides with the East. The direction of a Cos2 is coincides 
    with the upward direction. So, the [1,0,0] vector is looking at the North. 
    The [0,0,1] vector is looking at the East. And the [0,1,0] vector is looking
    upward. Positive values of Cos2 means upward direction of the vector, while 
    negative values of Dip angle means downward direction of a vector.
    
    '''
    __slots__ = ['dir', 'dip'] 
       
    
    def __init__(self, dir_or_pair, dip = None):
        if dip == None:
            self.dir, self.dip = self.dirdip(dir_or_pair[0], dir_or_pair[1])
        else:
            self.dir, self.dip = self.dirdip(dir_or_pair, dip)      
    
    def __len__(self):
        return 2
    
    def __getitem__(self, key):
        if key == 0:
            return self.dir
        elif key == 1:
            return self.dip
        else:
            raise IndexError("Invalid subscript "+str(key)+" to a plane")
            
    def __setitem__(self, key, value):
        if key == 0:
            self.dir, self.dip = self.dirdip(value, self.dip)
        elif key == 1:
            self.dir, self.dip = self.dirdip(self.dir, value)
        else:
            raise IndexError("Invalid subscript "+str(key)+" to a plane")
    
    # String representaion (for debugging)
    def __repr__(self):
        return 'plane(%06.2f, %05.2f)' % (self.dir, self.dip)

    def dirdip(self, dr, dp):
        dr = dr % 360
        dp = (dp + 90) % 360 - 90
        if dp > 90: 
            dr = (dr + 180) % 360
            dp = 180 - dp
        return dr, dp
    
    # Cos-sin operations
    def get_cos3(self):
        return sin(radians(self.dir)) * cos(radians(self.dip))
    def __setcos3(self, val):
        pass
    cos3 = property(get_cos3, __setcos3, None, "gets or sets the guide cos of the plane")
    
    def get_cos2(self):
        return -1 * sin(radians(self.dip))
    def __setcos2(self, val):
        pass
    cos2 = property(get_cos2, __setcos2, None, "gets or sets the guide cos of the plane")
    
    def get_cos1(self):
        return cos(radians(self.dir)) * cos(radians(self.dip))
    def __setcos1(self, val):
        pass
    cos1 = property(get_cos1, __setcos1, None, "gets or sets the guide cos of the plane")    
    
    
    
    def return_angle_between(self, other):
        '''
        Returns angle betwwen two 3d vectors defined by cos1-cos2-cos3. If you 
        want to find an angle between to planes it's better to call "normal"
        method before call this one
        '''
        if (other.dir == self.dir) and (other.dip == self.dip):
            return 0
        if abs(other.dir-self.dir)==180 and (other.dip==-1*self.dip):
            return 180
        angle = (self.cos1 * other.cos1 + \
                               self.cos2 * other.cos2 + self.cos3 * other.cos3)        
        return degrees(acos(angle))
        
    def normal(self):
        if not self.dip < 0:
            return Plane(self.dir, self.dip - 90)
        else:
            return Plane(self.dir, self.dip + 90)
    
        
    def define_by_plane_cos(self, cos1, cos2, cos3):
        '''
        Definy self by cos1-cos2-cos3 given for vector coincides with dipping 
        dirrection of a plane. This method raise Error while Dip angle of a 
        plane is equal is equal to +-1. If you want to define some plane by 
        cos1-cos2-cos3 it is better to use definy_by_normal_cos method.
        '''
        if cos1**2 + cos2**2 + cos3**2 < 0.98:
            raise ValueError ("cos1**2 + cos2**2 + cos3**3 must be equal to 1")
        
        dp = degrees(asin(cos2*(-1)))
        if dp > 90: dp = -1 * dp
        
        if abs(cos2) < 1:
            if cos1==0:
                if cos3>0:
                    dr = 90
                else:
                    dr = 270
            else:
                angl = degrees(atan(cos3/cos1))
                if cos1<0:
                    dr = 180 + angl
                else:
                    dr = 360 + angl
        else:
            dr = 0
        
        self.dir, self.dip = self.dirdip(dr, dp)
        self.dir = dr % 360     
        
        
    def define_by_normal_cos(self, cos1, cos2, cos3):
        '''
        Define self by cos1-cos2-cos given for vector coincides with normal of
        a plane. When horizontal plane is definded, Dir Angle of this plane is 
        set to 000
        '''
        if cos1**2 + cos2**2 + cos3**2 < 0.98:
            raise ValueError ("cos1**2 + cos2**2 + cos3**3 must be equal to 1")
        
        dp = degrees(acos(cos2))
        if dp > 90: dp = -1 * dp
            
        if abs(cos2) < 1:            
            if cos1 == 0:
                angl = 90
            else:
                angl = degrees(atan(cos3/cos1))
            if cos1<0:
                dr = 180 + angl
            else:
                dr = 360 + angl
        else:
            dr = 0
        
        self.dir, self.dip = self.dirdip(dr, dp)
        self.dir = dr % 360        

    def rotated(self, angel):
        '''
        Returns a vector rotated on 'angel' agains normal of self plane.
        '''        
        
        Dir_incr = degrees(atan(tan(radians(angel))/cos(radians(self.dip))))
        Dip_ = degrees(atan(cos(radians(Dir_incr))*tan(radians(self.dip))))
        if abs(angel) == 90:
            Dip_ = 0
            Dir_incr = angel
        return Plane(self.dir + Dir_incr, Dip_)
    
    def get_perpendicular_between(self, other):
        result = Plane(0, 0)
        angle = self.return_angle_between(other)
        sin_ = sin(radians(angle))
        result.define_by_normal_cos((self.cos2 * other.cos3 - self.cos3 * other.cos2)/sin_, \
                (self.cos3 * other.cos1 - self.cos1 * other.cos3) / sin_, \
                (self.cos1 * other.cos2 - self.cos2 * other.cos1) / sin_)
        
        return result.normal()

    def rotate_cos2(self, angle):
        angle = radians(angle)
        cosa = cos(angle); sina = sin(angle)
        cos1 = self.cos1*cosa - self.cos3*sina
        cos2 = self.cos2
        cos3 = self.cos1*sina + self.cos3*cosa
        self.define_by_plane_cos(cos1,cos2,cos3)
    
    def rotate_cos1(self, angle):
        angle = radians(angle)
        cosa = cos(angle); sina = sin(angle)
        cos1 = self.cos1
        cos2 = self.cos3*sina + self.cos2*cosa
        cos3 = self.cos3*cosa - self.cos2*sina
        self.define_by_plane_cos(cos1,cos2,cos3)

    def rotate_cos3(self, angle):
        angle = radians(angle)
        cosa = cos(angle); sina = sin(angle)
        cos1 = self.cos1*cosa + self.cos2*sina
        cos2 = -1*self.cos1*sina + self.cos2*cosa
        cos3 = self.cos3
        self.define_by_plane_cos(cos1,cos2,cos3) 

def make_plane_look_upward(plane):
    if plane.dip < 0:
        plane.dip *= -1
        plane.dir = (plane.dir + 180) % 360


def plane2xy(plane):
    make_plane_look_upward(plane)
    r = tan(radians(90 - plane.dip) / 2)
    x, y =  -1*r*sin(radians(plane.dir)), r*cos(radians(plane.dir))
    return   x,y
