import numpy as np

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return (vector / np.linalg.norm(vector))

def YPR_transform2(x,y,z,pt1,pt2,angle):
    acos=np.cos(np.radians(angle))
    asin=np.sin(np.radians(angle))
    a=pt1[0]
    b=pt1[1]
    c=pt1[2]
    d=pt2[0]
    e=pt2[1]
    f=pt2[2]
    unvec=unit_vector((a-d,b-e,c-f))
    u=unvec[0]
    v=unvec[1]
    w=unvec[2]

    #TRANSFORMATION MATRIX
    xpos=(a*(v**2+w**2) - u*(b*v + c*w - u*x - v*y - w*z))*(1-acos) + x*acos + (-c*v + b*w - w*y + v*z)*asin
    ypos=(b*(u**2+w**2) - v*(a*u + c*w - u*x - v*y - w*z))*(1-acos) + y*acos + (c*u - a*w + w*x - u*z)*asin
    zpos=(c*(u**2+v**2) - w*(a*u + b*v - u*x - v*y - w*z))*(1-acos) + z*acos + (-b*u + a*v - v*x + u*y)*asin

    return(xpos,ypos,zpos)


#ROTATE O
xo,yo,zo=YPR_transform2(15.260012,6.517810,16.639084,(15.769703,7.170020,15.728123),(15.250550,8.437755,15.245168),90)
print(xo,yo,zo)

#ROTATE H
xh,yh,zh=YPR_transform2(16.689814,6.816198,15.209056,(15.769703,7.170020,15.728123),(15.250550,8.437755,15.245168),90)
print(xh,yh,zh)
