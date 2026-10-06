import math
import numpy as np

# Format a list of floats into a string for output.
def vals_to_str(vals):
  text = f"  " + f"  ".join(f"{x:12.4f}" for x in vals) + f"\n"
  return text

# For given angle-of-attack,sideslip (alpha,beta),
# output a Tecplot file to plot vector components.
if __name__ == "__main__":
    filename_tec = 'aoa_vectors.tec'
    #
    pi = math.pi
    rad_per_deg = pi / 180.
    #
    alpha =    15 * rad_per_deg
    beta  =    10 * rad_per_deg
    Vmag  =  -1.5 # neg so vector points upstream
    xmin  = 0. ; xmax = 0.5
    ymin  = 0. ; ymax = 0.3
    zmin  = 0. ; zmax = 0.5
    zero  = 0.
    #
    u = Vmag * math.cos(beta) * math.cos(alpha)
    v = Vmag * math.cos(beta) * math.sin(alpha)
    w = Vmag * math.sin(beta)
    #
    vecXunit = np.array([1.,0.,0.])
    vecYunit = np.array([0.,1.,0.])
    vecZunit = np.array([0.,0.,1.])
    vecXaxis = vecXunit * xmax
    vecYaxis = vecYunit * ymax
    vecZaxis = vecZunit * zmax
    #
    vecZero  = np.array([0.,0.,0.])
    vecV     = np.array([u,v,w])
    vecVx    = np.dot(vecV,vecXunit) * vecXunit
    vecVy    = np.dot(vecV,vecYunit) * vecYunit
    vecVz    = np.dot(vecV,vecZunit) * vecZunit
    vecVxy   = vecVx + vecVy
    vecVxz   = vecVx + vecVz
    vecVyz   = vecVy + vecVz
    #
    try:
        file = open(filename_tec, 'w', encoding="ascii")
    except:
        print(f"Error: Unable to write file {filename_tec}.")
    with file:
        file.write("Variables=\"x\",\"y\",\"z\"\n")
        #
        file.write(f"Zone T=\"XY-plane\" I=2, J=2, F=POINT\n")
        file.write( vals_to_str([xmin,ymin,zero]) )
        file.write( vals_to_str([xmax,ymin,zero]) )
        file.write( vals_to_str([xmin,ymax,zero]) )
        file.write( vals_to_str([xmax,ymax,zero]) )
        #
        file.write(f"Zone T=\"XZ-plane\" I=2, J=2, F=POINT\n")
        file.write( vals_to_str([xmin,zero,zmin]) )
        file.write( vals_to_str([xmax,zero,zmin]) )
        file.write( vals_to_str([xmin,zero,zmax]) )
        file.write( vals_to_str([xmax,zero,zmax]) )
        #
        file.write(f"Zone T=\"YZ-plane\" I=2, J=2, F=POINT\n")
        file.write( vals_to_str([zero,ymin,zmin]) )
        file.write( vals_to_str([zero,ymax,zmin]) )
        file.write( vals_to_str([zero,ymin,zmax]) )
        file.write( vals_to_str([zero,ymax,zmax]) )
        #
        file.write(f"Zone T=\"V projections on Xmax\" I=2, J=2, F=POINT\n")
        file.write( vals_to_str([u,zero,zero]) )
        file.write( vals_to_str([u,v   ,zero]) )
        file.write( vals_to_str([u,zero,w   ]) )
        file.write( vals_to_str([u,v   ,w   ]) )
        #
        file.write(f"Zone T=\"X-axis\" F=POINT\n")
        file.write( vals_to_str(vecZero ) )
        file.write( vals_to_str(vecXaxis) )
        #
        file.write(f"Zone T=\"Y-axis\" F=POINT\n")
        file.write( vals_to_str(vecZero ) )
        file.write( vals_to_str(vecYaxis) )
        #
        file.write(f"Zone T=\"Z-axis\" F=POINT\n")
        file.write( vals_to_str(vecZero ) )
        file.write( vals_to_str(vecZaxis) )
        #
        file.write(f"Zone T=\"V\" F=POINT\n")
        file.write( vals_to_str(vecZero ) )
        file.write( vals_to_str(vecV    ) )
        #
        file.write(f"Zone T=\"Vxy\" F=POINT\n")
        file.write( vals_to_str(vecZero ) )
        file.write( vals_to_str(vecVxy  ) )
        #
        file.write(f"Zone T=\"Vxz\" F=POINT\n")
        file.write( vals_to_str(vecZero ) )
        file.write( vals_to_str(vecVxz  ) )
        #
        file.write(f"Zone T=\"Vyz\" F=POINT\n")
        file.write( vals_to_str(vecZero ) )
        file.write( vals_to_str(vecVyz  ) )
        #
        file.write(f"Zone T=\"Vx\" F=POINT\n")
        file.write( vals_to_str(vecZero ) )
        file.write( vals_to_str(vecVx   ) )
        #
        file.write(f"Zone T=\"Vy\" F=POINT\n")
        file.write( vals_to_str(vecZero ) )
        file.write( vals_to_str(vecVy   ) )
        #
        file.write(f"Zone T=\"Vz\" F=POINT\n")
        file.write( vals_to_str(vecZero ) )
        file.write( vals_to_str(vecVz   ) )
        #
        nmax = 21
        rad = Vmag * 0.8
        file.write(f"Zone T=\"Alpha\" F=POINT\n")
        for n in range(1,nmax+1):
            tmpa = alpha + (0.-alpha) * (n-1) / (nmax-1)
            tmpb = zero
            tmpu = rad * math.cos(tmpb) * math.cos(tmpa)
            tmpv = rad * math.cos(tmpb) * math.sin(tmpa)
            tmpw = rad * math.sin(tmpb)
            tmpV = np.array([tmpu,tmpv,tmpw])
            file.write( vals_to_str(tmpV) )
        #
        file.write(f"Zone T=\"Beta\" F=POINT\n")
        for n in range(1,nmax+1):
            tmpa = alpha
            tmpb = beta  + (0.-beta) * (n-1) / (nmax-1)
            tmpu = rad * math.cos(tmpb) * math.cos(tmpa)
            tmpv = rad * math.cos(tmpb) * math.sin(tmpa)
            tmpw = rad * math.sin(tmpb)
            tmpV = np.array([tmpu,tmpv,tmpw])
            file.write( vals_to_str(tmpV) )
        #
        file.write(f"Zone T=\"V in XY-plane\" I=2, J=2, F=POINT\n")
        file.write( vals_to_str([u   ,v   ,zero]) )
        file.write( vals_to_str([xmax,v   ,zero]) )
        file.write( vals_to_str([u   ,zero,zero]) )
        file.write( vals_to_str([xmax,zero,zero]) )


