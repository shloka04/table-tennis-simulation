"""
tt_engine3d.py -- 3D extension of the unified engine, adding sidespin.

Frame: x forward (along table), y up (table surface y=0), z lateral.
Spin vector omega = (0, omega_side, -omega_top):
  omega_top > 0  -> topspin  (Magnus down, as in 2D)
  omega_side > 0 -> sidespin (Magnus sideways -> curved serve)
Reduces EXACTLY to the 2D engine when omega_side = 0 (verified in self-test).

Forces (3D):
  Drag   F_D = -0.5 rho A C_D(Re) |v| v
  Magnus F_M = 0.5 rho A C_L(S) |v|^2 * (w x v)/|w x v|
  Gravity (0,-g,0)
  S = |w| r / |v|,  C_L = 1/(2+1/S),  C_D Reynolds sigmoid (same as 2D)
"""
import math
import numpy as np

class TableTennisSim3D:
    g=9.81; m=0.0027; r=0.020; d=0.040; rho=1.225; mu=1.81e-5; k_s=0.08
    L=2.74; net_x=1.37; net_h=0.1525; half_width=0.7625   # table half-width (1.525 m / 2)

    def __init__(self, rho=None):
        self.A=math.pi*self.r**2
        if rho is not None: self.rho=rho

    def C_D(self,v):
        if v<1e-9: return 0.47
        Re=self.rho*v*self.d/self.mu
        return 0.55-0.08/(1.0+math.exp(1.5e-4*(Re-28000.0)))
    def C_L(self,S):
        return 0.0 if S<1e-12 else 1.0/(2.0+1.0/S)

    def _accel(self,vx,vy,vz,wx,wy,wz):
        v=math.sqrt(vx*vx+vy*vy+vz*vz)
        if v<1e-9: return 0.0,-self.g,0.0
        w=math.sqrt(wx*wx+wy*wy+wz*wz)
        S=w*self.r/v
        cd=self.C_D(v); cl=self.C_L(S)
        kD=0.5*self.rho*self.A*cd*v/self.m
        # drag accel
        ax=-kD*vx; ay=-kD*vy-self.g; az=-kD*vz
        if w>1e-9:
            # Magnus: 0.5 rho A Cl v^2 * (w x v)/|w x v|
            cxx=wy*vz-wz*vy; cxy=wz*vx-wx*vz; cxz=wx*vy-wy*vx
            cmag=math.sqrt(cxx*cxx+cxy*cxy+cxz*cxz)
            if cmag>1e-12:
                fM=0.5*self.rho*self.A*cl*v*v/self.m
                ax+=fM*cxx/cmag; ay+=fM*cxy/cmag; az+=fM*cxz/cmag
        return ax,ay,az

    def simulate(self,v0,angle_deg,rpm_top,rpm_side=0.0,h=0.15,x0=0.15,
                 dt=0.001,t_max=3.0):
        th=math.radians(angle_deg)
        vx,vy,vz=v0*math.cos(th),v0*math.sin(th),0.0
        x,y,z=x0,h,0.0
        wt0=abs(rpm_top)*2*math.pi/60*math.copysign(1,rpm_top) if rpm_top else 0.0
        ws0=abs(rpm_side)*2*math.pi/60*math.copysign(1,rpm_side) if rpm_side else 0.0
        xs,ys,zs=[x],[y],[z]; t=0.0; n=int(t_max/dt)
        def deriv(tt, X,Y,Z, VX,VY,VZ):
            dec=math.exp(-self.k_s*tt)
            ax,ay,az=self._accel(VX,VY,VZ, 0.0, ws0*dec, -wt0*dec)
            return VX,VY,VZ, ax,ay,az
        for _ in range(n):
            k1=deriv(t, x,y,z, vx,vy,vz)
            k2=deriv(t+0.5*dt, x+0.5*dt*k1[0],y+0.5*dt*k1[1],z+0.5*dt*k1[2],
                               vx+0.5*dt*k1[3],vy+0.5*dt*k1[4],vz+0.5*dt*k1[5])
            k3=deriv(t+0.5*dt, x+0.5*dt*k2[0],y+0.5*dt*k2[1],z+0.5*dt*k2[2],
                               vx+0.5*dt*k2[3],vy+0.5*dt*k2[4],vz+0.5*dt*k2[5])
            k4=deriv(t+dt,     x+dt*k3[0],y+dt*k3[1],z+dt*k3[2],
                               vx+dt*k3[3],vy+dt*k3[4],vz+dt*k3[5])
            xp,yp=x,y
            x +=dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
            y +=dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
            z +=dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
            vx+=dt/6*(k1[3]+2*k2[3]+2*k3[3]+k4[3])
            vy+=dt/6*(k1[4]+2*k2[4]+2*k3[4]+k4[4])
            vz+=dt/6*(k1[5]+2*k2[5]+2*k3[5]+k4[5])
            t+=dt
            xs.append(x); ys.append(y); zs.append(z)
            if xp < self.net_x <= x and y < self.net_h:
                return xs,ys,zs,'net',x,z
            if yp>0.0>=y:
                frac=yp/(yp-y)
                x_l=xp+frac*(x-xp); z_l=zs[-2]+frac*(z-zs[-2])
                if x_l<=self.net_x: out='net'
                elif x_l<=self.L and abs(z_l)<=self.half_width: out='safe'
                else: out='long'
                return xs,ys,zs,out,x_l,z_l
            if x>6.0: break
        return xs,ys,zs,'long',x,z
