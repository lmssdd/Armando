import wx
import os
from numpy import *

def det3x3(mat):
    res = 0.
    res += mat[0][0] * mat[1][1] * mat[2][2]
    res += mat[0][1] * mat[1][2] * mat[2][0]
    res += mat[0][2] * mat[1][0] * mat[2][1]
    res -= mat[0][2] * mat[1][1] * mat[2][0]
    res -= mat[0][0] * mat[1][2] * mat[2][1]
    res -= mat[0][1] * mat[1][0] * mat[2][2]
    return res

class VectorEntryDialog(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(250, 210))
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_vx = wx.TextCtrl(self, -1, "0", size=(125, -1))
        self.tc_vy = wx.TextCtrl(self, -1, "0", size=(125, -1))
        self.tc_vz = wx.TextCtrl(self, -1, "0", size=(125, -1))
        
        sizer.Add(self.tc_vx, 0, \
            wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        sizer.Add(self.tc_vy, 0, \
            wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        sizer.Add(self.tc_vz, 0, \
            wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        
        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        
        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetHelpText("The OK button completes the dialog")
        btn.SetDefault()
        btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_CANCEL)
        btn.SetHelpText("The Cancel button cancels the dialog. (Cool, huh?)")
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        
        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        
        self.SetSizer(sizer)
        sizer.Fit(self)
    
    def GetValue(self):
        out = []
        out.append(self.tc_vx.GetValue())
        out.append(self.tc_vy.GetValue())
        out.append(self.tc_vz.GetValue())
        return out
    
    def SetValue(self, fields):
        self.tc_vx.SetValue(fields[0])
        self.tc_vy.SetValue(fields[1])
        self.tc_vz.SetValue(fields[2])

 
class FemData():
    """
    This class stores the data realtive to a Finite Element Grid,
    used to define the geometry of the model
    """
    def __init__(self):
        self.node = []
        self.tetra = []
        self.tria = []
    
    def Clean(self):
        self.node = []
        self.tetra = []
        self.tria = []
    
    def ReadFile(self, fileName):
        """Open and process the finite element grid"""
        
        self.node = []
        self.tetra = []
        self.tria = []
        
        Status = True
        
        try:
            fp = file(fileName, 'r')
            
            line = fp.readline()
            while line:
                if line.startswith('GRID'):
                    words = line.split(',')
                    for w in words:
                        w = w.strip()
                        if (len(w) == 0):
                            words.remove(w)
                    
                    coords = []
                    coords.append(float(words[-3]))
                    coords.append(float(words[-2]))
                    coords.append(float(words[-1]))
                    
                    self.node.append(coords)
            
                if line.startswith('CTETRA'):
                    words = line.split(',')
                    for w in words:
                        w = w.strip()
                        if (len(w) == 0):
                            words.remove(w)
                    
                    connectivity = []
                    connectivity.append(int(words[-4]))
                    connectivity.append(int(words[-3]))
                    connectivity.append(int(words[-2]))
                    connectivity.append(int(words[-1]))
                    material = int(words[-5])
                    
                    self.tetra.append(connectivity)
                
                if line.startswith('CTRIA'):
                    words = line.split(',')
                    for w in words:
                        w = w.strip()
                        if (len(w) == 0):
                            words.remove(w)
                    
                    connectivity = []
                    connectivity.append(int(words[3]))
                    connectivity.append(int(words[4]))
                    connectivity.append(int(words[5]))
                    material = int(words[2])
                    
                    self.tria.append(connectivity)
                
##                if line.startswith('CTRIA3'):
##                    words = line.split(' ')
##                    for w in words:
##                        w = w.strip()
##                        if (len(w) == 0):
##                            words.remove(w)
##                    print words
##                    connectivity = []
##                    connectivity.append(int(words[3]))
##                    connectivity.append(int(words[4]))
##                    connectivity.append(int(words[5]))
##                    material = int(words[2])
##                    
##                    self.tria.append(connectivity)
##                
                line = fp.readline()
            
            fp.close()
            
        except IOError:
            pass
            Status = False
        
        return Status
    
    def MeshTetra(self, d):
        """Reads a FEM grid and puts the particles inside it"""
        
        particles = []
        
        # Meshing inside the tetrahedrons 
        xmin =  1.0E30
        ymin =  1.0E30
        zmin =  1.0E30
        xmax = -1.0E30
        ymax = -1.0E30
        zmax = -1.0E30
        
        for coords in self.node:
            if coords[0] < xmin :
                xmin = coords[0]
            if coords[1] < ymin :
                ymin = coords[1]
            if coords[2] < zmin :
                zmin = coords[2]
            if coords[0] > xmax :
                xmax = coords[0]
            if coords[1] > ymax :
                ymax = coords[1]
            if coords[2] > zmax :
                zmax = coords[2]
        
        print "Creating cartesian grid \n"
        print "Limits: \n"
        print xmin, ymin, zmin
        print xmax, ymax, zmax
        
        nx = int((xmax - xmin) / d)
        ny = int((ymax - ymin) / d)
        nz = int((zmax - zmin) / d)
        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny
        dz = (zmax - zmin) / nz
        
        print "(nx,ny,nz) = (%d,%d,%d) \n" \
            % (nx, ny, nz)
        
        particleGrid = zeros((nx,ny,nz), dtype = int)
        
        # http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
        
        for el in self.tetra:
            n1 = self.node[el[0] -1]
            n2 = self.node[el[1] -1]
            n3 = self.node[el[2] -1]
            n4 = self.node[el[3] -1]
            
            elm = array([n1,n2,n3,n4])
            elxmin = min(elm[:,0])
            elymin = min(elm[:,1])
            elzmin = min(elm[:,2])
            elxmax = max(elm[:,0])
            elymax = max(elm[:,1])
            elzmax = max(elm[:,2])
            
            imin = max(int((elxmin - xmin) / dx), 0)
            jmin = max(int((elymin - ymin) / dy), 0)
            kmin = max(int((elzmin - zmin) / dz), 0)
            imax = min(int((elxmax - xmin) / dx), nx -1)
            jmax = min(int((elymax - ymin) / dy), ny -1)
            kmax = min(int((elzmax - zmin) / dz), nz -1)
            
            for k in range(kmin,kmax +1):
                for j in range(jmin,jmax +1):
                    for i in range(imin,imax +1):
                        xp = xmin + dx*(i +0.5)
                        yp = ymin + dy*(j +0.5)
                        zp = zmin + dz*(k +0.5)
                        np = [xp, yp, zp]
                        
                        d0 = det3x3([n1,n2,n3]) - det3x3([n1,n2,n4]) \
                            + det3x3([n1,n3,n4]) - det3x3([n2,n3,n4])
                        
                        d1 = det3x3([np,n2,n3]) - det3x3([np,n2,n4]) \
                            + det3x3([np,n3,n4]) - det3x3([n2,n3,n4])
                        
                        d2 = det3x3([n1,np,n3]) - det3x3([n1,np,n4]) \
                            + det3x3([n1,n3,n4]) - det3x3([np,n3,n4])
                        
                        d3 = det3x3([n1,n2,np]) - det3x3([n1,n2,n4]) \
                            + det3x3([n1,np,n4]) - det3x3([n2,np,n4])
                        
                        d4 = det3x3([n1,n2,n3]) - det3x3([n1,n2,np]) \
                            + det3x3([n1,n3,np]) - det3x3([n2,n3,np])
                        
                        if d0 < 0:
                            d0 *= -1.
                            d1 *= -1.
                            d2 *= -1.
                            d3 *= -1.
                            d4 *= -1.
                        
                        eps = 0.001 * d0
                        
                        if ((d1 >= -eps) and (d2 >= -eps) and \
                            (d3 >= -eps) and (d4 >= -eps)):
                            particleGrid[i,j,k] = 1
                
        
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if particleGrid[i,j,k] <> 0:
                        xp = xmin + dx*(i +0.5)
                        yp = ymin + dy*(j +0.5)
                        zp = zmin + dz*(k +0.5)
                        if (particleGrid[i,j,k] != 0):
                            particles.append([xp, yp, zp])
        
        return particles
    
    
    def MeshTria(self):
        """Reads a FEM grid and puts the particles inside it"""
        
        pos = []
        nor = []
        
        for el in self.tria:
            n1 = self.node[el[0] -1]
            n2 = self.node[el[1] -1]
            n3 = self.node[el[2] -1]
            
            elm = array([n1,n2,n3])
            
            xp = sum(elm[:,0]) / 3.0
            yp = sum(elm[:,1]) / 3.0
            zp = sum(elm[:,2]) / 3.0
            
            np = cross((elm[1,:] -elm[0,:]), (elm[2,:] -elm[1,:]))
            (nx,ny,nz) = np / linalg.norm(np)
            
            pos.append([xp, yp, zp])
            nor.append([nx, ny, nz])
        
        return pos, nor
    

class SphData():
    """
    This class stores the data relative to particles
    """
    def __init__(self):
        self.X = []
        self.V = []
        self.material = []
        self.h = []
        self.mass = []
        self.density = []
        self.pressure = []
        self.energy = []
        self.np = 0
        self.nv = 0

    def Clean(self):
        self.X = []
        self.V = []
        self.material = []
        self.h = []
        self.mass = []
        self.density = []
        self.pressure = []
        self.energy = []
        self.np = 0
        self.nv = 0
        
        return None
    
    def AddParticle(self, position, velocity,\
        material, h, mass, density, pressure, energy):
        
        self.X.append(position)
        self.V.append(velocity)
        self.material.append(material)
        self.h.append(h)
        self.mass.append(mass)
        self.density.append(density)
        self.pressure.append(pressure)
        self.energy.append(energy)
        if material > 0 :
            self.np += 1
        else:
            self.nv += 1
        
        return None
    
    def GetParticle(self, index):
        if (index >= 0 and index <= self.np +self.nv):
            return (self.X[index], self.V[index], \
                self.material[index], self.h[index], self.mass[index],\
                 self.density[index], self.pressure[index], self.energy[index])
        else:
            return None

    def GetParticleNumber(self):
        
        return (self.np, self.nv)

    def SetVectors(self, position=[], velocity=[], material=[], h=[],\
        mass=[], density=[], pressure=[], energy=[]):
        """Set the particles data"""
        
        self.X = position
        self.V = velocity
        self.material = material
        self.h = h
        self.mass = mass
        self.density = density
        self.pressure = pressure
        self.energy = energy
        
        return None
        
    def GetVectors(self):
        """Get the particles data"""
        
        return (self.X, self.V, self.material, self.h, self.mass,\
            self.density, self.pressure, self.energy)
        
    def ReadFile(self, dirName):
        """Open and process the particle data files"""
        
        self.X = []
        self.V = []
        self.material = []
        self.h = []
        self.mass = []
        self.density = []
        self.pressure = []
        self.energy = []
        self.np = 0
        self.nv = 0
        
        Error = False
        
        try:
            fileName = dirName + '/f_xv.dat'
            fp_xv = file(fileName, 'r')
            
            fileName = dirName + '/f_state.dat'
            fp_state = file(fileName, 'r')
            
            fileName = dirName + '/f_other.dat'
            fp_other = file(fileName, 'r')
            
            line = fp_xv.readline()
            np = int(line)
            
            for p in range(np):
                line = fp_xv.readline()
                tokens = line.split()
                
                xp = []
                xp.append(float(tokens[1]))
                xp.append(float(tokens[2]))
                xp.append(float(tokens[3]))
                
                vp = []
                vp.append(float(tokens[4]))
                vp.append(float(tokens[5]))
                vp.append(float(tokens[6]))
                
                line = fp_other.readline()
                tokens = line.split()
                
                matp = int(tokens[1])
                hp = float(tokens[2])
                
                line = fp_state.readline()
                tokens = line.split()
                
                massp = float(tokens[1])
                rhop = float(tokens[2])
                pp = float(tokens[3])
                ep = float(tokens[4])
                
                self.AddParticle(xp, vp, matp, hp, massp, rhop, pp, ep)
            
            fp_xv.close()
            fp_state.close()
            fp_other.close()
            
        except IOError:
            Error = True
            pass
        
        try:
            fileName = dirName + '/v_xv.dat'
            fp_xv = file(fileName, 'r')
            
            fileName = dirName + '/v_state.dat'
            fp_state = file(fileName, 'r')
            
            fileName = dirName + '/v_other.dat'
            fp_other = file(fileName, 'r')
            
            line = fp_xv.readline()
            nv = int(line)
            
            for p in range(nv):
                line = fp_xv.readline()
                tokens = line.split()
                
                xp = []
                xp.append(float(tokens[1]))
                xp.append(float(tokens[2]))
                xp.append(float(tokens[3]))
                
                vp = []
                vp.append(float(tokens[4]))
                vp.append(float(tokens[5]))
                vp.append(float(tokens[6]))
                
                line = fp_other.readline()
                tokens = line.split()
                
                matp = int(tokens[1])
                hp = float(tokens[2])
                
                line = fp_state.readline()
                tokens = line.split()
                
                massp = float(tokens[1])
                rhop = float(tokens[2])
                pp = float(tokens[3])
                ep = float(tokens[4])
                
                self.AddParticle(xp, vp, matp, hp, massp, rhop, pp, ep)
            
            fp_xv.close()
            fp_state.close()
            fp_other.close()
            
        except IOError:
            Error = True
            pass
        
        return Error
    
    def SaveFile(self,dirName):
        """Saves the particle data files"""
        
        Error = False
        
        np, nv = self.GetParticleNumber()
        
        try:
            fileName = dirName + '/f_xv.dat'
            fp_xv = file(fileName, 'w')
            
            fileName = dirName + '/f_state.dat'
            fp_state = file(fileName, 'w')
            
            fileName = dirName + '/f_other.dat'
            fp_other = file(fileName, 'w')
            
            fp_xv.write('%8i \n' % np)
            
            for p in range(np):
                xp, vp, matp, hp, massp, rhop, pp, ep = \
                    self.GetParticle(p)
                
                fp_xv.write('%8i %+14.8e %+14.8e %+14.8e ' \
                    % (p +1, xp[0], xp[1], xp[2]))
                fp_xv.write('%+14.8e %+14.8e %+14.8e \n' \
                    % (vp[0], vp[1], vp[2]))
                
                fp_other.write('%8i %5i %+14.8e \n' % (p +1, matp, hp))
            
                fp_state.write('%8i %+14.8e %+14.8e %+14.8e %+14.8e \n' \
                    % (p +1, massp, rhop, pp, ep))
                
            fp_xv.close()
            fp_state.close()
            fp_other.close()
            
        except IOError:
            Error = True
            pass
        
        try:
            fileName = dirName + '/v_xv.dat'
            fp_xv = file(fileName, 'w')
            
            fileName = dirName + '/v_state.dat'
            fp_state = file(fileName, 'w')
            
            fileName = dirName + '/v_other.dat'
            fp_other = file(fileName, 'w')
            
            fp_xv.write('%8i \n' % nv)
            
            for p in range(np, np +nv):
                xp, vp, matp, hp, massp, rhop, pp, ep = \
                    self.GetParticle(p)
                
                fp_xv.write('%8i %+14.8e %+14.8e %+14.8e ' \
                    % (p +1, xp[0], xp[1], xp[2]))
                fp_xv.write('%+14.8e %+14.8e %+14.8e \n' \
                    % (vp[0], vp[1], vp[2]))
                
                fp_other.write('%8i %5i %+14.8e \n' % (p +1, matp, hp))
            
                fp_state.write('%8i %+14.8e %+14.8e %+14.8e %+14.8e \n' \
                    % (p +1, massp, rhop, pp, ep))
            
            fp_xv.close()
            fp_state.close()
            fp_other.close()
            
        except IOError:
            Error = True
            pass
        
        return Error
        

class SimulationData():
    """
    This class stores the data relative to the simulation parameters
    """
    def __init__(self):
        self.pa_sph = 0
        self.nnps = 0
        self.sle = 0
        self.skf = 0
        self.dens_method = 0
        self.avg_vel = False
        self.v_part = False
        self.visc = False
        self.art_visc = False
        self.art_heat = False
        self.ext_heat = False
        self.gravity = False
        
        self.steps = 0
        self.timestep = 0.0
        self.refresh = 0
        self.savestep = 0
        
        self.trfluka = [0.0, 0.0, 0.0]
        self.sffluka = 0.0
        self.effluka = 0.0
        self.npfluka = 0.0
        self.dtfluka = 0.0
        
        self.matn = 1
        self.matprop = 16 * [0.0]
    
    
    def SetOptions(self, pa_sph=0, nnps=0, sle=0, skf=0, dens_method=0, 
        avg_vel=False, v_part=False, visc=False, art_visc=False, 
        art_heat=False, ext_heat=False, gravity=False):
        """Set the parameters for the project options"""
        
        self.pa_sph = pa_sph
        self.nnps = nnps
        self.sle = sle
        self.skf = skf
        self.dens_method = dens_method
        self.avg_vel = avg_vel
        self.v_part = v_part
        self.visc = visc
        self.art_visc = art_visc
        self.art_heat = art_heat
        self.ext_heat = ext_heat
        self.gravity = gravity
    
    def GetOptions(self):
        """Get the parameters for the project options"""
        
        return (self.pa_sph, self.nnps, self.sle, self.skf, \
            self.dens_method, self.avg_vel, self.v_part, self.visc, \
            self.art_visc, self.art_heat, self.ext_heat, self.gravity)

    def SetTiming(self, steps=0, timestep=0.0, refresh=1, savestep=1):
        """Set the parameters for the project timing"""
        
        self.steps = steps
        self.timestep = timestep
        self.refresh = refresh
        self.savestep = savestep
        
    def GetTiming(self):
        """Get the parameters for the project timing"""
        
        return self.steps, self.timestep, self.refresh, self.savestep
    
    def SetFluka(self, trfluka=[0.0, 0.0, 0.0], sffluka=0.0, \
        effluka=0.0, npfluka=0.0, dtfluka=0.0):
        """Set the parameters for the fluka input"""
        
        self.trfluka = trfluka
        self.sffluka = sffluka
        self.effluka = effluka
        self.npfluka = npfluka
        self.dtfluka = dtfluka
        
    def GetFluka(self):
        """get the parameters for the fluka input"""
        return (self.trfluka, self.sffluka, self.effluka, \
            self.npfluka, self.dtfluka)
    
    def SetMaterial(self, matn=1, matprop=16*[0.0]):
        """Set the parameters for the material"""
        
        self.matn = matn
        self.matprop = matprop
    
    def GetMaterial(self):
        """get the parameters for the material"""
        return (self.matn, self.matprop)
    
    def GetDensity(self):
        """get the material density"""
        return self.matprop[0]
    
    def ReadFile(self, dirName):
        """Open and process a project file."""
        
        pa_sph, nnps, sle, skf, dens_method, avg_vel, v_part, visc, art_visc, \
            art_heat, ext_heat, gravity = self.GetOptions()
        steps, timestep, refresh, savestep = self.GetTiming()
        trfluka, sffluka, effluka, npfluka, dtfluka = self.GetFluka()
        matn, matprop = self.GetMaterial()
        
        Error = False
        
        try:
            fileName = dirName + '/input.sph'
            fp = file(fileName, 'r')

            line = fp.readline()
            while len(line) > 0:
                tok = line.split()

                if tok[0] == 'PA_SPH':
                    pa_sph = int(tok[1])
                elif tok[0] == 'NNPS':
                    nnps = int(tok[1])
                elif tok[0] == 'SLE':
                    sle = int(tok[1])
                elif tok[0] == 'SKF':
                    skf = int(tok[1])
                elif tok[0] == 'DENS_METHOD':
                    dens_method = int(tok[1])
                elif tok[0] == 'AVG_VEL':
                    avg_vel = bool(int(tok[1]))
                elif tok[0] == 'V_PART':
                    v_part = bool(int(tok[1]))
                elif tok[0] == 'VISC':
                    visc = bool(int(tok[1]))
                elif tok[0] == 'ART_VISC':
                    art_visc = bool(int(tok[1]))
                elif tok[0] == 'ART_HEAT':
                    art_heat = bool(int(tok[1]))
                elif tok[0] == 'EXT_HEAT':
                    ext_heat = bool(int(tok[1]))
                elif tok[0] == 'GRAVITY':
                    gravity = bool(int(tok[1]))
                
                elif tok[0] == 'STEPS':
                    steps = int(tok[1])
                elif tok[0] == 'TIMESTEP':
                    timestep = float(tok[1])
                elif tok[0] == 'REFRESH':
                    refresh = int(tok[1])
                elif tok[0] == 'SAVESTEP':
                    savestep = int(tok[1])

                elif tok[0] == 'TRFLUKA':
                    trfluka = [float(tok[1]), float(tok[2]), float(tok[3])]
                elif tok[0] == 'SFFLUKA':
                    sffluka = float(tok[1])
                elif tok[0] == 'EFFLUKA':
                    effluka = float(tok[1])
                elif tok[0] == 'NPFLUKA':
                    npfluka = float(tok[1])
                elif tok[0] == 'DTFLUKA':
                    dtfluka = float(tok[1])
                
                elif tok[0] == 'MAT':
                    matn = int(tok[1])
                    ip1 = int(tok[2])
                    ip2 = int(tok[3])
                    for i in range(ip2 -ip1 +1):
                        matprop[ip1 -1 +i] = float(tok[4 +i])
                
                line = fp.readline()
            
            fp.close()
            
        except IOError:
            Error = True
            pass
            
        self.SetOptions(pa_sph, nnps, sle, skf, dens_method, avg_vel, \
             v_part, visc, art_visc, art_heat, ext_heat, gravity)
        self.SetTiming(steps, timestep, refresh, savestep)
        self.SetFluka(trfluka, sffluka, effluka, npfluka, dtfluka)
        self.SetMaterial(matn, matprop)
        
        return Error
    
    def SaveFile(self,dirName):
        """Saves a project file."""
        
        pa_sph, nnps, sle, skf, dens_method, avg_vel, v_part, visc, art_visc, \
            art_heat, ext_heat, gravity = self.GetOptions()
        steps, timestep, refresh, savestep = self.GetTiming()
        trfluka, sffluka, effluka, npfluka, dtfluka = self.GetFluka()
        matn, matprop = self.GetMaterial()
        
        Error = False
        
        try:
            fileName = dirName + '/input.sph'
            fp = file(fileName, 'w')
            
            fp.write('PA_SPH %i \n' % pa_sph)
            fp.write('NNPS %i \n' % nnps)
            fp.write('SLE %i \n' % sle)
            fp.write('SKF %i \n' % skf)
            fp.write('DENS_METHOD %i \n' % dens_method)
            
            fp.write('AVG_VEL %i \n' % avg_vel)
            fp.write('V_PART %i \n' % v_part)
            fp.write('VISC %i \n' % visc)
            fp.write('ART_VISC %i \n' % art_visc)
            fp.write('ART_HEAT %i \n' % art_heat)
            fp.write('EXT_HEAT %i \n' % ext_heat)
            fp.write('GRAVITY %i \n' % gravity)
            
            fp.write('TIMESTEP %e \n' % timestep)
            fp.write('STEPS %i \n' % steps)
            fp.write('REFRESH %i \n' % refresh)
            fp.write('SAVESTEP %i \n' % savestep)
            
            fp.write('TRFLUKA %10.4e %10.4e %10.4e \n' % \
                (trfluka[0], trfluka[1], trfluka[2]))
            fp.write('SFFLUKA %10.4e \n' % sffluka)
            fp.write('EFFLUKA %10.4e \n' % effluka)
            fp.write('NPFLUKA %10.4e \n' % npfluka)
            fp.write('DTFLUKA %10.4e \n' % dtfluka)
            
            if (matn > 0) and (matn <= 10):
                fp.write('MAT %d %d %d '% (matn, 1, 4))
                for i in range(0,4):
                    fp.write('%10.4e '% matprop[i])
                fp.write('\n')
            
            elif (matn > 10) and (matn <= 20):
                fp.write('MAT %d %d %d '% (matn, 1, 4))
                for i in range(0,4):
                    fp.write('%10.4e '% matprop[i])
                fp.write('\n')
                fp.write('MAT %d %d %d '% (matn, 5, 8))
                for i in range(4,8):
                    fp.write('%10.4e '% matprop[i])
                fp.write('\n')
                fp.write('MAT %d %d %d '% (matn, 9, 9))
                for i in range(8,9):
                    fp.write('%10.4e '% matprop[i])
                fp.write('\n')
            
            elif (matn > 20) and (matn <= 30):
                fp.write('MAT %d %d %d '% (matn, 1, 4))
                for i in range(0,4):
                    fp.write('%10.4e '% matprop[i])
                fp.write('\n')
                fp.write('MAT %d %d %d '% (matn, 5, 6))
                for i in range(4,6):
                    fp.write('%10.4e '% matprop[i])
                fp.write('\n')
            
            elif (matn > 30) and (matn <= 40):
                fp.write('MAT %d %d %d '% (matn, 1, 4))
                for i in range(0,4):
                    fp.write('%10.4e '% matprop[i])
                fp.write('\n')
                fp.write('MAT %d %d %d '% (matn, 5, 8))
                for i in range(4,8):
                    fp.write('%10.4e '% matprop[i])
                fp.write('\n')
                fp.write('MAT %d %d %d '% (matn, 9, 11))
                for i in range(8,11):
                    fp.write('%10.4e '% matprop[i])
                fp.write('\n')
            
            fp.close()
            
        except IOError:
            Error = True
            pass
           
        return Error


class ParticlesList(wx.ListCtrl):
    def __init__(self, parent, names, values):
        wx.ListCtrl.__init__(
            self, parent, -1, (25, 50), (425, 300), 
            style=wx.LC_REPORT|wx.LC_VIRTUAL|wx.LC_HRULES|wx.LC_VRULES)
        
        self.ClearAll()
        
        self.names = names
        self.values = values
        
        self.InsertColumn(0, "n.")
        self.SetColumnWidth(0, 100)
        for c in range(len(self.names)):
            self.InsertColumn(c +1, self.names[c])
            self.SetColumnWidth(c +1, 100)
            
        self.SetItemCount(len(self.values))
        
        self.attr1 = wx.ListItemAttr()
        self.attr1.SetBackgroundColour("yellow")

        self.attr2 = wx.ListItemAttr()
        self.attr2.SetBackgroundColour("light blue")
        
    
    def OnGetItemText(self, item, col):
        if col == 0:
            line = "%d" % (item +1)
        else:
            row = self.values[item]
            line = "%+10.4e" % row[col -1]
        return line

    def OnGetItemAttr(self, item):
        if item % 3 == 1:
            return self.attr1
        elif item % 3 == 2:
            return self.attr2
        else:
            return None


class ParticlesPanel(wx.Panel):
    def __init__(self, parent, data):
        
        wx.Panel.__init__(self, parent)
        
        self.data = data
        
        self.b_update = wx.Button(self, -1, "Update", \
            (25, 25), (200, -1))
        self.Bind(wx.EVT_BUTTON, self.OnUpdate, self.b_update)
        
        self.list = ParticlesList(self,['X','Y','Z'],[])
        
        self.Syncronize()
    
    def Syncronize(self):
        self.list.names = ['X','Y','Z']
        self.list.values = self.data.X
        self.list.SetItemCount(len(self.data.X))

    def OnUpdate(self, evt):
        self.list.names = ['X','Y','Z']
        self.list.values = self.data.X
        self.list.SetItemCount(len(self.data.X))
    

class OptionsPanel(wx.Panel):
    """
    This is OptionsPanel Class.  It just shows a few controls on a wxPanel,
    and has a simple menu.
    """
  
    def __init__(self, parent, data):
        self.data = data
        
        wx.Panel.__init__(self, parent)
        # Now create the Panel to put the other controls on.
        
        skfList = ['cubic spline kernel by W4 - Spline (Monaghan 1985)',
            'Gauss kernel   (Gingold and Monaghan 1981)',
            'Quintic kernel (Morris 1997)Keep unchanged']
        wx.StaticText(self, -1, "Smoothing kernel function", \
            (25, 25), (200, -1))
        self.ch_skf = wx.Choice(self, -1, (250, 25), (200, -1), choices = skfList)
        self.Bind(wx.EVT_CHOICE, self.OnChoice, self.ch_skf)
        
        sleList = ['Keep unchanged','h = fac * (m/rho)^(1/dim)',
            'dh/dt = (-1/dim)*(h/rho)*(drho/dt)',
            'Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )']
        wx.StaticText(self, -1, "Smoothing length evolution algorithm", \
            (25, 50), (200, -1))
        self.ch_sle = wx.Choice(self, -1, (250, 50), (200, -1), choices = sleList)
        self.Bind(wx.EVT_CHOICE, self.OnChoice, self.ch_sle)
        
        nnpsList = ['Simplest and direct searching','Sorting grid linked list',
            'Tree algorith']
        wx.StaticText(self, -1, "Nearest Neighbour searching method", \
            (25, 75), (200, -1))
        self.ch_nnps = wx.Choice(self, -1, (250, 75), (200, -1), choices = nnpsList)
        self.Bind(wx.EVT_CHOICE, self.OnChoice, self.ch_nnps)
        
        paList = ['(p(i)+p(j))/(rho(i)*rho(j))', 
            '(p(i)/rho(i)**2+p(j)/rho(j)**2)']
        wx.StaticText(self, -1, "Algorithm for particle approximation", \
            (25, 100), (200, -1))
        self.ch_pa = wx.Choice(self, -1, (250, 100), (200, -1), choices = paList)
        self.Bind(wx.EVT_CHOICE, self.OnChoice, self.ch_pa)
        
        densList = ['continuity equation','summation density',
            'normalized sum. density']
        wx.StaticText(self, -1, "Density approximation method", \
            (25, 125), (200, -1))
        self.ch_dens = wx.Choice(self, -1, (250, 125), (200, -1), choices = densList)
        self.Bind(wx.EVT_CHOICE, self.OnChoice, self.ch_dens)
        
        self.cb_avg_vel = wx.CheckBox(self, -1, 
            'Average velocity', (25, 150), (200, -1), style=wx.ALIGN_LEFT)
        self.Bind(wx.EVT_CHECKBOX, self.OnCheckBox, self.cb_avg_vel)
        
        self.cb_v_part = wx.CheckBox(self, -1, 
            'Boundary particles', (25, 175), (200, -1), style=wx.ALIGN_LEFT)
        self.Bind(wx.EVT_CHECKBOX, self.OnCheckBox, self.cb_v_part)
        
        self.cb_visc = wx.CheckBox(self, -1, 
            'Viscosity', (25, 200), (200, -1), style=wx.ALIGN_LEFT)
        self.Bind(wx.EVT_CHECKBOX, self.OnCheckBox, self.cb_visc)
        
        self.cb_art_visc = wx.CheckBox(self, -1, 
            'Artificial viscosity', (25, 225), (200, -1), style=wx.ALIGN_LEFT)
        self.Bind(wx.EVT_CHECKBOX, self.OnCheckBox, self.cb_art_visc)
        
        self.cb_art_heat = wx.CheckBox(self, -1, 
            'Artificial heat', (25, 250), (200, -1), style=wx.ALIGN_LEFT)
        self.Bind(wx.EVT_CHECKBOX, self.OnCheckBox, self.cb_art_heat)
        
        self.cb_ext_heat = wx.CheckBox(self, -1, 
            'External heat', (25, 275), (200, -1), style=wx.ALIGN_LEFT)
        self.Bind(wx.EVT_CHECKBOX, self.OnCheckBox, self.cb_ext_heat)
        
        self.cb_gravity = wx.CheckBox(self, -1, 
            'Gravity', (25, 300), (200, -1), style=wx.ALIGN_LEFT)
        self.Bind(wx.EVT_CHECKBOX, self.OnCheckBox, self.cb_gravity)
        
        self.Syncronize()
        
    
    def Syncronize(self):
        self.ch_skf.Select(self.data.skf)
        self.ch_sle.Select(self.data.sle)
        self.ch_nnps.Select(self.data.nnps)
        self.ch_pa.Select(self.data.pa_sph)
        self.ch_dens.Select(self.data.dens_method)
        self.cb_avg_vel.SetValue(self.data.avg_vel)
        self.cb_v_part.SetValue(self.data.v_part)
        self.cb_visc.SetValue(self.data.visc)
        self.cb_art_visc.SetValue(self.data.art_visc)
        self.cb_art_heat.SetValue(self.data.art_heat)
        self.cb_ext_heat.SetValue(self.data.ext_heat)
        self.cb_gravity.SetValue(self.data.gravity)
    
    def OnChoice(self, evt):
        ch = evt.GetEventObject()
        sn = ch.GetSelection()
        
        if ch == self.ch_skf :
            self.data.skf = sn
        
        if ch == self.ch_sle :
            self.data.sle = sn
        
        if ch == self.ch_nnps :
            self.data.nnps = sn
        
        if ch == self.ch_pa :
            self.data.pa_sph = sn
        
        if ch == self.ch_dens :
            self.data.dens_method = sn
        
        self.Syncronize()
    
    
    def OnCheckBox(self, evt):
        """Event handler for the check box"""
        cb = evt.GetEventObject()
        sn = cb.GetValue()
        
        if cb == self.cb_avg_vel :
            self.data.avg_vel = sn
        
        if cb == self.cb_v_part :
            self.data.v_part = sn
        
        if cb == self.cb_visc :
            self.data.visc = sn
        
        if cb == self.cb_art_visc :
            self.data.art_visc = sn
        
        if cb == self.cb_art_heat :
            self.data.art_heat = sn
        
        if cb == self.cb_ext_heat :
            self.data.ext_heat = sn
        
        if cb == self.cb_gravity :
            self.data.gravity = sn
        
        self.Syncronize()
    

class TimingPanel(wx.Panel):
    """
    This is TimingPanel Class.  It just shows a few controls on a wxPanel,
    and has a simple menu.
    """
  
    def __init__(self, parent, data):
        self.data = data
        
        wx.Panel.__init__(self, parent)
        # Now create the Panel to put the other controls on.
        
        wx.StaticText(self, -1, "Time step", (25, 25), (200, -1))
        self.tc_timestep = wx.TextCtrl(self, -1, "", (250, 25), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_timestep)
        
        wx.StaticText(self, -1, "Steps", (25, 50), (200, -1))
        self.tc_steps = wx.TextCtrl(self, -1, "", (250, 50), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_steps)
        
        wx.StaticText(self, -1, "Final time", (25, 75), (200, -1))
        self.tc_finaltime = wx.TextCtrl(self, -1, "", (250, 75), (200, -1),
            style=wx.TE_READONLY)
        
        wx.StaticText(self, -1, "Refresh steps", (25, 100), (200, -1))
        self.tc_refresh = wx.TextCtrl(self, -1, "", (250, 100), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_refresh)
        
        wx.StaticText(self, -1, "Save steps", (25, 125), (200, -1))
        self.tc_savestep = wx.TextCtrl(self, -1, "",  (250, 125), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_savestep)
        
        self.Syncronize()
    
    def Syncronize(self):
        self.tc_timestep.SetValue('%10.4e' % self.data.timestep)
        self.tc_steps.SetValue('%d' % self.data.steps)
        self.tc_finaltime.SetValue('%10.4e' 
            % (self.data.timestep * self.data.steps))
        self.tc_refresh.SetValue('%d' % self.data.refresh)
        self.tc_savestep.SetValue('%d' % self.data.savestep)
    
    
    def OnEvtText(self, evt):
        """Event handler for the textControl"""
        
        et = evt.GetEventObject()
        sn = et.GetValue()
        
        try:
            if et == self.tc_timestep :
                self.data.timestep = float(sn)
                
            if et == self.tc_steps :
                self.data.steps = int(sn)
            
            if et == self.tc_refresh :
                self.data.refresh = float(sn)
            
            if et == self.tc_savestep :
                self.data.savestep = float(sn)
            
        except:
            pass
        
        self.Syncronize()
    
    
class FlukaPanel(wx.Panel):
    """
    This is FlukaPanel Class.  It just shows a few controls on a wxPanel,
    and has a simple menu.
    """
  
    def __init__(self, parent, data):
        self.data = data
        
        wx.Panel.__init__(self, parent)
        # Now create the Panel to put the other controls on.
        
        wx.StaticText(self, -1, "Size factor", (25, 25), (200, -1))
        self.tc_sffluka = wx.TextCtrl(self, -1, "", (250, 25), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_sffluka)
        
        wx.StaticText(self, -1, "Energy factor", (25, 50), (200, -1))
        self.tc_effluka = wx.TextCtrl(self, -1, "", (250, 50), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_effluka)
        
        wx.StaticText(self, -1, "Number of particles", (25, 75), (200, -1))
        self.tc_npfluka = wx.TextCtrl(self, -1, "", (250, 75), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_npfluka)
        
        wx.StaticText(self, -1, "Deposition time", (25, 100), (200, -1))
        self.tc_dtfluka = wx.TextCtrl(self, -1, "", (250, 100), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_dtfluka)
        
        wx.StaticText(self, -1, "Translation X", (25, 125), (200, -1))
        self.tc_tXfluka = wx.TextCtrl(self, -1, "", (250, 125), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_tXfluka)
        
        wx.StaticText(self, -1, "Translation Y", (25, 150), (200, -1))
        self.tc_tYfluka = wx.TextCtrl(self, -1, "", (250, 150), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_tYfluka)
        
        wx.StaticText(self, -1, "Translation Z", (25, 175), (200, -1))
        self.tc_tZfluka = wx.TextCtrl(self, -1, "", (250, 175), (200, -1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_tZfluka)
        
        self.Syncronize()
        
    def Syncronize(self):
        self.tc_sffluka.SetValue('%10.4e' % self.data.sffluka)
        self.tc_effluka.SetValue('%10.4e' % self.data.effluka)
        self.tc_npfluka.SetValue('%10.4e' % self.data.npfluka)
        self.tc_dtfluka.SetValue('%10.4e' % self.data.dtfluka)
        self.tc_tXfluka.SetValue('%10.4e' % self.data.trfluka[0])
        self.tc_tYfluka.SetValue('%10.4e' % self.data.trfluka[1])
        self.tc_tZfluka.SetValue('%10.4e' % self.data.trfluka[2])
    
    
    def OnEvtText(self, evt):
        """Event handler for the textControl"""
        
        et = evt.GetEventObject()
        sn = et.GetValue()
        
        try:
            if et == self.tc_sffluka :
                self.data.sffluka = float(sn)
                
            if et == self.tc_effluka :
                self.data.effluka = float(sn)
                
            if et == self.tc_npfluka :
                self.data.npfluka = float(sn)
                
            if et == self.tc_dtfluka :
                self.data.dtfluka = float(sn)
            
            if et == self.tc_tXfluka :
                self.data.trfluka[0] = float(sn)
            
            if et == self.tc_tYfluka :
                self.data.trfluka[1] = float(sn)
            
            if et == self.tc_tZfluka :
                self.data.trfluka[2] = float(sn)
            
        except:
            pass
        
        self.Syncronize()
    


class MaterialPanel(wx.Panel):
    """
    This is MaterialPanel Class.  It just shows a few controls on a wxPanel,
    and has a simple menu.
    """
  
    def __init__(self, parent, data):
        self.data = data
        
        self.EOSname = ['Ideal gas', 'Polyynomial', 'Shock', 'Puff']
        
        self.propNameGas = 4 * ['']
        self.propNameGas[0] = 'Density'
        self.propNameGas[1] = 'Gamma'
        self.propNameGas[2] = 'Pressure shift'
        self.propNameGas[3] = 'Viscosity'
        
        self.propNamePoly = 9 * ['']
        self.propNamePoly[0] = 'Density'
        self.propNamePoly[1] = 'A1'
        self.propNamePoly[2] = 'A2'
        self.propNamePoly[3] = 'A3'
        self.propNamePoly[4] = 'B0'
        self.propNamePoly[5] = 'B1'
        self.propNamePoly[6] = 'T1'
        self.propNamePoly[7] = 'Min. Pressure'
        self.propNamePoly[8] = 'Visosity'
        
        self.propNameShock = 6 * ['']
        self.propNameShock[0] = 'Density'
        self.propNameShock[1] = 'C0'
        self.propNameShock[2] = 'Gamma0'
        self.propNameShock[3] = 'S0'
        self.propNameShock[4] = 'Min. Pressure'
        self.propNameShock[5] = 'Viscosity'
        
        self.propNamePuff = 11 * ['']
        self.propNamePuff[0] = 'Density'
        self.propNamePuff[1] = 'A1'
        self.propNamePuff[2] = 'A2'
        self.propNamePuff[3] = 'A3'
        self.propNamePuff[4] = 'Gamma0'
        self.propNamePuff[5] = 'H0'
        self.propNamePuff[6] = 'ES'
        self.propNamePuff[7] = 'T1'
        self.propNamePuff[8] = 'T2'
        self.propNamePuff[9] = 'Min. Pressure'
        self.propNamePuff[10] = 'Visosity'
        
        
        wx.Panel.__init__(self, parent)
        # Now create the Panel to put the other controls on.
        
        wx.StaticText(self, -1, "Material number", (25, 25), (200, -1))
        self.sc_material = wx.SpinCtrl(self, -1, "", (250, 25), (200, -1))
        self.sc_material.SetRange(1,40)
        self.Bind(wx.EVT_SPINCTRL, self.OnEvtSpinCtrl, self.sc_material)
        
        
        wx.StaticText(self, -1, "Equation Of State:", (25, 50), (200, -1))
        self.st_type = wx.StaticText(self, -1, "", (250, 50), (200, -1))
        
        self.st_prop = []
        self.tc_prop = []
        
        for i in range(16):
            self.st_prop.append(wx.StaticText(self, -1, '', 
                (25, 75 +i*25), (200 ,-1)))
            self.tc_prop.append(wx.TextCtrl(self, -1, '', 
                (250, 75 +i*25), (200 ,-1), style=wx.TE_PROCESS_ENTER))
            self.Bind(wx.EVT_TEXT_ENTER, self.OnEvtText, self.tc_prop[i])
        
        self.Syncronize()
        
    
    def Syncronize(self):
        
        self.sc_material.SetValue(self.data.matn)
        
        if ((self.data.matn > 0) and (self.data.matn <= 10)):
            propName = self.propNameGas
            self.st_type.SetLabel(self.EOSname[0])
        
        if ((self.data.matn > 10) and (self.data.matn <= 20)):
            propName = self.propNamePoly
            self.st_type.SetLabel(self.EOSname[1])

        if ((self.data.matn > 20) and (self.data.matn <= 30)):
            propName = self.propNameShock
            self.st_type.SetLabel(self.EOSname[2])
        
        if ((self.data.matn > 30) and (self.data.matn <= 40)):
            propName = self.propNamePuff
            self.st_type.SetLabel(self.EOSname[3])
            
        for i in range(len(propName)):
            self.st_prop[i].SetLabel(propName[i])
            self.tc_prop[i].SetValue('%10.4e' % self.data.matprop[i])
            self.st_prop[i].Show(True)
            self.tc_prop[i].Show(True)
        
        for i in range(len(propName), 16):
            self.st_prop[i].Show(False)
            self.tc_prop[i].Show(False)
        
    
    def OnEvtSpinCtrl(self, evt):
        """Event handler for the spinControl"""
        
        try:
            self.data.matn = self.sc_material.GetValue()
        except:
            pass
        
        self.Syncronize()
    
    def OnEvtText(self, evt):
        """Event handler for the textControl"""
        
        et = evt.GetEventObject()
        sn = et.GetValue()
        
        for i in range(16):
            try:
                if et == self.tc_prop[i]:
                    self.data.matprop[i] = float(sn)
            except:
                pass
        
        self.Syncronize()
    


class MeshingPanel(wx.Panel):
    """
    This is MeshingPanel Class.  It just shows a few controls on a wxPanel,
    and has a simple menu.
    """
  
    def __init__(self, parent, fem, sph, simulation):
        self.fem = fem
        self.sph = sph
        self.simulation = simulation
        
        wx.Panel.__init__(self, parent)
        # Now create the Panel to put the other controls on.
        
        self.b_fem = wx.Button(self, -1, "Read Finite Element Model", \
            (25, 25), (200, -1))
        self.Bind(wx.EVT_BUTTON, self.OnFEButton, self.b_fem)
        
        self.st_fem = wx.StaticText(self, -1, "", (25, 50), (450, -1))
        
        wx.StaticText(self, -1, "Smoothing length h", (25, 75), (200, -1))
        self.tc_h = wx.TextCtrl(self, -1, "", (250, 75), (200, -1), \
            style=wx.TE_PROCESS_ENTER)
        
        wx.StaticText(self, -1, "Particles in h", (25, 100), (200, -1))
        self.tc_k = wx.TextCtrl(self, -1, "", (250, 100), (200, -1), \
            style=wx.TE_PROCESS_ENTER)
        
        self.b_mesh = wx.Button(self, -1, "Mesh Finite Element Model", \
            (25, 125), (200, -1))
        self.Bind(wx.EVT_BUTTON, self.OnMeshButton, self.b_mesh)
        
        self.st_sph = wx.StaticText(self, -1, "", (25, 150), (450, -1))
        
        self.Syncronize()
    
    def Syncronize(self):
        self.st_fem.SetLabel(
            "Nodes: %6d \t Tetra Elements: %6d \t Tria Elements: %6d" 
            % (len(self.fem.node), len(self.fem.tetra), len(self.fem.tria)))
        self.st_sph.SetLabel(
            "Real particles: %6d \t Virtual particles: %6d" 
            % (self.sph.np, self.sph.nv))
    
    def OnFEButton(self, evt):
        """Event handler for the button to read the Finite Element Model"""
        
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultDir=os.getcwd(), 
            defaultFile="",
            wildcard="NASTRAN file (*.nas)|*.nas|" "All files (*.*)|*.*",
            style=wx.OPEN | wx.CHANGE_DIR
            )
        
        if dlg.ShowModal() == wx.ID_OK:
            fileName = dlg.GetPath()
            
            self.fem.ReadFile(fileName)
            
            self.Syncronize()
        
        dlg.Destroy()
        
        
    def OnMeshButton(self, evt):
        """Event handler for the button to read the Finite Element Model"""
        
        h = float(self.tc_h.GetValue())
        k = float(self.tc_k.GetValue())
        d = h / k
        
        self.sph.Clean()
        
        Xr = self.fem.MeshTetra(d)
        
        for xp in Xr:
            vp = [0.0, 0.0, 0.0]
            matp = self.simulation.matn
            hp = h
            rhop = self.simulation.matprop[0]
            massp = rhop * d**3
            pp = 0.0
            ep = 0.0
            
            self.sph.AddParticle(xp, vp, matp, hp, massp, rhop, pp, ep)
        
        
        Xv, Nv = self.fem.MeshTria()
        
        for (xp, vp) in zip(Xv, Nv):
            matp = -1 * self.simulation.matn
            hp = h
            rhop = self.simulation.matprop[0]
            massp = rhop * d**3
            pp = 0.0
            ep = 0.0
            
            self.sph.AddParticle(xp, vp, matp, hp, massp, rhop, pp, ep)
        
        
        self.Syncronize()
    

class MainFrame(wx.Frame):
    """
    This is MyFrame.  It just shows a few controls on a wxPanel,
    and has a simple menu.
    """
    
    def __init__(self, parent, title):
        
        self.fem = FemData()
        self.sph = SphData()
        self.simulation = SimulationData()
        
        wx.Frame.__init__(self, parent, -1, title,
                          pos=wx.DefaultPosition, size=(500, 500),
                          style=wx.DEFAULT_FRAME_STYLE^wx.RESIZE_BORDER)

        # Set up menu bar for the program.
        self.mainmenu = wx.MenuBar()                  # Create menu bar.

        # Make the 'File' menu.
        menu = wx.Menu()

        item = menu.Append(wx.ID_NEW, '&New', 'New project')  # Append a new menu
        self.Bind(wx.EVT_MENU, self.OnProjectNew, item)  # Create and assign a menu event.

        item = menu.Append(wx.ID_OPEN, '&Open', 'Open project')
        self.Bind(wx.EVT_MENU, self.OnProjectOpen, item)

        item = menu.Append(wx.ID_SAVE, '&Save', 'Save project')
        self.Bind(wx.EVT_MENU, self.OnProjectSave, item)

        item = menu.Append(wx.ID_EXIT, 'E&xit', 'Exit program')
        self.Bind(wx.EVT_MENU, self.OnProjectExit, item)

        self.mainmenu.Append(menu, '&File')  # Add the project menu to the menu bar.

        # Make the 'Particles' menu.
        menu = wx.Menu()
        
        item = menu.Append(201, 'Set &material', 'Set the material number')  # Append a new menu
        self.Bind(wx.EVT_MENU, self.OnSetMaterial, item)  # Create and assign a menu event.
        
        item = menu.Append(202, 'Set &mass', 'Set the particle mass')  # Append a new menu
        self.Bind(wx.EVT_MENU, self.OnSetMass, item)  # Create and assign a menu event.
        
        item = menu.Append(203, 'Set &h', 'Set the smoothing lenght')  # Append a new menu
        self.Bind(wx.EVT_MENU, self.OnSetH, item)  # Create and assign a menu event.
        
        item = menu.Append(204, 'Set &velocity', 'Set the initial veolcity')  # Append a new menu
        self.Bind(wx.EVT_MENU, self.OnSetVelocity, item)  # Create and assign a menu event.
        
        item = menu.Append(205, 'Set &density', 'Set the initial density')  # Append a new menu
        self.Bind(wx.EVT_MENU, self.OnSetDensity, item)  # Create and assign a menu event.
        
        item = menu.Append(206, 'Set &pressure', 'Set the initial pressure')  # Append a new menu
        self.Bind(wx.EVT_MENU, self.OnSetPressure, item)  # Create and assign a menu event.
        
        item = menu.Append(207, 'Set &energy', 'Set the initial energy')  # Append a new menu
        self.Bind(wx.EVT_MENU, self.OnSetEnergy, item)  # Create and assign a menu event.
        
        self.mainmenu.Append(menu, '&Particles')  # Add the project menu to the menu bar.
        
        
        # Attach the menu bar to the window.
        self.SetMenuBar(self.mainmenu)
    
        self.CreateStatusBar()
        
        # Here we create a panel and a notebook on the panel
        panelMain = wx.Panel(self)
        notebook = wx.Notebook(panelMain)
        
        # create the page windows as children of the notebook
        self.options = OptionsPanel(notebook, self.simulation)
        self.timing = TimingPanel(notebook, self.simulation)
        self.material = MaterialPanel(notebook, self.simulation)
        self.fluka = FlukaPanel(notebook, self.simulation)
        self.meshing = MeshingPanel(notebook, self.fem, self.sph, 
            self.simulation)
        self.particles = ParticlesPanel(notebook, self.sph)
        
        # add the pages to the notebook with the label to show on the tab
        notebook.AddPage(self.options, "Options")
        notebook.AddPage(self.timing, "Timing")
        notebook.AddPage(self.fluka, "Fluka")
        notebook.AddPage(self.material, "Material")
        notebook.AddPage(self.meshing, "Meshing")
        notebook.AddPage(self.particles, "Particles")
        
        # finally, put the notebook in a sizer for the panel to manage
        # the layout
        sizerMain = wx.BoxSizer()
        sizerMain.Add(notebook, 1, wx.EXPAND)
        panelMain.SetSizer(sizerMain)
        
        # Some global state variables.
        self.projectDirName = ''
        self.projectDirty = False
        

    # ----------------------------------------------------------------------------------------
    # Some nice little handlers.
    # ----------------------------------------------------------------------------------------

    def Syncronize(self):
        self.options.Syncronize()
        self.timing.Syncronize()
        self.fluka.Syncronize()
        self.material.Syncronize()
        self.meshing.Syncronize()
        self.particles.Syncronize()

    
    # ----------------------------------------------------------------------------------------
    # Event handlers from here on out.
    # ----------------------------------------------------------------------------------------
    
    def OnProjectNew(self, event):
        """Create a new project with defaults."""
        self.fem.Clean()
        self.sph.Clean()
        self.simulation.SetOptions()
        self.simulation.SetTiming()
        self.simulation.SetFluka()
        self.simulation.SetMaterial()
        
        self.Syncronize()
    
    def OnProjectOpen(self, event):
        """Open a file."""
        dlg = wx.DirDialog(self, "Choose a directory:",
            style=wx.DD_DEFAULT_STYLE
            #| wx.DD_DIR_MUST_EXIST
            #| wx.DD_CHANGE_DIR
            )
        
        if dlg.ShowModal() == wx.ID_OK:
            self.projectDirName = dlg.GetPath()
            
            if (self.sph.ReadFile(dlg.GetPath())):
                err = wx.MessageDialog(self, 
                'There was an error reading the particle data file', 
                'Error!', wx.OK)
                err.ShowModal()
                err.Destroy()
                    
            if (self.simulation.ReadFile(dlg.GetPath())):
                err = wx.MessageDialog(self, 
                'There was an error reading the simulation data file', 
                'Error!', wx.OK)
                err.ShowModal()
                err.Destroy()
                
            self.Syncronize()
            dlg.Destroy()
    
    def OnProjectSave(self, event):
        dlg = wx.DirDialog(self, "Choose a directory:",
            style=wx.DD_DEFAULT_STYLE
            #| wx.DD_DIR_MUST_EXIST
            #| wx.DD_CHANGE_DIR
            )
        
        if dlg.ShowModal() == wx.ID_OK:
            if (self.sph.SaveFile(dlg.GetPath())):
                err = wx.MessageDialog(self, 
                'There was an error writing the particle data file', 
                'Error!', wx.OK)
                err.ShowModal()
                err.Destroy()
            
            if (self.simulation.SaveFile(dlg.GetPath())):
                err = wx.MessageDialog(self, 
                'There was an error writing the simulation data file', 
                'Error!', wx.OK)
                err.ShowModal()
                err.Destroy()
            
            dlg.Destroy()
    
    def OnProjectExit(self, event):
        """Quit the program."""
        self.Close()
    
    def OnSetH(self, event):
        """Set the smoothing lenght"""
        
        X, V, material, h, mass, density, pressure, energy \
            = self.sph.GetVectors()
        np, nv = self.sph.GetParticleNumber()
        
        val = 0.0
        
        dlg = wx.TextEntryDialog(
                self, 'Enter the smoothing length',
                'Smoothing length h', 'SPH')
        
        dlg.SetValue("")
        
        if dlg.ShowModal() == wx.ID_OK:
            try:
                val = float(dlg.GetValue())
            except:
                dlg.SetValue("")
        
        dlg.Destroy()
        if val != 0.0:
            h = (np + nv) * [val]
            
            self.sph.SetVectors(X, V, material, h, \
                                mass, density, pressure, energy)
            
    
    def OnSetMaterial(self, event):
        """Set the particle material"""
        
        X, V, material, h, mass, density, pressure, energy \
            = self.sph.GetVectors()
        np, nv = self.sph.GetParticleNumber()
        
        val = 0
        
        dlg = wx.TextEntryDialog(
                self, 'Enter the particle material number',
                'Material number', 'SPH')
        
        dlg.SetValue("")
        
        if dlg.ShowModal() == wx.ID_OK:
            try:
                val = int(dlg.GetValue())
            except:
                dlg.SetValue("")
        
        dlg.Destroy()
        if val != 0:
            material = (np + nv) * [val]
            
            self.sph.SetVectors(X, V, material, h, \
                                mass, density, pressure, energy)
            
    
    def OnSetPressure(self, event):
        """Set the particle pressure"""
        
        X, V, material, h, mass, density, pressure, energy \
            = self.sph.GetVectors()
        np, nv = self.sph.GetParticleNumber()
        
        val = 0.0
        
        dlg = wx.TextEntryDialog(
                self, 'Enter the pressure',
                'Initial pressure', 'SPH')
        
        dlg.SetValue("")
        
        if dlg.ShowModal() == wx.ID_OK:
            try:
                val = float(dlg.GetValue())
            except:
                dlg.SetValue("")
        
        dlg.Destroy()
        if val != 0.0:
            pressure = (np + nv) * [val]
            
            self.sph.SetVectors(X, V, material, h, \
                                mass, density, pressure, energy)
            
        
    def OnSetEnergy(self, event):
        """Set the particle energy"""
        
        X, V, material, h, mass, density, pressure, energy \
            = self.sph.GetVectors()
        np, nv = self.sph.GetParticleNumber()
        
        val = 0.0
        
        dlg = wx.TextEntryDialog(
                self, 'Enter the energy',
                'Initial energy', 'SPH')
        
        dlg.SetValue("")
        
        if dlg.ShowModal() == wx.ID_OK:
            try:
                val = float(dlg.GetValue())
            except:
                dlg.SetValue("")
        
        dlg.Destroy()
        if val != 0.0:
            energy = (np + nv) * [val]
            
            self.sph.SetVectors(X, V, material, h, \
                                mass, density, pressure, energy)
            
    
    def OnSetMass(self, event):
        """Set the particle mass"""
        
        X, V, material, h, mass, density, pressure, energy \
            = self.sph.GetVectors()
        np, nv = self.sph.GetParticleNumber()
        
        val = 0.0
        
        dlg = wx.TextEntryDialog(
                self, 'Enter the mass',
                'Particle mass', 'SPH')
        
        dlg.SetValue("")
        
        if dlg.ShowModal() == wx.ID_OK:
            try:
                val = float(dlg.GetValue())
            except:
                dlg.SetValue("")
        
        dlg.Destroy()
        if val != 0.0:
            mass = (np + nv) * [val]
            
            self.sph.SetVectors(X, V, material, h, \
                                mass, density, pressure, energy)
            
    
    def OnSetDensity(self, event):
        """Set the particle density"""
        
        X, V, material, h, mass, density, pressure, energy \
            = self.sph.GetVectors()
        np, nv = self.sph.GetParticleNumber()
        
        val = 0.0
        
        dlg = wx.TextEntryDialog(
                self, 'Enter the density',
                'Initial density', 'SPH')
        
        dlg.SetValue("")
        
        if dlg.ShowModal() == wx.ID_OK:
            try:
                val = float(dlg.GetValue())
            except:
                dlg.SetValue("")
        
        dlg.Destroy()
        if val != 0.0:
            density = (np + nv) * [val]
            
            self.sph.SetVectors(X, V, material, h, \
                                mass, density, pressure, energy)
            
    
    def OnSetVelocity(self, event):
        """Set the particle velocity"""
        
        X, V, material, h, mass, density, pressure, energy \
            = self.sph.GetVectors()
        np, nv = self.sph.GetParticleNumber()
        
        val = []
        
        dlg = VectorEntryDialog(self, -1, 'Enter the velocity')
        
        dlg.SetValue(["", "", ""])
        
        if dlg.ShowModal() == wx.ID_OK:
            try:
                res = dlg.GetValue()
                val = [float(res[0]), float(res[1]), float(res[2])]
            except:
                dlg.SetValue(["", "", ""])
        
        dlg.Destroy()
        if val != []:
            V = (np + nv) * [val]
            
            self.sph.SetVectors(X, V, material, h, \
                                mass, density, pressure, energy)
            
    
class App(wx.App):
    def OnInit(self):
        frame = MainFrame(None, "Armando")
#        self.SetTopWindow(frame)
        frame.Show(True)
        return True
        
if __name__ == '__main__':
    app = App(0)
    app.MainLoop()
#app = MyApp(redirect=True)
#app.MainLoop()
