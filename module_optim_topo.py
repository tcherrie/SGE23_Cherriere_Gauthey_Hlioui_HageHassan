# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 11:30:20 2023

@author: cherriere
"""

# Module optim topo

from ngsolve import *
from ngsolve.webgui import Draw
from netgen.geom2d import CSG2d, Circle, Rectangle
from copy import copy
import numpy as np
import matplotlib.pyplot as plt

mu0 = 4e-7*np.pi

###############################################################################
# Maillage disque laminé

def meshLamDisk(Nbar, h = 1):
    geo = CSG2d()
    R=1
    x = R*2*(np.append(np.insert(np.arange(0.5,Nbar+0.5),0,0),Nbar)/Nbar-0.5)
    
    circle1 = Circle( center=(0,0), radius=R, bc="left_up" ) * Rectangle( pmin=(-R,0), pmax=(0,R))
    circle2 = Circle( center=(0,0), radius=R, bc="left_bot" ) * Rectangle( pmin=(-R,-R), pmax=(0,0))
    circle3 = Circle( center=(0,0), radius=R, bc="right_bot" ) * Rectangle( pmin=(0,-R), pmax=(R,0))
    circle4 = Circle( center=(0,0), radius=R, bc="right_up" ) * Rectangle( pmin=(0,0), pmax=(R,R))
    
    materials = ["iron","air"]
    
    for i in range(len(x)-1):
        geo.Add(Rectangle( pmin=(x[i],-R), pmax=(x[i+1],R), mat = materials[i%2] ) * (circle1 + circle2 + circle3 +circle4))

    # Attention à fixer la taille du maillage sinon le volume change à cause des elts grossiers
    m = geo.GenerateMesh(maxh=h)
    return Mesh(m)

###############################################################################
# Elements finis

def solvePrimal_linear(nu, hmoy ,mesh):
    # le champ 1 est vertical, le champ 2 est horizontal
    # on impose les forces magnétomotrices

    fespace_H1 = H1(mesh, order=1)
    fespace_H1.FreeDofs()[0] = False
    a = fespace_H1.TrialFunction()
    a_star = fespace_H1.TestFunction()
    K = BilinearForm(fespace_H1, symmetric=True)
    K +=  grad(a_star)* nu *grad(a) *dx
    
    # imposition des circulations sur le côté gauche et droit
    l1 = LinearForm(fespace_H1)
    l1 += -a_star* hmoy *sqrt(1-y*y)*ds(definedon=mesh.Boundaries("right_bot|right_up"))
    l1 += a_star* hmoy *sqrt(1-y*y)*ds(definedon=mesh.Boundaries("left_bot|left_up"))
    
    # imposition des circulations sur le haut et le bas
    l2 = LinearForm(fespace_H1)
    l2 += -a_star* hmoy * sqrt(1-x*x)* ds(definedon=mesh.Boundaries("right_bot|left_bot"))
    l2 += a_star* hmoy * sqrt(1-x*x)*ds(definedon=mesh.Boundaries("right_up|left_up"))

    K.Assemble()
    Kdec = K.mat.Inverse(inverse="sparsecholesky")
    
    a1 = GridFunction(fespace_H1)  # solution
    a1.vec.data =     Kdec * l1.Assemble().vec
    a2 = GridFunction(fespace_H1)  # solution
    a2.vec.data =     Kdec * l2.Assemble().vec
    
    return a1, a2


def solveDual_linear(mu, bmoy ,mesh):
    # le champ 1 est vertical, le champ 2 est horizontal
    # on impose les flux
        
    fespace_H1 = H1(mesh, order=1)
    fespace_H1.FreeDofs()[0] = False
    phi = fespace_H1.TrialFunction()
    phi_star = fespace_H1.TestFunction()
    K = BilinearForm(fespace_H1, symmetric=True)
    K +=  grad(phi_star)* mu *grad(phi) *dx
    
    # imposition du flux sur le côté haut et bas
    l1 = LinearForm(fespace_H1)
    l1 += -phi_star* bmoy * sqrt(1-x*x)* ds(definedon=mesh.Boundaries("right_bot|left_bot"))
    l1 += phi_star* bmoy * sqrt(1-x*x)*ds(definedon=mesh.Boundaries("right_up|left_up"))
    
    
    # imposition du flux sur le côté gauche et droit
    l2 = LinearForm(fespace_H1)
    l2 += phi_star* bmoy *sqrt(1-y*y)*ds(definedon=mesh.Boundaries("right_bot|right_up"))
    l2 += -phi_star* bmoy *sqrt(1-y*y)*ds(definedon=mesh.Boundaries("left_bot|left_up"))

    K.Assemble()
    Kdec = K.mat.Inverse(inverse="sparsecholesky")
    
    phi1 = GridFunction(fespace_H1)  # solution
    phi1.vec.data =     Kdec * l1.Assemble().vec
    phi2 = GridFunction(fespace_H1)  # solution
    phi2.vec.data =     Kdec * l2.Assemble().vec
    
    return phi1, phi2


###############################################################################
# Champs matériaux

def mu_defaut(mur,mesh):
    return mesh.MaterialCF({ "iron" : mur*mu0 }, default=mu0)

def rho_defaut(mesh):
    return mesh.MaterialCF({ "iron" : 1 }, default=0)

# Homogénéisation

def R(th):
    return CoefficientFunction( ( (cos(th), -sin(th)), (sin(th), cos(th)) ), dims = (2,2) )

def tR(th):
    return CoefficientFunction( ( (cos(th), sin(th)), (-sin(th), cos(th)) ), dims = (2,2) )

def muD(rh): 
    return mu0 * (rh * (mur - 1) + 1)

def muQ(rh):
    return mu0 * mur / (mur * (1 - rh) + rh)

def mu_xy(rh):
    return CoefficientFunction( ( (muD(rh),0),(0,muQ(rh)) ), dims = (2,2) )

def mu_star(rh,th):
    return R(th) * mu_xy(rh) * tR(th)



###############################################################################


# Fonction objectif

def compliance(potentiel, coeff, mesh):
    # équivalent à l'énergie / coénergie en linéaire
    return Integrate(grad(potentiel) * (coeff * grad(potentiel)), mesh)/2


def gradient_compliance(phi1, phi2, rh, th):
    
    mustar = mu_star(rh,th)
    Mm = mu_xy(rh)
    
    dmdrho = CoefficientFunction( ( ((mur-1)*mu0,0),(0,mur*mu0*(mur-1)/(mur+rh-mur*rh)**2 ) ), dims = (2,2) )
    dmudrho = R(th) * (dmdrho * tR(th))
    dmudtheta =  R(th+np.pi/2) * ( Mm * tR(th)) + R(th) * ( Mm * tR(th+np.pi/2))
    
    dJ_drho =  grad(phi1) * (dmudrho*grad(phi1) ) - grad(phi2) * (dmudrho*grad(phi2) )
    dJ_dtheta = grad(phi1) * (dmudtheta*grad(phi1) ) - grad(phi2) * (dmudtheta*grad(phi2) )

    return -dJ_drho, -dJ_dtheta

###############################################################################



###############################################################################

# Utilitaire

def a2b(a):
    return CoefficientFunction((grad(a)[1],-grad(a)[0]))


###############################################################################

# Affichage

def DrawVecField(vec,mesh):
    Draw(vec, mesh, vectors = { "grid_size":20})
    
def DrawMuStar(rh,th):
    Draw(rh*CoefficientFunction((cos(th),sin(th))),mesh,vectors = { "grid_size":20},min=0,max=1);