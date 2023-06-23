import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from monolayer_shared.analysis.ball_pivot import bp,BP_Params
from ovito.io import import_file
from ovito.modifiers import *
import quaternionic
from monolayer_shared.util.np import periodic_difference,normalized
import freud
import sys
import pandas as pd
import glob
import os
from tqdm.auto import tqdm

def calc_border_radius(box,pos,border):
    p = pos[border]
    r = periodic_difference(p-p[0],box)
    c = r.mean(axis=0)
    rc = r-c
    return np.linalg.norm(rc,axis=-1).max()

def calc_border_nor(nor,border):
    return normalized(nor[border].mean(axis=0))

def perim(box,pos,border):
    def dist(a,b):
        return np.linalg.norm(periodic_difference(a-b,box),axis=-1)
    if len(border)==0:
        return 0
    else:
        assert len(border)>2
        p = pos[border]
        return dist(p[1:],p[:-1]).sum()+dist(p[-1],p[0])

def surface(movie, frame):
    # Load movie
    pipeline = import_file(movie)
    pipeline.modifiers.append(SelectTypeModifier(operate_on = "particles", property = "Particle Type", types = {1,2,3}))
    pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=1.5, sort_by_size=True, compute_com=True))
    # Load frame
    data = pipeline.compute(frame)
    positions = np.array(data.particles['Position'])
    quaternions = np.array(data.particles['Orientation'])
    types = np.array(data.particles['Particle Type'])
    quaternions = quaternionic.array(data.particles['Orientation'])
    clusters = np.array(data.particles['Cluster'])
    Nog = len(quaternions)
    # Work only with vesicle
    if len(clusters[clusters == 1]) < 0.75*len(clusters):
        positions = positions[clusters <= 2]
        quaternions = quaternions[clusters <= 2]
        nvesicles = 2
    else:
        positions = positions[clusters == 1]
        quaternions = quaternions[clusters == 1]
        nvesicles = 1
    Nves = len(quaternions)
    # Transform quaternions to normals
    qux = quaternionic.array([0,1,0,0])
    normqs = quaternions*qux*np.conjugate(quaternions)
    norms = np.array(normqs)[:,1:]
    # Transform coordinates to match Miguel's code
    box = np.array([60,60,60])
    pos = positions+30
    nor = norms.copy()
    # Compute surface triangles
    tri_ind,tri_lb,tri_lb_n,borders,borders_ok = bp(box,pos,nor,
                                                    BP_Params(
                                                        radii=[1.0],
                                                        do_fill_holes=True,
                                                        do_remove=True,
                                                        seed_normal_alignment_max_angle=np.pi/2,
                                                        seed_normal_face_normal_max_angle=np.pi/2,
                                                        dihedral_max_angle=np.pi/2,
                                                     )
    )
    # Define triangle properties
    tri_pos = np.zeros((len(tri_ind),3))
    tri_nor = np.zeros((len(tri_ind),3))
    for i in range(len(tri_pos)):
    # for i in range(1):
        tri = pos[tri_ind[i]]
        trin = nor[tri_ind[i]]
        tri_pos[i] = tri.mean(0)
        tri_nor[i] = trin.mean(0)
        tri_nor[i] = tri_nor[i]/np.linalg.norm(tri_nor[i])
    # Define parameters of the ball-pivot method
    bead_radius = 1.0
    rho_radius = 2.0
    # Compute perimeters of the holes in the surface
    l2 = []
    b2 = []
    for border in borders:
    #         if len(border) < 6:
    #             continue
        if calc_border_radius(box,pos,border)<rho_radius:
            continue
        nz = calc_border_nor(nor,border)[-1]
        if np.abs(nz)>.75:
            lb = -1 if nz>0 else 1
        else:
            lb = 0
        perimeter = perim(box,pos,border)
        l2.append((lb,perimeter))
        b2.append(border)
    f_lb_per = np.array(l2)
    # print(f_lb_per)
    # print(b2)
    # print(pos[b2[0]])
    # print(pos[b2[0]]-30)
    # print(np.array(pos[b2[0]]-30).mean(0))
    # print('\n\n')
    # pcom = np.array(pos[b2[0]]-30).mean(0)
    # Compute radius and area of pore
    if len(f_lb_per) > 0:
        order = np.argsort(f_lb_per[:,1])[::-1]
        # print(order)
        f_lb_per = f_lb_per[order]
        # print(f_lb_per)
        perimeter = f_lb_per[0,1]
        radius = perimeter/(2*np.pi)
        area = np.pi*radius*radius
        pcom = np.array([pos[borders[order[0]]].mean(axis = 0)])
        # print(pos[borders[order[0]]])
        # print(pos[borders[order[0]]]-30)
        # print(pcom)
        # pcom = np.array([(pos[borders[order[0]]]-30).mean(axis = 0)])
        pcom = np.array([np.array(pos[b2[0]]-30).mean(0)])
        # print(pcom)
    else:
        perimeter = 0
        radius = 0
        area = 0
        pcom = np.array([[np.nan,np.nan,np.nan]])
        # print(pcom)
    # vcom = np.array([pos.mean(axis = 0)])
    vcom = np.array([(pos-30).mean(axis = 0)])
    pnor = pcom-vcom
    pnor = pnor/np.linalg.norm(pnor)
    return tri_pos, tri_nor, pcom, pnor, radius, nvesicles

def iotype(positions, types, interest_type, tri_pos, tri_nor, pcom, pnor):
    # Define points and query_points
    points = tri_pos - 30
    query_points = np.array([[0,0,0],[55,0,0]])
    query_points = positions[types == interest_type]
    # Define box and query
    box = freud.box.Box.cube(100)
    aq = freud.locality.AABBQuery(box, points)
    # Compute dot products
    dps = []
    # Compute list of neighbours
    bonds = np.array(aq.query(query_points, dict(num_neighbors=6)).toNeighborList())
    # Run through it
    for i in range(len(query_points)):
        lbonds = bonds[bonds[:,0] == i]
        assert len(lbonds) == 6
        dotps = []
        for lb in lbonds:
            dr = points[lb[1]]-query_points[lb[0]]
            dr = dr/np.linalg.norm(dr)
            dotp = np.dot(dr,tri_nor[lb[1]])
            dotps.append(dotp)
        dotpsm = np.mean(dotps)
        dps.append(dotpsm)
    dps = np.array(dps)
    # Reanalyse those particles close to the pore
    porepop = np.intersect1d(np.where(dps > -0.5)[0],np.where(dps < 0.5)[0])
    pore_dots = np.zeros(len(porepop))
    for i in range(len(porepop)):
        dr = query_points[porepop[i]]-(pcom[0,:]-30)
        dr = dr/np.linalg.norm(dr)
        dotp = np.dot(dr,pnor[0,:])
        pore_dots[i] = dotp
    # Define numbers
    Ninside = len(dps[dps >= 0.5])
    Noutside = len(dps[dps <= -0.5])
    Npore = len(porepop)
    Npout = len(pore_dots[pore_dots > 0])
    Npin = len(pore_dots[pore_dots <= 0])
    NAI = Ninside + Npin
    NAO = Noutside + Npout
    return NAI, NAO

def iodrops(drops, tri_pos, tri_nor, pcom, pnor):
    # Define points and query_points
    points = tri_pos - 30
    query_points = []
    for i in range(len(drops)):
        query_points.append([drops[i][1][0],drops[i][1][1],drops[i][1][2]])
    query_points = np.array(query_points)
    # Define box and query
    box = freud.box.Box.cube(100)
    aq = freud.locality.AABBQuery(box, points)
    # Compute dot products
    dps = []
    # Compute list of neighbours
    bonds = np.array(aq.query(query_points, dict(num_neighbors=6)).toNeighborList())
    # Run through it
    for i in range(len(query_points)):
        lbonds = bonds[bonds[:,0] == i]
        assert len(lbonds) == 6
        dotps = []
        for lb in lbonds:
            dr = points[lb[1]]-query_points[lb[0]]
            dr = dr/np.linalg.norm(dr)
            dotp = np.dot(dr,tri_nor[lb[1]])
            dotps.append(dotp)
        dotpsm = np.mean(dotps)
        dps.append(dotpsm)
    dps = np.array(dps)
    # Reanalyse those particles close to the pore
    porepop = np.intersect1d(np.where(dps > -0.5)[0],np.where(dps < 0.5)[0])
    pore_dots = np.zeros(len(porepop))
    for i in range(len(porepop)):
        dr = query_points[porepop[i]]-(pcom[0,:]-30)
        dr = dr/np.linalg.norm(dr)
        dotp = np.dot(dr,pnor[0,:])
        pore_dots[i] = dotp
    # Define numbers
    Ninside = len(dps[dps >= 0.5])
    Noutside = len(dps[dps <= -0.5])
    Npore = len(porepop)
    Npout = len(pore_dots[pore_dots > 0])
    Npin = len(pore_dots[pore_dots <= 0])
    NAI = Ninside + Npin
    NAO = Noutside + Npout
    return NAI, NAO

def insideout(movie, frame):
    # Triangles
    tri_pos, tri_nor, pcom, pnor, radius, nvesicles = surface(movie, frame)
    # Load movie
    pipeline = import_file(movie)
    pipeline.modifiers.append(SelectTypeModifier(operate_on = "particles", property = "Particle Type", types = {4}))
    pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=1.5, sort_by_size=True, compute_com=True, compute_gyration = True))
    # Load frame
    data = pipeline.compute(frame)
    positions = np.array(data.particles['Position'])
    types = np.array(data.particles['Particle Type'])
    clusters = np.array(data.particles['Cluster'])
    cluster_table = data.tables['clusters']
    comdrops = np.array(cluster_table['Center of Mass'][...])
    rogdrops = np.array(cluster_table['Radius of Gyration'][...])
    # OG numbers for each type
    N1OG = len(types[types == 1])
    N2OG = len(types[types == 2])
    N3OG = len(types[types == 3])
    # Delete droplets (critical size = 100)
    ccid = 0
    size = 1000
    csize = 100
    droplets = []
    while size > csize:
        ccid += 1
        size = len(clusters[clusters == ccid])
        comd = comdrops[ccid-1]
        rogd = rogdrops[ccid-1]
        if size > csize:
            droplets.append([size, comd, rogd])
    positions = positions[clusters >= ccid]
    types = types[clusters >= ccid]
    clusters = clusters[clusters >= ccid]
    # Free numbers for each type
    N1 = len(types[types == 1])
    N2 = len(types[types == 2])
    N3 = len(types[types == 3])
    # Droplet numbers for each type
    D1 = N1OG-N1
    D2 = N2OG-N2
    D3 = N3OG-N3
    # In/Out numbers for each type
    NI1, NO1 = iotype(positions, types, 1, tri_pos, tri_nor, pcom, pnor)
    NI2, NO2 = iotype(positions, types, 2, tri_pos, tri_nor, pcom, pnor)
    NI3, NO3 = iotype(positions, types, 3, tri_pos, tri_nor, pcom, pnor)
    Inside = np.array([NI1, NI2, NI3])
    Outside = np.array([NO1, NO2, NO3])
    # In/Out numbers of droplets
    if len(droplets) > 0:
        DI, DO = iodrops(droplets, tri_pos, tri_nor, pcom, pnor)
    else:
        DI = 0
        DO = 0
    Droplet = np.array([D1, D2, D3, len(droplets), DI, DO])
    return Inside, Outside, Droplet, radius, droplets, nvesicles

def droplets_and_pore_positions(movie, frame):
    # Triangles
    tri_pos, tri_nor, pcom, pnor, radius, nvesicles = surface(movie, frame)
    # print(pcom, pnor, radius)
    # input()
    # Load movie
    pipeline = import_file(movie)
    pipeline.modifiers.append(SelectTypeModifier(operate_on = "particles", property = "Particle Type", types = {4}))
    pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=1.5, sort_by_size=True, compute_com=True, compute_gyration = True))
    # Load frame
    data = pipeline.compute(frame)
    positions = np.array(data.particles['Position'])
    types = np.array(data.particles['Particle Type'])
    clusters = np.array(data.particles['Cluster'])
    cluster_table = data.tables['clusters']
    comdrops = np.array(cluster_table['Center of Mass'][...])
    rogdrops = np.array(cluster_table['Radius of Gyration'][...])
    # OG numbers for each type
    N1OG = len(types[types == 1])
    N2OG = len(types[types == 2])
    N3OG = len(types[types == 3])
    # Delete droplets (critical size = 100)
    ccid = 0
    size = 1000
    csize = 100
    droplets = []
    while size > csize:
        ccid += 1
        size = len(clusters[clusters == ccid])
        comd = comdrops[ccid-1]
        rogd = rogdrops[ccid-1]
        if size > csize:
            droplets.append([size, comd, rogd])
    pore_props = [pcom, pnor, radius]
    return pore_props, droplets

# def all_frames(movie):
#     pipeline = import_file(movie)
#     nframes = pipeline.source.num_frames
#     frames = np.arange(nframes+1)
#     I = np.zeros((len(frames),3))
#     O = np.zeros((len(frames),3))
#     D = np.zeros((len(frames),6))
#     R = np.zeros(len(frames))
#     N = np.zeros(len(frames))
#     for i, frame in enumerate(tqdm(frames)):
#         try:
#             I[i], O[i], D[i], R[i], drops, N[i] = insideout(movie, frame)
#         except ZeroDivisionError:
#             I[i], O[i], D[i], R[i], drops, N[i] = np.nan, np.nan, np.nan, np.nan, [], np.nan
#     return I, O, D, R, N, frames

def all_frames(movie):
    pipeline = import_file(movie)
    nframes = pipeline.source.num_frames
    # print(nframes)
    frames = np.arange(nframes+1)
    all_info = []
    # frames = np.arange(0,nframes+1,10)
    # frames = np.arange(11)
    # for i, frame in enumerate(frames):
    for i, frame in enumerate(tqdm(frames)):
        try:
            pore_props, droplets = droplets_and_pore_positions(movie, frame)
        except ZeroDivisionError:
            pore_props, droplets = [], []
        all_info.append([frame, pore_props, droplets])
        # print(i+1, frame)
        # print(pore_props)
        # print(droplets)
        # input()
    return all_info

def write_file(sim, all_info):
    f = open('%s_FrDist.txt'%(sim.split('/output.xyz')[0]), 'w')
    f.write('Frame PoreX PoreY PoreZ PorenX PorenY PorenZ PoreR Droplets[N,X,Y,Z,R]\n')
    for i in range(len(all_info)):
        frameinfo = all_info[i]
        if len(frameinfo[1]) > 0:
            f.write('%d\t%f %f %f %f %f %f %f\t'%(frameinfo[0], frameinfo[1][0][0][0], frameinfo[1][0][0][1], frameinfo[1][0][0][2], frameinfo[1][1][0][0], frameinfo[1][1][0][1], frameinfo[1][1][0][2], frameinfo[1][2]))
        else:
            f.write('%d\t%f %f %f %f %f %f %f\t'%(frameinfo[0], np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan))
        dropinfo = frameinfo[2]
        ntabs = 3
        for j in range(len(dropinfo)):
            f.write('%d %f %f %f %f\t'%(dropinfo[j][0], dropinfo[j][1][0], dropinfo[j][1][1], dropinfo[j][1][2], dropinfo[j][2]))
            ntabs += 1
        while ntabs < 52:
            f.write('\t')
            ntabs += 1
        f.write('\n')
    f.close()

bpath = '/nfs/scistore15/saricgrp/cvanhill/sims/SG_pore_condensation/vesicle/NR3.33'
# bpath = '/Users/christian/PhD/sims/mine/SG_pore_condensation/vesicle/IST/NR3.33'
p = int(sys.argv[1])
e = float(sys.argv[2])
m = float(sys.argv[3])
frac = float(sys.argv[4])

sims = glob.glob('%s/P%d/e%.1f_m%.2f_%.2ffrac_yes/sd*/output.xyz'%(bpath,p,e,m,frac))
nsims = len(sims)
# print(nsims)

for i, sim in enumerate(sims):
    if os.access('%s_FrDist.txt'%(sim.split('/output.xyz')[0]), os.F_OK):
        continue
    # if os.access('%s_FrData.txt'%(sim.split('/output.xyz')[0]), os.F_OK):
    #     continue
    # if os.access('%s/P%d/e%.1f_m%.2f_%.2ffrac_yes/FrData/%s_FrData.txt'%(bpath,p,e,m,frac,sim.split('/output.xyz')[0].split('_yes/')[1]), os.F_OK):
    #     continue
    print('  replica %d/%d'%(i+1,nsims))
    # I, O, D, R, N, frames = all_frames(sim)
    all_info = all_frames(sim)
    try:
        write_file(sim, all_info)
    except IndexError:
        print("Writing the file failed due to an IndexError: list index out of range")
        print(all_info)
        continue
    # r = os.system('cat %s_FrDist.txt'%(sim.split('/output.xyz')[0]))
    # input()
    # results = pd.DataFrame(index = frames, columns = ['P1_In','P1_Out','P1_Drop','P2_In','P2_Out','P2_Drop','P3_In','P3_Out','P3_Drop','Ndrops','InDrops','Nves','Prad'])
    # results['P1_In'] = I[:,0]
    # results['P1_Out'] = O[:,0]
    # results['P1_Drop'] = D[:,0]
    # results['P2_In'] = I[:,1]
    # results['P2_Out'] = O[:,1]
    # results['P2_Drop'] = D[:,1]
    # results['P3_In'] = I[:,2]
    # results['P3_Out'] = O[:,2]
    # results['P3_Drop'] = D[:,2]
    # results['Ndrops'] = D[:,3]
    # results['InDrops'] = D[:,4]
    # results['Nves'] = N
    # results['Prad'] = R
    # results.to_csv('%s_FrData.txt'%(sim.split('/output.xyz')[0]))

# if not os.access('%s/P%d/e%.1f_m%.2f_%.2ffrac_yes/FrData/%s_FrData.txt'%(bpath,p,e,m,frac,sim.split('/output.xyz')[0].split('_yes/')[1]), os.F_OK):
#     r = os.system('mkdir %s/P%d/e%.1f_m%.2f_%.2ffrac_yes/FrData/'%(bpath,p,e,m,frac))
# if len(glob.glob('%s/P%d/e%.1f_m%.2f_%.2ffrac_yes/sd*_FrData.txt'%(bpath,p,e,m,frac))) > 0:
#     r = os.system('mv %s/P%d/e%.1f_m%.2f_%.2ffrac_yes/sd*_FrData.txt %s/P%d/e%.1f_m%.2f_%.2ffrac_yes/FrData/'%(bpath,p,e,m,frac,bpath,p,e,m,frac))
